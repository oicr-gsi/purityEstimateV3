#!/usr/bin/env bash
# make_test_wdl.sh -- write a copy of the workflow that loads modules from a tree the
# scheduler does not know about yet.
#
#     scripts/make_test_wdl.sh <module_path>... -- <out.wdl> [src.wdl]
#
# Cromwell loads the runtime `modules` attribute in the qsub prolog, before the task script
# runs, so a module outside MODULEPATH fails there and the task dies with no output at all --
# which the engine reports as a synthetic rc=79. Setting MODULEPATH in the submitting shell
# does not help: Cromwell is a server and inherits the environment it was started with.
#
# The copy this writes clears the runtime attribute and loads the same modules inside each
# task instead, where `module use` can point at the tree first. Each task keeps its OWN module
# list, so a task that only wants samtools does not gain the rest.
#
# Delete the copy once the modules are published; nothing in the workflow itself needs undoing.

set -euo pipefail

tmp=""
trap 'rm -f "${tmp}"' EXIT

usage() {
    echo "usage: $(basename "$0") <module_path>... -- <out.wdl> [src.wdl]" >&2
    exit 2
}

main() {
    local -a paths=()
    while [ "$#" -gt 0 ] && [ "$1" != "--" ]; do paths+=("$1"); shift; done
    [ "${1:-}" = "--" ] || usage
    shift
    [ "$#" -ge 1 ] || usage
    [ "${#paths[@]}" -ge 1 ] || usage

    local here src out
    here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    out="$1"
    src="${2:-${here}/../purityEstimateV3.wdl}"

    [ -f "$src" ] || { echo "no workflow at $src" >&2; exit 1; }
    if [ "$(readlink -f "$src")" = "$(readlink -f "$out")" ]; then
        echo "refusing to write over the source workflow: $src" >&2
        exit 1
    fi

    local use_lines=""
    local p
    for p in "${paths[@]}"; do use_lines="${use_lines}        module use ${p}\n"; done

    tmp="$(mktemp "${out}.XXXXXX")"

    # A task's `modules` default is declared before its command block, so remembering the
    # most recent one is enough to emit the right load for the block that follows.
    awk -v uses="$use_lines" '
        # A task without a modules input must not inherit the previous task s; one that
        # declares modules with no default takes the workflow-level value, which is how the
        # head jobs get theirs.
        /^task / { mods = ""; intask = 1 }
        /^ *String +modules += +"/ {
            mods = $0
            sub(/^[^"]*"/, "", mods)
            sub(/".*$/, "", mods)
            if (!intask) wf_mods = mods
            sub(/= *"[^"]*"/, "= \"\"")
            print
            next
        }
        intask && /^ *String +modules *$/ { mods = wf_mods; print; next }
        { print }
        /^    command <<</ { inblock = 1; next }
        inblock && /^ *set -euo pipefail$/ {
            if (mods != "") {
                printf "%s", uses
                print "        module load " mods
            }
            inblock = 0
        }
    ' "$src" > "$tmp"

    local blocks loads
    blocks=$(grep -cE '^ *String +modules += +""' "$tmp" || true)
    loads=$(grep -c '^        module load ' "$tmp" || true)
    echo "wrote ${out}: ${blocks} module attribute(s) cleared, ${loads} in-task load(s)"
    mv "$tmp" "$out"
}

main "$@"
