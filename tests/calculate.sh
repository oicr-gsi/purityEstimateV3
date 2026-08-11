#!/bin/bash
set -euo pipefail

# Metrics for a purityEstimateV3 run. Deliberately structural, not numeric: WISP's
# purity estimates are floating point and shift with tool versions, so pinning them
# would make every upgrade a test failure. What IS worth pinning is the SHAPE of the
# output, which has already caught real differences -- an MRD-negative sample makes
# WISP omit somatic_vaf.png, so its tarball holds 7 files where a positive holds 8.

cd "$1" || exit 1

echo "### provisioned files"
# `! -type d`, NOT `-type f`: Vidarr provisions outputs as SYMLINKS, so -type f matches
# nothing and this section silently records an empty list, which would then never detect a
# missing output. Verified against a real run.
find . -maxdepth 1 ! -type d -printf '%f\n' | sort

echo "### archive members"
for archive in *.tar.gz; do
    [ -e "$archive" ] || continue
    echo "-- $archive"
    # Sort: tar member order follows directory order and is not stable across runs.
    tar -tzf "$archive" | sed 's|/$||' | sort
done

echo "### wisp summary columns"
for summary in *.wisp_summary.tsv; do
    [ -e "$summary" ] || continue
    echo "-- $summary"
    # Header only. Column names catch a WISP schema change (LOH columns appearing, a
    # renamed field); the values underneath are expected to drift.
    head -1 "$summary" | tr '\t' '\n'
    echo "rows: $(( $(wc -l < "$summary") - 1 ))"
done
