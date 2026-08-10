#!/bin/bash
set -euo pipefail

# Plain diff, not the usual `diff <(sort) <(sort)`: calculate.sh already sorts within
# each section, so its output is deterministic. Sorting the whole file again would
# scramble the section headers into the listings and turn a readable structural diff
# into an unreadable set difference.
diff -s "$1" "$2"
