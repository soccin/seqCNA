#!/bin/bash

#
# getRelease.sh -- print the release of this checkout, the
# "- **Version:** vX.Y.Z" line of VERSION.md
#
# The release is the client-facing version: it is set by hand in
# VERSION.md when a release is made, with a git tag of the same name,
# never derived from git. PIPE.sh records it as RELEASE= in
# 00.RUNSTATUS.txt; delivery.sh puts it in the email and the METHODS.md
# link. From ATAC-seq's bin/getRelease.sh.
#

set -euo pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"
VFILE=$(dirname "$SDIR")/VERSION.md

RELEASE=$(sed -n 's/^- \*\*Version:\*\* *//p' "$VFILE" 2>/dev/null || true)
if [[ ! "$RELEASE" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "getRelease.sh: no '- **Version:** vX.Y.Z' line in $VFILE" >&2
    exit 1
fi

echo "$RELEASE"
