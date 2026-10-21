#!/bin/bash
#
# Render the client-facing docs (output.md, methods.md) to PDF and HTML
# with pandoc; delivery.sh delivers the PDF and HTML, not the .md. Run
# after editing either .md and commit the four renders with it. Needs
# pandoc and xelatex with the DejaVu fonts.
#
# The level 1 heading of each .md becomes the document title
# (--shift-heading-level-by=-1). The PDF page header gives the title and
# the release of VERSION.md. The HTML uses pandoc's built-in style, so it
# loads nothing from the web.
#

set -euo pipefail

cd "$(dirname "$0")"

RELEASE=$(../bin/getRelease.sh)

#
# render MD TITLE [PANDOC_OPTS ...]
#
render() {
    local MD=$1
    local TITLE=$2
    local BASE=${MD%.md}
    shift 2

    pandoc $MD -o $BASE.pdf "$@" \
        --shift-heading-level-by=-1 \
        --pdf-engine=xelatex \
        -V sansfont="DejaVu Sans" \
        -V monofont="DejaVu Sans Mono" \
        -V geometry:margin=.75in \
        -V fontsize=12pt \
        -V header-includes="\\usepackage{fancyhdr} \\pagestyle{fancy} \\fancyhead[L]{$TITLE} \\fancyhead[R]{$RELEASE}"

    pandoc $MD -o $BASE.html "$@" \
        --shift-heading-level-by=-1 \
        --standalone
}

render output.md "seqCNA Output" --number-sections
render methods.md "seqCNA Methods"
