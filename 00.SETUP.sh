#!/bin/bash
#
# 00.SETUP.sh -- build the repo-local R library once per checkout
#
# Usage, from anywhere:
#    seqCNA/00.SETUP.sh
#
# R 4.5.x is taken from PATH (on IRIS: module load R/4.5.1_nds). The
# common packages DNAcopy (Bioconductor), pctGCdata (CRAN) and RJSONIO
# (CRAN) must already be in R's library path, shared install or your
# personal library; this script checks and stops with install hints if
# any is missing. Only seqDNAcopy, the fork in the submodule, is built
# here, into Rlib/ (gitignored). Safe to rerun; it rebuilds.
#

set -eu

SDIR=$( cd "$( dirname "$0" )" && pwd )
cd $SDIR

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

command -v Rscript >/dev/null 2>&1 || fatal "Rscript not on PATH (module load R/4.5.1_nds)"

RVERSION=$(Rscript --vanilla -e 'cat(as.character(getRversion()))')
case $RVERSION in
    4.5.*) ;;
    *) fatal "R $RVERSION on PATH; need 4.5.x (module load R/4.5.1_nds)" ;;
esac
echo "R $RVERSION: $(command -v Rscript)"

if [ ! -e seqDNAcopy/DESCRIPTION ]; then
    fatal "seqDNAcopy/ submodule is empty; run: git submodule update --init"
fi

#
# Common packages: present already?
#

Rscript --vanilla -e '
pkgs <- c("DNAcopy", "pctGCdata", "RJSONIO")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) {
    cat("\n    FATAL ERROR: R packages not installed:", missing, "\n\n")
    cat("    Install into your personal library with:\n")
    if ("DNAcopy" %in% missing)
        cat("        BiocManager::install(\"DNAcopy\")\n")
    if ("pctGCdata" %in% missing)
        cat("        install.packages(\"pctGCdata\")\n")
    if ("RJSONIO" %in% missing)
        cat("        install.packages(\"RJSONIO\")\n")
    cat("\n")
    quit(status=1)
}
for (p in pkgs) cat(sprintf("%-10s %s\n", p, as.character(packageVersion(p))))
'

#
# seqDNAcopy from the submodule into Rlib/
#

mkdir -p Rlib
R CMD INSTALL --clean --preclean -l Rlib seqDNAcopy

Rscript --vanilla -e '
library(seqDNAcopy, lib.loc="Rlib")
cat("seqDNAcopy", as.character(packageVersion("seqDNAcopy", lib.loc="Rlib")), "loads\n")
'

#
# Record what was built against
#

Rscript --vanilla -e '
.libPaths(c("Rlib", .libPaths()))
pkgs <- c("seqDNAcopy", "DNAcopy", "pctGCdata", "RJSONIO", "Rsamtools", "IRanges", "Rcpp")
cat("R", as.character(getRversion()), R.home(), "\n")
for (p in pkgs)
    cat(p, as.character(packageVersion(p)), dirname(system.file(package=p)), "\n")
cat("seqDNAcopy git:", system("git -C seqDNAcopy describe --always --dirty", intern=TRUE), "\n")
cat("built:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
' > Rlib/VERSIONS.txt

cat Rlib/VERSIONS.txt
echo
echo "Rlib/ ready"
