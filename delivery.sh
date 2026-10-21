#!/bin/bash

#
# delivery.sh -- deliver a finished seqCNA run to the BIC delivery area
#
# usage:
#   seqCNA/delivery.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]
#   seqCNA/delivery.sh [-n] -d|--default
#
# Run from the analysis directory, after the run has COMPLETED
# (bin/checkRun.sh exits 0). It
#
#   - checks the deliverables in <ProjNo>/ and that none of the tables
#     names a normal (pairing.tsv);
#   - writes project.yaml for the bicdelivery permissions database,
#     from the run and the first README.txt found going up from here
#     to Proj_<no> (bin/readme2yaml.R);
#   - copies the deliverables and the PDF and HTML of docs/methods.md
#     and docs/output.md (docs/MAKE_DOCS.sh) of this checkout to
#     DELIVERY/r_NNN/seqcna. The .md files, the PNGs (titled with the
#     pair ID) and pairing.tsv are not delivered;
#   - writes the delivery email to deliveryEmail_<no>_r_NNN_<yymmdd>.txt,
#     with the release of VERSION.md of this checkout;
#   - registers the run with the permissions database.
#
#   -d, --default   deliver to DROOT/<pi>/<invest>/Proj_<no>, taken
#                   from the path .../Users/<Pi>/<Invest>/Proj_<no>/...
#                   of the analysis directory (lowercased but for the
#                   project)
#   -n, --dry-run   print the commands that write to the delivery area
#                   or the database instead of running them. The email
#                   and project.yaml are still written here.
#
# Without a trailing r_NNN the run folder is the one after the highest
# r_NNN in the delivery folder, else r_001. An existing r_NNN can be
# given to add seqcna/ to it; nothing in it is deleted.
#
# DROOT is /data1/core002/res/bic/results. Every write to the delivery
# area runs as the service account svc_core001_bot01 (dzdo). From
# ATAC-seq's deliverResults.sh.
#

set -euo pipefail

#
# Physical paths: the service account copies from here and cannot
# enter a home directory that a symlink in the path may pass through.
#
SDIR="$( cd "$( dirname "$0" )" && pwd -P )"
CWD=$(pwd -P)

BOT_USER=svc_core001_bot01
DROOT=/data1/core002/res/bic/results
BIC_DELIVERY=$HOME/Code/BIC/Delivery/Version2j
DFOLDER=seqcna

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

WARNINGS=()
warn() {
    echo "    WARNING: $*"
    WARNINGS+=("$*")
}

#
# Default delivery path from .../Users/Aa/BBB/Proj_nnnnn/...
# -> aa/bbb/Proj_nnnnn
#
DPATH=""
if [[ "$CWD" == *"/Users/"* ]]; then
    rest="${CWD#*/Users/}"
    IFS=/ read -r u1 u2 proj _ <<< "$rest"
    if [[ -n "$u1" && -n "$u2" && "$proj" == Proj_* ]]; then
        DPATH="${u1,,}/${u2,,}/$proj"
    fi
fi

usage() {
    echo
    echo "usage: delivery.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]"
    echo "       delivery.sh [-n] -d|--default"
    echo
    echo "    -n, --dry-run   print the delivery and database commands, do not run them"
    echo
    echo "    default: ${DPATH:+$DROOT/}${DPATH:-(none)}"
    echo
}

ODIR=""
USE_DEFAULT=0
DRYRUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--default)
            USE_DEFAULT=1
            ;;
        -n|--dry-run)
            DRYRUN=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [[ -n "$ODIR" ]]; then
                echo "Only one delivery folder may be given"
                usage
                exit 1
            fi
            ODIR="$1"
            ;;
    esac
    shift
done

if [[ "$USE_DEFAULT" -eq 1 && -n "$ODIR" ]]; then
    echo "Give a delivery folder or -d, not both"
    usage
    exit 1
fi

if [[ -z "$ODIR" && "$USE_DEFAULT" -eq 0 ]]; then
    usage
    exit 1
fi

asbot() {
    if [[ "$DRYRUN" -eq 1 ]]; then
        echo "DRYRUN: dzdo -u $BOT_USER ${*@Q}"
    else
        dzdo -u "$BOT_USER" "$@"
    fi
}

##
# The run

STATUSF=00.RUNSTATUS.txt
if [ ! -e "$STATUSF" ]; then
    fatal "no $STATUSF; run delivery.sh from the analysis directory"
fi

if ! "$SDIR/bin/checkRun.sh"; then
    fatal "the run has not COMPLETED (bin/checkRun.sh above); nothing delivered"
fi
echo

getKey() {
    sed -n "s/^$1=//p" "$STATUSF" | tail -1
}

PROJ=$(getKey PROJNO)
GENOME_TAG=$(getKey GENOME_TAG)
for key in PROJ GENOME_TAG; do
    if [[ -z "${!key}" ]]; then
        fatal "no ${key/PROJ/PROJNO}= in $STATUSF"
    fi
done
if [[ "$PROJ" != Proj_* ]]; then
    fatal "PROJNO=$PROJ in $STATUSF is not Proj_<no>"
fi
PROJNO=${PROJ#Proj_}

#
# The assembly names init_impact_project_permissions.py accepts
#
case "$GENOME_TAG" in
    b37|b37_dmp)  ASSEMBLY=b37 ;;
    hg19)         ASSEMBLY=hg19 ;;
    b37+mm10)     ASSEMBLY=hybrid ;;
    mm10|GRC_m38) ASSEMBLY=mm10 ;;
    *)            fatal "no database assembly for GENOME_TAG=$GENOME_TAG" ;;
esac

if [[ "$USE_DEFAULT" -eq 1 ]]; then
    if [[ -z "$DPATH" ]]; then
        fatal "no default delivery path for $CWD (not .../Users/<Pi>/<Invest>/Proj_<no>/...)"
    fi
    if [[ "${DPATH##*/}" != "$PROJ" ]]; then
        fatal "the path names ${DPATH##*/} but the run is $PROJ; give the delivery folder"
    fi
    ODIR="$DROOT/$DPATH"
fi

#
# The email and its methods.md and output.md links name the release of
# VERSION.md of this checkout; it changes only when a release is made,
# never with git.
#
if ! VERSION=$("$SDIR/bin/getRelease.sh"); then
    fatal "no release in $SDIR/VERSION.md"
fi

#
# The docs go out as PDF and HTML, rendered from the .md by
# docs/MAKE_DOCS.sh and committed with them.
#
DOCS=()
for doc in methods output; do
    for ext in pdf html; do
        DOCS+=("$SDIR/docs/$doc.$ext")
    done
done
for f in "${DOCS[@]}"; do
    if [ ! -s "$f" ]; then
        fatal "missing $f (docs/MAKE_DOCS.sh)"
    fi
done

##
# Deliverables: the IGV .seg and the PDF, then either the three call
# tables or getGeneCalls' README.txt (nothing significant)

PDIR=$CWD/$PROJ
SEG=$PDIR/${PROJ}___IGV.seg
PDF=$PDIR/${PROJ}___seqSeg.pdf
TABLES=("$PDIR/GeneMatrix.csv" "$PDIR/SegmentMatrix.csv" "$PDIR/${PROJ}___GeneTable.xlsx")
NOCALLS=$PDIR/README.txt

for f in "$SEG" "$PDF"; do
    if [ ! -s "$f" ]; then
        fatal "missing or empty $f"
    fi
done

if [ -e "$NOCALLS" ]; then
    for f in "${TABLES[@]}"; do
        if [ -e "$f" ]; then
            fatal "both $NOCALLS (no significant calls) and $f"
        fi
    done
    HAS_CALLS=0
    FILES=("$SEG" "$PDF" "$NOCALLS")
else
    for f in "${TABLES[@]}"; do
        if [ ! -s "$f" ]; then
            fatal "missing or empty $f (and no $NOCALLS)"
        fi
    done
    HAS_CALLS=1
    FILES=("$SEG" "$PDF" "${TABLES[@]}")
fi

#
# No table may name a normal: the pipeline replaced each pair ID
# <tumor>__<normal> by the tumor ID (redactDeliverables.R). Check for
# __<normal> of every normal in pairing.tsv; the xlsx is a zip of XML.
# The PDF is built from the redacted plots and cannot be searched.
#
PAIRING=$CWD/pairing.tsv
if [ ! -s "$PAIRING" ]; then
    fatal "missing $PAIRING; was this run made before the normals were redacted?"
fi
mapfile -t NORMALS < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="NORMAL") c=i; next} c{print $c}' "$PAIRING" | sort -u)
if [[ ${#NORMALS[@]} -eq 0 ]]; then
    fatal "no NORMAL column or no rows in $PAIRING"
fi
PATTERNS=$(printf '__%s\n' "${NORMALS[@]}")

LEAKS=()
for f in "${FILES[@]}"; do
    case "$f" in
        *.pdf)
            continue
            ;;
        *.xlsx)
            if ! xml=$(unzip -p "$f"); then
                fatal "cannot read $f"
            fi
            if grep -qF "$PATTERNS" <<< "$xml"; then
                LEAKS+=("$f")
            fi
            ;;
        *)
            if grep -qF "$PATTERNS" "$f"; then
                LEAKS+=("$f")
            fi
            ;;
    esac
done
if [[ ${#LEAKS[@]} -gt 0 ]]; then
    fatal "these files name a normal (__<normal> of pairing.tsv): ${LEAKS[*]}"
fi

##
# Run folder: a trailing r_NNN as given, else the next one in ODIR

ODIR=$(realpath -m "$ODIR")
if [[ "$ODIR" =~ /(r_[0-9]+)$ ]]; then
    CURRDIR=${BASH_REMATCH[1]}
    ODIR=${ODIR%/*}
else
    CURRDIR=""
fi

#
# The deepest existing folder of ODIR must be readable by us; else
# the run folders in it could exist unseen, and r_001 would be picked
# or files added to a run folder without the note below.
#
seen=$ODIR
while [[ -n "$seen" && ! -e "$seen" ]]; do
    seen=${seen%/*}
done
seen=${seen:-/}
if [[ ! -r "$seen" || ! -x "$seen" ]]; then
    fatal "cannot read $seen, so cannot tell what is in $ODIR"
fi

if [[ -z "$CURRDIR" ]]; then
    last=0
    for d in "$ODIR"/r_[0-9]*; do
        n=${d##*/r_}
        if [[ -d "$d" && "$n" =~ ^[0-9]+$ ]] && ((10#$n > last)); then
            last=$((10#$n))
        fi
    done
    CURRDIR=$(printf 'r_%03d' $((last + 1)))
fi
DEST=$ODIR/$CURRDIR

echo "PROJECT=$PROJNO"
echo "GENOME_TAG=$GENOME_TAG (database: $ASSEMBLY)"
echo "VERSION=$VERSION"
echo "ODIR=$ODIR"
echo "CURRDIR=$CURRDIR"
if [[ -e "$DEST/$DFOLDER" ]]; then
    echo "NOTE: $DEST/$DFOLDER exists; files are added to it"
elif [[ -e "$DEST" ]]; then
    echo "NOTE: $DEST exists; $DFOLDER is added to it"
fi
echo

##
# project.yaml first: a missing or unreadable README.txt stops the
# delivery before anything is copied.

"$SDIR/bin/rscript.sh" "$SDIR/bin/readme2yaml.R" "$PROJNO" "$ODIR" "$CURRDIR" "$ASSEMBLY"
echo

#
# The database takes pi and invest from README.txt; the delivery folder
# is usually named after the same two users.
#
YAML_PI=$(sed -n 's/^pi: //p' project.yaml)
YAML_INVEST=$(sed -n 's/^invest: //p' project.yaml)
if [[ "$ODIR" != "$DROOT/$YAML_PI/$YAML_INVEST/$PROJ" ]]; then
    warn "$ODIR is not $DROOT/<pi>/<invest>/$PROJ for pi=$YAML_PI invest=$YAML_INVEST of README.txt"
fi

##
# Copy

asbot mkdir -p "$DEST/$DFOLDER"
asbot rsync -aP "${FILES[@]}" "${DOCS[@]}" "$DEST/$DFOLDER/"

##
# Email

EMAIL=deliveryEmail_${PROJNO}_${CURRDIR}_$(date +%y%m%d).txt

#
# A value as the replacement of sed s/// (git tags may contain / and &)
#
sedValue() {
    local v=${1//\\/\\\\}
    v=${v//\//\\/}
    printf '%s' "${v//&/\\&}"
}

fillTemplate() {
    sed -e "s/{PROJNO}/$(sedValue "$PROJNO")/g" \
        -e "s/{RUNDIR}/$(sedValue "$CURRDIR")/g" \
        -e "s/{VERSION}/$(sedValue "$VERSION")/g" \
        "$SDIR/$1"
}

echo
echo "========================================================================="
echo
{
    fillTemplate tmpldeliveryEmail_01.txt
    if [[ "$HAS_CALLS" -eq 1 ]]; then
        fillTemplate tmpldeliveryEmail_02.txt
    else
        fillTemplate tmpldeliveryEmail_03.txt
    fi
    fillTemplate tmpldeliveryEmail_99.txt
} | tee "$EMAIL"
echo "========================================================================="
echo

##
# Permissions database

PERMS=("$BIC_DELIVERY/authorization_db/init_impact_project_permissions.py" -p project.yaml)

#
# Run with set +eu: the module scripts reference unset variables, and
# a failure here must say what was already done.
#
registerRun() {
    set +eu
    if [ "$(type -t module)" != "function" ]; then
        . /etc/profile.d/modules.sh || return 1
    fi
    module purge || return 1
    module load python/3.8.0 || return 1
    local python38
    python38=$(which python3.8) || return 1
    module load py-python-ldap/3.4.2 || return 1
    "$python38" "${PERMS[@]}"
}

if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRYRUN: module load python/3.8.0 py-python-ldap/3.4.2"
    echo "DRYRUN: python3.8 ${PERMS[*]@Q}"
else
    rc=0
    registerRun || rc=$?
    set -eu
    if [[ $rc -ne 0 ]]; then
        echo
        echo "    FATAL ERROR: the files are delivered to $DEST/$DFOLDER and the"
        echo "    email is in $EMAIL, but registering the run with the permissions"
        echo "    database failed (exit $rc). Fix the cause and run, here:"
        echo
        echo "        module load python/3.8.0 py-python-ldap/3.4.2"
        echo "        python3.8 ${PERMS[*]@Q}"
        echo
        exit 1
    fi
fi

echo
if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRY RUN: would deliver to $DEST/$DFOLDER; nothing was written there or to the database"
else
    echo "Delivered to $DEST/$DFOLDER"
fi
echo "Email: $EMAIL"
for w in "${WARNINGS[@]}"; do
    echo "    WARNING: $w"
done
echo
