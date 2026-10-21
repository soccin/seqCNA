#!/bin/bash

#
# seqCNA control job for IRIS/Slurm. Run it from the project directory
# (the one holding pipeline/, tumors, normals). Every stage is fanned out
# with ~/bin/bsub and the control job blocks between stages with
# bSync/bCheck. The control job itself can run for hours on full-size
# BAMs, so the directives below put it on the long partition with
# --qos=priority. Flags given to sbatch on the command line override
# them.
#
# First, with no tumors or normals file yet, run it directly:
#
#    ./seqCNA/PIPE.sh TARGETTED
#
# It writes samples, one line per BAM, and stops (see "First step"
# below). Make tumors and normals from those lines, then:
#
# CMD:
#    sbatch ./seqCNA/PIPE.sh TARGETTED      # IMPACT, exome: BINSIZE=100
#    sbatch ./seqCNA/PIPE.sh sWGS           # shallow WGS:   BINSIZE=auto
#
# Then, at any time, from the same directory:
#
#    ./seqCNA/bin/checkRun.sh
#
# which exits 0 if the run worked, 1 if it failed, 2 if it is still
# running and 3 if sacct failed and the state cannot be read now.
#
# One run per project directory. A run started in a directory that
# already holds one (see checkPreviousRun) stops at once and changes
# nothing; rerun in a new directory. The exception is a run that failed
# in CHECK, before any job was submitted: fix the input and submit again
# in the same directory.
#
# sbatch here is the ~/bin wrapper. Slurm runs a copy of the batch script
# from its spool directory, so $0 cannot locate this checkout; the wrapper
# exports SBATCH_SCRIPT_DIR and SDIR is taken from that. Running PIPE.sh
# directly (no sbatch) falls back to $0.
#
# The full run can also be started directly, in a shell that outlives it
# (an interactive job, tmux, nohup):
#
#    ./seqCNA/PIPE.sh TARGETTED
#
# Then this script, not a Slurm job, is the control job: the run records
# hold CTRL_JOBID=none with its host and pid, and checkRun.sh looks for
# the pid there. That holds inside an interactive job too, where
# SLURM_JOB_ID is set but names the session, not this run (see RUN_MODE).
#
# Stages (job names qSeqCNA_<STAGE>_<pid of this script>):
#
#   CHECK   control job: inputs, genome, tumors/normals, project id
#   FIX     one job per BAM: fixChromosomeNames.sh      SLURM.01.FIX
#   COUNT   one job per tumor x normal: getPairedCounts  SLURM.02.COUNT
#   SEG     one job per pair: seqSegment                 SLURM.03.SEG
#   BEST    control job: selectBestMatch (pairing.tsv), outAll/, prune out/
#   POST    one job: postProcess.sh                      SLURM.04.POST
#           (deliverables name samples by tumor id only)
#   STAGING control job: sweep of every job, deliverables, COMPLETED
#
# Ported from the ChIP-seq control job (v1.0.0); see docs/SLURM_PORT.md.
#
#SBATCH -p cmobic_cpu
#SBATCH --qos=priority
#SBATCH -t 3-00:00:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -J CTRL.seqCNA
#SBATCH -o SLURM.CTRL/%j.out

set -eu

#
# RUN_MODE: sbatch if this process is the batch script of Slurm job
# $SLURM_JOB_ID, direct otherwise. SLURM_JOB_ID alone cannot tell: a run
# started by hand inside an interactive job inherits the session's id.
# Slurm runs a batch script as <SlurmdSpoolDir>/job<id>/slurm_script, so
# $0 can (checked on IRIS, Slurm 25.11). In direct mode SBATCH_SCRIPT_DIR
# is ignored; one inherited from the session would point elsewhere.
#
if [[ -n ${SLURM_JOB_ID:-} && $0 == */job$SLURM_JOB_ID/slurm_script ]]; then
    RUN_MODE=sbatch
    SDIR=${SBATCH_SCRIPT_DIR:-$( cd "$( dirname "$0" )" && pwd )}
else
    RUN_MODE=direct
    SDIR=$( cd "$( dirname "$0" )" && pwd )
fi

if [ ! -e "$SDIR/getPairedCounts" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo "    submit with ~/bin/sbatch, which exports SBATCH_SCRIPT_DIR"
    echo
    exit 1
fi

source $SDIR/bin/loadTools.sh
source $SDIR/bin/slurmTools.sh

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long --tags --dirty)

#
# The client-facing release, set by hand in VERSION.md; VERSION (git)
# is the exact code.
#
if ! SCRIPT_RELEASE=$($SDIR/bin/getRelease.sh); then
    echo
    echo "    FATAL ERROR: no release in $SDIR/VERSION.md"
    echo
    exit 1
fi

#
# Tool versions, recorded in the run status. R and the locked packages are
# checked against Rlib.lock.tsv in CHECK (bin/rscript.sh --check);
# samtools and ImageMagick are only recorded.
#
R_VERSION=$(Rscript --vanilla -e 'cat(as.character(getRversion()))')
RLIB_LOCK_MD5=$(md5sum <$SDIR/Rlib.lock.tsv | cut -d' ' -f1)
SAMTOOLS_VERSION=$(samtools --version | sed -n '1s/^samtools //p')
MAGICK_VERSION=$(magick -version | sed -n '1s/^Version: ImageMagick \([^ ]*\).*/\1/p')
PIPENAME="seqCNA"

ORIG_CMD="$SDIR/PIPE.sh $(printf '%q ' "$@")"

##
# Process command args

function usage {
    cat <<-END_USAGE

	usage: sbatch seqCNA/PIPE.sh [sWGS|TARGETTED]
	version=$SCRIPT_VERSION

	    sWGS        Shallow Whole Genome Seq (BINSIZE=auto)
	    TARGETTED   Targeted assays, IMPACT/EXOME/... (BINSIZE=100)

	Run from the project directory, which must hold:

	    pipeline/<alignments|bam>/*.bam   the BAMs, each with its index
	    tumors, normals                   one pattern per line, matched
	                                      (fgrep -f) against the relabeled
	                                      BAM paths, bamRelabel/<base>/
	                                      <base>__relabel.bam

	Without tumors or normals it writes samples, one line per BAM, and
	stops; copy those lines into tumors and normals.

	The genome is read from the BAM headers (b37 or mm10 families) and
	picks the gene annotation assay: Exome or Mouse-Twist. For other
	assays run postProcess.sh by hand afterwards.

	END_USAGE
    exit 1
}

if [ "$#" -ne 1 ]; then
    usage
fi

case $1 in
    -h|--help)
        usage
        ;;
esac

ASSAY_TYPE=$(echo $1 | tr '[a-z]' '[A-Z]')
case $ASSAY_TYPE in
    SWGS)
        BINSIZE=auto
        ;;
    TARGETTED)
        BINSIZE=100
        ;;
    *)
        echo
        echo "    FATAL ERROR: unknown assay type [$1]"
        usage
        ;;
esac

#
# findBams
#
# Set BAMDIR to pipeline/<alignments|bam> and BAMS to its BAMs; either is
# empty if not found.
#
findBams() {
    BAMDIR=$(ls -d pipeline/* 2>/dev/null | egrep "alignments|bam" | head -1) || true
    BAMS=""
    if [ -n "$BAMDIR" ] && [ -d "$BAMDIR" ]; then
        BAMS=$(ls $BAMDIR/*.bam 2>/dev/null) || true
    fi
}

#
# relabeledBam BAM
#
# Path of the relabeled copy of BAM written by the FIX stage; the base
# name is derived as in fixChromosomeNames.sh. tumors and normals are
# matched against these paths.
#
relabeledBam() {
    local base
    base=$(basename $1 | sed 's/.bam//')
    echo bamRelabel/$base/${base}__relabel.bam
}

##############################################################################
# First step: samples
#
# samples has one line per BAM, as the LSF version wrote it: the BAM base
# name with any Proj_..._s_ prefix removed and __ appended. The trailing
# __ matches the __relabel suffix of the relabeled BAM, so s_X__ does not
# also select s_X1. Copy the lines into tumors and normals.
#
# With no tumors or normals file, stop here after writing samples. No run
# record is written, so the run can then be submitted in this directory.
##############################################################################

findBams

if [ -n "$BAMS" ] && [ ! -e samples ]; then
    echo "$BAMS" \
        | xargs -n 1 basename \
        | sed 's/.*_recal_s_/s_/' \
        | sed 's/^Proj_.*s_/s_/' \
        | sed 's/.bam/__/' \
        >samples
    echo "Wrote samples ($(wc -l <samples) BAMs in $BAMDIR)"
fi

if [ ! -s tumors ] || [ ! -s normals ]; then
    if [ -z "$BAMS" ]; then
        cat <<-END_NOBAMS

		    FATAL ERROR: no BAMs found: need pipeline/<alignments|bam>/*.bam
		    Link the pipeline output directory as pipeline/

		END_NOBAMS
        exit 1
    fi
    cat <<-END_SAMPLES

	    Need to define the tumor and normal samples.

	    samples lists the BAMs in $BAMDIR, one per line. Put the
	    tumor lines in tumors and the normal lines in normals, then
	    submit the run:

	        sbatch seqCNA/PIPE.sh $1

	END_SAMPLES
    exit 0
fi

TAG=q$PIPENAME

echo SDIR=$SDIR
echo VERSION=$SCRIPT_VERSION
echo RELEASE=$SCRIPT_RELEASE
echo ASSAY_TYPE=$ASSAY_TYPE BINSIZE=$BINSIZE

#
# Every sacct query in bCheck is bounded by this timestamp. Job names
# embed $$ and pids recycle, so without it a name from an earlier run
# could be counted against this one.
#
export SEQCNA_START=$(date -d '-2 min' +%Y-%m-%dT%H:%M:%S)
RUNRE="^${TAG}_.*_$$\$"

#
# Run records, read by bin/checkRun.sh together with sacct.
#
#   00.RUNSTATUS.txt     KEY=VALUE. STATUS is RUNNING from here on. It is
#                        set to COMPLETED only as the last action of this
#                        script, after the deliverables check, and to
#                        FAILED (with the STAGE it failed in) by the EXIT
#                        trap. If the control job is killed outright
#                        (SIGKILL on its memory cap, node failure) the
#                        file stays at RUNNING and checkRun.sh finds the
#                        control job dead in sacct.
#
#   SLURM.CTRL/jobs.tsv  one line per stage job: job id, stage, sample
#                        and log path, written at submission.
#
#   out/pipeline_info/   pairs.tsv (PAIR, TUMOR_BAM, NORMAL_BAM, ODIR,
#                        TUMOR, NORMAL), read by the COUNT and SEG stages
#                        and by selectBestMatch; cmd.sh.log, the
#                        provenance record.
#
#   pairing.tsv          the best normal of each tumor (TUMOR, NORMAL,
#                        PAIR, rms.logr), written by selectBestMatch. It
#                        names the normals, so it stays out of $PROJNO/.
#
RUNSTATUS=00.RUNSTATUS.txt
JOBS=SLURM.CTRL/jobs.tsv
INFODIR=out/pipeline_info
PAIRS_TSV=$INFODIR/pairs.tsv

CTRL_JOBID=none
CTRL_LOG=none
if [ "$RUN_MODE" == "sbatch" ]; then
    CTRL_JOBID=$SLURM_JOB_ID
    CTRL_LOG=$(scontrol show job $CTRL_JOBID 2>/dev/null | sed -n 's/^ *StdOut=//p') || true
    CTRL_LOG=${CTRL_LOG:-SLURM.CTRL/$CTRL_JOBID.out}
fi
echo RUN_MODE=$RUN_MODE CTRL_JOBID=$CTRL_JOBID HOST=$(hostname) PID=$$

RUN_STATUS=RUNNING
RUN_STAGE=SETUP
RUN_RC=""
RUN_MESSAGE=""
RUN_STARTED=$(date '+%Y-%m-%d %H:%M:%S')
RUN_FINISHED=""
RUN_DONE=0
BAMDIR=""
GENOME_TAG=""
GENOME=""
ASSAY=""
PROJNO=""
SAMPLES=""
NPAIRS=""
BSYNC_PID=""

writeRunStatus() {
    cat >$RUNSTATUS.tmp <<EOS
STATUS=$RUN_STATUS
STAGE=$RUN_STAGE
RC=$RUN_RC
MESSAGE=$RUN_MESSAGE
STARTED=$RUN_STARTED
FINISHED=$RUN_FINISHED
VERSION=$SCRIPT_VERSION
RELEASE=$SCRIPT_RELEASE
R_VERSION=$R_VERSION
RLIB_LOCK_MD5=$RLIB_LOCK_MD5
SAMTOOLS_VERSION=$SAMTOOLS_VERSION
MAGICK_VERSION=$MAGICK_VERSION
RUN_MODE=$RUN_MODE
CTRL_JOBID=$CTRL_JOBID
CTRL_LOG=$CTRL_LOG
HOST=$(hostname)
PID=$$
SEQCNA_START=$SEQCNA_START
ASSAY_TYPE=$ASSAY_TYPE
BINSIZE=$BINSIZE
BAMDIR=$BAMDIR
GENOME_TAG=$GENOME_TAG
GENOME=$GENOME
ASSAY=$ASSAY
PROJNO=$PROJNO
NPAIRS=$NPAIRS
SAMPLES=$SAMPLES
JOBS=$JOBS
EOS
    mv $RUNSTATUS.tmp $RUNSTATUS
}

setStage() {
    RUN_STAGE=$1
    writeRunStatus
    echo
    echo "==== STAGE $RUN_STAGE $(date '+%Y-%m-%d %H:%M:%S')"
}

#
# fatal MESSAGE
#
# Stop the run with MESSAGE in the status file (through the EXIT trap).
#
fatal() {
    RUN_MESSAGE="$*"
    echo
    echo "    FATAL ERROR: $RUN_MESSAGE"
    echo
    exit 1
}

#
# One run per project directory; PIPE.sh is not built to be rerun. A
# second run would overwrite the first run's status file and job
# manifest, and the stages write into out/, outAll/ and bamRelabel/, so
# the earlier run's output would be mixed into this one. Stop before
# anything is written if the directory holds any of these. LSF and
# LSF.FIX are the log directories of LSF-era runs. SLURM.CTRL itself is
# not a sign: Slurm creates it for this job's log. A control job that
# Slurm requeues after a node failure also stops here, on the records of
# its own first attempt. The check runs before the EXIT trap is
# installed, so the earlier run's records are left as they are; the
# reason is in this job's log.
#
# The exception is a run that failed in CHECK (a typo in tumors, a
# missing index): it submitted nothing and wrote only $RUNSTATUS and the
# header of $JOBS, so a new run replaces those records.
#
failedInCheck() {
    local f

    for f in out outAll bamRelabel LSF LSF.FIX; do
        if [ -e "$f" ]; then
            return 1
        fi
    done

    if [ ! -f "$RUNSTATUS" ] \
        || ! grep -qx "STATUS=FAILED" $RUNSTATUS \
        || ! grep -qx "STAGE=CHECK" $RUNSTATUS; then
        return 1
    fi

    if [ -e "$JOBS" ] && grep -qv '^#' $JOBS; then
        return 1
    fi
}

checkPreviousRun() {
    local found=""
    local f

    for f in $RUNSTATUS $JOBS out outAll bamRelabel LSF LSF.FIX; do
        if [ -e "$f" ]; then
            found="$found $f"
        fi
    done

    if [ -n "$found" ] && failedInCheck; then
        echo
        echo "    The previous run here failed in CHECK and submitted nothing:"
        sed -n 's/^MESSAGE=/        /p' $RUNSTATUS
        echo "    Replacing its records."
        echo
        return
    fi

    if [ -n "$found" ]; then
        echo
        echo "    FATAL: this directory already holds a run:$found"
        echo "    PIPE.sh cannot be rerun in the same directory; use a new one."
        echo
        exit 1
    fi
}

checkPreviousRun

mkdir -p SLURM.CTRL
printf '#JOBID\tSTAGE\tSAMPLE\tLOG\n' >$JOBS
writeRunStatus

#
# Any exit before RUN_DONE=1 is a failure: an error under set -e, a failed
# bCheck, a failed check (fatal), or scancel / time limit (SIGTERM).
# Record it in the status file and cancel the rest of this run's stage
# jobs so they do not keep running after the control job is gone. SIGTERM
# reaches the trap promptly only while the control job is in bSync (see
# slurmTools.sh), which is where it spends nearly all of its time. SIGHUP
# is a direct run losing its terminal; under nohup it is ignored and the
# trap is not set.
#
onSignal() {
    RUN_MESSAGE="control job got SIG$1 (scancel, time limit, node shutdown or closed terminal)"
    exit $2
}
trap 'onSignal TERM 143' TERM
trap 'onSignal INT 130' INT
trap 'onSignal HUP 129' HUP

onExit() {
    local rc=$?
    if [ -n "$BSYNC_PID" ]; then
        kill $BSYNC_PID 2>/dev/null || true
    fi
    if [ "$RUN_DONE" != "1" ]; then
        #
        # A failed bCheck leaves the reason in the stage job logs; copy it
        # into MESSAGE.
        #
        local reason=""
        reason=$(failedJobsMessage "$RUN_STAGE" "$JOBS" 2>/dev/null) || true
        if [ -n "$reason" ]; then
            RUN_MESSAGE="$reason${RUN_MESSAGE:+; $RUN_MESSAGE}"
        fi
        RUN_STATUS=FAILED
        RUN_RC=$rc
        RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
        writeRunStatus
        echo
        echo "    PIPE.sh FAILED in stage $RUN_STAGE (rc=$rc); cancelling queued jobs of this run"
        echo
        bKill "$RUNRE"
    fi
}
trap onExit EXIT

##############################################################################
# Up-front checks. Nothing is submitted until all of them pass.
##############################################################################

setStage CHECK

#
# R and the packages in Rlib.lock.tsv must be at the locked versions.
# Every R script checks again when it starts (bin/rscript.sh).
#
$SDIR/bin/rscript.sh --check \
    || fatal "R or a package does not match Rlib.lock.tsv (details above)"

#
# BAMs: pipeline/<alignments|bam>/*.bam, each with an index.
#
findBams
if [ -z "$BAMDIR" ] || [ ! -d "$BAMDIR" ]; then
    fatal "no pipeline/alignments or pipeline/bam directory; link the pipeline output as pipeline/"
fi

if [ -z "$BAMS" ]; then
    fatal "no BAM files in $BAMDIR"
fi

for BAM in $BAMS; do
    if [ ! -r "$BAM" ]; then
        fatal "cannot read BAM [$BAM]"
    fi
    if [ ! -e "${BAM%.bam}.bai" ] && [ ! -e "$BAM.bai" ]; then
        fatal "no index (.bai) for [$BAM]"
    fi
done

#
# tumors and normals (present, see First step): one pattern per line,
# matched against the relabeled BAM paths, as tumorBams and normalBams
# are chosen after FIX. Lines copied from samples end in __ and match
# only those paths, not the original BAM paths.
#
RELABELED=$(for BAM in $BAMS; do relabeledBam $BAM; done)
for F in tumors normals; do
    if [ -z "$(echo "$RELABELED" | fgrep -f $F)" ]; then
        fatal "no BAM in $BAMDIR matches any pattern in $F"
    fi
done

#
# Genome. Every BAM is identified by GenomeData/getGenomeBuildBAM.sh and
# all must agree. The tag of the original BAM picks the annotation
# assay; the GENOME= argument of getPairedCounts is taken after the FIX
# stage from the relabeled BAM (see below), as seqCNA.sh did.
#
FIRST_TAG=""
for BAM in $BAMS; do
    BUILD=$($SDIR/GenomeData/getGenomeBuildBAM.sh $BAM) \
        || fatal "cannot identify the genome of [$BAM]"
    echo "    genome [$BUILD] $BAM"
    if [ -z "$FIRST_TAG" ]; then
        FIRST_TAG=$BUILD
    elif [ "$BUILD" != "$FIRST_TAG" ]; then
        fatal "BAMs on different genomes: [$FIRST_TAG] and [$BUILD] ($BAM)"
    fi
done
GENOME_TAG=$FIRST_TAG

if [[ $GENOME_TAG =~ b37|hg19 ]]; then
    ASSAY=Exome
elif [[ $GENOME_TAG =~ ^mm10|GRC_m38 ]]; then
    ASSAY=Mouse-Twist
else
    fatal "unknown genome [$GENOME_TAG]; supported: b37, b37_dmp, b37+mm10, mm10, GRC_m38"
fi
echo GENOME_TAG=$GENOME_TAG ASSAY=$ASSAY

#
# Project id: from pipeline/*request.txt, else from the path.
#
if ls pipeline/*request.txt >/dev/null 2>&1; then
    PROJNO=$(echo $(ls pipeline/*request.txt) | perl -ne 'm|/(Proj_.*)_request|; print $1')
else
    PROJNO=$($SDIR/extractProjectIDFromPath.py $(pwd -P))
fi
if [ -z "$PROJNO" ]; then
    fatal "cannot determine the project id"
fi
echo PROJNO=$PROJNO

writeRunStatus

#
# Provenance, written once the run has all its inputs settled.
#
writeCmdLog() {
    local CMD_LOG=$INFODIR/cmd.sh.log
    local GURL
    GURL=$(git --git-dir=$SDIR/.git --work-tree=$SDIR config --get remote.origin.url) || GURL=""

    cat <<-END_VERSION >$CMD_LOG
	DATE: $(date)
	SDIR: $SDIR
	GURL: $GURL
	GTAG: $SCRIPT_VERSION
	RELEASE: $SCRIPT_RELEASE
	PWD: $PWD
	RUN_MODE: $RUN_MODE
	CTRL_JOBID: $CTRL_JOBID

	ASSAY_TYPE: $ASSAY_TYPE
	BINSIZE: $BINSIZE
	GENOME_TAG: $GENOME_TAG
	GENOME: $GENOME
	ASSAY: $ASSAY
	PROJNO: $PROJNO
	BAMDIR: $BAMDIR
	SAMPLES: $SAMPLES
	NPAIRS: $NPAIRS

	TOOLS:
	  $(samtools --version | head -1 || true)
	  $(R --version | head -1 || true)
	  $(magick -version | head -1 || true)
	  $(sed -n 's/^seqDNAcopy /seqDNAcopy /p' $SDIR/Rlib/VERSIONS.txt | head -1 || true)
	  Rlib.lock.tsv md5 $RLIB_LOCK_MD5

	Script: $ORIG_CMD
	END_VERSION
}

##############################################################################
# Stages
##############################################################################

#
# Scrub the Slurm variables of the control job itself. bsub submits with
# the default --export=ALL, so they would otherwise ride along into every
# stage job. SLURM_CONF stays; the client tools need it. The SBATCH_*
# input variables go too (SDIR has already been taken from
# SBATCH_SCRIPT_DIR): one left in the submitting shell would apply to
# every stage job, and SBATCH_QOS=priority there gets every SHORT job
# rejected. Stage partition, walltime and qos come only from seqSub.
#
for V in ${!SLURM_@}; do
    if [ "$V" != "SLURM_CONF" ]; then
        unset $V
    fi
done
for V in ${!SBATCH_@}; do
    unset $V
done

#
# Walltime classes. Every stage job is SHORT or LONG.
#
#   SHORT  -W 1:59:00, no qos. ~/bin/bsub sends anything under two hours
#          to cmobic_short,cpushort. cpushort allows only qos=normal and
#          EnforcePartLimits=ALL is set, so a SHORT job must not carry
#          a qos or it is rejected ("Invalid qos specification").
#   LONG   -W 12:00:00 and qos priority. bsub sends it to cmobic_cpu;
#          SBATCH_QOS=priority is set on that one bsub call only.
#
# SHORT whenever possible: jobs on cmobic_cpu can wait hours to start.
# COUNT reads two full BAMs and its run time grows with their size; its
# class is chosen per job by runClass. The other stages are always
# SHORT. -M is the total memory for the job and a hard cgroup limit.
#
SHORT_W=1:59:00
LONG_W=12:00:00

#
# Minutes per GB of input (tumor + normal BAM). Measured 2026-10-20 in
# the foreground on an IRIS login node: 176 s for 9.0 GB of exome BAMs
# (0.33 min/GB, peak 2.0 GB RSS); 1 min/GB is a 3x margin. See
# docs/SLURM_PORT.md.
#
RATE_COUNT=1

#
# A job is SHORT if its estimated run time is at most SHORT_MAX_MIN, half
# the SHORT walltime, to leave room for slower nodes and a loaded shared
# filesystem. SEQCNA_SHORT_MAX_MIN overrides it (0 makes every estimated
# stage LONG).
#
SHORT_MAX_MIN=${SEQCNA_SHORT_MAX_MIN:-60}

#
# runClass MIN_PER_GB FILE [FILE ...]
#
# Print SHORT or LONG for a job whose input is FILE..., and log the
# estimate on stderr.
#
runClass() {
    local rate=$1
    local bytes=0
    local f size
    shift

    for f in "$@"; do
        size=$(stat -L -c %s "$f") || return 1
        bytes=$((bytes + size))
    done

    awk -v b=$bytes -v r=$rate -v m=$SHORT_MAX_MIN 'BEGIN {
        est = b / 1e9 * r
        class = (est <= m) ? "SHORT" : "LONG"
        printf "    runClass %.1f GB x %s min/GB = %.0f min -> %s\n", \
            b / 1e9, r, est, class >"/dev/stderr"
        print class
    }'
}

#
# seqSub STAGE SAMPLE LOGDIR CLASS "BSUB_OPTS" CMD [ARG ...]
#
# Submit one stage job as ${TAG}_STAGE_$$ in walltime class CLASS (SHORT
# or LONG) through bin/runStage.sh, which appends #SEQCNA_EXIT=<rc> to
# the job log, and record it in $JOBS.
#
seqSub() {
    local stage=$1
    local sample=$2
    local logdir=$3
    local class=$4
    local opts=$5
    local out jobid runtime
    local qos=()
    shift 5

    case $class in
        SHORT)
            runtime=$SHORT_W
            ;;
        LONG)
            runtime=$LONG_W
            qos=(SBATCH_QOS=priority)
            ;;
        *)
            fatal "bad class [$class] for stage $stage sample $sample"
            ;;
    esac

    if ! out=$(env "${qos[@]}" bsub -W $runtime -o $logdir/ \
                   -J ${TAG}_${stage}_$$ $opts \
                   $SDIR/bin/runStage.sh "$@"); then
        echo "$out"
        fatal "bsub failed for stage $stage sample $sample"
    fi
    echo "$out"

    jobid=$(echo "$out" | sed -n 's/^Job <\([0-9][0-9]*\)> is submitted.*/\1/p')
    if [ -z "$jobid" ]; then
        fatal "no job id in bsub output for stage $stage sample $sample"
    fi

    printf '%s\t%s\t%s\t%s\n' "$jobid" "$stage" "$sample" "$logdir/$jobid.out" >>$JOBS
}

#
# waitStage STAGE
#
# Block until every ${TAG}_STAGE_$$ job is done, then abort the run if any
# of them did not reach COMPLETED. The number of jobs submitted, from
# $JOBS, lets bCheck wait for every accounting record.
#
waitStage() {
    local njobs
    setStage $1
    njobs=$(awk -F'\t' -v s="$1" '$2 == s' $JOBS | wc -l)
    bSync ${TAG}_$1_$$
    bCheck ${TAG}_$1_$$ $njobs
}

#
# Deliverables. Everything has to be non-empty.
#
MISSING=""

checkFile() {
    if [ ! -s "$1" ]; then
        MISSING="$MISSING $1"
    fi
}

failIfMissing() {
    if [ -n "$MISSING" ]; then
        echo
        echo "    ERROR missing deliverables:"
        for F in $MISSING; do
            echo "        $F"
        done
        echo
        fatal "missing deliverables:$MISSING"
    fi
}

#
# finishRun MESSAGE
#
# The sweep of every stage job of the run, then COMPLETED.
#
finishRun() {
    local failed
    setStage STAGING
    failed=$(bCheckAll "$RUNRE")
    if [ "$failed" != "" ]; then
        echo -e "\n\n\nFailed Slurm jobs\n\n"
        echo "$failed"
        echo -e "\n\n"
        fatal "stage jobs not COMPLETED:$(echo "$failed" | awk '{printf " %s %s", $1, $2}')"
    fi

    RUN_STATUS=COMPLETED
    RUN_STAGE=DONE
    RUN_RC=0
    RUN_MESSAGE="$1"
    RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
    writeRunStatus
    RUN_DONE=1

    echo
    echo "==== seqCNA run COMPLETED $RUN_FINISHED"
    echo
}

setStage SUBMIT
mkdir -p out $INFODIR

#
# Stage 1: FIX, one job per BAM. Writes bamRelabel/<base>/<base>__relabel.bam.
#
for BAM in $BAMS; do
    seqSub FIX $(basename $BAM .bam) SLURM.01.FIX SHORT "-n 1 -M 4G" \
        $SDIR/fixChromosomeNames.sh $BAM
done

waitStage FIX

#
# Tumor and normal BAMs are the relabeled ones. The GENOME= argument of
# getPairedCounts comes from the relabeled tumor BAM, whose tag can
# differ from the original's (UCSC mm10 becomes mm10-relabel). This is
# the case table of the retired seqCNA.sh; keep it consistent with
# GenomeData/getGenomeBuildBAM.sh and the ASSAY choice above.
#
ls bamRelabel/*/*.bam | fgrep -f tumors >tumorBams
ls bamRelabel/*/*.bam | fgrep -f normals >normalBams

nTumBams=$(wc -l <tumorBams)
nNormBams=$(wc -l <normalBams)
if [ "$nTumBams" == "0" ] || [ "$nNormBams" == "0" ]; then
    fatal "no tumor ($nTumBams) or normal ($nNormBams) BAMs after FIX"
fi
echo "tumors=$nTumBams normals=$nNormBams"

RELABEL_TAG=$($SDIR/GenomeData/getGenomeBuildBAM.sh $(head -1 tumorBams)) \
    || fatal "cannot identify the genome of the relabeled BAM [$(head -1 tumorBams)]"

case $RELABEL_TAG in
    mm10-relabel|GRC_m38)
        GENOME="mm10"
        ;;
    b37+mm10|b37|b37_dmp)
        GENOME="hg19"
        ;;
    *)
        fatal "unknown or unusable genome after relabel [$RELABEL_TAG]"
        ;;
esac
echo GENOME=$GENOME
writeRunStatus

#
# Pairs: every tumor x normal. The pair id is <tumorId>__<normalId> from
# the @RG SM tags (getSampleTag.sh); its directory is
# out/<first 2 chars of tumorId>/<tumorId>/<pairId>/.
#
# A sample may be listed in both tumors and normals (e.g. every sample
# as a candidate normal for the others). Its self-pair is skipped; the
# LSF version ran it and selectBestMatch dropped it. Every tumor needs at
# least one normal other than itself.
#
printf 'PAIR\tTUMOR_BAM\tNORMAL_BAM\tODIR\tTUMOR\tNORMAL\n' >$PAIRS_TSV
for tumor in $(cat tumorBams); do
    tumorId=$($SDIR/getSampleTag.sh $tumor)
    SAMPLES="$SAMPLES $tumorId"
    scatter=${tumorId:0:2}
    nTumorPairs=0
    for normal in $(cat normalBams); do
        normalId=$($SDIR/getSampleTag.sh $normal)
        if [ "$tumorId" == "$normalId" ]; then
            echo "    skipping self-pair [$tumorId]: $tumor $normal"
            continue
        fi
        sID=${tumorId}__${normalId}
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sID" "$tumor" "$normal" \
            "out/$scatter/$tumorId/$sID" "$tumorId" "$normalId" >>$PAIRS_TSV
        nTumorPairs=$((nTumorPairs + 1))
    done
    if [ "$nTumorPairs" == "0" ]; then
        fatal "tumor [$tumorId] has no normal other than itself: $tumor"
    fi
done
SAMPLES=${SAMPLES# }
NPAIRS=$(awk 'NR > 1' $PAIRS_TSV | wc -l)
echo SAMPLES=$SAMPLES
echo NPAIRS=$NPAIRS
writeRunStatus
writeCmdLog

#
# Stage 2: COUNT, one job per pair: binned, GC-corrected counts.
#
while IFS=$'\t' read -r sID tumor normal oDir _; do
    CLASS=$(runClass $RATE_COUNT $tumor $normal) || fatal "cannot size [$tumor $normal]"
    seqSub COUNT $sID SLURM.02.COUNT $CLASS "-n 2 -M 8G" \
        $SDIR/getPairedCounts GENOME=$GENOME NORMAL=$normal TUMOR=$tumor \
            ODIR=$oDir SAMPLEID=$sID
done < <(awk 'NR > 1' $PAIRS_TSV)

waitStage COUNT

#
# Stage 3: SEG, one job per pair: segmentation, clustering, plot. With
# TUMORID, seqSegment also writes the plot titled with the tumor id only
# (<pair>_seqSeg.redacted.png), from which POST builds the PDF.
#
while IFS=$'\t' read -r sID tumor normal oDir tumorId _; do
    COUNTS=$oDir/${sID}_Counts.rda
    if [ ! -s "$COUNTS" ]; then
        fatal "no COUNT output [$COUNTS]"
    fi
    seqSub SEG $sID SLURM.03.SEG SHORT "-n 2 -M 8G" \
        $SDIR/seqSegment BINSIZE=$BINSIZE ODIR=$oDir COUNTS=$COUNTS \
            TUMORID=$tumorId
done < <(awk 'NR > 1' $PAIRS_TSV)

waitStage SEG

#
# Stage 4: BEST, in the control job. For each tumor keep the normal with
# the lowest rms.logr: selectBestMatch writes bestMatches____out and
# pairing.tsv, one best pair per tumor, with the tumor and normal ids
# from pairs.tsv. POST uses pairing.tsv to replace the pair ids in the
# deliverables by tumor ids. The full set of pairs is hard-linked into
# outAll/, then the non-best pair directories are removed from out/. The match is on whole lines: as a
# substring match, the best pair ..._s_N1 would also keep ..._s_N10.
#
setStage BEST

$SDIR/selectBestMatch out
if [ ! -s bestMatches____out ]; then
    fatal "selectBestMatch wrote no best matches (bestMatches____out)"
fi
nBest=$(wc -l <bestMatches____out)
if [ "$nBest" != "$nTumBams" ]; then
    fatal "selectBestMatch found $nBest best matches for $nTumBams tumors"
fi
if [ ! -s pairing.tsv ]; then
    fatal "selectBestMatch wrote no pairing.tsv"
fi

mkdir outAll
rsync -a --link-dest=../out out/ outAll

PRUNE=$(ls -d out/*/*/* | grep -vxFf bestMatches____out || true)
if [ -n "$PRUNE" ]; then
    echo "$PRUNE" | xargs -t rm -rf
fi

#
# Stage 5: POST, one job: deliverables directory $PROJNO/ with the IGV
# seg, gene calls and the PDF of the plots.
#
seqSub POST all SLURM.04.POST SHORT "-n 1 -M 8G" \
    $SDIR/postProcess.sh $ASSAY $PROJNO

waitStage POST

setStage DELIVERABLES

for oDir in $(cat bestMatches____out); do
    sID=$(basename $oDir)
    checkFile $oDir/${sID}_seqSeg.seg
    checkFile $oDir/${sID}_seqSeg.out
    checkFile $oDir/${sID}_seqSeg.png
    checkFile $oDir/${sID}_seqSeg.redacted.png
done
checkFile pairing.tsv
checkFile $PROJNO/${PROJNO}___IGV.seg
checkFile $PROJNO/${PROJNO}___seqSeg.pdf
#
# getGeneCalls writes README.txt instead of the matrices when no gene
# passes the significance cuts.
#
if [ ! -s $PROJNO/README.txt ]; then
    checkFile $PROJNO/SegmentMatrix.csv
    checkFile $PROJNO/GeneMatrix.csv
    checkFile $PROJNO/${PROJNO}___GeneTable.xlsx
fi
checkFile $INFODIR/cmd.sh.log

failIfMissing

finishRun ""

cat <<EOS

PIPE.sh done

  Deliverables in $PROJNO/
  Check the run with

      $SDIR/bin/checkRun.sh

EOS
