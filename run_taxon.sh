#!/bin/bash

# Design for a single taxon.

# Read the input
taxid=$1
segment=$2
refaccs=$3

# Set variables for design
PREP_MEMOIZE_DIR="/ebs/dgd-analysis/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=0.20
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="3"
ARG_PP="0.99"
ARG_MAXPRIMERSATSITE="10"
ARG_MAXTARGETLENGTH="250"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="10"

# Set a predictive model
# This is from commit da10963 of my adapt-seq-design repo
PREDICTIVE_MODEL="/home/hayden/adapt-seq-design/models/predictor_exp-and-pos_regress-on-active/model-8f534a8c"

# Make the memoize directory
mkdir -p $PREP_MEMOIZE_DIR

# Set --prep-influenza if needed
if [ "$taxid" == "11320" ] || [ "$taxid" == "11520" ] || [ "$taxid" == "11552" ]; then
    ARG_COVERBYYEARSTART=2015
    ARG_COVERBYYEARDECAY=0.95
    # Use --gp-over-all-seqs to improve runtime on diverse datasets
    ARGS_INFLUENZA="--prep-influenza --cover-by-year-decay $ARG_COVERBYYEARSTART $ARG_COVERBYYEARDECAY --gp-over-all-seqs"
else
    ARGS_INFLUENZA=""
fi

# Adjust some arguments for influenza A virus, for a reasonable runtime
if [ "$taxid" == "11320" ]; then
    ARG_PM="2"
    ARG_MAXPRIMERSATSITE="5"
fi

# Adjust some arguments for dengue, for a reasonable runtime
if [ "$taxid" == "12637" ]; then
    ARG_PM="2"
    ARG_MAXPRIMERSATSITE="5"
fi

# Adjust some arguments for Norwalk virus and Rhinovirus C, which
# have no designs with the above parameters (no suitable targets)
if [ "$taxid" == "11983" ] || [ "$taxid" == "463676" ]; then
    ARG_MAXPRIMERSATSITE="20"
    ARG_MAXTARGETLENGTH="500"
fi

# Set tmp directory
export TMPDIR="/ebs/tmpfs/tmp"

# Activate adapt conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate adapt-with-tf

# Determine an output directory and create it
timestamp=$(date +"%s")
segmentnospace=${segment// /-}
outdir="out/designs/${taxid}_${segmentnospace}/${timestamp}"
mkdir -p $outdir

# Run the design
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-args $taxid "$segment" $refaccs $outdir/design.tsv -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD $ARGS_INFLUENZA --write-input-seqs $outdir/input-sequences.txt --verbose &> $outdir/design.out

# gzip the stdout/stderr
gzip $outdir/design.out
