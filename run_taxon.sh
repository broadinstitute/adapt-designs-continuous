#!/bin/bash

# Design for a single taxon.

# Read the input
taxid=$1
segment=$2
refaccs=$3

# Set variables for design
PREP_MEMOIZE_DIR="/tmp/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=0.15
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="3"
ARG_PP="0.95"
ARG_MAXPRIMERSATSITE="5"
ARG_MAXTARGETLENGTH="500"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="15"

# Make the memoize directory
mkdir -p /tmp/prep-memoize-dir

# Set --prep-influenza if needed
if [ "$taxid" == "11320" ] || [ "$taxid" == "11520" ]; then
    ARG_COVERBYYEARSTART=2015
    ARG_COVERBYYEARDECAY=0.95
    ARGS_INFLUENZA="--prep-influenza --cover-by-year-decay $ARG_COVERBYYEARSTART $ARG_COVERBYYEARDECAY"
else
    ARGS_INFLUENZA=""
fi

# Activate dgd conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate dgd

# Determine an output directory and create it
timestamp=$(date +"%s")
outdir="out/designs/${taxid}_${segment}/${timestamp}"
mkdir -p $outdir

# Run the design
design.py complete-targets auto-from-args $taxid $segment $refaccs $outdir/design.tsv -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD $ARGS_INFLUENZA --write-input-seqs $outdir/input-sequences.txt --verbose &> $outdir/design.out

# gzip the stdout/stderr
gzip $outdir/design.out
