#!/bin/bash

# Design for a single taxon.

# Load run script from adapt-designs repo
source adapt-designs/scripts/run-adapt/run_common.sh

# Read the input
taxalist=$1
taxid=$2
segment=$3
specific=$4
obj=$5

# Create a temporary file containing taxonomies to use for specificity
# First get the family of this taxid from the taxalist, and then
# pull out all taxids (except this one) with that family
family=$(cat $taxalist | awk -F'\t' -v taxid="$taxid" -v segment="$segment" '$4==taxid && $5==segment {print $1}')
specificity_taxa=$(mktemp)
cat $taxalist | awk -F'\t' -v taxid="$taxid" -v family="$family" '$1==family && $4!=taxid {print $4"\t"$5}' > $specificity_taxa

# Find the RefSeq accessions
refseqs=$(cat $taxalist | awk -F'\t' -v taxid="$taxid" -v segment="$segment" '$4==taxid && $5==segment {print $6}')

# Change some parameters from their default values
CLUSTER_THRESHOLD="0.30"
ARG_PM="3"
ARG_MAXPRIMERSATSITE="10"

# Set --prep-influenza if needed
if [ "$taxid" == "11320" ] || [ "$taxid" == "11520" ] || [ "$taxid" == "11552" ]; then
    ARGS_INFLUENZA="--prep-influenza"
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
if [ -d "/ebs/tmpfs/tmp" ]; then
    export TMPDIR="/ebs/tmpfs/tmp"
fi

# Determine an output directory and create it
# Note this overrides the OUT_DIR set when sourcing run_common.sh
timestamp=$(date +"%s")
segmentnospace=${segment// /-}
experiment="${specific}_${obj}"
export OUT_DIR="out/designs/$experiment/${taxid}_${segmentnospace}/${timestamp}"
mkdir -p $OUT_DIR

# Make arguments depending on what is set for 'specific' and 'obj'
if [[ "$specific" == "specific" ]]; then
    ARG_SPECIFICITY="--id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --specific-against-taxa $specificity_taxa"
elif [[ "$specific" == "nonspecific" ]]; then
    ARG_SPECIFICITY=""
else
    echo "FATAL: Unknown value for 'specific': $specific"
    exit 1
fi
if [[ "$obj" == "max-activity" ]]; then
    ARG_OBJ="--obj maximize-activity --soft-guide-constraint $ARG_SOFTGUIDECONSTRAINT --hard-guide-constraint $ARG_HARDGUIDECONSTRAINT --penalty-strength $ARG_PENALTYSTRENGTH --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM"
elif [[ "$obj" == "min-guides" ]]; then
    ARG_OBJ="--obj minimize-guides -gm $ARG_GM -gp $ARG_GP --require-flanking3 H"
else
    echo "FATAL: Unknown value for 'obj': $obj"
    exit 1
fi

# To make sure runs for different "experiments" use the same input data,
# specify --use-accessions
ARG_USE_ACCESSIONS="--use-accessions taxonomies/accession-lists/all-vertebrate.20200501.tsv"

# Run the design
run-adapt design complete-targets auto-from-args $taxid "$segment" $refseqs $OUT_DIR/design.tsv -gl $ARG_GL -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --primer-gc-content-bounds $ARG_PRIMER_GC_LO $ARG_PRIMER_GC_HI --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH $ARG_OBJ --obj-fn-weights $ARG_OBJFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS $ARG_SPECIFICITY --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --cluster-threshold $CLUSTER_THRESHOLD $ARGS_INFLUENZA $ARG_USE_ACCESSIONS --verbose

# Remove tmp files
rm $specificity_taxa
