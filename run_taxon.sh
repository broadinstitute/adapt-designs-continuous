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

# Adjust some amplicon criteria for a short virus that
# cannot satisfy the default criteria
if [ "$taxid" == "365327" ]; then
    ARG_PRIMER_GC_LO="0.2"
    ARG_PRIMER_GC_HI="0.8"
fi

# Adjust specificity arguments for some viruses; these are ones that,
# with the default specificity arguments (--id-m and --id-frac) produced
# no designs when $specific is "specific" but *did* produce designs
# with $specific is "nonspecific"; likely the specificity arguments are
# too strict, so loosen them
# A few of them don't work with the adjusted (relaxed) arguments, so
# for these relaxe them even further
# These are: Porcine parvovivrus 2; Mamastrovirus 5; Respiratory syncytial virus; Ungulate bocaparvovirus 2; Ungulate bocaparvovirus 3; Panine gammaherpesvirus 1; Pongine gammaherpesvirus 2; Bovine rhinovirus 1; Rattus norvegicus polymoavirus 1; Sewage associated gemycircularvirus 3; Gemycircularvirus HV-GcV1; Bos taurus papillomavirus 13; Aichivirus F; Sewage derived gemykibivirus 1; Hepacivirus B; Goose paramyxivirus SF02; Bovine associated cyclovirus 1; Chiropteran bocaparvovirus 4; Puma feline foamy virus; Chicken associated huchismacovirus 1; Human associated huchismacovirus 1; Pestivirus K; Betaarterivirus suid 1; Fort Sherman orthobunyavirus; Equid gammaherpesvirus 7; Mopeia Lassa virus reassortant 29; Hippotragine gammaherpesvirus 1; Deltapapillomavirus 4; Finkel-Biskis-Jinkins murine sarcoma virus; Hepatitis GB virus B; Simian T-cell lymphotropic virus 6; Squirrel fibroma virus; Lyssavirus Ozernoe; Israel turkey meningoencephalomyelitis virus; Senegalese sole Iberrian betanodavirus; unidentified human coronavirus; Porcine circovirus type 1/2a; Human erythrovirus V9; Turkey parvovirus 1078; Bat mastadenovirus; Cervid alphaherpesvirus 1; Cyclovirus PKgoat21/PAK/2009
if [ "$taxid" == "1126383" ] || [ "$taxid" == "1239569" ] || [ "$taxid" == "12814" ] || [ "$taxid" == "1511874" ] || [ "$taxid" == "1511875" ] || [ "$taxid" == "159602" ] || [ "$taxid" == "159603" ] || [ "$taxid" == "1606765" ] || [ "$taxid" == "1679933" ] || [ "$taxid" == "1843761" ] || [ "$taxid" == "1862824" ] || [ "$taxid" == "1887213" ] || [ "$taxid" == "1986959" ] || [ "$taxid" == "2004967" ] || [ "$taxid" == "2008762" ] || [ "$taxid" == "204987" ] || [ "$taxid" == "2050927" ] || [ "$taxid" == "2169773" ] || [ "$taxid" == "2169894" ] || [ "$taxid" == "2169932" ] || [ "$taxid" == "2169934" ] || [ "$taxid" == "2170090" ] || [ "$taxid" == "2499685" ] || [ "$taxid" == "2560477" ] || [ "$taxid" == "291612" ] || [ "$taxid" == "300180" ] || [ "$taxid" == "333341" ] || [ "$taxid" == "337052" ] || [ "$taxid" == "353765" ] || [ "$taxid" == "39113" ] || [ "$taxid" == "481147" ] || [ "$taxid" == "538970" ] || [ "$taxid" == "642022" ] || [ "$taxid" == "64291" ] || [ "$taxid" == "683176" ] || [ "$taxid" == "694448" ] || [ "$taxid" == "720569" ] || [ "$taxid" == "72197" ] || [ "$taxid" == "740933" ] || [ "$taxid" == "740971" ] || [ "$taxid" == "79891" ] || [ "$taxid" == "942033" ]; then
    ARG_IDM="3"
    ARG_IDFRAC="0.2"
fi
if [ "$taxid" == "1679933" ] || [ "$taxid" == "1843761" ] || [ "$taxid" == "2004967" ] || [ "$taxid" == "1986959" ] || [ "$taxid" == "159602" ] || [ "$taxid" == "159603" ] || [ "$taxid" == "1606765" ] || [ "$taxid" == "1862824" ] || [ "$taxid" == "2008762" ] || [ "$taxid" == "2050927" ] || [ "$taxid" == "2169773" ] || [ "$taxid" == "353765" ] || [ "$taxid" == "538970" ] || [ "$taxid" == "642022" ] || [ "$taxid" == "740971" ] || [ "$taxid" == "942033" ]; then
    # Further relax arguments
    ARG_IDM="2"
    ARG_IDFRAC="0.4"
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
