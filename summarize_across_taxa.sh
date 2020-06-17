#!/bin/bash

# Summarize designs for all taxonomies.


JACCARD_THRES=0.5
TAXONOMIES_FILE="taxonomies/all-vertebrate.tsv"
IGNORE_BEFORE=1573900000
EXPERIMENTS=("nonspecific_max-activity" "nonspecific_min-guides" "specific_max-activity" "specific_min-guides")

# Determine an output directory and create it
timestamp=$(date +"%s")
summary_outdir="out/summaries/${timestamp}"
mkdir -p $summary_outdir

# Create blank files for summary
echo -n "" > $summary_outdir/summary.tsv
echo -n "" > $summary_outdir/summary.one-segment-per-taxid.tsv
echo -n "" > $summary_outdir/summary.taxa-not-designed.tsv

function summarize {
    experiment=$1
    if [ ! -d "out/designs/$experiment" ]; then
        echo "FATAL: unknown experiment '$experiment'"
        exit 1
    fi

    while read -r taxonomy; do
        family=$(echo "$taxonomy" | awk -F'\t' '{print $1}')
        genus=$(echo "$taxonomy" | awk -F'\t' '{print $2}')
        species=$(echo "$taxonomy" | awk -F'\t' '{print $3}')
        taxid=$(echo "$taxonomy" | awk -F'\t' '{print $4}')
        segment=$(echo "$taxonomy" | awk -F'\t' '{print $5}')
        refaccs=$(echo "$taxonomy" | awk -F'\t' '{print $6}')

        segmentnospace=${segment// /-}
        taxonomy_outdir="out/designs/$experiment/${taxid}_${segmentnospace}"
        if [ ! -d $taxonomy_outdir ]; then
            # Designs have not been started for this taxonomy
            continue
        fi

        # Summarize output, but redirect stderr to /dev/null
        summary=$(python summarize_design_on_taxon.py $taxonomy_outdir $JACCARD_THRES --ignore-before $IGNORE_BEFORE 2> /dev/null)
        if [ -z "$summary" ]; then
            # An error occurred; likely no designs could be found for this
            # taxonomy (e.g., none completed)
            continue
        fi
        while read -r summary_line; do
            timestamp=$(echo "$summary_line" | awk -F'\t' '{print $3}')

            echo -e "$family\t$genus\t$species\t$taxid\t$segment\t$experiment\t$summary_line" >> $summary_outdir/summary.tsv
        done <<< "$summary"
    done < <(tail -n +2 $TAXONOMIES_FILE)

    # Now rewrite summary.tsv, choosing just one segment per taxid for this experiment
    while read -r taxid; do
        # Previously this picked the segment with the best (smallest) total cluster cost; however,
        # this does not work as easily with the max-activity objective
        # Now, pick the segment with the best objective value for its largest cluster, for this experiment
        if [[ $experiment == *"min"* ]]; then
            segment_with_best_obj=$(cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" -v experiment="$experiment" '$4==taxid && $6==experiment && $7=="cluster" {print $5"\t"$11}' | sort -k2g | head -n 1 | awk -F'\t' '{print $1}')
        elif [[ $experiment == *"max"* ]]; then
            # Same as above, but pick highest objective value: `sort -k2gr`
            segment_with_best_obj=$(cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" -v experiment="$experiment" '$4==taxid && $6==experiment && $7=="cluster" {print $5"\t"$11}' | sort -k2gr | head -n 1 | awk -F'\t' '{print $1}')
        else
            echo "FATAL: unknown if objective is min or max"
            exit 1
        fi

        # Copy summary.tsv for this taxid and segment and experiment
        cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" -v segment="$segment_with_best_obj" -v experiment="$experiment" '$4==taxid && $5==segment && $6==experiment {print $0}' >> $summary_outdir/summary.one-segment-per-taxid.tsv
    done < <(cat $summary_outdir/summary.tsv | awk -F'\t' -v experiment="$experiment" '$6==experiment {print $4}' | sort | uniq)

    # Write a list of the taxa that did not yield designs for this experiment
    while read -r taxid; do
        cat $TAXONOMIES_FILE | awk -F'\t' -v taxid="$taxid" -v experiment="$experiment" '$4==taxid {print $1"\t"$2"\t"$3"\t"$4"\t"experiment}' | sort | uniq >> $summary_outdir/summary.taxa-not-designed.tsv
    done < <(comm -23 <(cat $TAXONOMIES_FILE | tail -n +2 | awk -F'\t' '{print $4}' | sort | uniq) <(cat $summary_outdir/summary.one-segment-per-taxid.tsv | awk -F'\t' -v experiment="$experiment" '$6==experiment {print $4}' | sort | uniq))
}

for experiment in "${EXPERIMENTS[@]}"; do
    summarize $experiment
done
