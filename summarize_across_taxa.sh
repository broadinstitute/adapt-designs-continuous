#!/bin/bash

# Summarize designs for all taxonomies.

JACCARD_THRES=0.5
TAXONOMIES_FILE="taxonomies/taxonomies.tsv"

# Determine an output directory and create it
timestamp=$(date +"%s")
summary_outdir="out/summaries/${timestamp}"
mkdir -p $summary_outdir

echo -n "" > $summary_outdir/summary.tsv
while read -r taxonomy; do
    family=$(echo "$taxonomy" | awk -F'\t' '{print $1}')
    genus=$(echo "$taxonomy" | awk -F'\t' '{print $2}')
    species=$(echo "$taxonomy" | awk -F'\t' '{print $3}')
    taxid=$(echo "$taxonomy" | awk -F'\t' '{print $4}')
    segment=$(echo "$taxonomy" | awk -F'\t' '{print $5}')
    refaccs=$(echo "$taxonomy" | awk -F'\t' '{print $6}')

    taxonomy_outdir="out/designs/${taxid}_${segment}"
    if [ ! -d $taxonomy_outdir ]; then
        # Designs have not been started for this taxonomy
        continue
    fi

    num_seqs=$(cat $taxonomy_outdir/input-sequences.txt | wc -l)

    summary=$(python summarize_design_on_taxon.py $taxonomy_outdir $JACCARD_THRES)
    if [ -z "$summary" ]; then
        # An error occurred; likely no designs could be found for this
        # taxonomy (e.g., none completed)
        continue
    fi
    while read -r summary_line; do
        echo -e "$family\t$genus\t$species\t$taxid\t$segment\t$num_seqs\t$summary_line" >> $summary_outdir/summary.tsv
    done <<< "$summary"
done < <(tail -n +2 $TAXONOMIES_FILE)
