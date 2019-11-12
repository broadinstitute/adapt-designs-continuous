#!/bin/bash

# Summarize designs for all taxonomies.

JACCARD_THRES=0.5
TAXONOMIES_FILE="taxonomies/all-human-associated-virus.tsv"

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

    segmentnospace=${segment// /-}
    taxonomy_outdir="out/designs/${taxid}_${segmentnospace}"
    if [ ! -d $taxonomy_outdir ]; then
        # Designs have not been started for this taxonomy
        continue
    fi

    # Summarize output, but redirect stderr to /dev/null
    summary=$(python summarize_design_on_taxon.py $taxonomy_outdir $JACCARD_THRES 2> /dev/null)
    if [ -z "$summary" ]; then
        # An error occurred; likely no designs could be found for this
        # taxonomy (e.g., none completed)
        continue
    fi
    while read -r summary_line; do
        timestamp=$(echo "$summary_line" | awk -F'\t' '{print $3}')
        num_seqs=$(cat ${taxonomy_outdir}/${timestamp}/input-sequences.txt | wc -l)

        echo -e "$family\t$genus\t$species\t$taxid\t$segment\t$num_seqs\t$summary_line" >> $summary_outdir/summary.tsv
    done <<< "$summary"
done < <(tail -n +2 $TAXONOMIES_FILE)

# Now rewirte summary.tsv, choosing just one segment per taxid
# Pick the one with the most number of sequences used for design, breaking ties
# arbitrarily
echo -n "" > $summary_outdir/summary.one-segment-per-taxid.tsv
while read -r taxid; do
    # Sort all segments by the number of sequences, and pick the one with the
    # most
    segment_with_most_num_seqs=$(cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" '$4==taxid {print $5"\t"$6}' | sort -k2nr | head -n 1 | awk '{print $1}')

    # Copy summary.tsv for this taxid and segment
    cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" -v segment="$segment_with_most_num_seqs" '$4==taxid && $5==segment {print $0}' >> $summary_outdir/summary.one-segment-per-taxid.tsv
done < <(tail -n +2 $TAXONOMIES_FILE | awk -F'\t' '{print $4}' | sort | uniq)
