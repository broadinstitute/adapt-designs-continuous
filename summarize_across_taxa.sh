#!/bin/bash

# Summarize designs for all taxonomies.

JACCARD_THRES=0.5
TAXONOMIES_FILE="taxonomies/all-viral-with-ge10-seqs.tsv"
IGNORE_BEFORE=1573900000

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
    summary=$(python summarize_design_on_taxon.py $taxonomy_outdir $JACCARD_THRES --ignore-before $IGNORE_BEFORE 2> /dev/null)
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
# Pick the one with the most number of input sequences (before
# curation), breaking ties arbitrarily
echo -n "" > $summary_outdir/summary.one-segment-per-taxid.tsv
while read -r taxid; do
    # Sort all segments by the number of sequences, and pick the one with the
    # most
    segment_with_most_num_seqs=$(cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" '$4==taxid && $7=="taxon" {print $5"\t"$11}' | sort -k2nr | head -n 1 | awk -F'\t' '{print $1}')

    # Copy summary.tsv for this taxid and segment
    cat $summary_outdir/summary.tsv | awk -F'\t' -v taxid="$taxid" -v segment="$segment_with_most_num_seqs" '$4==taxid && $5==segment {print $0}' >> $summary_outdir/summary.one-segment-per-taxid.tsv
done < <(tail -n +2 $TAXONOMIES_FILE | awk -F'\t' '{print $4}' | sort | uniq)

# Write a list of the taxa that did not yield designs
echo -n "" > $summary_outdir/summary.taxa-not-designed.tsv
while read -r taxid; do
    cat $TAXONOMIES_FILE | awk -F'\t' -v taxid="$taxid" '$4==taxid {print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq >> $summary_outdir/summary.taxa-not-designed.tsv
done < <(comm -23 <(cat $TAXONOMIES_FILE | tail -n +2 | awk -F'\t' '{print $4}' | sort | uniq) <(cat $summary_outdir/summary.one-segment-per-taxid.tsv | awk -F'\t' '{print $4}' | sort | uniq))
