#!/bin/bash

# Design for all taxonomies.

NJOBS=8
TAXONOMIES_FILE="taxonomies/taxonomies.tsv"

# Create commands for each taxonomy
commands="/tmp/commands-run-taxa"
echo -n "" > $commands
while read -r taxonomy; do
    taxid=$(echo "$taxonomy" | awk -F'\t' '{print $4}')
    segment=$(echo "$taxonomy" | awk -F'\t' '{print $5}')
    refaccs=$(echo "$taxonomy" | awk -F'\t' '{print $6}')
    echo "./run_taxon.sh $taxid $segment $refaccs" >> $commands
done < <(tail -n +2 $TAXONOMIES_FILE)

# Run commands in parallel
parallel --jobs $NJOBS --no-notice < $commands

rm $commands
