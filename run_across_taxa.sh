#!/bin/bash

# Design for all taxonomies.

NJOBS=32
TAXONOMIES_FILE="taxonomies/all-vertebrate.tsv"
MEMORY_LIMIT="100000000"    # kilobytes

# Limit memory usage for each process
ulimit -m $MEMORY_LIMIT
ulimit -v $MEMORY_LIMIT

# Create commands for each taxonomy
commands=$(mktemp)
echo -n "" > $commands
while read -r taxonomy; do
    taxid=$(echo "$taxonomy" | cut -f4)
    segment=$(echo "$taxonomy" | cut -f5)

    # Write for each of 4 "experiments"
    for specific in "specific" "nonspecific"; do
        for obj in "max-activity" "min-guides"; do
            echo "./run_taxon.sh $TAXONOMIES_FILE $taxid '$segment' $specific $obj" >> $commands
        done
    done
done < <(tail -n +2 $TAXONOMIES_FILE)

# Run commands in parallel
# Restart jobs 3 times on a delay of 60 sec
parallel --jobs $NJOBS --retries 3 --delay 60 --no-notice < $commands

rm $commands
