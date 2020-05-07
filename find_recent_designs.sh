#!/bin/bash

# Find designs produced recently and write them to stdout.

# Define a timestamp such that everything after that is
# considered recent; everything within the last 7 days
# is considered recent
recent=$(date +%s -d "7 days ago")

# Iterate over design.last-complete-timestamp files
for f in $(find out/designs -name "design.last-complete-timestamp"); do
    ts=$(cat "$f")
    if (( ts > recent )); then
        # Timestamp is recent
        # Check for a design.tsv.0 file, indicating a design successfully completed
        # (but there may not be any design options in it)
        f_tsv=$(echo "$f" | sed 's/design.last-complete-timestamp/design.tsv.0/')
        if [ -f $f_tsv ]; then
            # Design was successful
            # Print this
            exp=$(echo "$f" | cut -d'/' -f3)
            taxid_and_segment=$(echo "$f" | cut -d'/' -f4)
            taxid=$(echo "$taxid_and_segment" | cut -d'_' -f1)
            segment=$(echo "$taxid_and_segment" | cut -d'_' -f2)
            echo -e "$exp\t$taxid\t$segment"
        fi
    fi
done
