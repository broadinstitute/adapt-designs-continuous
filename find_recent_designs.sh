#!/bin/bash

# Find designs produced recently and write them to stdout.

# If set to true, require that the job have produced at least one
# design in order to be considered completed (if false, a job
# is completed if it successfully ran, even if it could not produce
# any designs)
REQUIRE_DESIGNS=true

# Define a timestamp such that everything after that is
# considered recent; everything within the last 150 days
# is considered recent
recent=$(date +%s -d "150 days ago")

# Write to a tmp file, and then find unique lines
tmpf=$(mktemp)

# Iterate over design.last-complete-timestamp files
for f in $(find out/designs -name "design.last-complete-timestamp"); do
    ts=$(cat "$f")
    if (( ts > recent )); then
        # Timestamp is recent
        # Check for a design.tsv.0 file, indicating a design successfully completed
        # (but there may not be any design options in it)
        f_tsv=$(echo "$f" | sed 's/design.last-complete-timestamp/design.tsv.0/')
        if [ -f $f_tsv ]; then
            # Design job was successful
            if [ "$REQUIRE_DESIGNS" = true ]; then
                if [[ $(wc -l < $f_tsv) -eq 1 ]]; then
                    # There is only one line in the .tsv file (the header); there
                    # are no actual designs
                    continue
                fi
            fi
            # Print this
            exp=$(echo "$f" | cut -d'/' -f3)
            taxid_and_segment=$(echo "$f" | cut -d'/' -f4)
            taxid=$(echo "$taxid_and_segment" | cut -d'_' -f1)
            segment=$(echo "$taxid_and_segment" | cut -d'_' -f2)
            echo -e "$exp\t$taxid\t$segment" >> $tmpf
        fi
    fi
done

# A design can appear multiple times, e.g., if run at
# multiple recent timestamps; filter to only output one
cat $tmpf | sort | uniq

rm $tmpf
