#!/bin/bash

# Design for all taxonomies.

NJOBS=8
TAXONOMIES_FILE="taxonomies/all-vertebrate.tsv"
RECENT_DESIGNS_FILE="recent-designs.tsv"
MEMORY_LIMIT="300000000"    # kilobytes

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
            if [ -f $RECENT_DESIGNS_FILE ]; then
                # Check if this was recently designed
                exp="${specific}_${obj}"
                if grep -Pq "^${exp}\t${taxid}\t${segment}$" $RECENT_DESIGNS_FILE; then
                    # Design exists as a recent design; skip it
                    continue
                fi
            fi

            echo "./run_taxon.sh $TAXONOMIES_FILE $taxid '$segment' $specific $obj" >> $commands
        done
    done
done < <(tail -n +2 $TAXONOMIES_FILE)

# Shuffle the commands so ones from the same family are not together; this
# way the more resource-intensive commands are spread out
commands_shuf=$(mktemp)
cat $commands | shuf > $commands_shuf

# Run commands in parallel
# Restart jobs 3 times on a delay of 180 sec
parallel --jobs $NJOBS --retries 1 --delay 60 --no-notice < $commands_shuf

rm $commands
rm $commands_shuf
