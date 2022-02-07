# adapt-designs-continuous

Scripts for running ADAPT routinely.

**For more information on ADAPT and on how to run it, please see the [ADAPT repository](https://github.com/broadinstitute/adapt) on GitHub.**

## Overview

`./summary_across_taxa.sh` parses a list of all vertebrate viruses (`taxonomies/all-vertebrate.tsv`) and runs ADAPT across all of them, placing output in a `out/` directory.
It uses `run_taxon.sh`, which runs ADAPT for a single virus.

## Vertebrate-associated designs

`out/all-vertebrate-2020.tsv.gz` is a compressed tar of designs output by ADAPT for 1,933 vertebrate-associated viral species, as reported and analyzed in the ADAPT paper.
These were designed in mid-2020.
The directory structure is `[DESIGN-TYPE]/[TAXID]_[SEGMENT]/[TIMESTAMP]/`.
The [adapt-analysis repository contains a file](https://github.com/broadinstitute/adapt-analysis/blob/main/adapt-benchmarking/panel-designs/all-vertebrate/data/design-results.tsv) with summary statistics on these designs; each row includes a species and the columns indicate the design (experiment) type, taxonomic ID, segment, timestamp, and more.
