# adapt-designs-continuous

Scripts for running ADAPT routinely.

**For more information on ADAPT and on how to run it, please see the [ADAPT repository](https://github.com/broadinstitute/adapt) on GitHub.**

## Overview

`./summary_across_taxa.sh` parses a list of all vertebrate viruses (`taxonomies/all-vertebrate.tsv`) and runs ADAPT across all of them, placing output in a `out/` directory.
It uses `run_taxon.sh`, which runs ADAPT for a single virus.
