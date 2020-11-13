Notes:
* NCBI's taxonomy identifier (https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) is useful
  for obtaining taxonomy IDs for species
* `all-human-associated-virus.tsv` was generated from a similar list used for CARMEN and recent CATCH designs -- all species with reported human infection
* `all-viral-with-ge10-seqs.tsv` was parsed from the NCBI viral accession (genome neighbor) list downloaded on Nov. 10, 2019 -- taking all species/segment combinations with >= 10 listed genomes. This used RefSeqs in the existing human-associated-viruses.tsv if corresponding species were already present in that, and I also added influenza viruses to the table. Note that this includes non-human-associated viruses as well.
* `all-vertebrate.tsv` was copied from /adapt-analysis/taxonomy-lists/2020-04/taxa.tsv and parsed as described in that directory. It contains all vertebrate-associated viruses (as of late April 2020).
