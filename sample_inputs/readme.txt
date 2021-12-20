sample_inputs/

This directory contains 2 sample datasets of aligned SARS-Cov-2 sequences, a reference genome to use as a root sequence, and their associated metadata.

All genome sequences and metadata were obtained from the GISAID database: https://www.gisaid.org

Files included:

jan31.fasta - Masked multiple sequence alignment of all SARS-CoV-2 sequences in GISAID whos collection date is prior to January 31st, 2020.

feb29.fasta - Masked multiple sequence alignment of all SARS-CoV-2 sequences in GISAID whos collection date is prior to Febraury 29th, 2020.

*For both of these datasets, we filtered out sequences with ambiguous collection dates (year only, or year and month only), as well as any sequences not from humans (Pangolin, etc.).

metadata.tsv - The file is a reduced version of the one available from GISAID. It contains sampling date and location metadata for every sequence in the jan31 and feb29 datasets.