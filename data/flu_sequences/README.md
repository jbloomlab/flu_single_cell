# Sequences for the influenza strains

This directory contains the sequences for the wildtype and synonymously barcoded A/WSN/1933 (H1N1) influenza strains used in the study.

There are maps for the 8 plasmids (pHW181-PB2, ...) used to generate the wildtype virus and the 8 plasmids (pHW181-PB2-syn, ...) used to generate the synonymously barcoded virus. These are the Genbank plasmid maps from the Bloom lab plasmid logs.

The Python script `make_fasta_and_gtf.py` then makes FASTA genome files and GTF files. These GTF files just annotate the entire viral RNAs with *gene_biotype* of *"vRNA"* for viral RNA. For the 3' end counting used by the Chromium, just annotating the viral RNAs is sufficient.

Running this Python script produces `flu-wsn.fasta`, `flu-wsn.gtf`, `flu-wsn-syn.fasta`, and `flu-wsn-syn.gtf`.
