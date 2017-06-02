# Single-cell sequencing of influenza-infected cells on Chromium 10X platform
This is designed to be a clean analysis of the infections with the "pure" influenza viruses that were sequenced on the Chromium 10X.

## Authors
Alistair Russell, Jesse Bloom, Cole Trapnell.
Alistair did all the experiments.
Alistair and Jesse did the analysis, Cole provided tips.

## Goal
The Chromium 10X platform was used to sequence A549 cells infected with influenza.
Questions of interest include:

  1. What is variability in viral transcription among cells?
   
  2. Are there viral or host patterns associated with different outcomes?

## Organization
The analysis consists of a set of iPython or RMarkdown notebooks.

The whole analysis can be run via the bash script [`run_analysis.bash`](run_analysis.bash).

This analysis consists of the following steps:

   1. The iPython notebook [`align_and_annotate.ipynb`][] aligns the reads, annotates the flu synonymous barcodes, and generates the cell-gene matrix.


## Input data
The BCL files that contain the deep sequencing data are on the Bloom lab `ngs` directory, and are linked to directly in [`align_and_annotate.ipynb`][].

Other input data are in `./data`:

  * `./data/flu_sequences/` contains the influenza genomes for both the wildtype A/WSN/1933 virus and the variants with synonymous mutations barcoding the 3' end of the mRNA, as taken from the Bloom lab reverse-genetics plasmids used to grow these viruses.

## Results and Conclusions
The read alignment, construction of the cell-gene matrices, and annotation of the synonymous barcodes at the 3' end of the flu mRNAs was successful.
The results are summarized in the iPython notebook [`align_and_annotate.ipynb`][].

## To-do
The next step is to analyze the data with `Monocle` by reading in the cell-gene matrices.


[`align_and_annotate.ipynb`]: align_and_annotate.ipynb
