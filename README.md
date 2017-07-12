# Single-cell sequencing of influenza-infected cells
This repository is an analysis of the transcriptional dynamics of influenza virus infection at the single-cell level.

Briefly, A549 cells were infected with A/WSN/1933 influenza virus at a low MOI.
The transcriptomes of the cells were then sequenced on the [Chromium 10X platform](https://www.10xgenomics.com/single-cell/).
The virus used was a mix of wildtype and virus with synonymous "barcodes" near the 3' end to help enable identification of co-infection.
The virus populations were also "pure" in the sense that they were passaged at low MOI to avoid accumulation of defective particles.
Multiple timepoints were collected and analyze.

## Authors
Alistair Russell, [Cole Trapnell](http://cole-trapnell-lab.github.io/), [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Organization of this repository
The analysis is performed by a set of Jupyter notebooks.

1. The Python Jupyter notebook [align_and_annotate.ipynb][] demultiplexes and aligns the reads, annotates the flu synonymous barcodes, and generates the cell-gene matrix. It requires installation of [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), which performs the demultiplexing and alignment. It also uses custom Python and bash scripts found in the `./scripts/` subdirectory, and requires installation of a few common Python modules. The notebook describes the software versions used. 

2. The R Jupyter notebook [monocle_analysis.ipynb][] analyzes the cell-gene matrix, making use of [Monocle][]. It generates most of the figures as well. The versions of R and associated packages that are used is described in the notebook.

The whole analysis can be run via the bash script [run_analysis.bash](run_analysis.bash).

## Input data
In addition to the notebooks / scripts themselves, the following input data is used:

1. The BCL files that contain the deep sequencing data are on the Bloom lab `ngs` directory, and are linked to directly in [align_and_annotate.ipynb][].

2. `./data/flu_sequences/` contains the influenza genomes for both the wildtype A/WSN/1933 virus and the variants with synonymous mutations barcoding the 3' end of the mRNA, as taken from the Bloom lab reverse-genetics plasmids used to grow these viruses.

## Results and Conclusions
All of the output from the analyses are written to the `./results/` subdirectory.

The Jupyter notebooks [align_and_annotate.ipynb][] [monocle_analysis.ipynb] contain a description of the results.


[align_and_annotate.ipynb]: align_and_annotate.ipynb
[monocle_analysis.ipynb]: monocle_analysis.ipynb
[Monocle]: http://cole-trapnell-lab.github.io/monocle-release/
