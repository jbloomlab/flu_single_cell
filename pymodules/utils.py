"""Python functions used in 10X single-cell analysis.

Written by Jesse Bloom."""


import math
import os
import json
import time
import subprocess
import csv
import scipy.io
import scipy.sparse
import pandas
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from IPython.display import Image, display


def mergeCellGeneMatrices(mergedgenes, mergedcells, mergedmatrix,
        samples, cellbarcodes, cellannotations, genes, matrices,
        removeprefix=None):
    '''Merges cell-gene matrices for a variety of ``cellranger`` samples.

    Creates file suitable for analysis with ``Monocle``.

    Does **not** perform any type of normalization, just merges the
    counts.

    Args:
        `mergedgenes`, `mergedcells`, `mergedmatrix` (str)
            Names of created file containing the gene names and
            cell names (as ``*.tsv``) and gene-cell matrix (in
            Market Exchange Format).
        `samples` (list)
            List of samples that we are merging.
        `cellbarcodes`, `cellannotations`, `genes`, `matrices` (dict)
            Each is a dictionary keyed by the each sample in `samples`.
            The values are the names of the files giving the cell barcodes,
            cell annotations, genes, and matrices for that sample. The
            annotations can be `None`. The annotations file should have a
            header, and the first column should give the cell barcode.
        `removeprefix` (str or `None`)
            If not `None`, should specify a prefix associated with each
            gene. This prefix is removed in `mergedgenes`.

    Result: 
        Creates the files `mergedgenes`, `mergedcells`, and `mergedmatrix`.
    '''
    assert len(samples) == len(set(samples)) >= 1
    genelist = annotations = None
    cells = []
    matrixlist = []
    for sample in samples:
        sample_m = scipy.io.mmread(matrices[sample])
        with open(genes[sample]) as f:
            sample_genes = f.readlines()
        if genelist is None:
            genelist = sample_genes
        elif genelist != sample_genes:
            raise ValueError("different gene sets for different matrices")
        with open(cellbarcodes[sample]) as f:
            sample_barcodes = ['{0}-{1}\t{1}'.format(bc.strip().split('-')[0],
                    sample) for bc in f.readlines()]
        assert sample_m.shape == (len(sample_genes), len(sample_barcodes))
        matrixlist.append(sample_m)
        if cellannotations:
            sample_annotations = pandas.read_csv(cellannotations[sample],
                    sep='\t')
            assert len(sample_annotations) == len(sample_barcodes)
            assert sample_annotations.columns.values[0] == 'cellbarcode'
            if annotations is None:
                annotations = ['CellBarcode', 'Sample'] + list(
                        sample_annotations.columns.values[1 : ])
            elif (list(sample_annotations.columns.values[1 : ]) !=
                    annotations[2 : ]):
                raise ValueError("Different annotations for different samples")
            sample_annotations['cellbarcode'] = sample_annotations.apply(lambda row:
                   '{0}-{1}'.format(row['cellbarcode'].strip().split('-')[0], sample), axis = 1)
            annotated_barcodes = []
            for (bc, sample) in map(lambda line: line.split('\t'),
                    sample_barcodes):
                bc_annotations = sample_annotations[sample_annotations[
                        'cellbarcode'] == bc]
                assert len(bc_annotations) == 1
                annotated_bc = [bc, sample] + list(bc_annotations.values[0][1 : ])
                assert len(annotated_bc) == len(annotations)
                annotated_barcodes.append('\t'.join(map(str, annotated_bc)))
            cells += annotated_barcodes
        else:
            annotations = ['CellBarcode', 'Sample']
            cells += sample_barcodes

    if removeprefix:
        assert all([2 == gene.count(removeprefix) for gene in genelist])
        genelist = [gene.replace(removeprefix, '') for gene in genelist]
    
    matrix = scipy.sparse.hstack(matrixlist)
    assert matrix.shape == (len(genelist), len(cells))

    scipy.io.mmwrite(mergedmatrix, matrix)

    with open(mergedgenes, 'w') as f:
        f.write('gene_long_name\tgene_short_name\n{0}'.format(
                ''.join(genelist)))

    with open(mergedcells, 'w') as f:
        f.write('{0}\n{1}'.format('\t'.join(annotations), '\n'.join(cells)))


def showPDF(pdfs, width=None):
    '''Displays images in *pdfs*. Multiple images displayed side-by-side.'''
    png = '_temp.png'
    if not isinstance(pdfs, list):
        pdfs = [pdfs]
    subprocess.check_call(['convert', '-density', '134', '-trim'] + pdfs + ['+append', png])
    time.sleep(0.5)
    display(Image(png, width=width))
    os.remove(png)


def tilePlots(inplots, tiledplot, ncols):
    """Tile plots in `inplots` in `ncols` columns to create `tiledplot`."""
    nrows = int(math.ceil(len(inplots) / float(ncols)))
    rowplots = []
    for irow in range(nrows):
        rowplot = '{0}/_{1}{2}'.format(os.path.dirname(tiledplot), irow,
                os.path.basename(tiledplot))
        subprocess.check_call(['convert', '+append', '-density', '250',
                ] + inplots[ncols * irow : ncols * irow + ncols] + [
                rowplot])
        rowplots.append(rowplot)
    subprocess.check_call(['convert', '-append', '-density', '250',
            ] + rowplots + [tiledplot])
    for rowplot in rowplots:
        os.remove(rowplot) 


def demuxSummary(summaryfile, plotfile, title):
    """Plot summary of ``cellranger mkfastq`` demultiplexing.

    Args:
        `summaryfile` (str)
            ``qc_summary.json`` file produced by ``cellranger mkdfastq``
        `plotfile` (str)
            Name of created PDF plot file.
        `title` (str)
            Title to put at top of plot.
    """
    assert os.path.isfile(summaryfile), "Can't find {0}".format(summaryfile)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf', ("non PDF "
            "file extension on plotfile: {0}".format(plotfile))
    with open(summaryfile) as f:
        qc = json.load(f)['sample_qc']
    # get data for all lanes, not individual ones
    qc_df = pandas.DataFrame.from_dict(dict([(sample, d['all']) for 
            (sample, d) in qc.items()]), orient='index')
    stats = ['number_reads', 'barcode_exact_match_ratio', 
            'read1_q30_base_ratio', 'read2_q30_base_ratio']
    (fig, axes) = plt.subplots(ncols=len(stats), nrows=1, figsize=(11, 3))
    for (stat, ax) in zip(stats, axes):
        qc_df[stat].plot.bar(ax=ax)
        ax.set_ylabel(stat.replace('_', ' '))
    plt.suptitle(title)
    plt.tight_layout(w_pad=3)
    plt.savefig(plotfile)


def countSummary(samples, metricfiles, plotfile, celltypes=['']):
    """Plot summary of ``cellranger count`` analysis.
    
    Args:
        `samples` (list)
            List of names of samples.
        `metricfiles` (list)
            List of ``metrics_summary.csv`` files corresponding to `samples`.
        `plotfile` (str)
            Name of created PDF plot.
        `celltypes` (list)
            List of cell types for which reads are counted.
            Leave at default value of `['']` if only one cell type.
    """
    assert os.path.splitext(plotfile)[1].lower() == '.pdf', "non PDF extension"
    assert len(samples) == len(metricfiles)

    countstats = {}
    for (sample, metricfile) in zip(samples, metricfiles):
        with open(metricfile) as f:
            metrics = list(csv.reader(f))
        metrics[0] = [x.replace('Confidently ', '') + ' (%)' if '%' in y else x
                for (x, y) in zip(metrics[0], metrics[1])]
        if not countstats:
            columns = metrics[0]
        assert columns == metrics[0], "Inconsistent column names"
        countstats[sample] = [float(x.replace('%', '').replace(',', '')) for x 
                in metrics[1]]

    countstats = pandas.DataFrame.from_dict(countstats, orient='index').reindex(
            index=samples)
    countstats.columns = columns

    stats = ['Number of Reads', 'Estimated Number of Cells']
    if len(celltypes) > 1: 
        stats += ['GEMs with >1 Cell', 'Fraction GEMs with >1 Cell (%)']
    celltypestats = ['Fraction Reads in Cells (%)', 'Median UMI Counts per Cell', 
            'Median Genes per Cell']
    if len(celltypes) > 1: 
        stats += ['GEMs with >1 Cell', 'Fraction GEMs with >1 Cell (%)']
        celltypestats += ['Estimated Number of Cell Partitions']
    nplots = len(stats) + len(celltypestats)
    if nplots == 5:
        ncols = 5
    else:
        ncols = 4
    nrows = int(math.ceil(nplots / float(ncols)))
    (fig, axes) = plt.subplots(ncols=ncols, nrows=nrows, figsize=(11, 3.1 * nrows))
    for (stat, ax) in zip(stats + celltypestats, axes.ravel()):
        addlegend = False
        if stat in celltypestats:
            df = countstats[['{0} {1}'.format(celltype, stat).strip() for
                    celltype in celltypes]].copy()
            df.columns = celltypes
            if len(celltypes) > 1:
                addlegend = True
            if stat in ['Estimated Number of Cell Partitions', 
                    'Fraction Reads in Cells (%)']:
                # Manually set to zero if there aren't reads for this cell type;
                # cellranger mis-estimates these states in this case. Use 50
                # as a cutoff for actually having reads for cell type.
                for celltype in celltypes:
                    rows = (50 > countstats[
                            '{0} Median UMI Counts per Cell'.format(
                            celltype).strip()])
                    df[celltype][rows.values] = 0
        else:
            df = countstats[stat]
        df.plot.bar(ax=ax)
        ax.set_ylabel(stat, fontsize=10)
        if addlegend:
            ax.set_ylim(0, 1.4 * df.values.max())
            ax.legend(borderaxespad=0, fontsize=10)
        elif ax.get_legend(): 
            ax.get_legend().set_visible(False)
    plt.tight_layout(w_pad=3)
    plt.savefig(plotfile)


def multiFlowCellMRO(sample, sample_demuxdirs, refgenome):
    """Creates `cellranger count` MRO file for multi-flow cell samples.

    The de-multiplexed reads for `sample` are listed in `sample_demuxdirs`.

    See here:
    https://support.10xgenomics.com/single-cell/software/pipelines/latest/advanced/multi-flowcell
    """
    return ('''
            @include "sc_rna_counter_cs.mro"
                                            
            call SC_RNA_COUNTER_CS(
                sample_id = "{0}",
                no_secondary_analysis = true,
                force_cells = null,
                sample_def = ['''.format(sample)
            + ','.join(['''
                    {{
                        "fastq_mode": "ILMN_BCL2FASTQ",
                        "gem_group": null,
                        "lanes": null,
                        "read_path": "{0}",
                        "sample_indices": ["any"],
                        "sample_names": ["{1}"]
                    }}'''.format(os.path.abspath(demuxdir), sample)
                    for demuxdir in sample_demuxdirs])
            + '''
                ],
                sample_desc = "",
                reference_path = "{0}",
                recovered_cells = null,
            )'''.format(os.path.abspath(refgenome)))
