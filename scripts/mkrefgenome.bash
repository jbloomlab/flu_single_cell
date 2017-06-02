#!/bin/bash 

# Make a cellranger reference genome for human cells with flu reads

# entire script exits if error on any command
set -e

ncpus=$(nproc --all)

echo "Creating reference genome"
ftp_site='ftp://ftp.ensembl.org/pub/release-87/'

# specify where get genomes by FTP
declare -A genomes
genomes['human']="${ftp_site}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# specify where get GTFs by FTP
declare -A gtfs
gtfs['human']="${ftp_site}gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz"

for species in "${!genomes[@]}"; do
    echo "Downloading $species genome..."
    wget -O - ${genomes[$species]} | gunzip -c > "${species}.fasta"
    echo "Downloading $species GTF..."
    wget -O - ${gtfs[$species]} | gunzip -c > temp.gtf
    echo "Filtering GTF for protein-coding genes..."
    cellranger mkgtf temp.gtf "${species}.gtf" --attribute=gene_biotype:protein_coding
    rm temp.gtf
done

echo "Concatenating human and flu genomes and GTF files..."
flugenome='./data/flu_sequences/flu-wsn.fasta'
flugtf='./data/flu_sequences/flu-wsn.gtf'
echo "Using flu genome and GTF in $flugenome and $flugtf"
if [ ! -f $flugenome ] || [ ! -f $flugtf ]; then
    echo "Cannot find flu genome or GTF files"
    exit 1
fi
cat human.fasta $flugenome > humanplusflu.fasta
cat human.gtf $flugtf > humanplusflu.gtf

echo "Using cellranger to make the reference genome..."
cellranger mkref --genome=humanplusflu --fasta=humanplusflu.fasta --genes=humanplusflu.gtf --nthreads=$ncpus
echo "Reference genome complete."

echo "Removing unneeded files."
for species in human humanplusflu; do
    rm "${species}.gtf"
    rm "${species}.fasta"
done
echo "Script complete."
