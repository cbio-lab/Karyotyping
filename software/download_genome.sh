#!/bin/bash

genome=$1
current_dir="$PWD"


mkdir -p ${current_dir}/${genome}_genome_dir
cd ${current_dir}/${genome}_genome_dir

echo "Downloading genome..."
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.gz ./

echo "Uncompressing files..."
gunzip ${genome}.fa.gz


# Mapping AluY sequence to genome
eval "$(conda shell.bash hook)"
conda activate hisat2

echo "Building genome index..."
hisat2-build ${genome}.fa ${genome}_idx

echo 'Mapping AluY sequence to genome...'
hisat2 --very-sensitive --non-deterministic --no-spliced-alignment -a --no-softclip -p 10 -x ${genome}_idx -c CGGTGGCTCAAGCCTGTAATCCCAGCACTTTG -S ${genome}_mapping.sam 2> mapping_stat.hisat
samtools view ${genome}_mapping.sam > ${genome}_mapping.txt
rm ${genome}_mapping.sam
conda deactivate
echo 'DONE!'
echo "-----------------------------------------"


# Visualisation coverage
echo "Visualisation..."
conda activate R_AluY_pkgs

mkdir -p ${current_dir}/plots
cd ${current_dir}/plots

Rscript ${current_dir}/src/AluY_distribution.R $genome $current_dir

zip -m chromosome_plots.zip *_chr*.pdf

conda deactivate
echo 'DONE!' 