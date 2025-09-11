#!/bin/bash

print_help() {
    cat <<EOF
Usage: $0 -g <ucsc_id> -s <sequence>

A script to visualize the location of annealed repeats in the genomes of different organisms.

Options:
  -g, --genome STR    Input UCSC genome ID (required)
  -s, --sequence STR    Input annealed sequence (required)
  -h, --help          Show this help message

EOF
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --genome|-g) genome="$2"; shift ;;
    --sequence|-s) sequence="$2"; shift ;;
    --help|-h) print_help ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

current_dir="$PWD"

mkdir -p ${current_dir}/${genome}_genome_dir
cd ${current_dir}/${genome}_genome_dir


# Downloading genome
if [ -f "${current_dir}/${genome}_genome_dir/${genome}.fa" ]; then
    echo "Genome already downloaded"
else
    echo "Downloading genome..."
    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.gz ./
    
    echo "Uncompressing files..."
    gunzip ${genome}.fa.gz
fi


# Mapping AluY sequence to genome
eval "$(conda shell.bash hook)"
conda activate hisat2

echo "Building genome index..."
hisat2-build -p 10 ${genome}.fa ${genome}_idx

echo 'Mapping AluY sequence to genome...'
hisat2 --very-sensitive --non-deterministic --no-spliced-alignment -a --no-softclip -p 10 -x ${genome}_idx -c $sequence -S ${genome}_mapping.sam 2> mapping_stat.hisat
samtools view ${genome}_mapping.sam > ${genome}_mapping.txt

rm ${genome}_mapping.sam
conda deactivate

# Log setup
LOG_FILE="params.log"

# Log all inputs
{
    echo "===== Parameter Log ====="
    echo "Timestamp: $(date)"
    echo "User: $(whoami)"
    # echo "Host: $(hostname)"
    # echo "PID: $$"
    echo "Genome: ${genome}"
    echo "Annealed sequence: ${sequence}"
    echo "========================"
} | tee "$LOG_FILE"

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