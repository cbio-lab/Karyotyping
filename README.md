# AluY-distribution

## Description
Pipeline for visualisation AluY repeats location in genomes of different organisms.

### 1. Download target genome from NCBI
Find genome ID.. run:

> wget /link/to/ncbi

### 2. Download genome info
Centomers position, chr ID, chr length

### 3. Map AluY to target genome
Activate `hisat2 env` and run:

> hisat2-build genome.fna hisat2_idx
> hisat2 --very-sensitive --non-deterministic --no-spliced-alignment -a --no-softclip -p 10 -x hisat2_idx -c target_seq -S mapping.sam 2> mapping_stat.hisat
> samtool view mapping.sam > view_map.txt

### 4. Visualisation
Run R-script
