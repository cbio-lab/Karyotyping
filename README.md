# AluY-distribution

## Description
Pipeline for visualisation AluY repeats location in genomes of different organisms.

### 1. Installation
Download src folder, which contains scripts and envirounments.
Click on ..

Install Miniconda (https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

For `Linux x86_64 (amd64)`:
> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

> bash Miniforge3-Linux-x86_64.sh

> conda init

If initialization failed, it can be done by calling `conda` with its full path, with something like
> ~/miniforge3/bin/conda init

Install the `hisat2` and `R_AluY_pkgs` envirounments from `.yml` files via command:
> conda env create -n hisat2 --file hisat2.yml

> conda env create -n R_AluY_pkgs --file R_AluY_pkgs.yml

### 2. Running
For running script you should find UCSC genome identificator (`ucsc_id`) of your genome of interest.

Go to https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=2524985107_dAb8wrCzhEaTXfQAKs3pAonzspgo, next `Tools` -> `Table Browser`

![](https://github.com/kanaevavera/AluY-distribution/blob/main/Ucsc_example.png)

Choose Clade, Genome and Assembly. `ucsc_id` is assebly ID after `/`, like `calJac4` in a picture.

From your current directory and `base` envirounment run `download_genome.sh`, that located in `src` folder. Also you need to provide to script `ucsc_id`:

> bash src/download_genome.sh ucsc_id

Plots for each chromosome will be located in `./plots/chromosome_plots.zip`. Main plot with all chromosomes in one picture will be in `./plots/ucsc_id_main.pdf`
