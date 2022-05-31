# HPV-genotyping
--David B. Stern, Ph.D.--

The pipeline was developed to process Illumina sequencing data generated from CD Genomic's [HPV Capture Kit](https://www.cd-genomics.com/diseasepanel/products_70.html) and designed to run on NIAID's HPC LOCUS.

## Usage

The [scripts](./scripts) directory contains bash and R scripts to process the data.

The [ref](./ref) directory contains a fasta file of papillomavirus reference genomes from [PaVE](https://pave.niaid.nih.gov/), including non-reference genomes, borrowed from [HPV-EM](https://github.com/jin-wash-u/HPV-EM) which has nicely reformatted names.
This reference fasta needs to be indexed by Bowtie 2

All scripts rely on a file called `files.txt` which contains the names of all the samples to be processed.

- clean_map_abundance.sh: LOCUS UGE array script to clean reads with bbduk, map reads to the reference genomes with Bowtie 2, and estimate the relative abundance of each genotype using [msamtools](https://github.com/arumugamlab/msamtools).

- stats.sh: Collects coverage and pairwise ID statistics from the bam files using awk and [NanoStat](https://github.com/wdecoster/nanostat), and generates coverage plots using `collect_stats.R`. Should be run from the directory with the bam files. Be sure to check paths for reference fasta, index, and `collect_stats.R`

- collect_stats.R: R script to generate table of statistics and coverage plots. Run automatically with `stats.sh`. Requires the tidyverse R package.

- merge_msamtools_stats.R: R script to merge the output of msamtools and `collect_stats.R` for each sample.
