#!/bin/bash -e

#$ -N hpv_map
#$ -cwd
#$ -V
#$ -j y
#$ -pe threaded 16
#$ -t 1-5


fname=$(sed "${SGE_TASK_ID}q;d" files.txt)

# load conda environment with software
eval "$(conda shell.bash hook)"
conda activate sterndb

#reads in directory called fastq
bbduk.sh -Xmx20g in1=fastq/${fname}_R1.fastq.gz \
                 in2=fastq/${fname}_R2.fastq.gz \
                 out1=fastq/${fname}_R1.clean.fq \
                 out2=fastq/${fname}_R2.clean.fq \
                 ref=/hpcdata/bcbb/stern/databases/adapters.fa \
                 ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 tbo minlen=75

bowtie2 --threads 16 -x ref/combined_pave_hpv.fa \
        -1 fastq/${fname}_R1.clean.fq \
        -2 fastq/${fname}_R2.clean.fq \
        --very-sensitive-local \
        --no-mixed \
        --no-discordant \
        -k 200 | \
        samtools view -bh -@ 16 | \
        samtools sort -n -@ 16 - > bams/${fname}.readsorted.bam #msamtools requires reads sorted by read name

# select all alignments where the read is at least 80bp long,
# the percent identity of alignment is >=95%,
# and at least 80% of the read is aligned.
msamtools filter -b -l 80 -p 95 -z 80 bams/${fname}.readsorted.bam \
        | msamtools profile --multi=proportional --label=bams/$fname --unit=rel -o bams/$fname.pave.profile.txt.gz -

#also generate bams sorted by reference fasta
samtools sort -@ 16 bams/${fname}.readsorted.bam > bams/${fname}.sorted.bam
