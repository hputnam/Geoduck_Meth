#!/bin/bash
## Job Name
#SBATCH --job-name=ct
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/062823-ct/

#grab deduplicated BAMS - example
#https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/032120-fds/EPI-230_S37_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.sorted.bam


/gscratch/srlab/programs/samtools-1.9/samtools merge \
Pg_merged.bam \
*.sorted.bam


perl /gscratch/srlab/programs/BS-Snper-master/BS-Snper.pl \
Pg_merged.bam \
--fa Panopea-generosa-v1.0.fa \
--output snp.candidate.out \
--methcg meth.cg \
--methchg meth.chg \
--methchh meth.chh \
--minhetfreq 0.1 \
--minhomfreq 0.85 \
--minquali 15 \
--mincover 10 \
--maxcover 1000 \
--minread2 2 \
--errorate 0.02 \
--mapvalue 20 \
>SNP.vcf 2>SNP.log

