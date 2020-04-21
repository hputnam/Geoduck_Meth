# Library Prep Methods
[https://hputnam.github.io/Putnam_Lab_Notebook/Geoduck_RRBS_Library_Prep/](https://hputnam.github.io/Putnam_Lab_Notebook/Geoduck_RRBS_Library_Prep/)


---


# Trimming

```
#!/bin/bash
## Job Name
#SBATCH --job-name=TrimGpgnrMeth
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-23:30:00
## Memory per node
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200320


%%bash

# Run TrimGalore using default adapter identification

reads_dir="/gscratch/scrubbed/strigg/analyses/20200320/FASTQS/"

find ${reads_dir}*_R1_001.fastq.gz | \
xargs basename -s _R1_001.fastq.gz | \
xargs -I{} /gscratch/srlab/strigg/bin/anaconda3/bin/trim_galore \
--cores 8 \
--output_dir /gscratch/scrubbed/strigg/analyses/20200320/TG_FASTQS \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/strigg/analyses/20200320/TG_FASTQS/FastQC \
--threads 28" \
--clip_R1 8 \
--clip_R2 8 \
--three_prime_clip_R1 8 \
--three_prime_clip_R2 8 \
--path_to_cutadapt /gscratch/srlab/strigg/bin/anaconda3/bin/cutadapt \
${reads_dir}{}_R1_001.fastq.gz \
${reads_dir}{}_R2_001.fastq.gz 


# Run multiqc for WGBS and MBD samples

/gscratch/srlab/strigg/bin/anaconda3/bin/multiqc \
/gscratch/scrubbed/strigg/analyses/20200320/TG_FASTQS/FastQC/.



```
---

**output from TrimGalore and MultiQC**  

[https://gannet.fish.washington.edu/metacarcinus/Pgenerosa/analyses/20200320/](https://gannet.fish.washington.edu/metacarcinus/Pgenerosa/analyses/20200320/)

---
# Genome Preparation

**File info**     

```
914M Dec 17 12:54 Panopea-generosa-v1.0.fa
b7b64f0ce79499d79a865348658d2e49  Panopea-generosa-v1.0.fa
mox location: /gscratch/srlab/sr320/data/geoduck/v01
```

---

# Bismark Alignment

```
#!/bin/bash
## Job Name
#SBATCH --job-name=fds
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=20-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/032120-fds/


#job following new trim by Shelly on 032020


# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/scrubbed/strigg/analyses/20200320/TG_FASTQS/"
genome_folder="/gscratch/srlab/sr320/data/geoduck/v01/"

source /gscratch/srlab/programs/scripts/paths.sh


#
#${bismark_dir}/bismark_genome_preparation \
#--verbose \
#--parallel 28 \
#--path_to_aligner ${bowtie2_dir} \
#${genome_folder}


find ${reads_dir}*_R1_001_val_1.fq.gz \
| xargs basename -s _R1_001_val_1.fq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
-1 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_R1_001_val_1.fq.gz \
-2 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_R2_001_val_2.fq.gz \




find *.bam | \
xargs basename -s .bam | \
xargs -I{} ${bismark_dir}/deduplicate_bismark \
--bam \
--paired \
{}.bam



${bismark_dir}/bismark_methylation_extractor \
--bedGraph --counts --scaffolds \
--multicore 14 \
--buffer_size 75% \
*deduplicated.bam



# Bismark processing report

${bismark_dir}/bismark2report

#Bismark summary report

${bismark_dir}/bismark2summary



# Sort files for methylkit and IGV

find *deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} ${samtools} \
sort --threads 28 {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ 16" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ 28 {}.sorted.bam





find *deduplicated.bismark.cov.gz \
| xargs basename -s _R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
| xargs -I{} ${bismark_dir}/coverage2cytosine \
--genome_folder ${genome_folder} \
-o {} \
--merge_CpG \
--zero_based \
{}_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz


#creating bedgraphs post merge

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x.bedgraph
done



for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x.bedgraph
done


#creating tab files with raw count for glms

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x.tab
done


for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x.tab
done
```
---

**All Output from Bismark Aligment**

[https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/032120-fds/](https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/032120-fds/)


---

# DMR analysis

```
# Directories and programs
working_dir="/Users/strigg/Desktop/20200327/"


########################
## Copy data locally ###
########################

# (base) ostrich:data strigg$ rsync --archive --progress --verbose /Volumes/web/seashell/bu-mox/scrubbed/032120-fds/*_5x.tab .
# building file list ... 
# 52 files to consider
# EPI-103_S27_L005_5x.tab
#    58567183 100%   38.53MB/s    0:00:01 (xfer#1, to-check=51/52)
# EPI-104_S28_L005_5x.tab
#    96730762 100%   32.31MB/s    0:00:02 (xfer#2, to-check=50/52)
# EPI-111_S29_L005_5x.tab
#    86734405 100%   27.61MB/s    0:00:02 (xfer#3, to-check=49/52)
# EPI-113_S30_L005_5x.tab
#    81726800 100%   25.86MB/s    0:00:03 (xfer#4, to-check=48/52)
# EPI-119_S31_L005_5x.tab
#    95243953 100%   36.64MB/s    0:00:02 (xfer#5, to-check=47/52)
# EPI-120_S32_L005_5x.tab
#    79782602 100%   30.58MB/s    0:00:02 (xfer#6, to-check=46/52)
# EPI-127_S33_L005_5x.tab
#    68596619 100%   29.72MB/s    0:00:02 (xfer#7, to-check=45/52)
# EPI-128_S34_L005_5x.tab
#    79254997 100%   34.45MB/s    0:00:02 (xfer#8, to-check=44/52)
# EPI-135_S35_L005_5x.tab
#    85616235 100%   31.71MB/s    0:00:02 (xfer#9, to-check=43/52)
# EPI-136_S36_L005_5x.tab
#    83193844 100%   29.04MB/s    0:00:02 (xfer#10, to-check=42/52)
# EPI-143_S37_L005_5x.tab
#    58629953 100%   23.17MB/s    0:00:02 (xfer#11, to-check=41/52)
# EPI-145_S38_L005_5x.tab
#    79088643 100%   31.85MB/s    0:00:02 (xfer#12, to-check=40/52)
# EPI-151_S2_L002_5x.tab
#    46073182 100%   29.00MB/s    0:00:01 (xfer#13, to-check=39/52)
# EPI-152_S3_L002_5x.tab
#    59585344 100%   27.58MB/s    0:00:02 (xfer#14, to-check=38/52)
# EPI-153_S4_L002_5x.tab
#    37898008 100%   35.36MB/s    0:00:01 (xfer#15, to-check=37/52)
# EPI-154_S5_L002_5x.tab
#    39451172 100%   36.53MB/s    0:00:01 (xfer#16, to-check=36/52)
# EPI-159_S6_L002_5x.tab
#    21036348 100%   32.94MB/s    0:00:00 (xfer#17, to-check=35/52)
# EPI-160_S7_L002_5x.tab
#    56408712 100%   25.99MB/s    0:00:02 (xfer#18, to-check=34/52)
# EPI-161_S8_L002_5x.tab
#    29480720 100%   33.27MB/s    0:00:00 (xfer#19, to-check=33/52)
# EPI-162_S9_L002_5x.tab
#    29542604 100%   17.51MB/s    0:00:01 (xfer#20, to-check=32/52)
# EPI-167_S10_L002_5x.tab
#    48017024 100%   25.33MB/s    0:00:01 (xfer#21, to-check=31/52)
# EPI-168_S11_L002_5x.tab
#    37374825 100%   19.96MB/s    0:00:01 (xfer#22, to-check=30/52)
# EPI-169_S12_L002_5x.tab
#    33823538 100%   19.26MB/s    0:00:01 (xfer#23, to-check=29/52)
# EPI-170_S13_L002_5x.tab
#    43871714 100%   22.64MB/s    0:00:01 (xfer#24, to-check=28/52)
# EPI-175_S14_L003_5x.tab
#    53554347 100%   23.11MB/s    0:00:02 (xfer#25, to-check=27/52)
# EPI-176_S15_L003_5x.tab
#    75152227 100%   33.48MB/s    0:00:02 (xfer#26, to-check=26/52)
# EPI-181_S16_L003_5x.tab
#    60226380 100%   35.04MB/s    0:00:01 (xfer#27, to-check=25/52)
# EPI-182_S17_L003_5x.tab
#    78826142 100%   28.82MB/s    0:00:02 (xfer#28, to-check=24/52)
# EPI-184_S18_L003_5x.tab
#    33143094 100%   21.46MB/s    0:00:01 (xfer#29, to-check=23/52)
# EPI-185_S19_L003_5x.tab
#     6999125 100%    9.07MB/s    0:00:00 (xfer#30, to-check=22/52)
# EPI-187_S20_L003_5x.tab
#    46817624 100%   23.47MB/s    0:00:01 (xfer#31, to-check=21/52)
# EPI-188_S21_L003_5x.tab
#    28751557 100%   17.01MB/s    0:00:01 (xfer#32, to-check=20/52)
# EPI-193_S22_L003_5x.tab
#    48662320 100%   23.77MB/s    0:00:01 (xfer#33, to-check=19/52)
# EPI-194_S23_L003_5x.tab
#    71100625 100%   24.24MB/s    0:00:02 (xfer#34, to-check=18/52)
# EPI-199_S24_L003_5x.tab
#    18494572 100%   13.81MB/s    0:00:01 (xfer#35, to-check=17/52)
# EPI-200_S25_L003_5x.tab
#    44112209 100%   29.71MB/s    0:00:01 (xfer#36, to-check=16/52)
# EPI-205_S26_L004_5x.tab
#    19259990 100%   19.54MB/s    0:00:00 (xfer#37, to-check=15/52)
# EPI-206_S27_L004_5x.tab
#    36301568 100%   18.52MB/s    0:00:01 (xfer#38, to-check=14/52)
# EPI-208_S28_L004_5x.tab
#    50451877 100%   22.74MB/s    0:00:02 (xfer#39, to-check=13/52)
# EPI-209_S29_L004_5x.tab
#    61728708 100%   35.19MB/s    0:00:01 (xfer#40, to-check=12/52)
# EPI-214_S30_L004_5x.tab
#    59565624 100%   26.34MB/s    0:00:02 (xfer#41, to-check=11/52)
# EPI-215_S31_L004_5x.tab
#    36880711 100%   32.72MB/s    0:00:01 (xfer#42, to-check=10/52)
# EPI-220_S32_L004_5x.tab
#    18214672 100%   32.41MB/s    0:00:00 (xfer#43, to-check=9/52)
# EPI-221_S33_L004_5x.tab
#    21247499 100%   18.51MB/s    0:00:01 (xfer#44, to-check=8/52)
# EPI-226_S34_L004_5x.tab
#    13707440 100%   27.87MB/s    0:00:00 (xfer#45, to-check=7/52)
# EPI-227_S35_L004_5x.tab
#    17050257 100%   17.54MB/s    0:00:00 (xfer#46, to-check=6/52)
# EPI-229_S36_L004_5x.tab
#    31275428 100%   17.10MB/s    0:00:01 (xfer#47, to-check=5/52)
# EPI-230_S37_L004_5x.tab
#    40107828 100%   21.56MB/s    0:00:01 (xfer#48, to-check=4/52)
# EPI-41_S38_L005_5x.tab
#    35461930 100%   19.72MB/s    0:00:01 (xfer#49, to-check=3/52)
# EPI-42_S39_L005_5x.tab
#    59467789 100%   25.04MB/s    0:00:02 (xfer#50, to-check=2/52)
# EPI-43_S40_L005_5x.tab
#    36347431 100%   27.58MB/s    0:00:01 (xfer#51, to-check=1/52)
# EPI-44_S41_L005_5x.tab
#    17482104 100%   23.48MB/s    0:00:00 (xfer#52, to-check=0/52)

# sent 2556436095 bytes  received 1164 bytes  38442665.55 bytes/sec
# total size is 2556120265  speedup is 1.00

#########################################
## convert 5x.tab files to allc format ##
#########################################

## 1. Scaffold
## 2. Position (subtract 1 from the end position to get same position as allc and cov files)
## 3. strand
## 4. context
## 5. number mC
## 6. total C (numer of mC ($5) + numer of unmC ($6))

%%bash

for f in merged_tab/*_5x.tab
do
awk -F"\t" '{print $1"\t"$3-1"\t""+""\t""CGA""\t"$5"\t"$5+$6"\t"1}' ${f} \
> $(basename ${f} _5x.tab)_5xmerg_allc.tsv
done

#########################################
#### run DMRfind for day 10 comparison ##
#########################################

methylpy DMRfind \
--allc-files EPI-103_S27_L005_5xmerg_allc.tsv \
EPI-104_S28_L005_5xmerg_allc.tsv \
EPI-111_S29_L005_5xmerg_allc.tsv \
EPI-113_S30_L005_5xmerg_allc.tsv \
EPI-119_S31_L005_5xmerg_allc.tsv \
EPI-120_S32_L005_5xmerg_allc.tsv \
EPI-135_S35_L005_5xmerg_allc.tsv \
EPI-136_S36_L005_5xmerg_allc.tsv \
EPI-127_S33_L005_5xmerg_allc.tsv \
EPI-128_S34_L005_5xmerg_allc.tsv \
EPI-143_S37_L005_5xmerg_allc.tsv \
EPI-145_S38_L005_5xmerg_allc.tsv \
--samples EPI-103 EPI-104 EPI-111 EPI-113 EPI-119 EPI-120 EPI-135 EPI-136 EPI-127 EPI-128 EPI-143 EPI-145 \
--mc-type "CGN" \
--num-procs 8 \
--output-prefix day10_AllpH_DMR250bp_cov5x \
--dmr-max-dist 250 \
--min-num-dms 3 

#########################################
## Run DMRfind for all ambient samples ##
#########################################

methylpy DMRfind \
--allc-files EPI-41_S38_L005_5xmerg_allc.tsv \
EPI-42_S39_L005_5xmerg_allc.tsv \
EPI-43_S40_L005_5xmerg_allc.tsv \
EPI-44_S41_L005_5xmerg_allc.tsv \
EPI-119_S31_L005_5xmerg_allc.tsv \
EPI-120_S32_L005_5xmerg_allc.tsv \
EPI-135_S35_L005_5xmerg_allc.tsv \
EPI-136_S36_L005_5xmerg_allc.tsv \
EPI-151_S2_L002_5xmerg_allc.tsv \
EPI-152_S3_L002_5xmerg_allc.tsv \
EPI-153_S4_L002_5xmerg_allc.tsv \
EPI-154_S5_L002_5xmerg_allc.tsv \
EPI-181_S16_L003_5xmerg_allc.tsv \
EPI-182_S17_L003_5xmerg_allc.tsv \
EPI-184_S18_L003_5xmerg_allc.tsv \
EPI-185_S19_L003_5xmerg_allc.tsv \
--samples EPI-41 EPI-42 EPI-43 EPI-44 EPI-119 EPI-120 EPI-135 EPI-136 EPI-151 EPI-152 EPI-153 EPI-154 EPI-181 EPI-182 EPI-184 EPI-185 \
--mc-type "CGN" \
--num-procs 8 \
--output-prefix amb_AllTimes_DMR250bp_cov5x \
--dmr-max-dist 250 \
--min-num-dms 3 

#########################################
#### Run DMRfind for Day 135 samples ####
#########################################

methylpy DMRfind \
--allc-files EPI-151_S2_L002_5xmerg_allc.tsv \
EPI-152_S3_L002_5xmerg_allc.tsv \
EPI-153_S4_L002_5xmerg_allc.tsv \
EPI-154_S5_L002_5xmerg_allc.tsv \
EPI-159_S6_L002_5xmerg_allc.tsv \
EPI-160_S7_L002_5xmerg_allc.tsv \
EPI-161_S8_L002_5xmerg_allc.tsv \
EPI-162_S9_L002_5xmerg_allc.tsv \
EPI-167_S10_L002_5xmerg_allc.tsv \
EPI-168_S11_L002_5xmerg_allc.tsv \
EPI-169_S12_L002_5xmerg_allc.tsv \
EPI-170_S13_L002_5xmerg_allc.tsv \
--samples EPI-151 EPI-152 EPI-153 EPI-154 EPI-159 EPI-160 EPI-161 EPI-162 EPI-167 EPI-168 EPI-169 EPI-170 \
--mc-type "CGN" \
--num-procs 8 \
--output-prefix day135_AllpH_DMR250bp_cov5x \
--dmr-max-dist 250 \
--min-num-dms 3 

###################################
# Run DMRfind for Day 145 samples
###################################

methylpy DMRfind \
--allc-files EPI-175_S14_L003_5xmerg_allc.tsv \
EPI-176_S15_L003_5xmerg_allc.tsv \
EPI-181_S16_L003_5xmerg_allc.tsv \
EPI-182_S17_L003_5xmerg_allc.tsv \
EPI-184_S18_L003_5xmerg_allc.tsv \
EPI-185_S19_L003_5xmerg_allc.tsv \
EPI-187_S20_L003_5xmerg_allc.tsv \
EPI-188_S21_L003_5xmerg_allc.tsv \
EPI-193_S22_L003_5xmerg_allc.tsv \
EPI-194_S23_L003_5xmerg_allc.tsv \
EPI-199_S24_L003_5xmerg_allc.tsv \
EPI-200_S25_L003_5xmerg_allc.tsv \
EPI-205_S26_L004_5xmerg_allc.tsv \
EPI-206_S27_L004_5xmerg_allc.tsv \
EPI-208_S28_L004_5xmerg_allc.tsv \
EPI-209_S29_L004_5xmerg_allc.tsv \
EPI-214_S30_L004_5xmerg_allc.tsv \
EPI-215_S31_L004_5xmerg_allc.tsv \
EPI-220_S32_L004_5xmerg_allc.tsv \
EPI-221_S33_L004_5xmerg_allc.tsv \
EPI-226_S34_L004_5xmerg_allc.tsv \
EPI-227_S35_L004_5xmerg_allc.tsv \
EPI-229_S36_L004_5xmerg_allc.tsv \
EPI-230_S37_L004_5xmerg_allc.tsv \
--samples EPI-175 EPI-176 EPI-181 EPI-182 EPI-184 EPI-185 EPI-187 EPI-188 EPI-193 EPI-194 EPI-199 EPI-200 EPI-205 EPI-206 EPI-208 EPI-209 EPI-214 EPI-215 EPI-220 EPI-221 EPI-226 EPI-227 EPI-229 EPI-230 \
--mc-type "CGN" \
--num-procs 8 \
--output-prefix day145_AllpH_DMR250bp_cov5x \
--dmr-max-dist 250 \
--min-num-dms 3 
```
---

**All Output from DMRfind**

[https://gannet.fish.washington.edu/metacarcinus/Pgenerosa/analyses/20200327/](https://gannet.fish.washington.edu/metacarcinus/Pgenerosa/analyses/20200327/)

---

### Group statistics and genomic feature analysis of regions

[https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/DMR_stats_and_anno.Rmd](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/DMR_stats_and_anno.Rmd)

---

# DMG analysis
