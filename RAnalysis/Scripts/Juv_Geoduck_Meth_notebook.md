## Obtain Geoduck genome file
#### http://gannet.fish.washington.edu/metacarcinus/Pgenerosa/GENOMES/v070/Pgenerosa_v070.fa

## Run Bismark Genome preparation
#### http://gannet.fish.washington.edu/metacarcinus/Pgenerosa/GENOMES/v070/Pgenerosa_v070.fa
/Applications/Bismark-0.20.0/bismark_genome_preparation ~/MyProjects/Geoduck_Meth/Genome/


## Genome-Wide Extraction
/Applications/Bismark-0.20.0/bismark_methylation_extractor --gzip -p --ignore_r2 2 --scaffolds --bedGraph --zero_based --no_overlap --multicore 20 --buffer_size 20G --cytosine_report --report --genome_folder ~/MyProjects/Geoduck_Meth/BSGenome ~/MyProjects/Geoduck_Meth/DeDup/EPI-41_S38_L005_dedup.sorted.bam

### Cytosine Methylation Report
http://gannet.fish.washington.edu/metacarcinus/Pgenerosa/20181011/BisMeExt/



### Extract regions with 5x seq coverage

#### Number of CpG sequenced at least 5x
cat EPI-41_S38_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-41.5x.cov.CpG
cat EPI-42_S39_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-42.5x.cov.CpG
cat EPI-43_S40_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-43.5x.cov.CpG
cat EPI-44_S41_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-44.5x.cov.CpG
cat EPI-103_S27_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-103.5x.cov.CpG
cat EPI-104_S28_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-104.5x.cov.CpG
cat EPI-111_S29_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-111.5x.cov.CpG
cat EPI-113_S30_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-113.5x.cov.CpG
cat EPI-119_S31_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-119.5x.cov.CpG
cat EPI-120_S32_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-120.5x.cov.CpG
cat EPI-127_S33_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-127.5x.cov.CpG
cat EPI-128_S34_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-128.5x.cov.CpG
cat EPI-135_S35_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-135.5x.cov.CpG
cat EPI-136_S36_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-136.5x.cov.CpG
cat EPI-143_S37_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-143.5x.cov.CpG
cat EPI-145_S38_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-145.5x.cov.CpG
cat EPI-151_S2_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-151.5x.cov.CpG
cat EPI-152_S3_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-152.5x.cov.CpG
cat EPI-153_S4_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-153.5x.cov.CpG
cat EPI-154_S5_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-154.5x.cov.CpG
cat EPI-159_S6_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz| awk '($5+$6 > 5)' > EPI-159.5x.cov.CpG
cat EPI-160_S7_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-160.5x.cov.CpG
cat EPI-161_S8_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-161.5x.cov.CpG
cat EPI-162_S9_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-162.5x.cov.CpG
cat EPI-167_S10_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-167.5x.cov.CpG
cat EPI-168_S11_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-168.5x.cov.CpG
cat EPI-169_S12_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-169.5x.cov.CpG
cat EPI-170_S13_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-170.5x.cov.CpG
cat EPI-175_S14_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-175.5x.cov.CpG
cat EPI-176_S15_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-176.5x.cov.CpG
cat EPI-181_S16_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-181.5x.cov.CpG
cat EPI-182_S17_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-182.5x.cov.CpG
cat EPI-184_S18_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-184.5x.cov.CpG
cat EPI-185_S19_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-185.5x.cov.CpG
cat EPI-187_S20_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-187.5x.cov.CpG
cat EPI-188_S21_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-188.5x.cov.CpG
cat EPI-193_S22_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-193.5x.cov.CpG
cat EPI-194_S23_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-194.5x.cov.CpG
cat EPI-199_S24_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-199.5x.cov.CpG
cat EPI-200_S25_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-200.5x.cov.CpG
cat EPI-205_S26_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-205.5x.cov.CpG
cat EPI-206_S27_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-206.5x.cov.CpG
cat EPI-208_S28_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-208.5x.cov.CpG
cat EPI-209_S29_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-209.5x.cov.CpG
cat EPI-214_S30_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-214.5x.cov.CpG
cat EPI-215_S31_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-215.5x.cov.CpG
cat EPI-220_S32_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-220.5x.cov.CpG
cat EPI-221_S33_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-221.5x.cov.CpG
cat EPI-226_S34_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-226.5x.cov.CpG
cat EPI-227_S35_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-227.5x.cov.CpG
cat EPI-229_S36_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-229.5x.cov.CpG
cat EPI-230_S37_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '($5+$6 > 5)' > EPI-230.5x.cov.CpG


### Identify regions common to all at 5X cov
multiIntersectBed -i a.bed b.bed c.bed
###### The basic concept of this approach is that it compares the intervals found in N sorted (-k1,1 -k2,2n for BED) BED/GFF/VCF files and reports whether 0 to N of those files are present at each interval.
multiIntersectBed -i EPI-41.5x.cov.CpG EPI-42.5x.cov.CpG EPI-43.5x.cov.CpG EPI-44.5x.cov.CpG EPI-103.5x.cov.CpG \
EPI-104.5x.cov.CpG EPI-111.5x.cov.CpG EPI-113.5x.cov.CpG EPI-119.5x.cov.CpG EPI-120.5x.cov.CpG EPI-127.5x.cov.CpG \
EPI-128.5x.cov.CpG EPI-135.5x.cov.CpG EPI-136.5x.cov.CpG EPI-143.5x.cov.CpG EPI-145.5x.cov.CpG EPI-151.5x.cov.CpG \
EPI-152.5x.cov.CpG EPI-153.5x.cov.CpG EPI-154.5x.cov.CpG EPI-159.5x.cov.CpG EPI-160.5x.cov.CpG EPI-161.5x.cov.CpG \
EPI-162.5x.cov.CpG EPI-167.5x.cov.CpG EPI-168.5x.cov.CpG EPI-169.5x.cov.CpG EPI-170.5x.cov.CpG EPI-175.5x.cov.CpG \
EPI-176.5x.cov.CpG EPI-181.5x.cov.CpG EPI-182.5x.cov.CpG EPI-184.5x.cov.CpG EPI-185.5x.cov.CpG EPI-187.5x.cov.CpG \
EPI-188.5x.cov.CpG EPI-193.5x.cov.CpG EPI-194.5x.cov.CpG EPI-199.5x.cov.CpG EPI-200.5x.cov.CpG EPI-205.5x.cov.CpG \
EPI-206.5x.cov.CpG EPI-208.5x.cov.CpG EPI-209.5x.cov.CpG EPI-214.5x.cov.CpG EPI-215.5x.cov.CpG EPI-220.5x.cov.CpG \
EPI-221.5x.cov.CpG EPI-226.5x.cov.CpG EPI-227.5x.cov.CpG EPI-229.5x.cov.CpG EPI-230.5x.cov.CpG > all.bed

cat all.bed | awk '$4 ==52' > filtered.52.all.bed


### Extract regions common to all at 5X cov
Done in R












































cat test.bed | awk '$5 ~/1,2,3/' > test.bed.txt