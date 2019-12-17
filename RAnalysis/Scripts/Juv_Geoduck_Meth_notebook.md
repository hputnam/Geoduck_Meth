## Geoduck data files 
https://osf.io/yem8n/

Home directory 
```Geoduck_Meth```

### Identify regions common to all 52 samples at 5X cov
multiIntersectBed -i a.bed b.bed c.bed
###### The basic concept of this approach is that it compares the intervals found in N sorted (-k1,1 -k2,2n for BED) BED/GFF/VCF files and reports whether 0 to N of those files are present at each interval.

```
multiIntersectBed -i /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/OSF/*_merge.cov.5x.cov.CpG > all.5x.bed

cat all.5x.bed | awk '$4 ==52' > filtered.52.5x.bed  
```

### Number of CpG sequenced at least 5x

wc -l filtered.52.5x.bed


* 21954 filtered.52.5x.bed


# Extracts CpG positions with 5x coverage from Gene regions for each sample

```
IntersectBed -a EPI-41_S38_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-41.5x.genes.bed

IntersectBed -a EPI-41.5x.genes.bed -b filtered.52.5x.bed > EPI-41.5x.genes_CpG.bed

IntersectBed -a EPI-42_S39_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-42.5x.genes.bed

IntersectBed -a EPI-42.5x.genes.bed -b filtered.52.5x.bed > EPI-42.5x.genes_CpG.bed

IntersectBed -a EPI-43_S40_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-43.5x.genes.bed

IntersectBed -a EPI-43.5x.genes.bed -b filtered.52.5x.bed > EPI-43.5x.genes_CpG.bed

IntersectBed -a EPI-44_S41_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-44.5x.genes.bed

IntersectBed -a EPI-44.5x.genes.bed -b filtered.52.5x.bed > EPI-44.5x.genes_CpG.bed

IntersectBed -a EPI-103_S27_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-103.5x.genes.bed

IntersectBed -a EPI-103.5x.genes.bed -b filtered.52.5x.bed > EPI-103.5x.genes_CpG.bed

IntersectBed -a EPI-104_S28_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-104.5x.genes.bed

IntersectBed -a EPI-104.5x.genes.bed -b filtered.52.5x.bed > EPI-104.5x.genes_CpG.bed

IntersectBed -a EPI-111_S29_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-111.5x.genes.bed

IntersectBed -a EPI-111.5x.genes.bed -b filtered.52.5x.bed > EPI-111.5x.genes_CpG.bed

IntersectBed -a EPI-113_S30_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-113.5x.genes.bed

IntersectBed -a EPI-113.5x.genes.bed -b filtered.52.5x.bed > EPI-113.5x.genes_CpG.bed

IntersectBed -a EPI-119_S31_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-119.5x.genes.bed

IntersectBed -a EPI-119.5x.genes.bed -b filtered.52.5x.bed > EPI-119.5x.genes_CpG.bed

IntersectBed -a EPI-120_S32_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-120.5x.genes.bed

IntersectBed -a EPI-120.5x.genes.bed -b filtered.52.5x.bed > EPI-120.5x.genes_CpG.bed

IntersectBed -a EPI-127_S33_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-127.5x.genes.bed

IntersectBed -a EPI-127.5x.genes.bed -b filtered.52.5x.bed > EPI-127.5x.genes_CpG.bed

IntersectBed -a EPI-128_S34_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-128.5x.genes.bed

IntersectBed -a EPI-128.5x.genes.bed -b filtered.52.5x.bed > EPI-128.5x.genes_CpG.bed

IntersectBed -a EPI-135_S35_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-135.5x.genes.bed

IntersectBed -a EPI-135.5x.genes.bed -b filtered.52.5x.bed > EPI-135.5x.genes_CpG.bed

IntersectBed -a EPI-136_S36_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-136.5x.genes.bed

IntersectBed -a EPI-136.5x.genes.bed -b filtered.52.5x.bed > EPI-136.5x.genes_CpG.bed

IntersectBed -a EPI-143_S37_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-143.5x.genes.bed

IntersectBed -a EPI-143.5x.genes.bed -b filtered.52.5x.bed > EPI-143.5x.genes_CpG.bed

IntersectBed -a EPI-145_S38_L005_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-145.5x.genes.bed

IntersectBed -a EPI-145.5x.genes.bed -b filtered.52.5x.bed > EPI-145.5x.genes_CpG.bed

IntersectBed -a EPI-151_S2_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-151.5x.genes.bed

IntersectBed -a EPI-151.5x.genes.bed -b filtered.52.5x.bed > EPI-151.5x.genes_CpG.bed

IntersectBed -a EPI-152_S3_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-152.5x.genes.bed

IntersectBed -a EPI-152.5x.genes.bed -b filtered.52.5x.bed > EPI-152.5x.genes_CpG.bed

IntersectBed -a EPI-153_S4_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-153.5x.genes.bed

IntersectBed -a EPI-153.5x.genes.bed -b filtered.52.5x.bed > EPI-153.5x.genes_CpG.bed

IntersectBed -a EPI-154_S5_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-154.5x.genes.bed

IntersectBed -a EPI-154.5x.genes.bed -b filtered.52.5x.bed > EPI-154.5x.genes_CpG.bed

IntersectBed -a EPI-159_S6_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-159.5x.genes.bed

IntersectBed -a EPI-159.5x.genes.bed -b filtered.52.5x.bed > EPI-159.5x.genes_CpG.bed

IntersectBed -a EPI-160_S7_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-160.5x.genes.bed

IntersectBed -a EPI-160.5x.genes.bed -b filtered.52.5x.bed > EPI-160.5x.genes_CpG.bed

IntersectBed -a EPI-161_S8_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-161.5x.genes.bed

IntersectBed -a EPI-161.5x.genes.bed -b filtered.52.5x.bed > EPI-161.5x.genes_CpG.bed

IntersectBed -a EPI-162_S9_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-162.5x.genes.bed

IntersectBed -a EPI-162.5x.genes.bed -b filtered.52.5x.bed > EPI-162.5x.genes_CpG.bed

IntersectBed -a EPI-167_S10_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-167.5x.genes.bed

IntersectBed -a EPI-167.5x.genes.bed -b filtered.52.5x.bed > EPI-167.5x.genes_CpG.bed

IntersectBed -a EPI-168_S11_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-168.5x.genes.bed

IntersectBed -a EPI-168.5x.genes.bed -b filtered.52.5x.bed > EPI-168.5x.genes_CpG.bed

IntersectBed -a EPI-169_S12_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-169.5x.genes.bed

IntersectBed -a EPI-169.5x.genes.bed -b filtered.52.5x.bed > EPI-169.5x.genes_CpG.bed

IntersectBed -a EPI-170_S13_L002_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-170.5x.genes.bed

IntersectBed -a EPI-170.5x.genes.bed -b filtered.52.5x.bed > EPI-170.5x.genes_CpG.bed

IntersectBed -a EPI-175_S14_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-175.5x.genes.bed

IntersectBed -a EPI-175.5x.genes.bed -b filtered.52.5x.bed > EPI-175.5x.genes_CpG.bed

IntersectBed -a EPI-176_S15_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-176.5x.genes.bed

IntersectBed -a EPI-176.5x.genes.bed -b filtered.52.5x.bed > EPI-176.5x.genes_CpG.bed

IntersectBed -a EPI-181_S16_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-181.5x.genes.bed

IntersectBed -a EPI-181.5x.genes.bed -b filtered.52.5x.bed > EPI-181.5x.genes_CpG.bed

IntersectBed -a EPI-182_S17_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-182.5x.genes.bed

IntersectBed -a EPI-182.5x.genes.bed -b filtered.52.5x.bed > EPI-182.5x.genes_CpG.bed

IntersectBed -a EPI-184_S18_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-184.5x.genes.bed

IntersectBed -a EPI-184.5x.genes.bed -b filtered.52.5x.bed > EPI-184.5x.genes_CpG.bed

IntersectBed -a EPI-185_S19_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-185.5x.genes.bed

IntersectBed -a EPI-185.5x.genes.bed -b filtered.52.5x.bed > EPI-185.5x.genes_CpG.bed

IntersectBed -a EPI-187_S20_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-187.5x.genes.bed

IntersectBed -a EPI-187.5x.genes.bed -b filtered.52.5x.bed > EPI-187.5x.genes_CpG.bed

IntersectBed -a EPI-188_S21_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-188.5x.genes.bed

IntersectBed -a EPI-188.5x.genes.bed -b filtered.52.5x.bed > EPI-188.5x.genes_CpG.bed

IntersectBed -a EPI-193_S22_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-193.5x.genes.bed

IntersectBed -a EPI-193.5x.genes.bed -b filtered.52.5x.bed > EPI-193.5x.genes_CpG.bed

IntersectBed -a EPI-194_S23_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-194.5x.genes.bed

IntersectBed -a EPI-194.5x.genes.bed -b filtered.52.5x.bed > EPI-194.5x.genes_CpG.bed

IntersectBed -a EPI-199_S24_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-199.5x.genes.bed

IntersectBed -a EPI-199.5x.genes.bed -b filtered.52.5x.bed > EPI-199.5x.genes_CpG.bed

IntersectBed -a EPI-200_S25_L003_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-200.5x.genes.bed

IntersectBed -a EPI-200.5x.genes.bed -b filtered.52.5x.bed > EPI-200.5x.genes_CpG.bed

IntersectBed -a EPI-205_S26_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-205.5x.genes.bed

IntersectBed -a EPI-205.5x.genes.bed -b filtered.52.5x.bed > EPI-205.5x.genes_CpG.bed

IntersectBed -a EPI-206_S27_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-206.5x.genes.bed

IntersectBed -a EPI-206.5x.genes.bed -b filtered.52.5x.bed > EPI-206.5x.genes_CpG.bed

IntersectBed -a EPI-208_S28_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-208.5x.genes.bed

IntersectBed -a EPI-208.5x.genes.bed -b filtered.52.5x.bed > EPI-208.5x.genes_CpG.bed

IntersectBed -a EPI-209_S29_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-209.5x.genes.bed

IntersectBed -a EPI-209.5x.genes.bed -b filtered.52.5x.bed > EPI-209.5x.genes_CpG.bed


IntersectBed -a EPI-214_S30_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-214.5x.genes.bed

IntersectBed -a EPI-214.5x.genes.bed -b filtered.52.5x.bed > EPI-214.5x.genes_CpG.bed

IntersectBed -a EPI-215_S31_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-215.5x.genes.bed

IntersectBed -a EPI-215.5x.genes.bed -b filtered.52.5x.bed > EPI-215.5x.genes_CpG.bed

IntersectBed -a EPI-220_S32_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-220.5x.genes.bed

IntersectBed -a EPI-220.5x.genes.bed -b filtered.52.5x.bed > EPI-220.5x.genes_CpG.bed

IntersectBed -a EPI-221_S33_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-221.5x.genes.bed

IntersectBed -a EPI-221.5x.genes.bed -b filtered.52.5x.bed > EPI-221.5x.genes_CpG.bed

IntersectBed -a EPI-226_S34_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-226.5x.genes.bed

IntersectBed -a EPI-226.5x.genes.bed -b filtered.52.5x.bed > EPI-226.5x.genes_CpG.bed

IntersectBed -a EPI-227_S35_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-227.5x.genes.bed

IntersectBed -a EPI-227.5x.genes.bed -b filtered.52.5x.bed > EPI-227.5x.genes_CpG.bed


IntersectBed -a EPI-229_S36_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-229.5x.genes.bed

IntersectBed -a EPI-229.5x.genes.bed -b filtered.52.5x.bed > EPI-229.5x.genes_CpG.bed

IntersectBed -a EPI-230_S37_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov_merge.cov.5x.cov.CpG -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3 > EPI-230.5x.genes.bed

IntersectBed -a EPI-230.5x.genes.bed -b filtered.52.5x.bed > EPI-230.5x.genes_CpG.bed


```

















### Identify regions common to all at 10X cov

multiIntersectBed -i /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/*_merge.cov.10x.cov.CpG > all.10x.bed

cat all.10x.bed | awk '$4 ==52' > filtered.52.10x.bed

Number of CpG sequenced at least 10x
```wc -l filtered.52.10x.bed```
9550 filtered.52.10x.bed


### Extract regions common to all at 10x cov
Done in R


## Calculate number of CpG in the genome

```grep -v '>' /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Pgenerosa_v074.fa | grep -o -i 'CG' | wc -l```

* number of CG's in the assembled genome
15,449,592 /2 =  7,724,796 

Number of CpG sequenced at least 10x
```wc -l filtered.52.10x.bed```
10584 filtered.52.all.10x.bed

### Compare regions common to all at 10x cov against gene regions

* Intersect the gene track with the 10x methylation track
intersectBed -a ~/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/filtered.52.all.10x.bed  -b /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Pgenerosa_v074.gene.gff > 10x.meth.genes.intersect

wc -l 10x.meth.genes.intersect
259 10x.meth.genes.intersect


### Compare regions common to all at 5x cov against gene regions

* Intersect the gene track with the 5x methylation track
intersectBed -a ~/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/filtered.52.5x.bed  -b ~/MyProjects/Geoduck_Meth/RAnalysis/Data/Genome/Pgenerosa_v074.gene.gff > 5x.meth.genes.intersect

wc -l 5x.meth.genes.intersect
675 5x.meth.genes.intersect







































cat test.bed | awk '$5 ~/1,2,3/' > test.bed.txt