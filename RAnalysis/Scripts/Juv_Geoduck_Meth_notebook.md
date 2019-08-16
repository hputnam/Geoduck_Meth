## Obtain Geoduck genome assembly file
#### 
wget http://owl.fish.washington.edu/halfshell/genomic-databank/Pgenerosa_v074.fa

## Obtain Geoduck genome annotation file

wget http://owl.fish.washington.edu/halfshell/genomic-databank/Pgenerosa_v074.gene.gff

wget https://gannet.fish.washington.edu/Atumefaciens/20190701_pgen_maker_v074_annotation/blastp_annotation/blastp.outfmt6

## Run Bismark Genome preparation
#### Taking the trimmed files @ /gscratch/srlab/strigg/data/Pgenr/FASTQS and aligning to version 74 of the genome.


## Genome-Wide Extraction
https://sr320.github.io/Geoduck-alignments-v74/

### Cytosine Methylation Report
wget -r \
--no-directories --no-parent \
-A deduplicated.bismark.cov.gz https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/0807/

wget -r \
--no-directories --no-parent \
-A deduplicated.bismark.cov.gz https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/0807-003/

wget -r \
--no-directories --no-parent \
-A deduplicated.bismark.cov.gz https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/0807-004/

wget -r \
--no-directories --no-parent \
-A deduplicated.bismark.cov.gz  https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/0809-005/



### Extract regions with 5x seq coverage

#### Number of CpG sequenced at least 5x


for file in /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/*.cov
do
cat $file | awk '($5+$6 > 5)' > ${file}.5x.cov.CpG
done


#### Number of CpG sequenced at least 10x

for file in /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/*.cov
do
cat $file | awk '($5+$6 > 10)' > ${file}.10x.cov.CpG
done




### Identify regions common to all at 5X cov
multiIntersectBed -i a.bed b.bed c.bed
###### The basic concept of this approach is that it compares the intervals found in N sorted (-k1,1 -k2,2n for BED) BED/GFF/VCF files and reports whether 0 to N of those files are present at each interval.


multiIntersectBed -i /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/*5x.cov.CpG > all.5x.bed
cat all.5x.bed | awk '$4 ==52' > filtered.52.5x.bed

Number of CpG sequenced at least 5x
```wc -l filtered.52.5x.bed```
20172 filtered.52.5x.bed

### Identify regions common to all at 10X cov
multiIntersectBed -i /Users/hputnam/MyProjects/Geoduck_Meth/RAnalysis/Data/Extracted/*10x.cov.CpG > all.5x.bed > all.10x.bed

cat all.10x.bed | awk '$4 ==52' > filtered.52.all.10x.bed
Number of CpG sequenced at least 10x
```wc -l filtered.52.all.10x.bed```
10584 filtered.52.all.10x.bed


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