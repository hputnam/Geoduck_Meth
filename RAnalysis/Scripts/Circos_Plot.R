library(tidyverse)
library(circlize) 

#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track
#The input data for circos.genomicInitialize() is also a data frame with at least three columns. The first column is genomic category (for cytoband data, it is chromosome name), and the next two columns are positions in each genomic category.

#Goal: Outer Coordinates, gene regions band
#Load Chromosomes
duck.genome <- data.frame(chr=c("Scaffold_01", "Scaffold_02", "Scaffold_03", "Scaffold_04", "Scaffold_05", "Scaffold_06",
                                "Scaffold_07", "Scaffold_08", "Scaffold_09", "Scaffold_10", "Scaffold_11", "Scaffold_12",
                                "Scaffold_13", "Scaffold_14", "Scaffold_15", "Scaffold_16", "Scaffold_17", "Scaffold_18"),
                          start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ), 
                          end=c(89643857,	69596280,	57743597, 65288255, 67248332, 61759565, 43120122, 61151155, 
                                38581958, 53961475, 51449921, 50438331, 44396874, 45393038, 47938513, 31980953, 34923512,	27737463),
                          lev1=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ), 
                          lev2=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ))
duck.genome$chr <- as.character(duck.genome$chr)
#duck.genome <- subset(duck.genome, chr=="Scaffold_01")
str(duck.genome)


gene <- read.table("RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.gene.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene <- gene[,c(1,4,5,6)]
gene$V4 <-as.numeric(gene$V4)
gene$V5 <- as.numeric(gene$V5)
gene$V6 <- as.numeric(0.5)
gene$V7 <- as.numeric(1)
str(gene)
#gene <- subset(gene, V1=="Scaffold_01")

pdf("RAnalysis/Output/CircosPlot.pdf")
circos.clear()
circos.initializeWithIdeogram(duck.genome, species = NULL, sort.chr = TRUE)
circos.initializeWithIdeogram(gene, species = NULL, sort.chr = TRUE)
dev.off()

#Goal: For a small region only
#Outer coordinate band, genes band, CpG band, and methylome band

CpG <- read.table("RAnalysis/Data/Genome/Panopea-generosa-v1.0.CpG.gff", sep = "\t", header = FALSE, stringsAsFactors = FALSE,nrows=10000)
CpG <- CpG[,c(1,4,5,6)]
CpG$V4 <-as.numeric(CpG$V4)
CpG$V5 <- as.numeric(CpG$V5)
CpG$V6 <- as.numeric(0.5)
CpG$V7 <- as.numeric(0.5)
str(CpG)
#CpG <- subset(CpG, V1=="Scaffold_01")

small.duck.genome <- data.frame(chr=c("Scaffold_01"),
                          start=c(1), 
                          end=c(100000),
                          lev1=c(1), 
                          lev2=c(1))

small.gene <- subset(gene, V1=="Scaffold_01")
small.gene <- small.gene[1:6,]
small.CpG <- CpG[1:2079,]

  
circos.clear()
circos.initializeWithIdeogram(small.duck.genome, species = NULL, sort.chr = TRUE)
circos.initializeWithIdeogram(small.gene, species = NULL, sort.chr = TRUE)
circos.genomicTrack(small.CpG, ylim=c(0,0.5), track.heigh=0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col="gray")
                    })











