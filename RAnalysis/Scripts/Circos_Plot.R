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
#circos.initializeWithIdeogram(duck.genome, species = NULL, sort.chr = TRUE)
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

Meth <- read.table("RAnalysis/Data/OSF/bedgraphs/5x/union_5x.bedgraph", sep = "\t", 
                   header = TRUE, stringsAsFactors = FALSE,nrows=863)
str(Meth)
tail(Meth)

Meth$X103 <-as.numeric(Meth$X103)
Meth$X104 <-as.numeric(Meth$X104)
Meth$X111 <-as.numeric(Meth$X111)
Meth$X113<-as.numeric(Meth$X113)
Meth$X119 <-as.numeric(Meth$X119)
Meth$X120 <-as.numeric(Meth$X120)
Meth$X127 <-as.numeric(Meth$X127)
Meth$X128 <-as.numeric(Meth$X128)
Meth$X135 <-as.numeric(Meth$X135)
Meth$X136 <-as.numeric(Meth$X136)
Meth$X143 <-as.numeric(Meth$X143)

Meth$X145 <-as.numeric(Meth$X145)
Meth$X151 <-as.numeric(Meth$X151)
Meth$X152 <-as.numeric(Meth$X152)
Meth$X153<-as.numeric(Meth$X153)
Meth$X154 <-as.numeric(Meth$X154)
Meth$X159 <-as.numeric(Meth$X159)
Meth$X160 <-as.numeric(Meth$X160)
Meth$X161 <-as.numeric(Meth$X161)
Meth$X162 <-as.numeric(Meth$X162)
Meth$X167 <-as.numeric(Meth$X167)
Meth$X168 <-as.numeric(Meth$X168)

Meth$X169 <-as.numeric(Meth$X169)
Meth$X170 <-as.numeric(Meth$X170)
Meth$X175 <-as.numeric(Meth$X175)
Meth$X176<-as.numeric(Meth$X176)
Meth$X181 <-as.numeric(Meth$X181)
Meth$X182 <-as.numeric(Meth$X182)
Meth$X184 <-as.numeric(Meth$X184)
Meth$X185 <-as.numeric(Meth$X185)
Meth$X187 <-as.numeric(Meth$X187)
Meth$X188 <-as.numeric(Meth$X188)
Meth$X193 <-as.numeric(Meth$X193)

Meth$X194 <-as.numeric(Meth$X194)
Meth$X199 <-as.numeric(Meth$X199)
Meth$X200 <-as.numeric(Meth$X200)
Meth$X205<-as.numeric(Meth$X205)
Meth$X206 <-as.numeric(Meth$X206)
Meth$X208 <-as.numeric(Meth$X208)
Meth$X209 <-as.numeric(Meth$X209)
Meth$X214 <-as.numeric(Meth$X214)
Meth$X215 <-as.numeric(Meth$X215)
Meth$X220 <-as.numeric(Meth$X220)
Meth$X221 <-as.numeric(Meth$X221)

Meth$X226 <-as.numeric(Meth$X226)
Meth$X227 <-as.numeric(Meth$X227)
Meth$X229 <-as.numeric(Meth$X229)
Meth$X230<-as.numeric(Meth$X230)
Meth$X41 <-as.numeric(Meth$X41)
Meth$X42 <-as.numeric(Meth$X42)
Meth$X43 <-as.numeric(Meth$X43)
Meth$X44 <-as.numeric(Meth$X44)

str(Meth)

smps <- colnames(Meth[4:55])
Meth[4:55] <- if_else(Meth[4:55] >= 0, 1, 0)
Meth[is.na(Meth)] <- 0
Meth$Mval <- rowMeans(Meth[4:55], na.rm=TRUE)

  #if_else(condition, true, false, 
head(Meth)
small.Meth <- Meth[,c(1,2,3,56)]
str(small.Meth)
max(small.Meth$Mval)

pdf("RAnalysis/Output/CircosPlot_Inset.pdf")
circos.clear()
circos.initializeWithIdeogram(small.duck.genome, species = NULL, sort.chr = TRUE)
circos.initializeWithIdeogram(small.gene, species = NULL, sort.chr = TRUE)
circos.genomicTrack(small.CpG, ylim=c(0,0.5), track.heigh=0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col="gray")
                    })
circos.genomicTrack(small.Meth, ylim=c(0,1), track.heigh=0.1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col="red")
                    })
dev.off()










