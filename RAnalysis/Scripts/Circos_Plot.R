library(tidyverse)

library(circlize) 
#https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track
#The input data for circos.genomicInitialize() is also a data frame with at least three columns. The first column is genomic category (for cytoband data, it is chromosome name), and the next two columns are positions in each genomic category.


#Load Chromosomes
duck.genome <- data.frame(chr=c("Scaffold_01", "Scaffold_02", "Scaffold_03", "Scaffold_04", "Scaffold_05", "Scaffold_06",
                                "Scaffold_07", "Scaffold_08", "Scaffold_09", "Scaffold_10", "Scaffold_11", "Scaffold_12",
                                "Scaffold_13", "Scaffold_14", "Scaffold_15", "Scaffold_16", "Scaffold_17", "Scaffold_18"),
                          start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ), 
                          end=c(100,	69596280,	57743597, 65288255, 67248332, 61759565, 43120122, 61151155, 
                                38581958, 53961475, 51449921, 50438331, 44396874, 45393038, 47938513, 31980953, 34923512,	27737463),
                          lev1=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ), 
                          lev2=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ))
duck.genome$chr <- as.character(duck.genome$chr)
duck.genome <- subset(duck.genome, chr=="Scaffold_01")
str(duck.genome)

CDS <- read.table("RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.CDS.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
CDS <- CDS[,c(1,4,5,6)]
CDS$V4 <-as.numeric(CDS$V4)
CDS$V5 <- as.numeric(CDS$V5)
CDS$V6 <- as.numeric(0.5)
str(CDS)

CDS <- subset(CDS, V1=="Scaffold_01")

gene <- read.table("RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.gene.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene <- gene[,c(1,4,5,6)]
gene$V4 <-as.numeric(gene$V4)
gene$V5 <- as.numeric(gene$V5)
gene$V6 <- as.numeric(0.5)
str(gene)

CpG <- read.table("RAnalysis/Data/Genome/Panopea-generosa-v1.0.CpG.gff", sep = "\t", header = FALSE, stringsAsFactors = FALSE,nrows=1435111)
tail(CpG)
CpG <- CpG[,c(1,4,5,6)]
CpG$V4 <-as.numeric(CpG$V4)
CpG$V5 <- as.numeric(CpG$V5)
CpG$V6 <- as.numeric(0.5)
str(CpG)






circos.clear()

circos.initializeWithIdeogram(duck.genome, species = NULL, sort.chr = TRUE)
circos.genomicTrack(gene, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="blue")
                    })
circos.genomicTrack(CDS, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="cyan")
                    })
circos.genomicTrack(CpG, stack=FALSE, ylim=c(0,1), track.height = 0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.2, col="black")
                    })







library(BioCircos)
#https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html

duck.genome <- toGRanges(data.frame(chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18" ), 
                                    start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ), 
                                    end=c(89643857,	69596280,	57743597, 65288255, 67248332, 61759565, 43120122, 61151155, 38581958, 53961475, 
                                          51449921, 50438331, 44396874, 45393038, 47938513, 31980953, 34923512,	27737463)))

duck.genome = list("Scaffold_01" = 89643857)
# ,
#                    "Scaffold_02"= 69596280,
#                    "Scaffold_03"= 57743597, 
#                    "Scaffold_04"= 65288255, 
#                    "Scaffold_05"= 67248332, 
#                    "Scaffold_06"= 61759565,
#                    "Scaffold_07"= 43120122, 
#                    "Scaffold_08"= 61151155, 
#                    "Scaffold_09"= 38581958, 
#                    "Scaffold_10"= 53961475, 
#                    "Scaffold_11"= 51449921, 
#                    "Scaffold_12"= 50438331,
#                    "Scaffold_13"= 44396874, 
#                    "Scaffold_14"= 45393038, 
#                    "Scaffold_15"= 47938513, 
#                    "Scaffold_16"= 31980953, 
#                    "Scaffold_17"= 34923512, 
#                    "Scaffold_18"= 27737463)



# Arcs coordinates
snvChr = gene$V1
snvStart = gene$V4
snvEnd = gene$V5
snvValues <- gene$V6

# Create gene track
tracks <- BioCircosCNVTrack('gene_track', as.character(snvChr), snvStart, snvEnd, snvValues, 
                           color = "red", range = c(0,6))
BioCircos( genome = duck.genome, genomeFillColor = c("gray", 'lightgray'),
          genomeTicksScale = 10e+5)

# Create gene track
tracks = BioCircosCNVTrack('gene_track', as.character(snvChr), snvStart, snvEnd, snvValues, 
                           color = "red", range = c(0,6))
# Add background
#tracks = tracks + BioCircosBackgroundTrack("arcs_background", colors = "#2222EE")

BioCircos(tracks, genomeFillColor = "gray", genomeTicksDisplay = F, genomeLabelDy = 0)



