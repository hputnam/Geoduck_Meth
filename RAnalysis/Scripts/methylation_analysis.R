#Title: Geoduck juvenile methylation
#Project: NOAA FFAR
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20191120
#See Readme file for details

rm(list=ls()) #clears workspace 

#install and load required libraries
library(plotrix) 
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(topGO)
library(methylGSA)
library(genefilter)


# Set Working Directory:
setwd("~/MyProjects/Geoduck_Meth/RAnalysis/") #set working

#load sample information
sample.info <- read.csv("Data/Sample.Info.csv", header=T, sep=",", na.string="NA", stringsAsFactors = F) #read in info
samp <- sample.info$Sample.ID # set sample info
samp <- gsub("[_]", "-", samp) #remove extra characters

#load genes gff
Genes <- read.csv("Data/Genome/Panopea-generosa-vv0.74.a3.gene.gff3", header=F, sep="\t", na.string="NA", skip=3, stringsAsFactors = F) #read in data file
Genes <- Genes[,c(1,4,5,9)] #select desired columns only
colnames(Genes) <- c("scaffold", "start", "stop", "gene") #rename columns
Genes$gene <- gsub(";.*","",Genes$gene) #remove extra characters
Genes$gene <- gsub("ID=","",Genes$gene) #remove extra characters

#Load annotation file
Annot <- read.csv("Data/Genome/Panopea-generosa-vv0.74.a3.5d66e5e537f02-blast_functional.tab", header=TRUE, sep="\t", na.string="NA", stringsAsFactors = F)
colnames(Annot) <- c("gene","start", "end", "score", "uniprot.id", "Match.ID", "m_start", "m_end", "E.value", "identity", "al")
Annot$gene <- gsub("20544-","",Annot$gene) #remove extra characters
Annot$gene  <- gsub(".m01","",Annot$gene) #remove extra characters
Uniprots <- Annot[,c(1,5)] #pull out uniprot information
colnames(Uniprots) <- c("gene", "uniprot.id") #rename columns
write.csv(Uniprots, "Data/Genome/Uniprot.IDS.csv") #save to file
GOs <- read.csv("Data/Genome/Uniprot_GOs.tab", header=TRUE, sep="\t", na.string="NA", stringsAsFactors = F) #read in gene ontology information
colnames(GOs) <- c("uniprot.id", "protein.name.short", "status", "protein.name.long", "gene.names", "organism", "length", "GO.IDs") #rename columns of gene ontology information

###### Load in methylation counts #####
# epi.41 <- read.csv("Data/Extracted/EPI-41.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.41 <- epi.41[,-3]
# colnames(epi.41) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.41$position <- as.numeric(epi.41$position)
# epi.41$Sample.ID <- "EPI_41"
# epi.41 <- full_join(epi.41, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# # epi.41$Sample.ID <- "EPI_41"
# # epi.41$merger <- paste0(epi.41$scaffold,"_", epi.41$position)
# # epi.41 <- merge(Data, epi.41, by="merger")
# # epi.41 <- epi.41[,6:14]
# # colnames(epi.41) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")
# 
# epi.42 <- read.csv("Data/Extracted/EPI-42.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.42 <- epi.42[,-3]
# colnames(epi.42) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.42$position <- as.numeric(epi.42$position)
# epi.42$Sample.ID <- "EPI_42"
# epi.42 <- full_join(epi.42, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.43 <- read.csv("Data/Extracted/EPI-43.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.43 <- epi.43[,-3]
# colnames(epi.43) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.43$position <- as.numeric(epi.43$position)
# epi.43$Sample.ID <- "EPI_43"
# epi.43 <- full_join(epi.43, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.44 <- read.csv("Data/Extracted/EPI-44.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.44 <- epi.44[,-3]
# colnames(epi.44) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.44$position <- as.numeric(epi.44$position)
# epi.44$Sample.ID <- "EPI_44"
# epi.44 <- full_join(epi.44, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.103 <- read.csv("Data/Extracted/EPI-103.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.103 <- epi.103[,-3]
# colnames(epi.103) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.103$position <- as.numeric(epi.103$position)
# epi.103$Sample.ID <- "EPI_103"
# epi.103 <- full_join(epi.103, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.104 <- read.csv("Data/Extracted/EPI-104.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.104 <- epi.104[,-3]
# colnames(epi.104) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.104$position <- as.numeric(epi.104$position)
# epi.104$Sample.ID <- "EPI_104"
# epi.104 <- full_join(epi.104, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.111 <- read.csv("Data/Extracted/EPI-111.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.111 <- epi.111[,-3]
# colnames(epi.111) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.111$position <- as.numeric(epi.111$position)
# epi.111$Sample.ID <- "EPI_111"
# epi.111 <- full_join(epi.111, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.113 <- read.csv("Data/Extracted/EPI-113.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.113 <- epi.113[,-3]
# colnames(epi.113) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.113$position <- as.numeric(epi.113$position)
# epi.113$Sample.ID <- "EPI_113"
# epi.113 <- full_join(epi.113, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.119 <- read.csv("Data/Extracted/EPI-119.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.119 <- epi.119[,-3]
# colnames(epi.119) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.119$position <- as.numeric(epi.119$position)
# epi.119$Sample.ID <- "EPI_119"
# epi.119 <- full_join(epi.119, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.120 <- read.csv("Data/Extracted/EPI-120.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.120 <- epi.120[,-3]
# colnames(epi.120) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.120$position <- as.numeric(epi.120$position)
# epi.120$Sample.ID <- "EPI_120"
# epi.120 <- full_join(epi.120, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.127 <- read.csv("Data/Extracted/EPI-127.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.127 <- epi.127[,-3]
# colnames(epi.127) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.127$position <- as.numeric(epi.127$position)
# epi.127$Sample.ID <- "EPI_127"
# epi.127 <- full_join(epi.127, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.128 <- read.csv("Data/Extracted/EPI-128.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.128 <- epi.128[,-3]
# colnames(epi.128) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.128$position <- as.numeric(epi.128$position)
# epi.128$Sample.ID <- "EPI_128"
# epi.128 <- full_join(epi.128, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.135 <- read.csv("Data/Extracted/EPI-135.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.135 <- epi.135[,-3]
# colnames(epi.135) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.135$position <- as.numeric(epi.135$position)
# epi.135$Sample.ID <- "EPI_135"
# epi.135 <- full_join(epi.135, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.136 <- read.csv("Data/Extracted/EPI-136.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.136 <- epi.136[,-3]
# colnames(epi.136) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.136$position <- as.numeric(epi.136$position)
# epi.136$Sample.ID <- "EPI_136"
# epi.136 <- full_join(epi.136, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.143 <- read.csv("Data/Extracted/EPI-143.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.143 <- epi.143[,-3]
# colnames(epi.143) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.143$position <- as.numeric(epi.143$position)
# epi.143$Sample.ID <- "EPI_143"
# epi.143 <- full_join(epi.143, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.145 <- read.csv("Data/Extracted/EPI-145.5x.genes_CpG.bedgraph", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.145 <- epi.145[,-3]
# colnames(epi.145) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.145$position <- as.numeric(epi.145$position)
# epi.145$Sample.ID <- "EPI_145"
# epi.145 <- full_join(epi.145, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.151 <- read.csv("Data/Extracted/EPI-151.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.151 <- epi.151[,-3]
# colnames(epi.151) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.151$position <- as.numeric(epi.151$position)
# epi.151$Sample.ID <- "EPI_151"
# epi.151 <- full_join(epi.151, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.152 <- read.csv("Data/Extracted/EPI-152.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.152 <- epi.152[,-3]
# colnames(epi.152) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.152$position <- as.numeric(epi.152$position)
# epi.152$Sample.ID <- "EPI_152"
# epi.152 <- full_join(epi.152, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.153 <- read.csv("Data/Extracted/EPI-153.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.153 <- epi.153[,-3]
# colnames(epi.153) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.153$position <- as.numeric(epi.153$position)
# epi.153$Sample.ID <- "EPI_153"
# epi.153 <- full_join(epi.153, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.154 <- read.csv("Data/Extracted/EPI-154.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.154 <- epi.154[,-3]
# colnames(epi.154) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.154$position <- as.numeric(epi.154$position)
# epi.154$Sample.ID <- "EPI_154"
# epi.154 <- full_join(epi.154, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.159 <- read.csv("Data/Extracted/EPI-159.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.159 <- epi.159[,-3]
# colnames(epi.159) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.159$position <- as.numeric(epi.159$position)
# epi.159$Sample.ID <- "EPI_159"
# epi.159 <- full_join(epi.159, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.160 <- read.csv("Data/Extracted/EPI-160.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.160 <- epi.160[,-3]
# colnames(epi.160) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.160$position <- as.numeric(epi.160$position)
# epi.160$Sample.ID <- "EPI_160"
# epi.160 <- full_join(epi.160, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.161 <- read.csv("Data/Extracted/EPI-161.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.161 <- epi.161[,-3]
# colnames(epi.161) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.161$position <- as.numeric(epi.161$position)
# epi.161$Sample.ID <- "EPI_161"
# epi.161 <- full_join(epi.161, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.162 <- read.csv("Data/Extracted/EPI-162.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.162 <- epi.162[,-3]
# colnames(epi.162) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.162$position <- as.numeric(epi.162$position)
# epi.162$Sample.ID <- "EPI_162"
# epi.162 <- full_join(epi.162, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.167 <- read.csv("Data/Extracted/EPI-167.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.167 <- epi.167[,-3]
# colnames(epi.167) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.167$position <- as.numeric(epi.167$position)
# epi.167$Sample.ID <- "EPI_167"
# epi.167 <- full_join(epi.167, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.168 <- read.csv("Data/Extracted/EPI-168.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.168 <- epi.168[,-3]
# colnames(epi.168) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.168$position <- as.numeric(epi.168$position)
# epi.168$Sample.ID <- "EPI_168"
# epi.168 <- full_join(epi.168, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.169 <- read.csv("Data/Extracted/EPI-169.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.169 <- epi.169[,-3]
# colnames(epi.169) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.169$position <- as.numeric(epi.169$position)
# epi.169$Sample.ID <- "EPI_169"
# epi.169 <- full_join(epi.169, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.170 <- read.csv("Data/Extracted/EPI-170.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.170 <- epi.170[,-3]
# colnames(epi.170) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.170$position <- as.numeric(epi.170$position)
# epi.170$Sample.ID <- "EPI_170"
# epi.170 <- full_join(epi.170, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.175 <- read.csv("Data/Extracted/EPI-175.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.175 <- epi.175[,-3]
# colnames(epi.175) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.175$position <- as.numeric(epi.175$position)
# epi.175$Sample.ID <- "EPI_175"
# epi.175 <- full_join(epi.175, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.176 <- read.csv("Data/Extracted/EPI-176.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.176 <- epi.176[,-3]
# colnames(epi.176) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.176$Sample.ID <- "EPI_176"
# epi.176 <- full_join(epi.176, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.181 <- read.csv("Data/Extracted/EPI-181.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.181 <- epi.181[,-3]
# colnames(epi.181) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.181$position <- as.numeric(epi.181$position)
# epi.181$Sample.ID <- "EPI_181"
# epi.181 <- full_join(epi.181, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.182 <- read.csv("Data/Extracted/EPI-182.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.182 <- epi.182[,-3]
# colnames(epi.182) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.182$position <- as.numeric(epi.182$position)
# epi.182$Sample.ID <- "EPI_182"
# epi.182 <- full_join(epi.182, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.184 <- read.csv("Data/Extracted/EPI-184.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.184 <- epi.184[,-3]
# colnames(epi.184) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.184$position <- as.numeric(epi.184$position)
# epi.184$Sample.ID <- "EPI_184"
# epi.184 <- full_join(epi.184, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.185 <- read.csv("Data/Extracted/EPI-185.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.185 <- epi.185[,-3]
# colnames(epi.185) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.185$position <- as.numeric(epi.185$position)
# epi.185$Sample.ID <- "EPI_185"
# epi.185 <- full_join(epi.185, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.187 <- read.csv("Data/Extracted/EPI-187.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.187 <- epi.187[,-3]
# colnames(epi.187) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.187$position <- as.numeric(epi.187$position)
# epi.187$Sample.ID <- "EPI_187"
# epi.187 <- full_join(epi.187, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.188 <- read.csv("Data/Extracted/EPI-188.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.188 <- epi.188[,-3]
# colnames(epi.188) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.188$position <- as.numeric(epi.188$position)
# epi.188$Sample.ID <- "EPI_188"
# epi.188 <- full_join(epi.188, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.193 <- read.csv("Data/Extracted/EPI-193.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.193 <- epi.193[,-3]
# colnames(epi.193) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.193$position <- as.numeric(epi.193$position)
# epi.193$Sample.ID <- "EPI_193"
# epi.193 <- full_join(epi.193, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.194 <- read.csv("Data/Extracted/EPI-194.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.194 <- epi.194[,-3]
# colnames(epi.194) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.194$position <- as.numeric(epi.194$position)
# epi.194$Sample.ID <- "EPI_194"
# epi.194 <- full_join(epi.194, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.199 <- read.csv("Data/Extracted/EPI-199.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.199 <- epi.199[,-3]
# colnames(epi.199) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.199$position <- as.numeric(epi.199$position)
# epi.199$Sample.ID <- "EPI_199"
# epi.199 <- full_join(epi.199, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.200 <- read.csv("Data/Extracted/EPI-200.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.200 <- epi.200[,-3]
# colnames(epi.200) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.200$position <- as.numeric(epi.200$position)
# epi.200$Sample.ID <- "EPI_200"
# epi.200 <- full_join(epi.200, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.205 <- read.csv("Data/Extracted/EPI-205.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.205 <- epi.205[,-3]
# colnames(epi.205) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.205$position <- as.numeric(epi.205$position)
# epi.205$Sample.ID <- "EPI_205"
# epi.205 <- full_join(epi.205, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.206 <- read.csv("Data/Extracted/EPI-206.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.206 <- epi.206[,-3]
# colnames(epi.206) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.206$position <- as.numeric(epi.206$position)
# epi.206$Sample.ID <- "EPI_206"
# epi.206 <- full_join(epi.206, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.208 <- read.csv("Data/Extracted/EPI-208.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.208 <- epi.208[,-3]
# colnames(epi.208) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.208$position <- as.numeric(epi.208$position)
# epi.208$Sample.ID <- "EPI_208"
# epi.208 <- full_join(epi.208, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.209 <- read.csv("Data/Extracted/EPI-209.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.209 <- epi.209[,-3]
# colnames(epi.209) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.209$position <- as.numeric(epi.209$position)
# epi.209$Sample.ID <- "EPI_209"
# epi.209 <- full_join(epi.209, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.214 <- read.csv("Data/Extracted/EPI-214.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.214 <- epi.214[,-3]
# colnames(epi.214) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.214$position <- as.numeric(epi.214$position)
# epi.214$Sample.ID <- "EPI_214"
# epi.214 <- full_join(epi.214, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.215 <- read.csv("Data/Extracted/EPI-215.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.215 <- epi.215[,-3]
# colnames(epi.215) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.215$position <- as.numeric(epi.215$position)
# epi.215$Sample.ID <- "EPI_215"
# epi.215 <- full_join(epi.215, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.220 <- read.csv("Data/Extracted/EPI-220.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.220 <- epi.220[,-3]
# colnames(epi.220) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.220$position <- as.numeric(epi.220$position)
# epi.220$Sample.ID <- "EPI_220"
# epi.220 <- full_join(epi.220, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.221 <- read.csv("Data/Extracted/EPI-221.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.221 <- epi.221[,-3]
# colnames(epi.221) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.221$position <- as.numeric(epi.221$position)
# epi.221$Sample.ID <- "EPI_221"
# epi.221 <- full_join(epi.221, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.226 <- read.csv("Data/Extracted/EPI-226.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.226 <- epi.226[,-3]
# colnames(epi.226) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.226$position <- as.numeric(epi.226$position)
# epi.226$Sample.ID <- "EPI_226"
# epi.226 <- full_join(epi.226, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.227 <- read.csv("Data/Extracted/EPI-227.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.227 <- epi.227[,-3]
# colnames(epi.227) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.227$position <- as.numeric(epi.227$position)
# epi.227$Sample.ID <- "EPI_227"
# epi.227 <- full_join(epi.227, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.229 <- read.csv("Data/Extracted/EPI-229.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.229 <- epi.229[,-3]
# colnames(epi.229) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.229$position <- as.numeric(epi.229$position)
# epi.229$Sample.ID <- "EPI_229"
# epi.229 <- full_join(epi.229, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
# 
# epi.230 <- read.csv("Data/Extracted/EPI-230.5x.genes_CpG.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
# epi.230 <- epi.230[,-3]
# colnames(epi.230) <- c("scaffold", "position","per.meth","meth","unmeth")
# epi.230$position <- as.numeric(epi.230$position)
# epi.230$Sample.ID <- "EPI_230"
# epi.230 <- full_join(epi.230, Genes, by="scaffold") %>% filter(position >= start & position <= stop)

##### Plotting mean percent methylation #####

##### Plot Ambient Treatments through time #####
time.amb.meth <- rbind(epi.41,epi.42,epi.43,epi.44,
                       epi.119,epi.120,epi.135,epi.136,
                       epi.151,epi.152,epi.153,epi.154,
                       epi.181,epi.182,epi.184,epi.185)
time.amb.meth <- merge(time.amb.meth,sample.info, by="Sample.ID")
time.amb.meth <- time.amb.meth[,c(1,2,3,5,6,9,12,14)]
colnames(time.amb.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "Treatment", "TimePoint")

#
# colnames(time.amb.meth) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")
# timepoint.amb.meth <- merge(time.amb.meth,sample.info, by="Sample.ID")
# time.amb.meth <- timepoint.amb.meth[,c(1,2,3,5,6,9,14)]
# colnames(time.amb.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment")
#
# sams <- aggregate(per.meth~Sample.ID, data=timepoint.amb.meth, FUN=mean)
# sams <- merge(sams, sample.info, by="Sample.ID")
# time.amb.means <- aggregate(per.meth~TimePoint, data=sams, FUN=mean)
# time.amb.se <- aggregate(per.meth~TimePoint, data=sams, FUN=std.error)
# time.amb.meth.per <- cbind(time.amb.means, time.amb.se$per.meth)
# colnames(time.amb.meth.per) <- c("TimePoint", "mean","se")
#
# Fig.time.amb.per.meth <- ggplot(time.amb.meth.per, aes(x=TimePoint, y=mean)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
#   geom_point(pch=16, size = 3, position = position_dodge(width = 0.05)) +
#   #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
#   #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
#   xlab("Time") +
#   ylab("Methylation Level") +
#   #ylim(19.5,20.5)  +
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(), #Set the plot background
#         legend.position='none',
#         axis.title.y = element_text(size = rel(1.8), angle = 90),
#         axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
#   ggtitle("Ambient Through Time") +
#   theme(plot.title = element_text(face = 'bold',
#                                   size = 12,
#                                   hjust = 0))
# Fig.time.amb.per.meth
#
# sams
# time.mod <- lm(per.meth ~TimePoint, data=sams)
# summary(aov(time.mod))

# Treatments after 10 days of exposure  #####
D10.trt.meth <- rbind(epi.103,epi.104,epi.111,epi.113,
                      epi.119,epi.120,epi.127,epi.128,
                      epi.135,epi.136,epi.143,epi.145)
D10.trt.meth <- merge(D10.trt.meth,sample.info, by="Sample.ID")
D10.trt.meth <- D10.trt.meth[,c(1,2,3,5,6,9,12)]
colnames(D10.trt.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment1")


#####
# D10.trt.meth$per.meth <- (D10.trt.meth$meth/(D10.trt.meth$meth+D10.trt.meth$unmeth))*100
# sams <- aggregate(per.meth~Sample.ID, data=D10.trt.meth, FUN=mean)
# sams <- merge(sams, sample.info, by="Sample.ID")
# D10.trt.means <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=mean)
# D10.trt.se <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=std.error)
# D10.meth.per <- cbind(D10.trt.means, D10.trt.se$per.meth)
# colnames(D10.meth.per) <- c("Initial.Treatment", "mean","se")
#
# #Exposure1 Plotting
# Fig.D10.per.meth <- ggplot(D10.meth.per, aes(x=Initial.Treatment, y=mean, group=Initial.Treatment)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
#   #geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +
#   geom_point(pch=c(16,15,17), size = 3, position = position_dodge(width = 0.05)) +
#   #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
#   #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
#   xlab("Initial Treatment") +
#   ylab("Methylation Level") +
#   #ylim(19.5,20.5)  +
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(), #Set the plot background
#         legend.position='none',
#         axis.title.y = element_text(size = rel(1.8), angle = 90),
#         axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
#   ggtitle("Exposure 1 Day 10") +
#   theme(plot.title = element_text(face = 'bold',
#                                   size = 12,
#                                   hjust = 0))
# Fig.D10.per.meth
#
# sams
# D10.mod <- lm(per.meth ~Initial.Treatment, data=sams)
# summary(aov(D10.mod))

#  D135 treatments after common garden #####
D135.trt.meth <- rbind(epi.151,epi.152,epi.153,epi.154,
                       epi.159,epi.160,epi.161,epi.162,
                       epi.167,epi.168,epi.169,epi.170)
D135.trt.meth <- merge(D135.trt.meth, sample.info, by="Sample.ID")
D135.trt.meth <- D135.trt.meth[,c(1,2,3,5,6,9,12)]
colnames(D135.trt.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment")


#####
# D135.trt.meth$per.meth <- (D135.trt.meth$meth/(D135.trt.meth$meth+D135.trt.meth$unmeth))*100
# sams <- aggregate(per.meth~Sample.ID, data=D135.trt.meth, FUN=mean)
# sams <- merge(sams, sample.info, by="Sample.ID")
# D135.trt.means <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=mean)
# D135.trt.se <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=std.error)
# D135.meth.per <- cbind(D135.trt.means, D135.trt.se$per.meth)
# colnames(D135.meth.per) <- c("Initial.Treatment", "mean","se")
#
# #Common Garden Plotting
# Fig.D135.per.meth <- ggplot(D135.meth.per, aes(x=Initial.Treatment, y=mean, group=Initial.Treatment)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
#   #geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +
#   geom_point(pch=c(16,15,17), size = 3, position = position_dodge(width = 0.05)) +
#   #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
#   #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
#   xlab("Initial Treatment") +
#   ylab("Methylation Level") +
#   #ylim(19.5,20.5)  +
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(), #Set the plot background
#         legend.position='none',
#         axis.title.y = element_text(size = rel(1.8), angle = 90),
#         axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
#   ggtitle("Post Common Garden") +
#   theme(plot.title = element_text(face = 'bold',
#                                   size = 12,
#                                   hjust = 0))
# Fig.D135.per.meth
#
# sams
# D135.mod <- lm(per.meth ~Initial.Treatment, data=sams)
# summary(aov(D135.mod))

# Secondary exposure D145 #####
D145.trt.meth <- rbind(epi.175,epi.176,epi.181,epi.182,
                  epi.184,epi.185,epi.187,epi.188,
                  epi.193,epi.194,epi.199,epi.200,
                  epi.205,epi.206,epi.208,epi.209,
                  epi.214,epi.215,epi.220,epi.221,
                  epi.226,epi.227,epi.229,epi.230)
D145.trt.meth <- merge(D145.trt.meth,sample.info, by="Sample.ID")
D145.trt.meth <- D145.trt.meth[,c(1,2,3,5,6,9,12,13)]
colnames(D145.trt.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment1", "treatment2")

#####
# sec.meth$per.meth <- (sec.meth$meth/(sec.meth$meth+sec.meth$unmeth))*100
# sams <- aggregate(per.meth~Sample.ID, data=sec.meth, FUN=mean)
# sams <- merge(sams, sample.info, by="Sample.ID")
# D145.trt.means <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=mean)
# D145.trt.se <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=std.error)
# D145.meth.per <- cbind(D145.trt.means, D145.trt.se$per.meth)
# colnames(D145.meth.per) <- c("Initial.Treatment", "Secondary.Treatment", "mean","se")
# D145.meth.per$Initial.Treatment <- factor(D145.meth.per$Initial.Treatment,levels=c("Ambient","Super.Low", "Low"))
#
# #Exposure2 Plotting
# Fig.D145.per.meth <- ggplot(D145.meth.per, aes(x=Secondary.Treatment, y=mean, group=Initial.Treatment)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
#   geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +
#   geom_point(aes(shape=Initial.Treatment), size = 3, position = position_dodge(width = 0.05)) +
#   #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
#   #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
#   xlab("Secondary Treatment") +
#   ylab("Methylation Level") +
#   #ylim(19.5,20.5) +
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(), #Set the plot background
#         legend.position='none',
#         axis.title.y = element_text(size = rel(1.8), angle = 90),
#         axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
#   ggtitle("Secondary Exposure") +
#   theme(plot.title = element_text(face = 'bold',
#                                   size = 12,
#                                   hjust = 0))
# Fig.D145.per.meth
#
# sams
# D145.mod <- lm(per.meth ~Initial.Treatment*Secondary.Treatment, data=sams)
# summary(aov(D145.mod))

# save plots

# Meth.Level <- arrangeGrob(Fig.D10.per.meth, Fig.D135.per.meth, Fig.D145.per.meth, ncol=2)
# ggsave(file="Output/Methylation_Level.pdf", Meth.Level, width =6, height = 6, units = c("in"))
#
# D10 <- arrangeGrob(Fig.D10.per.meth, ncol=1)
# ggsave(file="Output/D10.trt.pdf", D10, width =6, height = 6, units = c("in"))
#
# D135 <- arrangeGrob(Fig.D135.per.meth, ncol=1)
# ggsave(file="Output/D135.trt.pdf", D135, width =6, height = 6, units = c("in"))
#
# D145 <- arrangeGrob(Fig.D145.per.meth, ncol=1)
# ggsave(file="Output/D145.trt.pdf", D145, width =6, height = 6, units = c("in"))

# Save data tables #####
write.table(time.amb.meth, 'Output/Time.Amb.tsv', sep='\t', row.names=FALSE)
write.table(D10.trt.meth, 'Output/D10.tsv', sep='\t', row.names=FALSE)
write.table(D135.trt.meth, 'Output/D135.tsv', sep='\t', row.names=FALSE)
write.table(D145.trt.meth, 'Output/D145.tsv', sep='\t', row.names=FALSE)


##### READ IN DATAFRAMES #####
time.amb.meth <- read.csv('Output/Time.Amb.tsv', sep='\t') #read in Day10 data for 5x coverage for CpGs in all samples
D10.trt.meth <- read.csv('Output/D10.tsv', sep='\t') #read in Day10 data for 5x coverage for CpGs in all samples
D135.trt.meth <- read.csv( 'Output/D135.tsv', sep='\t') #read in Day135 data for 5x coverage for CpGs in all samples
D145.trt.meth <- read.csv('Output/D145.tsv', sep='\t') #read in Day145 data for 5x coverage for CpGs in all samples

##### Statistical comparisons of methylation #####
# Ambient through Time #####
time.amb.meth.annot <- merge(time.amb.meth, Annot, by="gene", all.x=TRUE) #merge the data with the annotation
time.amb.meth.annot <- merge(time.amb.meth.annot, GOs, by="uniprot.id", all.x=TRUE) #merge the annotated data with GO terms from uniprot
time.amb.meth.annot[is.na(time.amb.meth.annot)] <- "unknown" #fill unknown for genes lacking annotation indicated by NA
time.amb.meth.annot[time.amb.meth.annot==""] <- "unknown" #fill unknown for genes lacking annotation indicated by empty space
time.amb.meth.annot.MWU <- as.data.frame(cbind(as.character(time.amb.meth.annot$gene), time.amb.meth.annot$GO.IDs)) #select the gene name and GO terms into new dataframe
colnames(time.amb.meth.annot.MWU) <- c("gene", "GO.IDs") #rename new data frame columns
nonredun.time.amb.meth.annot <- time.amb.meth.annot.MWU[!duplicated(time.amb.meth.annot.MWU$gene),]
colnames(nonredun.time.amb.meth.annot) <- c("gene", "GO.IDs") #rename columns of nonredundant data

sub_meth_table  <- time.amb.meth
sub_meth_table$group <- paste0(time.amb.meth$Sample.ID, sub_meth_table$gene)

#filter for genes with >5 methylated positions
min.filt <- count(sub_meth_table, vars = c( group))
newdata <- min.filt[ which(min.filt$n > 4), ]
sub_meth_table <- sub_meth_table[sub_meth_table$group %in% newdata$vars,]

# create data frame to stored results
results <- data.frame()
gs <- unique(sub_meth_table$gene)

#first subset the unique dataframes and second run the GLMs
for(i in 1:length(sub_meth_table$gene)){
  
  #subset the dataframe gene by gene
  sub_meth_table1 <- subset(sub_meth_table, gene ==gs[i])
  
  # fit glm position model
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ TimePoint, 
             data=sub_meth_table1, family=binomial)
  a <- anova(fit, test="Chisq")
  
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,6],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   #pval.position = a$`Pr(>Chi)`[3], #uncomment if you want to include position of CpG within a gene
                   #pval.treatment_x_position = a$`Pr(>Chi)`[4], #uncomment if you want to include position of CpG within a gene interaction with treatment
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.TimePoint <- p.adjust(results$pval.treatment, method='BH')
#results$adj.pval.position <- p.adjust(results$pval.position, method='BH') #uncomment if you want to include position of CpG within a gene
#results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH') #uncomment if you want to include position of CpG within a gene interaction with treatment

write.table(results, 'Output/Time_Amb_HP_GLM.tsv', sep='\t', row.names=FALSE)

Time.Amb.sig <-results
Time.Amb.sig <- Time.Amb.sig[,c(1,3)]
Time.Amb.sig <- Time.Amb.sig[order(Time.Amb.sig$adj.pval.TimePoint),]
Time.Amb.sig <- Time.Amb.sig[which(Time.Amb.sig$adj.pval.TimePoint<0.05), ]

write.table(Time.Amb.sig, 'Output/Time.Amb_HP_GLM_SortedResults.tsv', sep='\t', row.names=FALSE)

sub_meth_table$per.meth <- (sub_meth_table$meth/(sub_meth_table$meth+sub_meth_table$unmeth))*100 #calculate percent methylation
sub_meth_table.Time.Amb <- sub_meth_table #store in day-specific dataframe

# Annotation of Time DMG 
Time.Amb.annot <- merge(Time.Amb.sig , Annot, by="gene", all.x=TRUE)
Time.Amb.annot <- merge(Time.Amb.annot, GOs, by="uniprot.id", all.x=TRUE)
Time.Amb.annot <- Time.Amb.annot[order(Time.Amb.annot$adj.pval.TimePoint),]
Time.Amb.annot <- Time.Amb.annot[!duplicated(Time.Amb.annot$gene),]
write.table(Time.Amb.annot, 'Output/Table_S1_Time_Amb_sig_annot.tsv', sep='\t', row.names=FALSE)

#####GO Enrichment of Time.Amb DMGs
#Day 10 GO Enrichment for Time.Amb DMGs with GO MWU approach https://github.com/z0on/GO_MWU #####
#identify MWU data
Time.Amb.MWU <- nonredun.time.amb.meth.annot
tail(Time.Amb.MWU$gene)
#write MWU data to file
write.table(Time.Amb.MWU, 'GO_MWU/Time.Amb.GO.MWU', sep='\t', row.names=FALSE)
#merge all genes with DMGs
Time.Amb.MWU <- merge(Time.Amb.MWU, Time.Amb.annot, by="gene", all=TRUE)
tail(Time.Amb.MWU$gene)
Time.Amb.MWU <- Time.Amb.MWU[,c(1,4)]
colnames(Time.Amb.MWU) <- c("gene", "DMG")
Time.Amb.MWU$DMG[Time.Amb.MWU$DMG < 0.05] <- 1
which(Time.Amb.MWU$DMG ==1)
Time.Amb.MWU[is.na(Time.Amb.MWU)] <- 0
which(Time.Amb.MWU$DMG ==1)

write.table(Time.Amb.MWU, 'GO_MWU/Time.Amb_interest_gsig.MWU', sep=',', row.names=FALSE)


# Day10 #####
D10.trt.meth.annot <- merge(D10.trt.meth, Annot, by="gene", all.x=TRUE) #merge the data with the annotation
D10.trt.meth.annot <- merge(D10.trt.meth.annot, GOs, by="uniprot.id", all.x=TRUE) #merge the annotated data with GO terms from uniprot
D10.trt.meth.annot[is.na(D10.trt.meth.annot)] <- "unknown" #fill unknown for genes lacking annotation indicated by NA
D10.trt.meth.annot[D10.trt.meth.annot==""] <- "unknown" #fill unknown for genes lacking annotation indicated by empty space
D10.trt.meth.annot.MWU <- as.data.frame(cbind(as.character(D10.trt.meth.annot$gene), D10.trt.meth.annot$GO.IDs)) #select the gene name and GO terms into new dataframe
colnames(D10.trt.meth.annot.MWU) <- c("gene", "GO.IDs") #rename new data frame columns
write.table(D10.trt.meth.annot.MWU, 'GO_MWU/D10.trt.meth.annot.MWU', sep='\t', row.names=FALSE) #save new dataframe for use in GO enrichement analysis
nonredun.D10.trt.meth.annot <- D10.trt.meth.annot.MWU[!duplicated(D10.trt.meth.annot.MWU$gene),]
colnames(nonredun.D10.trt.meth.annot) <- c("gene", "GO.IDs") #rename columns of nonredundant data

meth_table  <- D10.trt.meth 
sub_meth_table <- meth_table[meth_table$Initial.Treatment == 'Ambient' | meth_table$Initial.Treatment == 'Low' | meth_table$Initial.Treatment == 'Super.Low', ]
sub_meth_table$group <- paste0(sub_meth_table$Sample.ID, sub_meth_table$gene)
  
#filter for genes with >5 methylated positions
min.filt <- count(sub_meth_table, vars = c( group))
newdata <- min.filt[ which(min.filt$n > 4), ]
sub_meth_table <- sub_meth_table[sub_meth_table$group %in% newdata$vars,]

# create data frame to stored results
results <- data.frame()
gs <- unique(sub_meth_table$gene)

#first subset the unique dataframes and second run the GLMs
for(i in 1:length(sub_meth_table$gene)){
  
  #subset the dataframe gene by gene
  sub_meth_table1 <- subset(sub_meth_table, gene ==gs[i])
  
  # fit glm position model
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ Initial.Treatment, 
             data=sub_meth_table1, family=binomial)
  a <- anova(fit, test="Chisq")
  
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,9],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   #pval.position = a$`Pr(>Chi)`[3], #uncomment if you want to include position of CpG within a gene
                   #pval.treatment_x_position = a$`Pr(>Chi)`[4], #uncomment if you want to include position of CpG within a gene interaction with treatment
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.Initial.Treatment <- p.adjust(results$pval.treatment, method='BH')
#results$adj.pval.position <- p.adjust(results$pval.position, method='BH') #uncomment if you want to include position of CpG within a gene
#results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH') #uncomment if you want to include position of CpG within a gene interaction with treatment

write.table(results, 'Output/D10_HP_GLM.tsv', sep='\t', row.names=FALSE)

D10.sig <-results
D10.sig <- D10.sig[,c(1,3)]
D10.sig <- D10.sig[order(D10.sig$adj.pval.Initial.Treatment),]
D10.sig <- D10.sig[which(D10.sig$adj.pval.Initial.Treatment<0.05), ]

write.table(D10.sig, 'Output/D10_HP_GLM_SortedResults.tsv', sep='\t', row.names=FALSE)

sub_meth_table$per.meth <- (sub_meth_table$meth/(sub_meth_table$meth+sub_meth_table$unmeth))*100 #calculate percent methylation
sub_meth_table.D10 <- sub_meth_table #store in day-specific dataframe

# Annotation of D10 DMG 
D10.sig.annot <- merge(D10.sig, Annot, by="gene", all.x=TRUE)
D10.sig.annot <- merge(D10.sig.annot, GOs, by="uniprot.id", all.x=TRUE)
D10.sig.annot <- D10.sig.annot[order(D10.sig.annot$adj.pval.Initial.Treatment),]
D10.sig.annot <- D10.sig.annot[!duplicated(D10.sig.annot$gene),]
write.table(D10.sig.annot, 'Output/Table_S2_D10_sig_annot.tsv', sep='\t', row.names=FALSE)

#####GO Enrichment of D10 DMGs
#Day 10 GO Enrichment for D10 DMGs with GO MWU approach https://github.com/z0on/GO_MWU #####
#identify MWU data
D10.MWU <- nonredun.D10.trt.meth.annot
#write MWU data to file
write.table(D10.MWU, 'GO_MWU/D10.GO.MWU', sep='\t', row.names=FALSE)
#merge all genes with DMGs
D10.MWU <- merge(D10.MWU, D10.sig.annot, by="gene", all=TRUE)
D10.MWU <- D10.MWU[,c(1,4)]
colnames(D10.MWU) <- c("gene", "DMG")
D10.MWU$DMG[D10.MWU$DMG < 0.05] <- 1
D10.MWU[is.na(D10.MWU)] <- 0

write.table(D10.MWU, 'GO_MWU/D10_interest_gsig.MWU', sep=',', row.names=FALSE)

##### Day135 #####
D135.trt.meth.annot <- merge(D135.trt.meth, Annot, by="gene", all.x=TRUE)
D135.trt.meth.annot <- merge(D135.trt.meth.annot, GOs, by="uniprot.id", all.x=TRUE)
D135.trt.meth.annot[is.na(D135.trt.meth.annot)] <- "unknown"
D135.trt.meth.annot[D135.trt.meth.annot==""] <- "unknown"
D135.trt.meth.annot.MWU <- as.data.frame(cbind(as.character(D135.trt.meth.annot$gene), D135.trt.meth.annot$GO.IDs))
colnames(D135.trt.meth.annot.MWU) <- c("gene", "GO.IDs")
write.table(D135.trt.meth.annot.MWU, 'GO_MWU/D135.trt.meth.annot.MWU', sep='\t', row.names=FALSE)
nonredun.D135.trt.meth.annot <- D135.trt.meth.annot.MWU[!duplicated(D135.trt.meth.annot.MWU$gene),]
colnames(nonredun.D135.trt.meth.annot) <- c("gene", "GO.IDs")

meth_table  <- read.csv("Output/D135.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment == 'Ambient' | meth_table$treatment == 'Low' | meth_table$treatment == 'Super.Low', ]

sub_meth_table$group <- paste0(sub_meth_table$Sample.ID, sub_meth_table$gene)

#filter for genes with >5 methylated positions
min.filt <- count(sub_meth_table, vars = c( group))
newdata <- min.filt[ which(min.filt$n > 4), ]
sub_meth_table <- sub_meth_table[sub_meth_table$group %in% newdata$vars,]

# HP Model
# create data frame to store results
results <- data.frame()
gs <- unique(sub_meth_table$gene)

#first subset the unique dataframes and second run the GLMs
for(i in 1:length(sub_meth_table$gene)){
  
  #subset the dataframe gene by gene
  sub_meth_table1 <- subset(sub_meth_table, gene ==gs[i])
  
  #if(length(unique(sub_meth_table1$position)) >1){
  
  # fit glm position model
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment, 
             data=sub_meth_table1, family=binomial)
  #} else {
  #  # fit glm model
  #  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
  #             data=sub_meth_table1, family=binomial)
  #}
  #s <- step(fit, trace=0)
  a <- anova(fit, test="Chisq")
  a
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,6],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   #pval.position = a$`Pr(>Chi)`[3],
                   #pval.treatment_x_position = a$`Pr(>Chi)`[4],
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.treatment <- p.adjust(results$pval.treatment, method='BH')
#results$adj.pval.position <- p.adjust(results$pval.position, method='BH')
#results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH')

write.table(results, 'Output/D135_HP_GLM.tsv', sep='\t', row.names=FALSE)

D135.sig <-results
D135.sig <- D135.sig[,c(1,3)]
D135.sig <- D135.sig[order(D135.sig$adj.pval.treatment),]
D135.sig <- D135.sig[which(D135.sig$adj.pval.treatment<0.05), ]

write.table(D135.sig, 'Output/D135_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)

sub_meth_table$per.meth <- (sub_meth_table$meth/(sub_meth_table$meth+sub_meth_table$unmeth))*100
sub_meth_table.D135 <- sub_meth_table

# Annotation of D135 DMG 
D135.sig.annot <- merge(D135.sig, Annot, by="gene", all.x=TRUE)
D135.sig.annot <- merge(D135.sig.annot, GOs, by="uniprot.id", all.x=TRUE)
D135.sig.annot <- D135.sig.annot[order(D135.sig.annot$adj.pval.treatment),]
D135.sig.annot <- D135.sig.annot[!duplicated(D135.sig.annot$gene),]
write.table(D135.sig.annot, 'Output/Table_S2_D135_sig_annot.tsv', sep='\t', row.names=FALSE)

#Day 135 GO MWU approach https://github.com/z0on/GO_MWU #####
#identify MWU data
D135.MWU <- nonredun.D135.trt.meth.annot
#write MWU data to file
write.table(D135.MWU, 'GO_MWU/D135.GO.MWU', sep='\t', row.names=FALSE)
#merge all genes with DMGs
D135.MWU <- merge(D135.MWU, D135.sig.annot, by="gene", all=TRUE)
D135.MWU <- D135.MWU[,c(1,4)]
colnames(D135.MWU) <- c("gene", "DMG")
D135.MWU$DMG[D135.MWU$DMG <0.05] <- 1
D135.MWU[is.na(D135.MWU)] <- 0
write.table(D135.MWU, 'GO_MWU/D135_interest_gsig.MWU', sep=',', row.names=FALSE)


# Day 145 Initial by Secondary Interaction Test #####
D145.trt.meth.annot <- merge(D145.trt.meth, Annot, by="gene", all.x=TRUE)
D145.trt.meth.annot <- merge(D145.trt.meth.annot, GOs, by="uniprot.id", all.x=TRUE)
D145.trt.meth.annot[is.na(D145.trt.meth.annot)] <- "unknown"
D145.trt.meth.annot[D145.trt.meth.annot==""] <- "unknown"
D145.trt.meth.annot.MWU <- as.data.frame(cbind(as.character(D145.trt.meth.annot$gene), D145.trt.meth.annot$GO.IDs))
colnames(D145.trt.meth.annot.MWU) <- c("gene", "GO.IDs")
write.table(D145.trt.meth.annot.MWU, 'GO_MWU/D145.trt.meth.annot.MWU', sep='\t', row.names=FALSE)
nonredun.D145.trt.meth.annot <- D145.trt.meth.annot.MWU[!duplicated(D145.trt.meth.annot.MWU$gene),]
colnames(nonredun.D145.trt.meth.annot) <- c("gene", "GO.IDs")

meth_table <- read.csv("Output/sec.meth.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment1 == "Ambient" | meth_table$treatment1 == "Low" | meth_table$treatment1 == "Super.Low"
                             | meth_table$treatment2 == "Ambient" | meth_table$treatment2 == "Low", ]
sub_meth_table$group <- paste0(sub_meth_table$Sample.ID, sub_meth_table$gene)

#filter for genes with >5 methylated positions
min.filt <- count(sub_meth_table, vars = c( group))
newdata <- min.filt[ which(min.filt$n > 4), ]
sub_meth_table <- sub_meth_table[sub_meth_table$group %in% newdata$vars,]

# create data frame to store results
results <- data.frame()
gs <- unique(sub_meth_table$gene)

#first subset the unique dataframes and second run the GLMs
for(i in 1:length(sub_meth_table$gene)){
    
    #subset the dataframe gene by gene
    sub_meth_table1 <- subset(sub_meth_table, gene ==gs[i])
    
    # fit glm position model
    fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
               data=sub_meth_table1, family=binomial)

    a <- anova(fit, test="LRT")
    
    # capture summary stats to data frame
    df <- data.frame(gene = sub_meth_table1[1,6],
                     pval.treatment1 = a$`Pr(>Chi)`[2],
                     pval.treatment2 = a$`Pr(>Chi)`[3],
                     #position = a$`Pr(>Chi)`[2], 
                     pval.treatment1_x_treatment2 = a$`Pr(>Chi)`[4],
                     #pval.treatment1_x_position = a$`Pr(>Chi)`[2],
                     #pval.treatment2_x_position = a$`Pr(>Chi)`[2],
                     #pval.treatment1_x_pval.treatment2_x_position = a$`Pr(>Chi)`[2],
                     stringsAsFactors = F)

    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
    
  }

results[is.na(results)] <- 0

results$adj.pval.treatment1 <- p.adjust(results$pval.treatment1, method='BH')
results$adj.pval.treatment2 <- p.adjust(results$pval.treatment2, method='BH')
results$adj.pval.treatment1_x_treatment2 <- p.adjust(results$pval.treatment1_x_treatment2, method='BH')
#results$adj.pval.treatment1_x_position <- p.adjust(results$pval.treatment1_x_position, method='BH')
#results$adj.pval.treatment2_x_position <- p.adjust(results$pval.treatment2_x_position, method='BH')
#results$adj.pval.treatment1_x_pval.treatment2_x_position <- p.adjust(results$pval.treatment1_x_pval.treatment2_x_position, method='BH')

write.table(results, 'Output/sec.meth_HP_GLM.all.tsv', sep='\t', row.names=FALSE)

D145.sig <-results
D145.sig <- D145.sig[,c(1,5:7)]
D145.sig <- D145.sig[order(results$adj.pval.treatment1_x_treatment2),]
D145.sig <- D145.sig[which(D145.sig$adj.pval.treatment1_x_treatment2<0.05), ]

write.table(D145.sig, 'Output/Interaction_D145_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)

sub_meth_table$per.meth <- (sub_meth_table$meth/(sub_meth_table$meth+sub_meth_table$unmeth))*100
sub_meth_table.D145 <- sub_meth_table

# Annotation of Secondary Exposure Interaction DMG 
D145.sig.annot <- merge(D145.sig, Annot, by="gene", all.x=TRUE)
D145.sig.annot <- merge(D145.sig.annot, GOs, by="uniprot.id", all.x=TRUE)
D145.sig.annot <- D145.sig.annot[order(D145.sig.annot$adj.pval.treatment1_x_treatment2),]
D145.sig.annot <- D145.sig.annot[!duplicated(D145.sig.annot$gene),]
write.table(D145.sig.annot, 'Output/Table_S2_D145_sig_annot.tsv', sep='\t', row.names=FALSE)

#Day 145 GO MWU approach https://github.com/z0on/GO_MWU #####
#identify MWU data
D145.MWU <- nonredun.D145.trt.meth.annot
#write MWU data to file
write.table(D145.MWU, 'GO_MWU/D145.GO.MWU', sep='\t', row.names=FALSE)
#merge all genes with DMGs
D145.MWU <- merge(D145.MWU, D145.sig.annot, by="gene", all=TRUE)
D145.MWU <- D145.MWU[,c(1,4)]
colnames(D145.MWU) <- c("gene", "DMG")
D145.MWU$DMG[D145.MWU$DMG <0.05] <- 1
D145.MWU[is.na(D145.MWU)] <- 0

write.table(D145.MWU, 'GO_MWU/D145_interest_gsig.MWU', sep=',', row.names=FALSE)

# for (i in 1:20){
#   sig.gene <- subset(sub_meth_table, gene ==sec.data.sig.int.sig.annot$gene[i])
#   means <- aggregate(per.meth ~ treatment1*treatment2*gene*position, data=sig.gene, FUN=mean)
#   
#   temp_plot <- ggplot(means, aes(x=treatment2, y=per.meth, fill=treatment1)) +
#     geom_violin(aes(fill = treatment1), trim = FALSE) +
#     geom_boxplot(width = 0.1, position = position_dodge(0.9), color="black")+
#     scale_fill_manual(values = c("blue", "purple", "red"))+
#     theme_bw()
#   
#   ggsave(temp_plot, file=paste0("Output/D145_Plots/plot_", sec.data.sig.int.sig.annot$protein.name[i],".png"), width = 14, height = 10, units = "cm")
# }

##### DMGs in all comparisons #####

ints_10_135 <- merge(D10.sig, D135.sig, by="gene")
ints_10_145 <- merge(D10.sig, D145.sig, by="gene")
ints_135_145 <- merge(D135.sig, D145.sig, by="gene")
ints_all <- merge(ints_10_135, D145.sig, by="gene")
ints_all.annot <- merge(ints_all, Annot, by="gene", all.x=TRUE)
ints_all.annot <- ints_all.annot[!duplicated(ints_all.annot$gene),]
write.table(ints_all.annot, 'Output/DMGs_all_timepoints_annotated.tsv', sep='\t', row.names=FALSE)

ints_all.annot

Consistent.DMG.GO <- merge(ints_all.annot, GOs, by="uniprot.id", all=TRUE)
Consistent.DMG.GO <- na.omit(Consistent.DMG.GO)

All.Ids <- as.character(Consistent.DMG.GO$GO.IDs)
myCollection <- GOCollection(All.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
#fl <- "http://current.geneontology.org/ontology/go-basic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
All.MF.slim <- goSlim(myCollection, slim, "MF")
All.BP.slim <- goSlim(myCollection, slim, "BP")
All.CC.slim <- goSlim(myCollection, slim, "CC")


#####
# for (i in 1:nrow(ints_all.annot)){
#   sig.gene <- subset(sub_meth_table.D10, gene ==ints_all.annot$gene[i])
#   means <- aggregate(per.meth ~ Initial.Treatment*gene, data=sig.gene, FUN=mean)
#   ses <- aggregate(per.meth ~ Initial.Treatment*gene, data=sig.gene, FUN=std.error)
#   means$ses <- ses$per.meth
#   
#   temp_plot <- ggplot(means, aes(x=Initial.Treatment, y=per.meth, fill=Initial.Treatment)) +
#     geom_bar(stat="identity", color="black", 
#              position=position_dodge()) +
#     geom_errorbar(aes(ymin=per.meth-ses, ymax=per.meth+ses), width=.2,
#                   position=position_dodge(.9)) +
#     xlab("Treatment") + #plot x axis label
#     ylab("Percent Methylation") + #plot y axis label
#     scale_fill_manual(values = c("blue", "purple", "red"))+
#     theme_bw()
#   
#   ggsave(temp_plot, file=paste0("Output/Consistent_plots/D10plot_", ints_all.annot$gene[i],"_",ints_all.annot$protein.name[i],".png"), width = 14, height = 10, units = "cm")
# }
# 
# for (i in 1:nrow(ints_all.annot)){
#   sig.gene <- subset(sub_meth_table.D135, gene ==ints_all.annot$gene[i])
#   means <- aggregate(per.meth ~ treatment*gene, data=sig.gene, FUN=mean)
#   ses <- aggregate(per.meth ~ treatment*gene, data=sig.gene, FUN=std.error)
#   means$ses <- ses$per.meth
#   
#   temp_plot <- ggplot(means, aes(x=treatment, y=per.meth, fill=treatment)) +
#     geom_bar(stat="identity", color="black", 
#              position=position_dodge()) +
#     geom_errorbar(aes(ymin=per.meth-ses, ymax=per.meth+ses), width=.2,
#                   position=position_dodge(.9)) +
#     xlab("Treatment") + #plot x axis label
#     ylab("Percent Methylation") + #plot y axis label
#     scale_fill_manual(values = c("blue", "purple", "red"))+
#     theme_bw()
#   
#   ggsave(temp_plot, file=paste0("Output/Consistent_plots/D135plot_", ints_all.annot$gene[i],"_",ints_all.annot$protein.name[i],".png"), width = 14, height = 10, units = "cm")
# }
# 
# for (i in 1:nrow(ints_all.annot)){
#   sig.gene <- subset(sub_meth_table.D145, gene ==ints_all.annot$gene[i])
#   means <- aggregate(per.meth ~ treatment1*treatment2*gene, data=sig.gene, FUN=mean)
#   ses <- aggregate(per.meth ~ treatment1*treatment2*gene, data=sig.gene, FUN=std.error)
#   means$ses <- ses$per.meth
#   
#   temp_plot <- ggplot(data=means, aes(x=treatment2, y=per.meth, group=treatment1, colour=treatment1, shape=treatment1)) + #plot data
#     geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
#     scale_colour_manual(values = c("blue", "purple", "red")) + #set line color
#     geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
#     scale_shape_manual(values=c(1,18,16)) + #set shapes
#     geom_errorbar(aes(ymin=per.meth-ses, ymax=per.meth+ses), #plot error bars
#                   width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
#     xlab("Secondary Treatment") + #plot x axis label
#     ylab("Percent Methylation") + #plot y axis label
#     theme_bw()
#   
#   ggsave(temp_plot, file=paste0("Output/Consistent_plots/D145plot_",ints_all.annot$gene[i],"_", ints_all.annot$protein.name[i],".png"), width = 14, height = 10, units = "cm")
# }


##### Heatmap of DMG in Ambient by time #####
Time.Amb.meth_table.means <- aggregate(per.meth ~ TimePoint*gene, data=sub_meth_table.Time.Amb, FUN=mean)

Time.Amb.meth_table.means <- Time.Amb.meth_table.means[Time.Amb.meth_table.means$gene %in% Time.Amb.sig$gene,]

Time.Amb <- Time.Amb.meth_table.means %>%
  pivot_wider(names_from = TimePoint, values_from = per.meth)

Time.Amb.mat <- Time.Amb[,2:5]
Time.Amb.mat <- Time.Amb.mat-rowMeans(Time.Amb.mat)



pheatmap(Time.Amb.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =TRUE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

pdf("Output/Time.Amb.DMG_relative_Meth_All.pdf")
pheatmap(Time.Amb.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

##### Heatmap of D10 DMGs Sig Treatment #####
D10.meth_table.means <- aggregate(per.meth ~ Initial.Treatment*gene, data=sub_meth_table.D10, FUN=mean)

D10.meth_table.means <- D10.meth_table.means[D10.meth_table.means$gene %in% D10.sig$gene,]

D10 <- D10.meth_table.means %>%
  pivot_wider(names_from = Initial.Treatment, values_from = per.meth)
#D10 <- D10[D10$gene %in% ints_all.annot$gene,]

D10.mat <- D10[,2:4]
D10.mat <- D10.mat-rowMeans(D10.mat)

pheatmap(D10.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =TRUE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

pdf("Output/D10.DMG_relative_Meth_All.pdf")
pheatmap(D10.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()
##### Heatmap of D135 DMGs Sig Treatment #####
D135.meth_table.means <- aggregate(per.meth ~ treatment*gene, data=sub_meth_table.D135, FUN=mean)

D135.meth_table.means <- D135.meth_table.means[D135.meth_table.means$gene %in% D135.sig$gene,]

D135 <- D135.meth_table.means %>%
  pivot_wider(names_from = treatment, values_from = per.meth)

#D135 <- D135[D135$gene %in% ints_all.annot$gene,]

D135.mat <- D135[,2:4]

D135.mat <- D135.mat-rowMeans(D135.mat)

pheatmap(D135.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =TRUE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

pdf("Output/D135.DMG_relative_Meth_All.pdf")
pheatmap(D135.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

##### Heatmap of Secondary Exposure Interaction DMG  DMGs Sig Treatment #####
sec.int.meth_table.means <- aggregate(per.meth ~ treatment1*treatment2*gene, data=sub_meth_table.D145, FUN=mean)

sec.int.meth_table.means <- sec.int.meth_table.means[sec.int.meth_table.means$gene %in% D145.sig$gene,]

D145 <- sec.int.meth_table.means %>%
  pivot_wider(names_from = c(treatment1,treatment2), values_from = per.meth)

#D145 <- D145[D145$gene %in% ints_all.annot$gene,]

D145.mat <- D145[,2:7]
D145.mat <- D145.mat-rowMeans(D145.mat)

col.order <- c("Ambient_Ambient", "Ambient_Low", "Low_Ambient", "Low_Low", "Super.Low_Ambient", "Super.Low_Low")
D145.mat  <- D145.mat [,col.order]

pdf("Output/D145.DMG_relative_Meth_All.pdf")
pheatmap(D145.mat, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=FALSE, gaps_col = c(2,4),
          show_colnames =TRUE) #plot heatmap of all DEG by group
dev.off()

#####GO Enrichment of consistent DMG

#GO MWU approach https://github.com/z0on/GO_MWU #####
#identify MWU data
ALL.MWU <- nonredun.D10.trt.meth.annot
#write MWU data to file
write.table(ALL.MWU, 'GO_MWU/ALL.GO.MWU', sep='\t', row.names=FALSE)
#merge all genes with DMGs
ALL.MWU.Data <- merge(ALL.MWU, ints_all.annot, by="gene", all=TRUE)
#ALL.MWU.Data <- ALL.MWU.Data[,c(1,2)]
#ALL.MWU.Data <- ALL.MWU.Data[-1,]
#colnames(ALL.MWU.Data) <- c("gene", "DMG")
KNOWN.ALL.MWU.Data <- subset(ALL.MWU.Data, GO.IDs!="unknown")

KNOWN.ALL.MWU.Data$DMG <- rowSums(KNOWN.ALL.MWU.Data[,3:7])
KNOWN.ALL.MWU.Data$DMG[KNOWN.ALL.MWU.Data$DMG >0] <- 1
KNOWN.ALL.MWU.Data[is.na(KNOWN.ALL.MWU.Data)] <- 0

KNOWN.ALL.MWU.Data <- KNOWN.ALL.MWU.Data[,c(1,18)]
write.table(KNOWN.ALL.MWU.Data, 'GO_MWU/ALL_interest_gsig.MWU', sep=',', row.names=FALSE)



