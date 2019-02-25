#Title: Geoduck juvenile methylation
#Project: NOAA FFAR
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20190217
#See Readme file for details

rm(list=ls()) #clears workspace 


library(plotrix) 
library(ggplot2)
library(gridExtra)
library(vegan)
library(pca3d)
library(pheatmap)
library(dplyr)
library(tidyverse)

# Set Working Directory:
setwd("~/MyProjects/Geoduck_Meth/RAnalysis/") #set working

sample.info <- read.csv("Data/Sample.Info.csv", header=T, sep=",", na.string="NA", stringsAsFactors = F)

Data <- read.csv("Data/Extracted/10x.meth.genes.intersect", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
Data <- Data[,1:3]
colnames(Data) <- c("scaffold", "start", "stop")
Data$location <- Data$start +1
Data$merger <- paste0(Data$scaffold,"_", Data$location)

Genes <- read.csv("Data/Genome/Pgenerosa_v070.a.makergene.gff", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
Genes <- Genes[,c(1,4,5,9)]
colnames(Genes) <- c("scaffold", "start", "stop", "gene")
Genes$start <-as.numeric(Genes$start)
Genes$stop <-as.numeric(Genes$stop)


for(i in 1:length(sample.info$Sample.ID)){
  epi <- read.csv(paste0("Data/Extracted/",sample.info$Sample.ID[i],".10x.cov.CpG"), header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
  epi <- epi[,-3]
  colnames(epi) <- c("scaffold", "position","per.meth","meth","unmeth")
  epi$position <- as.numeric(epi$position)
  epi <- full_join(epi, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
  epi$Sample.ID <- sample.info$Sample.ID[i]
  epi$merger <- paste0(epi$scaffold,"_", epi$position)
  epi <- merge(Data, epi, by="merger")
  epi <- epi[,6:14]
  colnames(epi) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")
  
}

  
  

##### Load in methylation counts #####
epi.41 <- read.csv("Data/Extracted/EPI-41.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.41 <- epi.41[,-3]
colnames(epi.41) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.41$position <- as.numeric(epi.41$position)
epi.41 <- full_join(epi.41, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.41$Sample.ID <- "EPI_41"
epi.41$merger <- paste0(epi.41$scaffold,"_", epi.41$position)
epi.41 <- merge(Data, epi.41, by="merger")
epi.41 <- epi.41[,6:14]
colnames(epi.41) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.42 <- read.csv("Data/Extracted/EPI-42.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.42 <- epi.42[,-3]
colnames(epi.42) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.42$position <- as.numeric(epi.42$position)
epi.42 <- full_join(epi.42, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.42$Sample.ID <- "EPI_42"
epi.42$merger <- paste0(epi.42$scaffold,"_", epi.42$position)
epi.42 <- merge(Data, epi.42, by="merger")
epi.42 <- epi.42[,6:14]
colnames(epi.42) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.43 <- read.csv("Data/Extracted/EPI-43.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.43 <- epi.43[,-3]
colnames(epi.43) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.43$position <- as.numeric(epi.43$position)
epi.43 <- full_join(epi.43, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.43$Sample.ID <- "EPI_43"
epi.43$merger <- paste0(epi.43$scaffold,"_", epi.43$position)
epi.43 <- merge(Data, epi.43, by="merger")
epi.43 <- epi.43[,6:14]
colnames(epi.43) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.44 <- read.csv("Data/Extracted/EPI-44.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.44 <- epi.44[,-3]
colnames(epi.44) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.44$position <- as.numeric(epi.44$position)
epi.44 <- full_join(epi.44, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.44$Sample.ID <- "EPI_44"
epi.44$merger <- paste0(epi.44$scaffold,"_", epi.44$position)
epi.44 <- merge(Data, epi.44, by="merger")
epi.44 <- epi.44[,6:14]
colnames(epi.44) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.103 <- read.csv("Data/Extracted/EPI-103.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.103 <- epi.103[,-3]
colnames(epi.103) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.103$position <- as.numeric(epi.103$position)
epi.103 <- full_join(epi.103, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.103$Sample.ID <- "EPI_103"
epi.103$merger <- paste0(epi.103$scaffold,"_", epi.103$position)
epi.103 <- merge(Data, epi.103, by="merger")
epi.103 <- epi.103[,6:14]
colnames(epi.103) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.104 <- read.csv("Data/Extracted/EPI-104.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.104 <- epi.104[,-3]
colnames(epi.104) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.104$position <- as.numeric(epi.104$position)
epi.104 <- full_join(epi.104, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.104$Sample.ID <- "EPI_104"
epi.104$merger <- paste0(epi.104$scaffold,"_", epi.104$position)
epi.104 <- merge(Data, epi.104, by="merger")
epi.104 <- epi.104[,6:14]
colnames(epi.104) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.111 <- read.csv("Data/Extracted/EPI-111.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.111 <- epi.111[,-3]
colnames(epi.111) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.111$position <- as.numeric(epi.111$position)
epi.111 <- full_join(epi.111, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.111$Sample.ID <- "EPI_111"
epi.111$merger <- paste0(epi.111$scaffold,"_", epi.111$position)
epi.111 <- merge(Data, epi.111, by="merger")
epi.111 <- epi.111[,6:14]
colnames(epi.111) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.113 <- read.csv("Data/Extracted/EPI-113.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.113 <- epi.113[,-3]
colnames(epi.113) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.113$position <- as.numeric(epi.113$position)
epi.113 <- full_join(epi.113, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.113$Sample.ID <- "EPI_113"
epi.113$merger <- paste0(epi.113$scaffold,"_", epi.113$position)
epi.113 <- merge(Data, epi.113, by="merger")
epi.113 <- epi.113[,6:14]
colnames(epi.113) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.119 <- read.csv("Data/Extracted/EPI-119.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.119 <- epi.119[,-3]
colnames(epi.119) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.119$position <- as.numeric(epi.119$position)
epi.119 <- full_join(epi.119, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.119$Sample.ID <- "EPI_119"
epi.119$merger <- paste0(epi.119$scaffold,"_", epi.119$position)
epi.119 <- merge(Data, epi.119, by="merger")
epi.119 <- epi.119[,6:14]
colnames(epi.119) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.120 <- read.csv("Data/Extracted/EPI-120.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.120 <- epi.120[,-3]
colnames(epi.120) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.120$position <- as.numeric(epi.120$position)
epi.120 <- full_join(epi.120, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.120$Sample.ID <- "EPI_120"
epi.120$merger <- paste0(epi.120$scaffold,"_", epi.120$position)
epi.120 <- merge(Data, epi.120, by="merger")
epi.120 <- epi.120[,6:14]
colnames(epi.120) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.127 <- read.csv("Data/Extracted/EPI-127.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.127 <- epi.127[,-3]
colnames(epi.127) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.127$position <- as.numeric(epi.127$position)
epi.127 <- full_join(epi.127, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.127$Sample.ID <- "EPI_127"
epi.127$merger <- paste0(epi.127$scaffold,"_", epi.127$position)
epi.127 <- merge(Data, epi.127, by="merger")
epi.127 <- epi.127[,6:14]
colnames(epi.127) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.128 <- read.csv("Data/Extracted/EPI-128.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.128 <- epi.128[,-3]
colnames(epi.128) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.128$position <- as.numeric(epi.128$position)
epi.128 <- full_join(epi.128, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.128$Sample.ID <- "EPI_128"
epi.128$merger <- paste0(epi.128$scaffold,"_", epi.128$position)
epi.128 <- merge(Data, epi.128, by="merger")
epi.128 <- epi.128[,6:14]
colnames(epi.128) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.135 <- read.csv("Data/Extracted/EPI-135.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.135 <- epi.135[,-3]
colnames(epi.135) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.135$position <- as.numeric(epi.135$position)
epi.135 <- full_join(epi.135, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.135$Sample.ID <- "EPI_135"
epi.135$merger <- paste0(epi.135$scaffold,"_", epi.135$position)
epi.135 <- merge(Data, epi.135, by="merger")
epi.135 <- epi.135[,6:14]
colnames(epi.135) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.136 <- read.csv("Data/Extracted/EPI-136.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.136 <- epi.136[,-3]
colnames(epi.136) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.136$position <- as.numeric(epi.136$position)
epi.136 <- full_join(epi.136, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.136$Sample.ID <- "EPI_136"
epi.136$merger <- paste0(epi.136$scaffold,"_", epi.136$position)
epi.136 <- merge(Data, epi.136, by="merger")
epi.136 <- epi.136[,6:14]
colnames(epi.136) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.143 <- read.csv("Data/Extracted/EPI-143.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.143 <- epi.143[,-3]
colnames(epi.143) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.143$position <- as.numeric(epi.143$position)
epi.143 <- full_join(epi.143, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.143$Sample.ID <- "EPI_143"
epi.143$merger <- paste0(epi.143$scaffold,"_", epi.143$position)
epi.143 <- merge(Data, epi.143, by="merger")
epi.143 <- epi.143[,6:14]
colnames(epi.143) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.145 <- read.csv("Data/Extracted/EPI-145.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.145 <- epi.145[,-3]
colnames(epi.145) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.145$position <- as.numeric(epi.145$position)
epi.145 <- full_join(epi.145, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.145$Sample.ID <- "EPI_145"
epi.145$merger <- paste0(epi.145$scaffold,"_", epi.145$position)
epi.145 <- merge(Data, epi.145, by="merger")
epi.145 <- epi.145[,6:14]
colnames(epi.145) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.151 <- read.csv("Data/Extracted/EPI-151.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.151 <- epi.151[,-3]
colnames(epi.151) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.151$position <- as.numeric(epi.151$position)
epi.151 <- full_join(epi.151, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.151$Sample.ID <- "EPI_151"
epi.151$merger <- paste0(epi.151$scaffold,"_", epi.151$position)
epi.151 <- merge(Data, epi.151, by="merger")
epi.151 <- epi.151[,6:14]
colnames(epi.151) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.152 <- read.csv("Data/Extracted/EPI-152.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.152 <- epi.152[,-3]
colnames(epi.152) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.152$position <- as.numeric(epi.152$position)
epi.152 <- full_join(epi.152, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.152$Sample.ID <- "EPI_152"
epi.152$merger <- paste0(epi.152$scaffold,"_", epi.152$position)
epi.152 <- merge(Data, epi.152, by="merger")
epi.152 <- epi.152[,6:14]
colnames(epi.152) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.153 <- read.csv("Data/Extracted/EPI-153.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.153 <- epi.153[,-3]
colnames(epi.153) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.153$position <- as.numeric(epi.153$position)
epi.153 <- full_join(epi.153, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.153$Sample.ID <- "EPI_153"
epi.153$merger <- paste0(epi.153$scaffold,"_", epi.153$position)
epi.153 <- merge(Data, epi.153, by="merger")
epi.153 <- epi.153[,6:14]
colnames(epi.153) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.154 <- read.csv("Data/Extracted/EPI-154.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.154 <- epi.154[,-3]
colnames(epi.154) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.154$position <- as.numeric(epi.154$position)
epi.154 <- full_join(epi.154, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.154$Sample.ID <- "EPI_154"
epi.154$merger <- paste0(epi.154$scaffold,"_", epi.154$position)
epi.154 <- merge(Data, epi.154, by="merger")
epi.154 <- epi.154[,6:14]
colnames(epi.154) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.159 <- read.csv("Data/Extracted/EPI-159.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.159 <- epi.159[,-3]
colnames(epi.159) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.159$position <- as.numeric(epi.159$position)
epi.159 <- full_join(epi.159, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.159$Sample.ID <- "EPI_159"
epi.159$merger <- paste0(epi.159$scaffold,"_", epi.159$position)
epi.159 <- merge(Data, epi.159, by="merger")
epi.159 <- epi.159[,6:14]
colnames(epi.159) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.160 <- read.csv("Data/Extracted/EPI-160.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.160 <- epi.160[,-3]
colnames(epi.160) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.160$position <- as.numeric(epi.160$position)
epi.160 <- full_join(epi.160, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.160$Sample.ID <- "EPI_160"
epi.160$merger <- paste0(epi.160$scaffold,"_", epi.160$position)
epi.160 <- merge(Data, epi.160, by="merger")
epi.160 <- epi.160[,6:14]
colnames(epi.160) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.161 <- read.csv("Data/Extracted/EPI-161.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.161 <- epi.161[,-3]
colnames(epi.161) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.161$position <- as.numeric(epi.161$position)
epi.161 <- full_join(epi.161, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.161$Sample.ID <- "EPI_161"
epi.161$merger <- paste0(epi.161$scaffold,"_", epi.161$position)
epi.161 <- merge(Data, epi.161, by="merger")
epi.161 <- epi.161[,6:14]
colnames(epi.161) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.162 <- read.csv("Data/Extracted/EPI-162.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.162 <- epi.162[,-3]
colnames(epi.162) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.162$position <- as.numeric(epi.162$position)
epi.162 <- full_join(epi.162, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.162$Sample.ID <- "EPI_162"
epi.162$merger <- paste0(epi.162$scaffold,"_", epi.162$position)
epi.162 <- merge(Data, epi.162, by="merger")
epi.162 <- epi.162[,6:14]
colnames(epi.162) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.167 <- read.csv("Data/Extracted/EPI-167.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.167 <- epi.167[,-3]
colnames(epi.167) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.167$position <- as.numeric(epi.167$position)
epi.167 <- full_join(epi.167, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.167$Sample.ID <- "EPI_167"
epi.167$merger <- paste0(epi.167$scaffold,"_", epi.167$position)
epi.167 <- merge(Data, epi.167, by="merger")
epi.167 <- epi.167[,6:14]
colnames(epi.167) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.168 <- read.csv("Data/Extracted/EPI-168.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.168 <- epi.168[,-3]
colnames(epi.168) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.168$position <- as.numeric(epi.168$position)
epi.168 <- full_join(epi.168, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.168$Sample.ID <- "EPI_168"
epi.168$merger <- paste0(epi.168$scaffold,"_", epi.168$position)
epi.168 <- merge(Data, epi.168, by="merger")
epi.168 <- epi.168[,6:14]
colnames(epi.168) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.169 <- read.csv("Data/Extracted/EPI-169.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.169 <- epi.169[,-3]
colnames(epi.169) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.169$position <- as.numeric(epi.169$position)
epi.169 <- full_join(epi.169, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.169$Sample.ID <- "EPI_169"
epi.169$merger <- paste0(epi.169$scaffold,"_", epi.169$position)
epi.169 <- merge(Data, epi.169, by="merger")
epi.169 <- epi.169[,6:14]
colnames(epi.169) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.170 <- read.csv("Data/Extracted/EPI-170.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.170 <- epi.170[,-3]
colnames(epi.170) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.170$position <- as.numeric(epi.170$position)
epi.170 <- full_join(epi.170, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.170$Sample.ID <- "EPI_170"
epi.170$merger <- paste0(epi.170$scaffold,"_", epi.170$position)
epi.170 <- merge(Data, epi.170, by="merger")
epi.170 <- epi.170[,6:14]
colnames(epi.170) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.175 <- read.csv("Data/Extracted/EPI-175.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.175 <- epi.175[,-3]
colnames(epi.175) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.175$position <- as.numeric(epi.175$position)
epi.175 <- full_join(epi.175, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.175$Sample.ID <- "EPI_175"
epi.175$merger <- paste0(epi.175$scaffold,"_", epi.175$position)
epi.175 <- merge(Data, epi.175, by="merger")
epi.175 <- epi.175[,6:14]
colnames(epi.175) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.176 <- read.csv("Data/Extracted/EPI-176.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.176 <- epi.176[,-3]
colnames(epi.176) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.176$position <- as.numeric(epi.176$position)
epi.176 <- full_join(epi.176, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.176$Sample.ID <- "EPI_176"
epi.176$merger <- paste0(epi.176$scaffold,"_", epi.176$position)
epi.176 <- merge(Data, epi.176, by="merger")
epi.176 <- epi.176[,6:14]
colnames(epi.176) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.181 <- read.csv("Data/Extracted/EPI-181.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.181 <- epi.181[,-3]
colnames(epi.181) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.181$position <- as.numeric(epi.181$position)
epi.181 <- full_join(epi.181, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.181$Sample.ID <- "EPI_181"
epi.181$merger <- paste0(epi.181$scaffold,"_", epi.181$position)
epi.181 <- merge(Data, epi.181, by="merger")
epi.181 <- epi.181[,6:14]
colnames(epi.181) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.182 <- read.csv("Data/Extracted/EPI-182.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.182 <- epi.182[,-3]
colnames(epi.182) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.182$position <- as.numeric(epi.182$position)
epi.182 <- full_join(epi.182, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.182$Sample.ID <- "EPI_182"
epi.182$merger <- paste0(epi.182$scaffold,"_", epi.182$position)
epi.182 <- merge(Data, epi.182, by="merger")
epi.182 <- epi.182[,6:14]
colnames(epi.182) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.184 <- read.csv("Data/Extracted/EPI-184.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.184 <- epi.184[,-3]
colnames(epi.184) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.184$position <- as.numeric(epi.184$position)
epi.184 <- full_join(epi.184, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.184$Sample.ID <- "EPI_184"
epi.184$merger <- paste0(epi.184$scaffold,"_", epi.184$position)
epi.184 <- merge(Data, epi.184, by="merger")
epi.184 <- epi.184[,6:14]
colnames(epi.184) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.185 <- read.csv("Data/Extracted/EPI-185.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.185 <- epi.185[,-3]
colnames(epi.185) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.185$position <- as.numeric(epi.185$position)
epi.185 <- full_join(epi.185, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.185$Sample.ID <- "EPI_185"
epi.185$merger <- paste0(epi.185$scaffold,"_", epi.185$position)
epi.185 <- merge(Data, epi.185, by="merger")
epi.185 <- epi.185[,6:14]
colnames(epi.185) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.187 <- read.csv("Data/Extracted/EPI-187.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.187 <- epi.187[,-3]
colnames(epi.187) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.187$position <- as.numeric(epi.187$position)
epi.187 <- full_join(epi.187, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.187$Sample.ID <- "EPI_187"
epi.187$merger <- paste0(epi.187$scaffold,"_", epi.187$position)
epi.187 <- merge(Data, epi.187, by="merger")
epi.187 <- epi.187[,6:14]
colnames(epi.187) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.188 <- read.csv("Data/Extracted/EPI-188.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.188 <- epi.188[,-3]
colnames(epi.188) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.188$position <- as.numeric(epi.188$position)
epi.188 <- full_join(epi.188, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.188$Sample.ID <- "EPI_188"
epi.188$merger <- paste0(epi.188$scaffold,"_", epi.188$position)
epi.188 <- merge(Data, epi.188, by="merger")
epi.188 <- epi.188[,6:14]
colnames(epi.188) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.193 <- read.csv("Data/Extracted/EPI-193.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.193 <- epi.193[,-3]
colnames(epi.193) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.193$position <- as.numeric(epi.193$position)
epi.193 <- full_join(epi.193, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.193$Sample.ID <- "EPI_193"
epi.193$merger <- paste0(epi.193$scaffold,"_", epi.193$position)
epi.193 <- merge(Data, epi.193, by="merger")
epi.193 <- epi.193[,6:14]
colnames(epi.193) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.194 <- read.csv("Data/Extracted/EPI-194.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.194 <- epi.194[,-3]
colnames(epi.194) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.194$position <- as.numeric(epi.194$position)
epi.194 <- full_join(epi.194, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.194$Sample.ID <- "EPI_194"
epi.194$merger <- paste0(epi.194$scaffold,"_", epi.194$position)
epi.194 <- merge(Data, epi.194, by="merger")
epi.194 <- epi.194[,6:14]
colnames(epi.194) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.199 <- read.csv("Data/Extracted/EPI-199.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.199 <- epi.199[,-3]
colnames(epi.199) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.199$position <- as.numeric(epi.199$position)
epi.199 <- full_join(epi.199, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.199$Sample.ID <- "EPI_199"
epi.199$merger <- paste0(epi.199$scaffold,"_", epi.199$position)
epi.199 <- merge(Data, epi.199, by="merger")
epi.199 <- epi.199[,6:14]
colnames(epi.199) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.200 <- read.csv("Data/Extracted/EPI-200.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.200 <- epi.200[,-3]
colnames(epi.200) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.200$position <- as.numeric(epi.200$position)
epi.200 <- full_join(epi.200, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.200$Sample.ID <- "EPI_200"
epi.200$merger <- paste0(epi.200$scaffold,"_", epi.200$position)
epi.200 <- merge(Data, epi.200, by="merger")
epi.200 <- epi.200[,6:14]
colnames(epi.200) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.205 <- read.csv("Data/Extracted/EPI-205.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.205 <- epi.205[,-3]
colnames(epi.205) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.205$position <- as.numeric(epi.205$position)
epi.205 <- full_join(epi.205, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.205$Sample.ID <- "EPI_205"
epi.205$merger <- paste0(epi.205$scaffold,"_", epi.205$position)
epi.205 <- merge(Data, epi.205, by="merger")
epi.205 <- epi.205[,6:14]
colnames(epi.205) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.206 <- read.csv("Data/Extracted/EPI-206.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.206 <- epi.206[,-3]
colnames(epi.206) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.206$position <- as.numeric(epi.206$position)
epi.206 <- full_join(epi.206, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.206$Sample.ID <- "EPI_206"
epi.206$merger <- paste0(epi.206$scaffold,"_", epi.206$position)
epi.206 <- merge(Data, epi.206, by="merger")
epi.206 <- epi.206[,6:14]
colnames(epi.206) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.208 <- read.csv("Data/Extracted/EPI-208.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.208 <- epi.208[,-3]
colnames(epi.208) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.208$position <- as.numeric(epi.208$position)
epi.208 <- full_join(epi.208, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.208$Sample.ID <- "EPI_208"
epi.208$merger <- paste0(epi.208$scaffold,"_", epi.208$position)
epi.208 <- merge(Data, epi.208, by="merger")
epi.208 <- epi.208[,6:14]
colnames(epi.208) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.209 <- read.csv("Data/Extracted/EPI-209.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.209 <- epi.209[,-3]
colnames(epi.209) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.209$position <- as.numeric(epi.209$position)
epi.209 <- full_join(epi.209, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.209$Sample.ID <- "EPI_209"
epi.209$merger <- paste0(epi.209$scaffold,"_", epi.209$position)
epi.209 <- merge(Data, epi.209, by="merger")
epi.209 <- epi.209[,6:14]
colnames(epi.209) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.214 <- read.csv("Data/Extracted/EPI-214.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.214 <- epi.214[,-3]
colnames(epi.214) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.214$position <- as.numeric(epi.214$position)
epi.214 <- full_join(epi.214, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.214$Sample.ID <- "EPI_214"
epi.214$merger <- paste0(epi.214$scaffold,"_", epi.214$position)
epi.214 <- merge(Data, epi.214, by="merger")
epi.214 <- epi.214[,6:14]
colnames(epi.214) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.215 <- read.csv("Data/Extracted/EPI-215.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.215 <- epi.215[,-3]
colnames(epi.215) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.215$position <- as.numeric(epi.215$position)
epi.215 <- full_join(epi.215, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.215$Sample.ID <- "EPI_215"
epi.215$merger <- paste0(epi.215$scaffold,"_", epi.215$position)
epi.215 <- merge(Data, epi.215, by="merger")
epi.215 <- epi.215[,6:14]
colnames(epi.215) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.220 <- read.csv("Data/Extracted/EPI-220.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.220 <- epi.220[,-3]
colnames(epi.220) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.220$position <- as.numeric(epi.220$position)
epi.220 <- full_join(epi.220, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.220$Sample.ID <- "EPI_220"
epi.220$merger <- paste0(epi.220$scaffold,"_", epi.220$position)
epi.220 <- merge(Data, epi.220, by="merger")
epi.220 <- epi.220[,6:14]
colnames(epi.220) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.221 <- read.csv("Data/Extracted/EPI-221.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.221 <- epi.221[,-3]
colnames(epi.221) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.221$position <- as.numeric(epi.221$position)
epi.221 <- full_join(epi.221, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.221$Sample.ID <- "EPI_221"
epi.221$merger <- paste0(epi.221$scaffold,"_", epi.221$position)
epi.221 <- merge(Data, epi.221, by="merger")
epi.221 <- epi.221[,6:14]
colnames(epi.221) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.226 <- read.csv("Data/Extracted/EPI-226.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.226 <- epi.226[,-3]
colnames(epi.226) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.226$position <- as.numeric(epi.226$position)
epi.226 <- full_join(epi.226, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.226$Sample.ID <- "EPI_226"
epi.226$merger <- paste0(epi.226$scaffold,"_", epi.226$position)
epi.226 <- merge(Data, epi.226, by="merger")
epi.226 <- epi.226[,6:14]
colnames(epi.226) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.227 <- read.csv("Data/Extracted/EPI-227.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.227 <- epi.227[,-3]
colnames(epi.227) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.227$position <- as.numeric(epi.227$position)
epi.227 <- full_join(epi.227, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.227$Sample.ID <- "EPI_227"
epi.227$merger <- paste0(epi.227$scaffold,"_", epi.227$position)
epi.227 <- merge(Data, epi.227, by="merger")
epi.227 <- epi.227[,6:14]
colnames(epi.227) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.229 <- read.csv("Data/Extracted/EPI-229.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.229 <- epi.229[,-3]
colnames(epi.229) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.229$position <- as.numeric(epi.229$position)
epi.229 <- full_join(epi.229, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.229$Sample.ID <- "EPI_229"
epi.229$merger <- paste0(epi.229$scaffold,"_", epi.229$position)
epi.229 <- merge(Data, epi.229, by="merger")
epi.229 <- epi.229[,6:14]
colnames(epi.229) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")

epi.230 <- read.csv("Data/Extracted/EPI-230.10x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
epi.230 <- epi.230[,-3]
colnames(epi.230) <- c("scaffold", "position","per.meth","meth","unmeth")
epi.230$position <- as.numeric(epi.230$position)
epi.230 <- full_join(epi.230, Genes, by="scaffold") %>% filter(position >= start & position <= stop)
epi.230$Sample.ID <- "EPI_230"
epi.230$merger <- paste0(epi.230$scaffold,"_", epi.230$position)
epi.230 <- merge(Data, epi.230, by="merger")
epi.230 <- epi.230[,6:14]
colnames(epi.230) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")


##### Statistical comparisons of methylation #####

# Compare Ambient Treatments through time #####
time.amb.meth <- rbind(epi.41,epi.42,epi.43,epi.44,
                       epi.119,epi.120,epi.135,epi.136,
                       epi.151,epi.152,epi.153,epi.154,
                       epi.181,epi.182,epi.184,epi.185)

colnames(time.amb.meth) <- c("scaffold", "position","per.meth","meth","unmeth", "start","stop", "gene", "Sample.ID")
timepoint.amb.meth <- merge(time.amb.meth,sample.info, by="Sample.ID")
time.amb.meth <- timepoint.amb.meth[,c(1,2,3,5,6,9,14)]
colnames(time.amb.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment")

sams <- aggregate(per.meth~Sample.ID, data=timepoint.amb.meth, FUN=mean)
sams <- merge(sams, sample.info, by="Sample.ID")
time.amb.means <- aggregate(per.meth~TimePoint, data=sams, FUN=mean)
time.amb.se <- aggregate(per.meth~TimePoint, data=sams, FUN=std.error)
time.amb.meth.per <- cbind(time.amb.means, time.amb.se$per.meth)
colnames(time.amb.meth.per) <- c("TimePoint", "mean","se")

Fig.time.amb.per.meth <- ggplot(time.amb.meth.per, aes(x=TimePoint, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_point(pch=16, size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Time") +
  ylab("Methylation Level") +
  #ylim(19.5,20.5)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none', 
        axis.title.y = element_text(size = rel(1.8), angle = 90),
        axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
  ggtitle("Ambient Through Time") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.time.amb.per.meth


#Compare treatments after 10 days of exposure  #####
D10.trt.meth <- rbind(epi.103,epi.104,epi.111,epi.113,
                      epi.119,epi.120,epi.127,epi.127,
                      epi.135,epi.136,epi.143,epi.145)
colnames(D10.trt.meth) <- c("scaffold", "position", "per.meth", "meth", "unmeth", "start", "stop", "gene", "Sample.ID")
D10.trt.meth <- merge(D10.trt.meth,sample.info, by="Sample.ID")
D10.trt.meth <- D10.trt.meth[,c(1,2,3,5,6,9,12)]
colnames(D10.trt.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment")

D10.trt.meth$per.meth <- (D10.trt.meth$meth/(D10.trt.meth$meth+D10.trt.meth$unmeth))*100
sams <- aggregate(per.meth~Sample.ID, data=D10.trt.meth, FUN=mean)
sams <- merge(sams, sample.info, by="Sample.ID")
D10.trt.means <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=mean)
D10.trt.se <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=std.error)
D10.meth.per <- cbind(D10.trt.means, D10.trt.se$per.meth)
colnames(D10.meth.per) <- c("Initial.Treatment", "mean","se")

#Exposure1 Plotting
Fig.D10.per.meth <- ggplot(D10.meth.per, aes(x=Initial.Treatment, y=mean, group=Initial.Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  #geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(pch=c(16,15,17), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Initial Treatment") +
  ylab("Methylation Level") +
  #ylim(19.5,20.5)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none', 
        axis.title.y = element_text(size = rel(1.8), angle = 90),
        axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
  ggtitle("Exposure 1 Day 10") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.D10.per.meth


#Compare treatments after common garden
D135.trt.meth <- rbind(epi.151,epi.152,epi.153,epi.154,
                       epi.159,epi.160,epi.161,epi.162,
                       epi.167,epi.168,epi.169,epi.170)
colnames(D135.trt.meth)  <- c("scaffold", "position", "per.meth", "meth", "unmeth", "start", "stop", "gene", "Sample.ID")
D135.trt.meth <- merge(D135.trt.meth,sample.info, by="Sample.ID")
D135.trt.meth <- D135.trt.meth[,c(1,2,3,5,6,9,12)]
colnames(D135.trt.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment")



D135.trt.meth$per.meth <- (D135.trt.meth$meth/(D135.trt.meth$meth+D135.trt.meth$unmeth))*100
sams <- aggregate(per.meth~Sample.ID, data=D135.trt.meth, FUN=mean)
sams <- merge(sams, sample.info, by="Sample.ID")
D135.trt.means <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=mean)
D135.trt.se <- aggregate(per.meth~Initial.Treatment, data=sams, FUN=std.error)
D135.meth.per <- cbind(D135.trt.means, D135.trt.se$per.meth)
colnames(D135.meth.per) <- c("Initial.Treatment", "mean","se")

#Common Garden Plotting
Fig.D135.per.meth <- ggplot(D135.meth.per, aes(x=Initial.Treatment, y=mean, group=Initial.Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  #geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(pch=c(16,15,17), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Initial Treatment") +
  ylab("Methylation Level") +
  #ylim(19.5,20.5)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none', 
        axis.title.y = element_text(size = rel(1.8), angle = 90),
        axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
  ggtitle("Post Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.D135.per.meth


#Compare secondary exposure
sec.meth <- rbind(epi.175,epi.176,epi.181,epi.182,
                  epi.184,epi.185,epi.187,epi.188,
                  epi.193,epi.194,epi.199,epi.200,
                  epi.205,epi.206,epi.208,epi.209,
                  epi.214,epi.215,epi.220,epi.221,
                  epi.226,epi.227,epi.229,epi.230)
colnames(sec.meth) <- c("scaffold", "position", "per.meth", "meth", "unmeth", "start", "stop", "gene", "Sample.ID")
sec.meth <- merge(sec.meth,sample.info, by="Sample.ID")
sec.meth <- sec.meth[,c(1,2,3,5,6,9,12,13)]
colnames(sec.meth) <-c("Sample.ID","scaffold", "position", "meth", "unmeth", "gene", "treatment1", "treatment2")

#####

sec.meth$per.meth <- (sec.meth$meth/(sec.meth$meth+sec.meth$unmeth))*100
sams <- aggregate(per.meth~Sample.ID, data=sec.meth, FUN=mean)
sams <- merge(sams, sample.info, by="Sample.ID")
D145.trt.means <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=mean)
D145.trt.se <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=std.error)
D145.meth.per <- cbind(D145.trt.means, D145.trt.se$per.meth)
colnames(D145.meth.per) <- c("Initial.Treatment", "Secondary.Treatment", "mean","se")
D145.meth.per$Initial.Treatment <- factor(D145.meth.per$Initial.Treatment,levels=c("Ambient","Super.Low", "Low"))

#Exposure2 Plotting
Fig.D145.per.meth <- ggplot(D145.meth.per, aes(x=Secondary.Treatment, y=mean, group=Initial.Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Initial.Treatment), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Secondary Treatment") +
  ylab("Methylation Level") +
  #ylim(19.5,20.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none', 
        axis.title.y = element_text(size = rel(1.8), angle = 90),
        axis.title.x = element_text(size = rel(1.8), angle = 00)) + #remove legend background
  ggtitle("Secondary Exposure") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.D145.per.meth




Meth.Level <- arrangeGrob(Fig.time.amb.per.meth,
                          Fig.D10.per.meth, Fig.D135.per.meth, Fig.D145.per.meth, ncol=2)
ggsave(file="Output/Methylation_Level.pdf", Meth.Level, width =6, height = 6, units = c("in"))

Time <- arrangeGrob(Fig.time.amb.per.meth, ncol=1)
ggsave(file="Output/Time.pdf", Time, width =6, height = 6, units = c("in"))

D10 <- arrangeGrob(Fig.D10.per.meth, ncol=1)
ggsave(file="Output/D10.trt.pdf", D10, width =6, height = 6, units = c("in"))

D135 <- arrangeGrob(Fig.D135.per.meth, ncol=1)
ggsave(file="Output/D135.trt.pdf", D135, width =6, height = 6, units = c("in"))

D145 <- arrangeGrob(Fig.D145.per.meth, ncol=1)
ggsave(file="Output/D145.trt.pdf", D145, width =6, height = 6, units = c("in"))

# Save data tables #
write.table(time.amb.meth, 'Output/Time.tsv', sep='\t', row.names=FALSE)
write.table(D10.trt.meth, 'Output/D10.tsv', sep='\t', row.names=FALSE)
write.table(D135.trt.meth, 'Output/D135.tsv', sep='\t', row.names=FALSE)
write.table(sec.meth, 'Output/sec.meth.tsv', sep='\t', row.names=FALSE)

##### One way tests #####

# TIME

meth_table  <- read.csv("Output/Time.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment == 'Day0' | meth_table$treatment == 'Day10' | meth_table$treatment == 'Day135' | meth_table$treatment == 'Day145', ]
#genes <- do_GLM_gene(sub_meth_table, min_cpg=2)
#write.table(genes, 'Output/Time_Liew.tsv', sep='\t', row.names=FALSE)


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
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment * position, 
             data=sub_meth_table1, family=binomial)
  #} else {
  #  # fit glm model
  #  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
  #             data=sub_meth_table1, family=binomial)
  #}
  #s <- step(fit, trace=0)
  a <- anova(fit, test="LRT")
  a
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,6],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   pval.position = a$`Pr(>Chi)`[3],
                   pval.treatment_x_position = a$`Pr(>Chi)`[4],
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.treatment <- p.adjust(results$pval.treatment, method='BH')
results$adj.pval.position <- p.adjust(results$pval.position, method='BH')
results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH')

write.table(results, 'Output/Time_HP_GLM.tsv', sep='\t', row.names=FALSE)

time.sig <- results[order(results$adj.pval.treatment),]
time.sig <- time.sig[ which(time.sig$adj.pval.treatment <0.05), ]
write.table(time.sig, 'Output/time_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)



# Day10

meth_table  <- read.csv("Output/D10.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment == 'Ambient' | meth_table$treatment == 'Low' | meth_table$treatment == 'Super.Low', ]
#genes <- do_GLM_gene(sub_meth_table, min_cpg=2)
#write.table(genes, 'Output/Time_Liew.tsv', sep='\t', row.names=FALSE)


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
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment * position, 
             data=sub_meth_table1, family=binomial)
  #} else {
  #  # fit glm model
  #  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
  #             data=sub_meth_table1, family=binomial)
  #}
  #s <- step(fit, trace=0)
  a <- anova(fit, test="LRT")
  a
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,6],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   pval.position = a$`Pr(>Chi)`[3],
                   pval.treatment_x_position = a$`Pr(>Chi)`[4],
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.treatment <- p.adjust(results$pval.treatment, method='BH')
results$adj.pval.position <- p.adjust(results$pval.position, method='BH')
results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH')

write.table(results, 'Output/D10_HP_GLM.tsv', sep='\t', row.names=FALSE)

D10.sig <- results[order(results$adj.pval.treatment),]
D10.sig <- D10.sig[ which(D10.sig$adj.pval.treatment <0.05), ]
write.table(D10.sig, 'Output/D10_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)


# Day135
meth_table  <- read.csv("Output/D135.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment == 'Ambient' | meth_table$treatment == 'Low' | meth_table$treatment == 'Super.Low', ]
#genes <- do_GLM_gene(sub_meth_table, min_cpg=2)
#write.table(genes, 'Output/Time_Liew.tsv', sep='\t', row.names=FALSE)


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
  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment * position, 
             data=sub_meth_table1, family=binomial)
  #} else {
  #  # fit glm model
  #  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
  #             data=sub_meth_table1, family=binomial)
  #}
  #s <- step(fit, trace=0)
  a <- anova(fit, test="LRT")
  a
  # capture summary stats to data frame
  df <- data.frame(gene = sub_meth_table1[1,6],
                   pval.treatment = a$`Pr(>Chi)`[2],
                   pval.position = a$`Pr(>Chi)`[3],
                   pval.treatment_x_position = a$`Pr(>Chi)`[4],
                   stringsAsFactors = F)
  
  # bind rows of temporary data frame to the results data frame
  results <- rbind(results, df)
  
}

results[is.na(results)] <- 0

results$adj.pval.treatment <- p.adjust(results$pval.treatment, method='BH')
results$adj.pval.position <- p.adjust(results$pval.position, method='BH')
results$adj.pval.treatment_x_position <- p.adjust(results$pval.treatment_x_position, method='BH')

write.table(results, 'Output/D135_HP_GLM.tsv', sep='\t', row.names=FALSE)

D135.sig <- results[order(results$adj.pval.treatment),]
D135.sig <- D135.sig[ which(D135.sig$adj.pval.treatment <0.05), ]
write.table(D135.sig, 'Output/D135_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)

##### Secondary Interaction Test #####

meth_table <- read.csv("Output/sec.meth.tsv", header=T, sep="\t", na.string="NA", stringsAsFactors = F)
sub_meth_table <- meth_table[meth_table$treatment1 == "Ambient" | meth_table$treatment1 == "Low" | meth_table$treatment1 == "Super.Low"
                             | meth_table$treatment2 == "Ambient" | meth_table$treatment2 == "Low", ]

# create data frame to store results
results <- data.frame()
gs <- unique(sub_meth_table$gene)

#first subset the unique dataframes and second run the GLMs
for(i in 1:length(sub_meth_table$gene)){
    
    #subset the dataframe gene by gene
    sub_meth_table1 <- subset(sub_meth_table, gene ==gs[i])
    
    #if(length(unique(sub_meth_table1$position)) >1){
      
    # fit glm position model
    fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2 * position, 
               data=sub_meth_table1, family=binomial)
    #} else {
    #  # fit glm model
    #  fit <- glm(matrix(c(meth, unmeth), ncol=2) ~ treatment1 * treatment2, 
    #             data=sub_meth_table1, family=binomial)
    #}
    #s <- step(fit, trace=0)
    a <- anova(fit, test="LRT")
    a
    # capture summary stats to data frame
    df <- data.frame(gene = sub_meth_table1[1,6],
                     pval.treatment1 = a$`Pr(>Chi)`[2],
                     pval.treatment2 = a$`Pr(>Chi)`[3],
                     position = a$`Pr(>Chi)`[4], 
                     pval.treatment1_x_treatment2 = a$`Pr(>Chi)`[5],
                     pval.treatment1_x_position = a$`Pr(>Chi)`[6],
                     pval.treatment2_x_position = a$`Pr(>Chi)`[7],
                     pval.treatment1_x_pval.treatment2_x_position = a$`Pr(>Chi)`[8],
                     stringsAsFactors = F)

    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
    
  }

results[is.na(results)] <- 0

results$adj.pval.treatment1 <- p.adjust(results$pval.treatment1, method='BH')
results$adj.pval.treatment2 <- p.adjust(results$pval.treatment2, method='BH')
results$adj.pval.treatment1_x_treatment2 <- p.adjust(results$pval.treatment1_x_treatment2, method='BH')
results$adj.pval.treatment1_x_position <- p.adjust(results$pval.treatment1_x_position, method='BH')
results$adj.pval.treatment2_x_position <- p.adjust(results$pval.treatment2_x_position, method='BH')
results$adj.pval.treatment1_x_pval.treatment2_x_position <- p.adjust(results$pval.treatment1_x_pval.treatment2_x_position, method='BH')

write.table(results, 'Output/sec.meth_HP_GLM.all.tsv', sep='\t', row.names=FALSE)

sec.data.sig.int <- results[order(results$adj.pval.treatment1_x_treatment2),]
sec.data.sig.int <- subset(sec.data.sig.int, by=adj.pval.treatment1_x_treatment2<0.05)
sec.data.sig.int <- sec.data.sig.int[ which(sec.data.sig.int$adj.pval.treatment1_x_treatment2 <0.05), ]
write.table(sec.data.sig.int, 'Output/sec.meth_HP_GLM_BH_0_05.tsv', sep='\t', row.names=FALSE)


for (i in 1:length(sec.data.sig.int$gene)){
inter.gene <- subset(sub_meth_table, gene ==sec.data.sig.int$gene[i])
means <- aggregate(per.meth ~ treatment1*treatment2*position, data=inter.gene, FUN=mean)

temp_plot <- ggplot(means, aes(x=treatment2, y=per.meth, fill=treatment1)) +
  geom_violin(aes(fill = treatment1), trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color="black")+
  scale_fill_manual(values = c("blue", "purple", "red"))+
  theme_bw()

ggsave(temp_plot, file=paste0("Output/plot_", i,".png"), width = 14, height = 10, units = "cm")
}


#heatmap of mean expression across all positions within a gene
inter.genes <- merge(sec.data.sig.int, sub_meth_table, by="gene")
all.means <- aggregate(per.meth ~ treatment1*treatment2*gene*Sample.ID, data=inter.genes, FUN=mean)

wide <- spread(all.means,  key=c("gene"), value=per.meth)
rownames(wide) <- wide$Sample.ID
annot <- wide[,1:2]
#annot$group <- paste0(annot$treatment1, "_", annot$treatment2)
mat <- wide[,c(4:length(wide))] #make an expression object
mat <- t(mat)
#mat <- mat - rowMeans(mat) #difference in expression compared to average across all samples
#ann_colors <- list(CO2 = c(Ambient="blue", Low="purple", Low="red")) #manually set colors
#jpeg(file="/Output/Heatmap.DMG.sec.meth.jpg") #save file
pheatmap(mat, annotation = annot,
         show_rownames =F, cluster_cols = FALSE,
         show_colnames =TRUE) #plot heatmap of all DEG by group
#dev.off()

