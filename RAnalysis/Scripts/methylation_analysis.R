#Title: Geoduck juvenile methylation
#Project: NOAA FFAR
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20180111
#See Readme file for details

rm(list=ls()) #clears workspace 


library(plotrix) 
library(ggplot2)
library(gridExtra)

# Set Working Directory:
setwd("~/MyProjects/Geoduck_Meth/RAnalysis/") #set working

sample.info <- read.csv("Data/Sample.Info.csv", header=T, sep=",", na.string="NA", stringsAsFactors = F)

Data <- read.csv("Data/Extracted/filtered.52.all.bed", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
Data <- Data[,1:3]
colnames(Data) <- c("scaffold", "start", "stop")
Data$location <- Data$start +1
Data$merger <- paste0(Data$scaffold,"_", Data$location)

##### Load in methylation counts #####
epi.41 <- read.csv("Data/Extracted/EPI-41.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.41) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.41$location <- as.numeric(epi.41$location)
epi.41$merger <- paste0(epi.41$scaffold,"_", epi.41$location)
epi.41 <- merge(Data, epi.41, by="merger")
EPI_41 <- epi.41[,c(1,5,10,11)]
EPI_41$Sample.ID <- "EPI_41"

epi.42 <- read.csv("Data/Extracted/EPI-42.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.42) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.42$location <- as.numeric(epi.42$location)
epi.42$merger <- paste0(epi.42$scaffold,"_", epi.42$location)
epi.42 <- merge(Data, epi.42, by="merger")
EPI_42 <- epi.42[,c(1,5,10,11)]
EPI_42$Sample.ID <- "EPI_42"

epi.43 <- read.csv("Data/Extracted/EPI-43.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.43) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.43$location <- as.numeric(epi.43$location)
epi.43$merger <- paste0(epi.43$scaffold,"_", epi.43$location)
epi.43 <- merge(Data, epi.43, by="merger")
EPI_43 <- epi.43[,c(1,5,10,11)]
EPI_43$Sample.ID <- "EPI_43"

epi.44 <- read.csv("Data/Extracted/EPI-44.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.44) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.44$location <- as.numeric(epi.44$location)
epi.44$merger <- paste0(epi.44$scaffold,"_", epi.44$location)
epi.44 <- merge(Data, epi.44, by="merger")
EPI_44 <- epi.44[,c(1,5,10,11)]
EPI_44$Sample.ID <- "EPI_44"

epi.103 <- read.csv("Data/Extracted/EPI-103.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.103) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.103$location <- as.numeric(epi.103$location)
epi.103$merger <- paste0(epi.103$scaffold,"_", epi.103$location)
epi.103 <- merge(Data, epi.103, by="merger")
EPI_103 <- epi.103[,c(1,5,10,11)]
EPI_103$Sample.ID <- "EPI_103"

epi.104 <- read.csv("Data/Extracted/EPI-104.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.104) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.104$location <- as.numeric(epi.104$location)
epi.104$merger <- paste0(epi.104$scaffold,"_", epi.104$location)
epi.104 <- merge(Data, epi.104, by="merger")
EPI_104 <- epi.104[,c(1,5,10,11)]
EPI_104$Sample.ID <- "EPI_104"

epi.111 <- read.csv("Data/Extracted/EPI-111.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.111) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.111$location <- as.numeric(epi.111$location)
epi.111$merger <- paste0(epi.111$scaffold,"_", epi.111$location)
epi.111 <- merge(Data, epi.111, by="merger")
EPI_111 <- epi.111[,c(1,5,10,11)]
EPI_111$Sample.ID <- "EPI_111"

epi.113 <- read.csv("Data/Extracted/EPI-113.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.113) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.113$location <- as.numeric(epi.113$location)
epi.113$merger <- paste0(epi.113$scaffold,"_", epi.113$location)
epi.113 <- merge(Data, epi.113, by="merger")
EPI_113 <- epi.113[,c(1,5,10,11)]
EPI_113$Sample.ID <- "EPI_113"

epi.119 <- read.csv("Data/Extracted/EPI-119.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.119) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.119$location <- as.numeric(epi.119$location)
epi.119$merger <- paste0(epi.119$scaffold,"_", epi.119$location)
epi.119 <- merge(Data, epi.119, by="merger")
EPI_119 <- epi.119[,c(1,5,10,11)]
EPI_119$Sample.ID <- "EPI_119"

epi.120 <- read.csv("Data/Extracted/EPI-120.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.120) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.120$location <- as.numeric(epi.120$location)
epi.120$merger <- paste0(epi.120$scaffold,"_", epi.120$location)
epi.120 <- merge(Data, epi.120, by="merger")
EPI_120 <- epi.120[,c(1,5,10,11)]
EPI_120$Sample.ID <- "EPI_120"

epi.127 <- read.csv("Data/Extracted/EPI-127.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.127) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.127$location <- as.numeric(epi.127$location)
epi.127$merger <- paste0(epi.127$scaffold,"_", epi.127$location)
epi.127 <- merge(Data, epi.127, by="merger")
EPI_127 <- epi.127[,c(1,5,10,11)]
EPI_127$Sample.ID <- "EPI_127"

epi.128 <- read.csv("Data/Extracted/EPI-128.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.128) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.128$location <- as.numeric(epi.128$location)
epi.128$merger <- paste0(epi.128$scaffold,"_", epi.128$location)
epi.128 <- merge(Data, epi.128, by="merger")
EPI_128 <- epi.128[,c(1,5,10,11)]
EPI_128$Sample.ID <- "EPI_128"

epi.135 <- read.csv("Data/Extracted/EPI-135.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.135) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.135$location <- as.numeric(epi.135$location)
epi.135$merger <- paste0(epi.135$scaffold,"_", epi.135$location)
epi.135 <- merge(Data, epi.135, by="merger")
EPI_135 <- epi.135[,c(1,5,10,11)]
EPI_135$Sample.ID <- "EPI_135"

epi.136 <- read.csv("Data/Extracted/EPI-136.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.136) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.136$location <- as.numeric(epi.136$location)
epi.136$merger <- paste0(epi.136$scaffold,"_", epi.136$location)
epi.136 <- merge(Data, epi.136, by="merger")
EPI_136 <- epi.136[,c(1,5,10,11)]
EPI_136$Sample.ID <- "EPI_136"

epi.143 <- read.csv("Data/Extracted/EPI-143.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.143) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.143$location <- as.numeric(epi.143$location)
epi.143$merger <- paste0(epi.143$scaffold,"_", epi.143$location)
epi.143 <- merge(Data, epi.143, by="merger")
EPI_143 <- epi.143[,c(1,5,10,11)]
EPI_143$Sample.ID <- "EPI_143"

epi.145 <- read.csv("Data/Extracted/EPI-145.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.145) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.145$location <- as.numeric(epi.145$location)
epi.145$merger <- paste0(epi.145$scaffold,"_", epi.145$location)
epi.145 <- merge(Data, epi.145, by="merger")
EPI_145 <- epi.145[,c(1,5,10,11)]
EPI_145$Sample.ID <- "EPI_145"

epi.151 <- read.csv("Data/Extracted/EPI-151.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.151) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.151$location <- as.numeric(epi.151$location)
epi.151$merger <- paste0(epi.151$scaffold,"_", epi.151$location)
epi.151 <- merge(Data, epi.151, by="merger")
EPI_151 <- epi.151[,c(1,5,10,11)]
EPI_151$Sample.ID <- "EPI_151"

epi.152 <- read.csv("Data/Extracted/EPI-152.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.152) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.152$location <- as.numeric(epi.152$location)
epi.152$merger <- paste0(epi.152$scaffold,"_", epi.152$location)
epi.152 <- merge(Data, epi.152, by="merger")
EPI_152 <- epi.152[,c(1,5,10,11)]
EPI_152$Sample.ID <- "EPI_152"

epi.153 <- read.csv("Data/Extracted/EPI-153.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.153) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.153$location <- as.numeric(epi.153$location)
epi.153$merger <- paste0(epi.153$scaffold,"_", epi.153$location)
epi.153 <- merge(Data, epi.153, by="merger")
EPI_153 <- epi.153[,c(1,5,10,11)]
EPI_153$Sample.ID <- "EPI_153"

epi.154 <- read.csv("Data/Extracted/EPI-154.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.154) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.154$location <- as.numeric(epi.154$location)
epi.154$merger <- paste0(epi.154$scaffold,"_", epi.154$location)
epi.154 <- merge(Data, epi.154, by="merger")
EPI_154 <- epi.154[,c(1,5,10,11)]
EPI_154$Sample.ID <- "EPI_154"

epi.159 <- read.csv("Data/Extracted/EPI-159.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.159) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.159$location <- as.numeric(epi.159$location)
epi.159$merger <- paste0(epi.159$scaffold,"_", epi.159$location)
epi.159 <- merge(Data, epi.159, by="merger")
EPI_159 <- epi.159[,c(1,5,10,11)]
EPI_159$Sample.ID <- "EPI_159"

epi.160 <- read.csv("Data/Extracted/EPI-160.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.160) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.160$location <- as.numeric(epi.160$location)
epi.160$merger <- paste0(epi.160$scaffold,"_", epi.160$location)
epi.160 <- merge(Data, epi.160, by="merger")
EPI_160 <- epi.160[,c(1,5,10,11)]
EPI_160$Sample.ID <- "EPI_160"

epi.161 <- read.csv("Data/Extracted/EPI-161.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.161) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.161$location <- as.numeric(epi.161$location)
epi.161$merger <- paste0(epi.161$scaffold,"_", epi.161$location)
epi.161 <- merge(Data, epi.161, by="merger")
EPI_161 <- epi.161[,c(1,5,10,11)]
EPI_161$Sample.ID <- "EPI_161"

epi.162 <- read.csv("Data/Extracted/EPI-162.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.162) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.162$location <- as.numeric(epi.162$location)
epi.162$merger <- paste0(epi.162$scaffold,"_", epi.162$location)
epi.162 <- merge(Data, epi.162, by="merger")
EPI_162 <- epi.162[,c(1,5,10,11)]
EPI_162$Sample.ID <- "EPI_162"

epi.167 <- read.csv("Data/Extracted/EPI-167.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.167) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.167$location <- as.numeric(epi.167$location)
epi.167$merger <- paste0(epi.167$scaffold,"_", epi.167$location)
epi.167 <- merge(Data, epi.167, by="merger")
EPI_167 <- epi.167[,c(1,5,10,11)]
EPI_167$Sample.ID <- "EPI_167"

epi.168 <- read.csv("Data/Extracted/EPI-168.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.168) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.168$location <- as.numeric(epi.168$location)
epi.168$merger <- paste0(epi.168$scaffold,"_", epi.168$location)
epi.168 <- merge(Data, epi.168, by="merger")
EPI_168 <- epi.168[,c(1,5,10,11)]
EPI_168$Sample.ID <- "EPI_168"

epi.169 <- read.csv("Data/Extracted/EPI-169.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.169) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.169$location <- as.numeric(epi.169$location)
epi.169$merger <- paste0(epi.169$scaffold,"_", epi.169$location)
epi.169 <- merge(Data, epi.169, by="merger")
EPI_169 <- epi.169[,c(1,5,10,11)]
EPI_169$Sample.ID <- "EPI_169"

epi.170 <- read.csv("Data/Extracted/EPI-170.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.170) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.170$location <- as.numeric(epi.170$location)
epi.170$merger <- paste0(epi.170$scaffold,"_", epi.170$location)
epi.170 <- merge(Data, epi.170, by="merger")
EPI_170 <- epi.170[,c(1,5,10,11)]
EPI_170$Sample.ID <- "EPI_170"

epi.175 <- read.csv("Data/Extracted/EPI-175.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.175) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.175$location <- as.numeric(epi.175$location)
epi.175$merger <- paste0(epi.175$scaffold,"_", epi.175$location)
epi.175 <- merge(Data, epi.175, by="merger")
EPI_175 <- epi.175[,c(1,5,10,11)]
EPI_175$Sample.ID <- "EPI_175"

epi.176 <- read.csv("Data/Extracted/EPI-176.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.176) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.176$location <- as.numeric(epi.176$location)
epi.176$merger <- paste0(epi.176$scaffold,"_", epi.176$location)
epi.176 <- merge(Data, epi.176, by="merger")
EPI_176 <- epi.176[,c(1,5,10,11)]
EPI_176$Sample.ID <- "EPI_176"

epi.181 <- read.csv("Data/Extracted/EPI-181.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.181) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.181$location <- as.numeric(epi.181$location)
epi.181$merger <- paste0(epi.181$scaffold,"_", epi.181$location)
epi.181 <- merge(Data, epi.181, by="merger")
EPI_181 <- epi.181[,c(1,5,10,11)]
EPI_181$Sample.ID <- "EPI_181"

epi.182 <- read.csv("Data/Extracted/EPI-182.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.182) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.182$location <- as.numeric(epi.182$location)
epi.182$merger <- paste0(epi.182$scaffold,"_", epi.182$location)
epi.182 <- merge(Data, epi.182, by="merger")
EPI_182 <- epi.182[,c(1,5,10,11)]
EPI_182$Sample.ID <- "EPI_182"

epi.184 <- read.csv("Data/Extracted/EPI-184.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.184) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.184$location <- as.numeric(epi.184$location)
epi.184$merger <- paste0(epi.184$scaffold,"_", epi.184$location)
epi.184 <- merge(Data, epi.184, by="merger")
EPI_184 <- epi.184[,c(1,5,10,11)]
EPI_184$Sample.ID <- "EPI_184"

epi.185 <- read.csv("Data/Extracted/EPI-185.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.185) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.185$location <- as.numeric(epi.185$location)
epi.185$merger <- paste0(epi.185$scaffold,"_", epi.185$location)
epi.185 <- merge(Data, epi.185, by="merger")
EPI_185 <- epi.185[,c(1,5,10,11)]
EPI_185$Sample.ID <- "EPI_185"

epi.187 <- read.csv("Data/Extracted/EPI-187.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.187) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.187$location <- as.numeric(epi.187$location)
epi.187$merger <- paste0(epi.187$scaffold,"_", epi.187$location)
epi.187 <- merge(Data, epi.187, by="merger")
EPI_187 <- epi.187[,c(1,5,10,11)]
EPI_187$Sample.ID <- "EPI_187"

epi.188 <- read.csv("Data/Extracted/EPI-188.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.188) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.188$location <- as.numeric(epi.188$location)
epi.188$merger <- paste0(epi.188$scaffold,"_", epi.188$location)
epi.188 <- merge(Data, epi.188, by="merger")
EPI_188 <- epi.188[,c(1,5,10,11)]
EPI_188$Sample.ID <- "EPI_188"

epi.193 <- read.csv("Data/Extracted/EPI-193.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.193) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.193$location <- as.numeric(epi.193$location)
epi.193$merger <- paste0(epi.193$scaffold,"_", epi.193$location)
epi.193 <- merge(Data, epi.193, by="merger")
EPI_193 <- epi.193[,c(1,5,10,11)]
EPI_193$Sample.ID <- "EPI_193"

epi.194 <- read.csv("Data/Extracted/EPI-194.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.194) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.194$location <- as.numeric(epi.194$location)
epi.194$merger <- paste0(epi.194$scaffold,"_", epi.194$location)
epi.194 <- merge(Data, epi.194, by="merger")
EPI_194 <- epi.194[,c(1,5,10,11)]
EPI_194$Sample.ID <- "EPI_194"

epi.199 <- read.csv("Data/Extracted/EPI-199.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.199) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.199$location <- as.numeric(epi.199$location)
epi.199$merger <- paste0(epi.199$scaffold,"_", epi.199$location)
epi.199 <- merge(Data, epi.199, by="merger")
EPI_199 <- epi.199[,c(1,5,10,11)]
EPI_199$Sample.ID <- "EPI_199"

epi.200 <- read.csv("Data/Extracted/EPI-200.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.200) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.200$location <- as.numeric(epi.200$location)
epi.200$merger <- paste0(epi.200$scaffold,"_", epi.200$location)
epi.200 <- merge(Data, epi.200, by="merger")
EPI_200 <- epi.200[,c(1,5,10,11)]
EPI_200$Sample.ID <- "EPI_200"

epi.205 <- read.csv("Data/Extracted/EPI-205.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.205) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.205$location <- as.numeric(epi.205$location)
epi.205$merger <- paste0(epi.205$scaffold,"_", epi.205$location)
epi.205 <- merge(Data, epi.205, by="merger")
EPI_205 <- epi.205[,c(1,5,10,11)]
EPI_205$Sample.ID <- "EPI_205"

epi.206 <- read.csv("Data/Extracted/EPI-206.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.206) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.206$location <- as.numeric(epi.206$location)
epi.206$merger <- paste0(epi.206$scaffold,"_", epi.206$location)
epi.206 <- merge(Data, epi.206, by="merger")
EPI_206 <- epi.206[,c(1,5,10,11)]
EPI_206$Sample.ID <- "EPI_206"

epi.208 <- read.csv("Data/Extracted/EPI-208.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.208) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.208$location <- as.numeric(epi.208$location)
epi.208$merger <- paste0(epi.208$scaffold,"_", epi.208$location)
epi.208 <- merge(Data, epi.208, by="merger")
EPI_208 <- epi.208[,c(1,5,10,11)]
EPI_208$Sample.ID <- "EPI_208"

epi.209 <- read.csv("Data/Extracted/EPI-209.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.209) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.209$location <- as.numeric(epi.209$location)
epi.209$merger <- paste0(epi.209$scaffold,"_", epi.209$location)
epi.209 <- merge(Data, epi.209, by="merger")
EPI_209 <- epi.209[,c(1,5,10,11)]
EPI_209$Sample.ID <- "EPI_209"

epi.214 <- read.csv("Data/Extracted/EPI-214.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.214) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.214$location <- as.numeric(epi.214$location)
epi.214$merger <- paste0(epi.214$scaffold,"_", epi.214$location)
epi.214 <- merge(Data, epi.214, by="merger")
EPI_214 <- epi.214[,c(1,5,10,11)]
EPI_214$Sample.ID <- "EPI_214"

epi.215 <- read.csv("Data/Extracted/EPI-215.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.215) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.215$location <- as.numeric(epi.215$location)
epi.215$merger <- paste0(epi.215$scaffold,"_", epi.215$location)
epi.215 <- merge(Data, epi.215, by="merger")
EPI_215 <- epi.215[,c(1,5,10,11)]
EPI_215$Sample.ID <- "EPI_215"

epi.220 <- read.csv("Data/Extracted/EPI-220.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.220) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.220$location <- as.numeric(epi.220$location)
epi.220$merger <- paste0(epi.220$scaffold,"_", epi.220$location)
epi.220 <- merge(Data, epi.220, by="merger")
EPI_220 <- epi.220[,c(1,5,10,11)]
EPI_220$Sample.ID <- "EPI_220"

epi.221 <- read.csv("Data/Extracted/EPI-221.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.221) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.221$location <- as.numeric(epi.221$location)
epi.221$merger <- paste0(epi.221$scaffold,"_", epi.221$location)
epi.221 <- merge(Data, epi.221, by="merger")
EPI_221 <- epi.221[,c(1,5,10,11)]
EPI_221$Sample.ID <- "EPI_221"

epi.226 <- read.csv("Data/Extracted/EPI-226.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.226) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.226$location <- as.numeric(epi.226$location)
epi.226$merger <- paste0(epi.226$scaffold,"_", epi.226$location)
epi.226 <- merge(Data, epi.226, by="merger")
EPI_226 <- epi.226[,c(1,5,10,11)]
EPI_226$Sample.ID <- "EPI_226"

epi.227 <- read.csv("Data/Extracted/EPI-227.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.227) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.227$location <- as.numeric(epi.227$location)
epi.227$merger <- paste0(epi.227$scaffold,"_", epi.227$location)
epi.227 <- merge(Data, epi.227, by="merger")
EPI_227 <- epi.227[,c(1,5,10,11)]
EPI_227$Sample.ID <- "EPI_227"

epi.229 <- read.csv("Data/Extracted/EPI-229.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.229) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.229$location <- as.numeric(epi.229$location)
epi.229$merger <- paste0(epi.229$scaffold,"_", epi.229$location)
epi.229 <- merge(Data, epi.229, by="merger")
EPI_229 <- epi.229[,c(1,5,10,11)]
EPI_229$Sample.ID <- "EPI_229"

epi.230 <- read.csv("Data/Extracted/EPI-230.5x.cov.CpG", header=F, sep="", na.string="NA", stringsAsFactors = F) #read in data file
colnames(epi.230) <- c("scaffold", "location", "stop","per.meth","c.meth","c.unmeth")
epi.230$location <- as.numeric(epi.230$location)
epi.230$merger <- paste0(epi.230$scaffold,"_", epi.230$location)
epi.230 <- merge(Data, epi.230, by="merger")
EPI_230 <- epi.230[,c(1,5,10,11)]
EPI_230$Sample.ID <- "EPI_230"


##### Statistical comparisons of methylation #####

# Compare Ambient Treatments through time
time.amb.meth <- rbind(EPI_41,EPI_42,EPI_43,EPI_44,
                       EPI_119,EPI_120,EPI_135,EPI_136,
                       EPI_151,EPI_152,EPI_153,EPI_154,
                       EPI_181,EPI_182,EPI_184,EPI_185)
colnames(time.amb.meth) <- c("scaffold", "location", "meth","unmeth", "Sample.ID")
time.amb.meth <- merge(time.amb.meth,sample.info, by="Sample.ID")
time.mod <- glm(cbind(meth, unmeth) ~ TimePoint, data=time.amb.meth, family = "binomial")
summary(time.mod)

time.amb.meth$per.meth <- (time.amb.meth$meth/(time.amb.meth$meth+time.amb.meth$unmeth))*100
sams <- aggregate(per.meth~Sample.ID, data=time.amb.meth, FUN=mean)
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
  ylim(20.7,21.7)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("Ambient Through Time") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.time.amb.per.meth


#Compare treatments after 10 days of exposure
D10.trt.meth <- rbind(EPI_103,EPI_104,EPI_111,EPI_113,
                       EPI_119,EPI_120,EPI_127,EPI_127,
                       EPI_135,EPI_136,EPI_143,EPI_145)
colnames(D10.trt.meth) <- c("scaffold", "location", "meth","unmeth", "Sample.ID")
D10.trt.meth <- merge(D10.trt.meth,sample.info, by="Sample.ID")
D10.trt.mod <- glm(cbind(meth, unmeth) ~ Initial.Treatment, data=D10.trt.meth, family = "binomial")
summary(D10.trt.mod)

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
  geom_point(aes(shape=Initial.Treatment), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Initial Treatment") +
  ylab("Methylation Level") +
  ylim(20.7,21.7)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("Exposure 1 Day 10") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.D10.per.meth


#Compare treatments after common garden
D135.trt.meth <- rbind(EPI_151,EPI_152,EPI_153,EPI_154,
                      EPI_159,EPI_160,EPI_161,EPI_162,
                      EPI_167,EPI_168,EPI_169,EPI_170)
colnames(D135.trt.meth) <- c("scaffold", "location", "meth","unmeth", "Sample.ID")
D135.trt.meth <- merge(D135.trt.meth,sample.info, by="Sample.ID")
D135.trt.mod <- glm(cbind(meth, unmeth) ~ Initial.Treatment, data=D135.trt.meth, family = "binomial")
summary(D135.trt.mod)


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
  geom_point(aes(shape=Initial.Treatment), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Initial Treatment") +
  ylab("Methylation Level") +
  ylim(20.7,21.7)  +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("Post Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.D135.per.meth


#Compare secondary exposure
sec.meth <- rbind(EPI_175,EPI_176,EPI_181,EPI_182,
                  EPI_184,EPI_185,EPI_187,EPI_188,
                  EPI_193,EPI_194,EPI_199,EPI_200,
                  EPI_205,EPI_206,EPI_208,EPI_209,
                  EPI_214,EPI_215,EPI_220,EPI_221,
                  EPI_226,EPI_227,EPI_229,EPI_230)
colnames(sec.meth) <- c("scaffold", "location", "meth","unmeth", "Sample.ID")
sec.meth <- merge(sec.meth,sample.info, by="Sample.ID")
sec.meth.mod <- glm(cbind(meth, unmeth) ~ as.factor(Initial.Treatment)*as.factor(Secondary.Treatment), data=sec.meth, family = "binomial")
summary(sec.meth.mod)

sec.meth$per.meth <- (sec.meth$meth/(sec.meth$meth+sec.meth$unmeth))*100
sams <- aggregate(per.meth~Sample.ID, data=sec.meth, FUN=mean)
sams <- merge(sams, sample.info, by="Sample.ID")
D145.trt.means <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=mean)
D145.trt.se <- aggregate(per.meth~Initial.Treatment*Secondary.Treatment, data=sams, FUN=std.error)
D145.meth.per <- cbind(D145.trt.means, D145.trt.se$per.meth)
colnames(D145.meth.per) <- c("Initial.Treatment", "Secondary.Treatment", "mean","se")

#Exposure2 Plotting
Fig.D145.per.meth <- ggplot(D145.meth.per, aes(x=Secondary.Treatment, y=mean, group=Initial.Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Initial.Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Initial.Treatment), size = 3, position = position_dodge(width = 0.05)) +
  #annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  xlab("Secondary Treatment") +
  ylab("Methylation Level") +
  ylim(20.7,21.7) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
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



# meth <- rbind(EPI_41,EPI_42,EPI_43,EPI_44,EPI_103,
#               EPI_104,EPI_111,EPI_113,EPI_119,EPI_120,
#               EPI_127,EPI_128,EPI_135,EPI_136,EPI_143,
#               EPI_145,EPI_151,EPI_152,EPI_153,EPI_154,
#               EPI_159,EPI_160,EPI_161,EPI_162,EPI_167,
#               EPI_168,EPI_169,EPI_170,EPI_175,EPI_176,
#               EPI_181,EPI_182,EPI_184,EPI_185,EPI_187,
#               EPI_188,EPI_193,EPI_194,EPI_199,EPI_200,
#               EPI_205,EPI_206,EPI_208,EPI_209,EPI_214,
#               EPI_215,EPI_220,EPI_221,EPI_226,EPI_227,
#               EPI_229,EPI_230)
# 
# 
# glm(methylated, non methylated ~ pH, data=data, family = "binomial")


