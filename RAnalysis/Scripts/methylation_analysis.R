#Title: Geoduck juvenile methylation
#Project: NOAA FFAR
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20191120
#See Readme file for details

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

rm(list=ls()) #clears workspace 

#load sample information
sample.info <- read.csv("RAnalysis/Data/Sample.Info.csv", header=T, sep=",", na.string="NA", stringsAsFactors = F) #read in info
samp <- sample.info$Sample.ID # set sample info
samp <- gsub("[_]", "-", samp) #remove extra characters

#load genes gff
Genes <- read.csv("RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.gene.gff3", header=F, sep="\t", na.string="NA", skip=3, stringsAsFactors = F) #read in data file
Genes <- Genes[,c(1,4,5,9)] #select desired columns only
colnames(Genes) <- c("scaffold", "start", "stop", "gene") #rename columns
Genes$gene <- gsub(";.*","",Genes$gene) #remove extra characters
Genes$gene <- gsub("ID=","",Genes$gene) #remove extra characters

#Load annotation file
Annot <- read.csv("RAnalysis/Data/Genome/Panopea-generosa-vv0.74.a4.5d951a9b74287-blast_functional.tab", header=TRUE, sep="\t", na.string="NA", stringsAsFactors = F, skip=12)
colnames(Annot) <- c("gene","start", "end", "score", "uniprot.id", "Match.ID", "m_start", "m_end", "E.value", "identity", "al")
Annot$gene <- gsub("21910-","",Annot$gene) #remove extra characters
Annot$gene  <- gsub(".m01","",Annot$gene) #remove extra characters
Uniprots <- Annot[,c(1,5)] #pull out uniprot information
colnames(Uniprots) <- c("gene", "uniprot.id") #rename columns
write.csv(Uniprots, "RAnalysis/Data/Genome/Uniprot.IDS.csv") #save to file
#GOs <- read.csv("RAnalysis/Data/Genome/Uniprot_GOs.tab", header=TRUE, sep="\t", na.string="NA", stringsAsFactors = F) #read in gene ontology information
#colnames(GOs) <- c("uniprot.id", "protein.name.short", "status", "protein.name.long", "gene.names", "organism", "length", "GO.IDs") #rename columns of gene ontology information

###### Load in methylation counts #####
meth.data <- list.files(path = "RAnalysis/Data/OSF/bedgraphs/10x", pattern = ".bedgraph$", full.names=TRUE) %>%
  set_names(.) %>% 
  map_dfr(read.csv,.id="Sample.ID", header=F, sep="", na.string="NA", stringsAsFactors = F) %>% 
  dplyr::select(-V3) %>%
  group_by(Sample.ID)
colnames(meth.data) <- c("Sample.ID", "scaffold", "position","per.meth")

meth.data$Sample.ID <- gsub("RAnalysis/Data/OSF/bedgraphs/10x/","",meth.data$Sample.ID) #remove extra characters
meth.data$Sample.ID <- gsub("_.*","",meth.data$Sample.ID) #remove extra characters 
meth.data$Sample.ID <- gsub("-","_",meth.data$Sample.ID) #remove extra characters

MD <- merge(meth.data,sample.info, by="Sample.ID")

Time <- subset(MD, Initial.Treatment == "Ambient" & Secondary.Treatment == "Ambient") 
D10 <-  subset(MD, TimePoint == "Day10" ) 
D135 <-  subset(MD, TimePoint == "Day135" ) 
D145 <-  subset(MD, TimePoint == "Day145" ) 

#not enough memory on laptop for this full join
meth.data.x <- full_join(D135, Genes, by="scaffold") %>% filter(position >= start & position <= stop)

# Save data tables #####
write.table(Time, 'Output/Time.Amb.tsv', sep='\t', row.names=FALSE)
write.table(D10, 'Output/D10.tsv', sep='\t', row.names=FALSE)
write.table(D135, 'Output/D135.tsv', sep='\t', row.names=FALSE)
write.table(D145, 'Output/D145.tsv', sep='\t', row.names=FALSE)


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



