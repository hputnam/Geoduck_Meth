# Geoduck_Meth

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Issues](https://github.com/hputnam/Geoduck_Meth/issues)
- [Citation](#citation)

## Overview
This repository includes data and analysis scripts to accompany:

**Title:** Environmental “memory” through dynamic DNA methylation provides acclimatization of geoduck clams to seawater pH

**Authors:** Putnam HM*, Trigg SA, White SJ, Spencer L,  Vadopalas B, Natarajan A, Hetzel J, Jaeger E, Soohoo J, Gallardo C, Goetz FW,  Roberts S

**Abstract:** Environmental conditioning in aquaculture occurs through the husbandry and handling of organisms throughout breeding and rearing. This approach is often non-targeted, but could be an efficient and effective avenue to generate gains in production output, if the mechanisms and transferability of environmental conditioning practices are thoroughly investigated. Here we generated a draft genome (XMB) of the economically significant aquaculture geoduck clam, Panopea generosa, as a community resource and genomic foundation for mechanistic studies of environmental conditioning, such as epigenetic priming.  We tested the sensitivity of early life stages and the potential for geoduck clams to display environmental hardening, or acclimatization to low pH through epigenetic “memory” a series of  repeat exposure experiments (23 days of exposure to either pH 7.9, 7.4, or 7.0, followed by 112 days of common garden at ambient pH, and 10 days of reciprocal exposure to either pH 7.9 or 7.4) and assessed physiological and epigenetic response. In ambient conditions geoduck grew through time and displayed differential methylation with age at the region (DMR) and gene (DMGs) levels. In comparison to ambient conditions (pH 7.9), juvenile geoduck displayed decreased shell size in the two low pH conditions (pH 7.4 and pH 7.0) and significant DMRs and DMGs between treatments. When replaced in the ambient ambient common garden conditions for several months, the initial exposure to low pH resulted in compensatory growth, such that the juveniles grew larger in the two low pH treatments compared to those initially in the ambient pH. Further, DMRs and DMGs corresponding to initial treatment were present after 112 days of common garden, indicating a molecular “memory” of initial exposure.  Re-exposure to ambient pH (7.9) and low pH (7.4) for another 10 days revealed the pre-exposed juveniles (initially exposed to pH 7.4 or 7.0) were also more resistant to low pH (pH 7.4) in the second exposure, in comparison to the naive, ambient (pH 7.9) initial treatment clams, as shown by a smaller change in shell size. Further the mean methylation change across the 434 DMGs with a significant interaction between initial and secondary treatments was 4.4x higher in naive clams, relative to conditioned clams. Genes with functions of X, Y, Z showed significant interaction of initial by secondary treatments, or intragenerational acclimatization modulated by DNA methylation. This suggests that acclimatization to OA can result in benefits to geoduck growth at a critical life stage, with exposure memory that is epigenetic in nature, thereby providing a mechanism for environmental hardening in aquaculture for industry benefit.

## Repo Contents: 

- ### [code](https://github.com/hputnam/Geoduck_Meth/tree/master/code):
	- [Geoduck\_Meth\_Pipeline.md](https://github.com/hputnam/Geoduck_Meth/blob/master/code/Geoduck_Meth_Pipeline.md): Bioinformatic analysis pipeline for processing raw bisulfite sequencing reads through DMR and DMG analysis. 
 
- ### [RAnalysis](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis): 
			
	- **[Data](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Data):**
		- [Avtech\_pH\_data.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Avtech_pH_data.csv): pH logger data
		- [Avtech\_temp\_data.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Avtech_temp_data.csv): temperature logger data
		- [Batch114.pdf](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Batch114.pdf): Batch information for CO2 in seawater reference material used in total alkalinity measurements
		- [Cell\_Counts\_Juvenile\_Geoduck\_Exp1.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Cell_Counts_Juvenile_Geoduck_Exp1.csv):  
		- [Cell\_Counts\_Juvenile\_Geoduck\_Exp2.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Cell_Counts_Juvenile_Geoduck_Exp2.csv): 
		- [Cumulative_TA_Output.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Cumulative_TA_Output.csv): Output from total alkalinity analysis
		- [Flow1\_Juvenile\_Geoduck.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Flow1_Juvenile_Geoduck.csv):  Flow rates for initial low pH exposure
		- [Flow2\_Juvenile\_Geoduck.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Flow2_Juvenile_Geoduck.csv):  Flow rates for secondary low pH exposure 
		- [Hobo\_Temperature\_Juvenile\_Geoduck_Exp2.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Hobo_Temperature_Juvenile_Geoduck_Exp2.csv):
		- [Hobo\_Temperature\_Juvenile\_Geoduck_OCG.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Hobo_Temperature_Juvenile_Geoduck_OCG.csv):  temperature logger data from common garden rearing
		- [pH\_Calibration_Files](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Data/pH_Calibration_Files):  14 files containing temperature and pH (mV) measurements of Tris standard used for pH probe calibration. Files are names by the date that measurements were recorded.
		- [Sample.Info.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Sample.Info.csv):  Sample meta data for bisulfite sequencing libraries
		- [Size\_Juvenile\_Geoduck.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Size_Juvenile_Geoduck.csv):  Size measurements of juvenile geoducks 
		- [SW\_Chem\_Juvenile\_Geoduck.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/SW_Chem_Juvenile_Geoduck.csv):  Table of seawater chemistry measurements (temperature, pH (mV), Salinity, and total alkalinity) from initial and secondary pH exposures. 
		- [TA_CRMs.csv](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/TA_CRMs.csv):  total alkalinity measurements for CO2 in seawater reference material
		- **CONFIRM THESE FILES NEED TO LIVE HERE:**
			- [Genome](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Data/Genome): Are these files are also on OSF? do they need to be in this repo?
				- [Panopea-generosa-genes-annotations.tab](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Genome/Panopea-generosa-genes-annotations.tab):
				- [Panopea-generosa-v1.0.a4.exon.gff3](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Data/Genome/Panopea-generosa-v1.0.a4.exon.gff3):
								
	- **[Figs](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs):** This directory contains main and supplemental manuscript figures
		- [Fig.4.jpg](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig.4.jpg):  Heatmaps showing relative percent methylation of significantly differentially methylated regions from comparison of (A) day 10 samples, (B) day 135 samples, and (C) day 145 samples.
		- [Fig.5.jpg](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig.5.jpg): DMR Features broken down by the percent of DMRs in each of the genomic regions indicated. Asterisks indicate a significant P value from a chi square test of proportions after FDR correction. 
		- [Fig.S1.pdf](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig.S1.pdf): plots of pH and temperature logger data from initial and secondary low pH exposure 
		- [Fig.S2.jpg](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig.S2.jpg):  Size and genome location of DMRs. A) Violin plots showing the distribution of DMR lengths (in base pairs) across different sample comparisons (day 10, all samples from day 10; day 135, all samples from day 135; day 145, all samples from day 145;, time, all ambient samples over time). Number of DMRs is listed across the top of the plot. B) Distribution of DMRs from different comparisons across genome scaffolds. Y axis is the proportion of DMRs in each scaffold out of the total number of DMRs in each comparison expressed as a percentage. C) Relative scaffold position of DMRs from different comparisons across genome scaffolds. 
		- [Fig.S3.jpg](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig.S3.jpg):  Characteristics of differentially methylated regions among animals reared in ambient conditions over time. (A) Heatmap showing relative percent methylation of significantly differentially methylated regions. (B) Distribution of DMRs across genomic features. Asterisks indicate a significant P value from a chi square test of proportions after FDR correction.
		- [Fig2_ExpDesign.pdf](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig2_ExpDesign.pdf): Figure 2. Experimental design.
		- [Fig2_ExpDesign.png](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Fig2_ExpDesign.png): Experimental design. **DO WE NEED BOTH exp. design figs?**
		- [FigX\_DMG\_heatmaps_treatmentonly.pdf](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/FigX_DMG_heatmaps_treatmentonly.pdf): Figure 6. Heatmaps displaying the mean of the differentially methylated genes by A) initial treatment, B) treatment history while in the common garden, and C) by initial and secondary treatments.  Numbers above the paired columns in C) represent the average change in methylation between the pairs.
		- [Geoduck_Size.png](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Figs/Geoduck_Size.png): Figure 3. Geoduck shell size (mean ± sem) relative to the mean size of that treatment at that time point for A) the initial exposure, B) the ambient common garden period, and C) the secondary reciprocal exposure. Results of posthoc analyses are indicated by dissimilar letters for significantly different treatments.
		
	- **[Output](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Output):** The directory containing the analysis results, figures, tables, and supplementary material.
		- **THIS DIRECTORY NEEDS AN OVERHAUL. Can any files be consolidated and/or removed if not needed?** 
	
	- **[Scripts](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Scripts):** R scripts for all statistical analysis
		- [CarbChem.R](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/CarbChem.R):  
		- [Circos_Plot.R](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/Circos_Plot.R):
		- [DMR\_stats\_and\_anno.Rmd](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/DMR_stats_and_anno.Rmd): Script used to identify DMRs across treatment groups by ANOVA, generate DMR size and genomic location distribution plots, and heatmaps. Performs chi square test for genomic feature analysis and generates stacked bar plots. 
		- [GM.Rmd](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/GM.Rmd):
		- [Seed_Exp.R](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/Seed_Exp.R):
		- **CONFIRM THESE ARE USED OR NEEDED IN MANUSCRIPT:**
			- [Characterizing-CpG-Methylation.ipynb](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/Characterizing-CpG-Methylation.ipynb):
			- [index.html](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/index.html):
			- [methyl_island_sliding_window.pl](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/methyl_island_sliding_window.pl)
			- [.ipynb_checkpoints](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Scripts/.ipynb_checkpoints)
			- [.Rhistory](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/.Rhistory)


- [Geoduck_Meth.RProj](https://github.com/hputnam/Geoduck_Meth/blob/master/Geoduck_Meth.Rproj):  **Is this RProj needed?**

## Citation
Putnam HM*, Trigg SA, White SJ, Spencer L,  Vadopalas B, Natarajan A, Hetzel J, Jaeger E, Soohoo J, Gallardo C, Goetz FW,  Roberts S. (2020) Environmental “memory” through dynamic DNA methylation provides acclimatization of geoduck clams to seawater pH.