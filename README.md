# Geoduck_Meth

## Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Issues](https://github.com/hputnam/Geoduck_Meth/issues)
- [Citation](#citation)

## Overview
This repository includes data and analysis scripts to accompany:

**Title:** Dynamic DNA methylation contributes to carryover effects and beneficial acclimatization in geoduck clams

**Authors:** Hollie M. Putnam1, Shelly A. Trigg, Samuel J. White2, Laura H. Spencer2, Brent Vadopalas2,3, Aparna Natarajan4, Jonathan Hetzel4, Erich Jaeger4, Jonathan Soohoo4, Cristian Gallardo-Escárate5,6, Frederick W. Goetz7, Steven B. Roberts2

**Abstract:** To investigate the role of DNA methylation in environmental acclimatization in the long-lived and commercially important geoduck clam Panopea generosa, we generated the first draft reference genome and surveyed the physiology and DNA methylomes of juveniles initially exposed to environments differing in seawater pH and then following depuration and reexposure. Juveniles were initially exposed to one of three conditions (pH 7.9, 7.4, or 7.0) for 23 days, followed by an ambient common garden (~4 months), then a second reciprocal exposure to pH 7.9 or pH 7.4 for 23 days. Within 10 days of the initial exposure to pH 7.4 and 7.0, juvenile clams showed decreased shell size relative to ambient pH (7.9), which was associated with differentially methylated regions (DMRs) and genes (DMGs). After several months under ambient common-garden conditions, juveniles initially exposed to low pH (pH 7.4 and 7.0) grew larger than those initially exposed to ambient. DMRs and DMGs reflected these phenotypic differences, demonstrating epigenetic carryover effects persisted months after initial exposure. Functional enrichment analysis of the DMGs revealed regulation of signal transduction through widespread changes in the Wnt signaling pathways that influence cell growth, proliferation, tissue and skeletal formation, and cytoskeletal change. Correspondingly, developmental processes, protein metabolism, and transport were also enriched within the DMGs. Upon secondary exposure to pH 7.4, naive juvenile clams (pH 7.9) were more sensitive to low pH, with more reduced growth compared to those initially exposed to lower pH (pH 7.4 and 7.0). Similarly, naive clam methylation had a nearly 2-fold greater magnitude of methylation change in differentially methylated genes after 10 days of secondary exposure. Collectively, this new genomic resource and coordinated phenotypic and methylomic response support that environmental hardening via epigenetic mechanisms can provide beneficial phenotypes in molluscs, with important implications for wild and aquaculture populations under climate change.


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
		- [Supplemental_Tables](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Output/Supplemental_Tables)
	
	- **[Scripts](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Scripts):** R scripts for all statistical analysis
		- [Circos_Plot.R](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/Circos_Plot.R):
		- [DMR\_stats\_and\_anno.Rmd](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/DMR_stats_and_anno.Rmd): Script used to identify DMRs across treatment groups by ANOVA, generate DMR size and genomic location distribution plots, and heatmaps. Performs chi square test for genomic feature analysis and generates stacked bar plots. 
		- [GM.Rmd](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/GM.Rmd): Script to perform differentially methylated gene analysis.
		- **CONFIRM THESE ARE USED OR NEEDED IN MANUSCRIPT:**
			- [Characterizing-CpG-Methylation.ipynb](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/Characterizing-CpG-Methylation.ipynb):
			- [index.html](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/index.html):
			- [methyl_island_sliding_window.pl](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/methyl_island_sliding_window.pl)
			- [.ipynb_checkpoints](https://github.com/hputnam/Geoduck_Meth/tree/master/RAnalysis/Scripts/.ipynb_checkpoints)
			- [.Rhistory](https://github.com/hputnam/Geoduck_Meth/blob/master/RAnalysis/Scripts/.Rhistory)


- [Geoduck_Meth.RProj](https://github.com/hputnam/Geoduck_Meth/blob/master/Geoduck_Meth.Rproj):  **Is this RProj needed?**

