# Capstone Project

## Microarray Data

### Materials & Methods

One dataset was chosen by setting the screening criteria for the species type as "Homo sapiens", from GEO database of National Center for Biotechnology Information (NCBI) ([https://www.ncbi.nlm.nih.gov/geo/),](https://www.ncbi.nlm.nih.gov/geo/),) with such keywords being searched as "**Renal Cell Carcinoma**", "miRNA", "mRNA". Then study type "Expression Profiling by Array" was selected. As a result, we included the dataset **GSE16441**, which is divided into two platforms. Platform GPL6480 (Agilent-014850 Whole Human Genome Microarray 4x44K G4112F), Samples GSM413237-GSM413270, contains the total RNA expression data from 17 RCC tumor samples and 17 corresponding non-tumor samples and was used to find differentially expressed lncRNA and mRNA. Platform GPL8659 (Agilent Human miRNA Microarray Rel12.0) Samples GSM413271-GSM413304, contains the miRNA microarray expression data from 17 RCC tumors and 17 corresponding non-tumor samples and was used to find differentially expressed microRNAs.

## DEGs Analysis

### Materials & Methods

EGs were analyzed by R package "Linear Models for Microarray Data (**limma**)" function for datasets and "**GEOquery**" R package was used to retrieve to GSE from GEO database. Let me note that the data was already LOWESS normalized and log2 transformed as indicated in the matrix file. For the GPL6480 platform (the total RNA one), **log Foldchange\>2 and adjust. p\<0.05** were regarded as threshold values for selecting DEGs and DELs. For the GPL8659 platform samples (the micro-RNA one), l**og Foldchange\>1 and adjust. p\<0.05** were regarded as threshold values for selecting DEGs and DElncRNAs. Statistical significance for the selection of this threshold was found, and those genes that were up- and down- regulated can also be selected for performing the subsequent analysis. R software was also used to draw volcano map of DELs, DEmiRs and DEmRNAs. This data would be used in the following ceRNA network construction and protein interaction network construction.

### Results

According to the cut-off criteria mentioned above, volcano plots were obtained. We finally found that in the samples of total RNA, 880 were differentially expressed, having 213 upregulated and 667 downregulated. Among these total RNA, we have 6 DElncRNAs, 5 downregulated and 1 upregulated. Also, in the miRNA sample we found 87 differentially expressed miRNAs, 55 being upregulated and 32 downregulated.

### Figures

![Figure 1: Volcano Plot for the DEmRNAs samples; Blue points are the ones selected](figures/figure1.png)

![Figure 2: Box Plot for the log2FC expression levels in the samples for the DEmRNAs samples](figures/figure2.png)

![Figure 3: Volcano Plot for the DEmiRs samples; Blue points are the ones selected](figures/figure3.png)

![Figure 4: Box Plot for the log2FC expression levels in the samples for the DEmiRs samples](figures/figure4.png)
