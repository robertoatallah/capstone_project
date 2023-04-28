---
---
---

# Capstone Project - Roberto Atallah

# Integrated analysis of lncRNA-miRNA-mRNA ceRNA network in human renal cell carcinoma

## Introduction

Renal cell carcinoma (RCC) is the most common type of kidney cancer, accounting for approximately 90% of all kidney cancer cases (Hsieh et al., 2017). Despite advances in diagnosis and treatment, RCC still poses a significant health threat, with a 5-year survival rate of less than 10% for advanced stage disease (Hsieh et al., 2017). In recent years, non-coding RNAs (ncRNAs) have emerged as important regulators of gene expression and have been implicated in the pathogenesis of various cancers, including RCC (Slack & Chinnaiyan, 2019). Among ncRNAs, long non-coding RNAs (lncRNAs) have garnered particular interest due to their involvement in numerous biological processes, such as cell differentiation, proliferation, and apoptosis (Kung et al., 2013). Furthermore, lncRNAs are the ncRNAs that are over 200 bp in length and have been shown to play a critical role in the development and progression of RCC, thus they can be targeted for the treatment of the cancer (Haiyan et al., 2021).

Recently, competing endogenous RNA (ceRNA) networks have been proposed as a novel regulatory mechanism involving lncRNAs, microRNAs (miRNAs), and messenger RNAs (mRNAs) (Ala, 2020). In ceRNA networks, lncRNAs can act as sponges for miRNAs, thereby regulating the expression of downstream target mRNAs. Dysregulation of ceRNA networks has been implicated in various diseases, including cancer, like for example in the research conducted by (Shen et al., 2022) where the dysfunction of the ceRNA network of the tumor suppressor miR-637 led to clear cell renal cell carcinoma (ccRCC) development and prognosis.

In this study, I aimed to comprehensively analyze the lncRNA-miRNA-mRNA ceRNA network in human RCC, by using bioinformatics analyses of public datasets. Firstly, this study applied Gene Expression Omnibus (GEO) database for obtaining the RCC-related lncRNA, mRNA and miRNA expression data. Therefore the DEGs (Differentially expressed genes), DEMs (Differentially expressed miRNAs) and DELs (Differentially expressed lncRNAs) were obtained. Then I used the DEGs to perform functional enrichment analysis to identify key biological pathways and processes involved in RCC. Cytoscape was further utilized to construct a lncRNA-miRNA-mRNA network, followed by the construction of a PPI network to finally come up with key genes expressed differently in RCC. The findings may provide insights into the molecular mechanisms underlying RCC and potential therapeutic targets for this deadly disease.

This is the flow diagram of the work that has been done.

![](figures/Figure13.png)

## Microarray Data

### Materials & Methods

One dataset was chosen by setting the screening criteria for the species type as "Homo sapiens", from GEO database of National Center for Biotechnology Information (NCBI) ([https://www.ncbi.nlm.nih.gov/geo/),](https://www.ncbi.nlm.nih.gov/geo/),) with such keywords being searched as "**Renal Cell Carcinoma**", "miRNA", "mRNA". Then study type "Expression Profiling by Array" was selected. As a result, I included the dataset **GSE16441**, which is divided into two platforms. Platform GPL6480 (Agilent-014850 Whole Human Genome Microarray 4x44K G4112F), Samples GSM413237-GSM413270, contains the total RNA expression data from 17 RCC tumor samples and 17 corresponding non-tumor samples and was used to find differentially expressed lncRNA and mRNA. Platform GPL8659 (Agilent Human miRNA Microarray Rel12.0) Samples GSM413271-GSM413304, contains the miRNA microarray expression data from 17 RCC tumors and 17 corresponding non-tumor samples and was used to find differentially expressed microRNAs.

## DEGs Analysis

### Materials & Methods

DEGs are analyzed by R package "Linear Models for Microarray Data (**limma**)" (Ritchie et al., 2015) function for datasets and "**GEOquery**" R package was used to retrieve to GSE from GEO database. Let me note that the data was already LOISS normalized and log2 transformed as indicated in the matrix file. For the GPL6480 platform (the total RNA one), **log Foldchange\>2 and adjust. p\<0.05** Ire regarded as threshold values for selecting DEGs and DELs. For the GPL8659 platform samples (the micro-RNA one), l**og Foldchange\>1 and adjust. p\<0.05** Ire regarded as threshold values for selecting DEGs and DElncRNAs. Statistical significance for the selection of this threshold was found, and those genes that Ire up- and down- regulated can also be selected for performing the subsequent analysis. R software was also used to draw volcano map of DELs, DEmiRs and DEmRNAs. This data would be used in the following ceRNA network construction and protein interaction network construction.

### Code

DEG

``` r
#   Differential expression analysis with limma
setwd("~/Documents/capstone/rscripts/")
library(GEOquery)
library(limma)
library(dplyr)

# load series and platform data from GEO

gset <- getGEO("GSE16441", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pData(gset)$data_processing[1] #Data has already been LOISS normalized and log2 transformed, it is indicated in file.
# make proper column names to match toptable

fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000000000011111111111111111"
sml <- strsplit(gsms, split="")[[1]]
#Data has already been LOISS normalized and log2 transformed, it is indicated in file.
# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("RCC","Normal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01) 

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title", "Gene.ID"))

#Take table and mutate it based on condition to extract DEGs
deg <- tT %>%
  mutate(condition = abs(logFC) > 2 & adj.P.Val < 0.05 & Gene.symbol != "")
# summarize test results 
sum(deg$condition == TRUE)
sum(deg$condition== TRUE & deg$logFC > 0)
sum(deg$condition== TRUE & deg$logFC < 0)

deg <- subset(deg, deg$condition == TRUE)

write.csv(deg, file="../files/DEG_table.csv", row.names=FALSE)
write.table(deg$Gene.ID, file= "../files/DEG_IDs.txt", row.names = FALSE, col.names = FALSE)
write.table(deg$Gene.symbol, file= "../files/DEG_symbols.txt", row.names = FALSE, col.names = FALSE)
#REGEX to take long non coding RNA from the data
write.table(grep("LINC.*", deg$Gene.symbol, value = TRUE), "../files/LNC.txt", row.names = FALSE, col.names = FALSE)
#finding top 20 DEGs significant and top 6 LNC
head_deg <- subset(head(deg, n=20), select= c("Gene.symbol","adj.P.Val","P.Value","logFC"))

lnc <- deg[grep("LINC.*", deg$Gene.symbol), ]
lnc_head <- subset(head(lnc, n=6), select= c("Gene.symbol","adj.P.Val","P.Value","logFC"))
#writing to file
write.table(head_deg, file="../files/DEG_table_head.csv", row.names=FALSE)
write.table(lnc_head, file="../files/LNC_table_head.csv", row.names=FALSE)
# volcano plot (log P-value vs log fold change)
library(ggplot2)

tT %>% 
  mutate(Significant = adj.P.Val < 0.05 & abs(logFC) > 2) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col=Significant)) + geom_point()

ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE16441", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxIx=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
```

DEM

``` r
#   Differential expression analysis with limma
setwd("~/Documents/capstone/rscripts/")
library(GEOquery)
library(limma)
library(dplyr)

# load series and platform data from GEO

gset <- getGEO("GSE16441", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8659", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000000000011111111111111111"
sml <- strsplit(gsms, split="")[[1]]
#Data has already been LOISS normalized and log2 transformed, it is indicated in file.
# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("RCC","Normal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01) 
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","miRNA_ID","SPOT_ID"))

dem <- tT %>%
  mutate(condition = abs(logFC) > 1 & adj.P.Val < 0.05)
# summarize test results 
sum(dem$condition == TRUE)
sum(dem$condition== TRUE & dem$logFC > 0)
sum(dem$condition== TRUE & dem$logFC < 0)
dem <- subset(dem, dem$condition == TRUE)
write.csv(dem, file="../files/dem_table.csv", row.names=FALSE)
write.table(dem$miRNA_ID, file= "../files/dem_symbols.txt", row.names = FALSE, col.names = FALSE)
#finding top 20 DEMs significant
head_dem <- subset(head(dem, n=20), select= c("miRNA_ID","adj.P.Val","P.Value","logFC"))
write.table(head_dem, file="../files/dem_table_head.csv", row.names=FALSE)

# summarize test results 


# volcano plot (log P-value vs log fold change)
library(ggplot2)

tT %>% 
  mutate(Significant = adj.P.Val < 0.05 & abs(logFC) > 1) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col=Significant)) + geom_point()

ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02","#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE16441", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxIx=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
```

### Results

According to the cut-off criteria mentioned above, volcano plots Ire obtained. I finally found that in the samples of total RNA, 880 differentially expressed, having 213 upregulated and 667 downregulated. Among these total RNA, I have 6 DElncRNAs, 5 downregulated and 1 upregulated. Also, in the miRNA sample I found 87 differentially expressed miRNAs, 55 being upregulated and 32 downregulated.

### Figures

![Figure 1: Volcano Plot for the DEmRNAs samples; Blue points are the ones selected](figures/figure1.png){width="339"}

![Figure 2: Box Plot for the log2FC expression levels in the samples for the DEmRNAs samples](figures/figure2.png){width="372"}

![Figure 3: Volcano Plot for the DEmiRs samples; Blue points are the ones selected](figures/figure3.png){width="337"}

![Figure 4: Box Plot for the log2FC expression levels in the samples for the DEmiRs samples](figures/figure4.png){width="371"}

### Tables

![Table 1: Top 20 significant DEGs](table_figures/table1.png){width="287"}

![Table 2: All DELncRNAs](table_figures/table2.png){width="286"}

![Table 3: Top 20 significant DEMiRs](table_figures/table3.png){width="285"}

## Gene Ontology and Pathway Analysis

### Materials & Methods

The Database for Annotation, Visualization and Integrated Discovery (DAVID, [http://david.abcc.ncifcrf.gov/)](http://david.abcc.ncifcrf.gov/)) (Dennis et al., 2003) is a resource that offers functional interpretation of plenty of genes derived from genomic researches. In present study, DAVID database was used to perform Gene Ontology (GO) (Ashburner et al., 2000) and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway analysis (Kanehisa et al., 2019). The ontology contains three hierarchies: biological process (BP), cellular component (CC) and molecular function (MF). Pathway analysis is a functional analysis that maps genes to KEGG pathways. Herein, I performed GO analysis and KEGG pathway analysis using only the differentially expressed genes that I obtained. The P value denoted the significance of the GO and pathway term enrichment in the DEGs. "P value \<0.05" was set as the cut-off criterion. All of that was performed with an R script.

### Results

I found significant results in GO terms and KEGG pathways which will be displayed below. Most of them are related to functions in the kidney that might be affected by Renal Cell Carcinoma, like excretion, membrane transports, mineral absorption, etc\... Also I see notice that 43 genes are related in pathways in cancer in the KEGG pathways, and others are related to other metabolic pathways in kidney.

### Table

![Table 4: Top 5 significant GO Terms by Category and KEGG pathways](table_figures/table4.png){alt="Table 4: Top 5 significant GO Terms by Category and KEGG pathways" width="462"}

### Figure

![Figure 5: Display of the GO Terms and KEGG pathways by Count](figures/figure5.png)

### Code

``` r
setwd("~/Documents/capstone/rscripts/")
library(dplyr)
#This table was downloaded from DAVID database, and taking top 5 of each category by P value
#PS the data of table is already ordered by P.value
chart<-read.table("../files/chart_4CBB66A493291681316993753.txt", sep="\t", header = TRUE)
cc <- chart[grep("^GOTERM_CC", chart$Category), ][1:5, ]
bp <- chart[grep("^GOTERM_BP", chart$Category), ][1:5, ]
mf <- chart[grep("^GOTERM_MF", chart$Category), ][1:5, ]
kegg <- chart[grep("^KEGG", chart$Category), ][1:5, ]

#Restructuring the tables
cc <- cc %>% 
  mutate(Category = "CC") %>%
  select(Category, everything())
  
cc<-cc[,c(1:3, 5)]


bp <- bp %>% 
  mutate(Category = "BP") %>%
  select(Category, everything())

bp<-bp[,c(1:3, 5)]


mf <- mf %>% 
  mutate(Category = "MF") %>%
  select(Category, everything())

mf<-mf[,c(1:3, 5)]


kegg <- kegg %>% 
  mutate(Category = "KEGG") %>%
  select(Category, everything())

kegg<-kegg[,c(1:3, 5)]
#Binding all rows together
enrichement<- rbind(cc, bp, mf, kegg)

temp <- enrichement
temp$Term <- substr(temp$Term, 1, 10)
library(ggplot2)
# create a bar plot of the enrichment data
ggplot(enrichement, aes(x = Count, y = reorder(Term, Count), fill = Category)) +
  geom_bar(stat = "identity") +
  xlab("Count") +
  ylab("Row") +
  ggtitle("Enrichement Data by Row Name")

write.csv(enrichement, "../files/enrichement.csv")
  
```

## Prediction of miRNA-mRNA and lncRNA-miRNA pairs & Establishment of the ceRNA network of lncRNA-miRNA-mRNA

### Materials & Methods

Interaction between lncRNAs and miRNAs were taken from **DIANA-LncBase v3, miRNA-lncRNA targets** (<https://diana.e-ce.uth.gr/lncbasev3/interactions>) (Karagkouni et al., 2019). It is a reference repository with experimentally supported miRNA targets on non-coding transcripts. It catalogues approximately \~500,000 entries, corresponding to \~240,000 unique tissue and cell-type specific miRNA-lncRNA interactions. I chose only miRNA with high confidence level.\
I then intersected the predicted miRNAs with the DEmiRs found previously and obtained final miRNAs that I then used in the **mirWalk database** ([http://mirwalk.umm.uniheidelberg.de/search_mirnas/](http://mirwalk.umm.uni-heidelberg.de/search_mirnas/)) (Sticht et al., 2018), which is a comprehensive database that provides predicted miRNA binding sites of known genes of human, rat, cow, mouse, and dog. The predicted target genes found in this database search were intersected with my DEmRNAs to find final DEmRNAs. The final ceRNA network was visualized using cytoscape.

### Results

From the six DElncRNAs, only four were found to interact with miRNAs by using the DIANA-LncBase database. These interactions which are a total 64 can be visualized in the figure 6. However 56 unique miRNA were predicted (1 miRNA may interact with many lncRNA), and by taking the intersection of these 56 miRNA with the DEmiRNAs already found before, we will find only three miRNAs interacting with three LncRNAs: LINC00982 interacts with hsa-miR-28-5p; hsa-miR-532-5p LINC00944 interacts with hsa-miR-532-5p and LINC00473 interacts with hsa-miR-140-5p.\
Then using mirWalk database, these three miRNAs have been found to interact with 12,248 target mRNAs. By doing the intersection of these predicted mRNAs (which are gene targets) with the 880 differentially expressed mRNAs already found before, and I found 208 genes common between the two, these were the final DEmRNAs. 267 miRNA-mRNA interactions were selected (because 1 mRNA can be targeted by many miRNA), and combined with the lncRNA-miRNA interaction, I obtained the ceRNA network in figure 7.

### Code

### Figures

![Figure 6: LncRNA-miRNA interactions predicted by DIANA lncBase. The triangle represents lncRNA and the quadrilateral represents miRNA.](figures/figure6.png){width="453"}

![Figure 7: lncRNA-miRNA-mRNA ceRNA network. The triangle, lncRNA, the quadrilateral is miRNA, the oval is mRNA (red = up regulation, blue = down regulation)](figures/figure7.png){width="703"}

## PPI network analysis and key gene identification

### Materials & Methods

The plug-in CytoHubba in Cytoscape software is a visualization software that obtains the dense relationship through the degree of connectivity. Those hub genes in ceRNA network were identified by CytoHubba. The present study adopted the Search Tool for the Retrieval of interacting Genes/Proteins (STRING; version 11.5) for retrieving protein interactions between the 208 DEGs found earlier. The PPI network was visualized using cytoscape then CytoHubba was used to identify the hub genes of PPI network. The key genes were determined by taking the intersection of hub genes in the ceRNA network and PPI network. They were used for follow-up analysis.

### Results

The top 30 DEGs of ceRNA network acquired based on degree of connectivity using cytohubba were visualized in figure 8. When uploading these DEGs identified in ceRNA network to STRING website, 188 edges were found and 208 nodes in constructing the PPI network (figure 9). The top 30 DEGs of PPI network obtained by cytoHubba were visualized in figure 10. There were three key genes (BSND, RGMA, THSD4) which were determined by Venn Diagram of ceRNA network and PPI network (figure 11). They key genes were then visualized with cytoscape (figure 12)

### Figures

![Figure 8: Relationship network diagram of hub genes of ceRNA network](figures/figure8.png){width="446"}

![Figure 9: PPI generated by STRING database](figures/figure9.png){width="702"}

![Figure 10: Hub genes from the STRING PPI network, top 30](figures/figure10.png){width="444"}

![Figure 11: Venn diagram of key genes](figures/figure11.png){width="424"}

![Figure 12: Key genes lncRNA-miRNA-mRNA interaction](figures/figure12.png){width="421"}

### Code

``` r
library(dplyr)
setwd("~/Documents/capstone/rscripts/")
deg <- read.csv("../files/DEG_table.csv")
dem <- read.csv("../files/dem_table.csv")
#Reading file of lncRNA miRNA interactions then taking predicted miRNA ids without duplicates
lncmiR <- read.csv("../files/lncRNA_miRNA.csv")

write.table(unique(lncmiR$miRNA), "../files/miRNA_predicted.txt", row.names = FALSE, col.names = FALSE)
miRNA_predicted <- read.table("../files/miRNA_predicted.txt")
#Finding intersections of miRNA in DEM file and predicted ones
miRNA_intersection <- intersect(dem$miRNA_ID, miRNA_predicted$V1)
write.table(miRNA_intersection, "../files/miRNA_intersection.txt", row.names = FALSE, col.names=FALSE)
#We got 3 miRNAs found in both DEM and predicted."hsa-miR-532-5p" "hsa-miR-140-5p"and "hsa-miR-28-5p"
#Now we want to find the genes that are targets to these three using 3 DB

#MIRWalk Target of miRNA prediction
mirwalk_pred <- read.csv(file = "../files/miRWalk_miRNA_Targets-3.csv")
mirwalk_pred <- mirwalk_pred[, c(1,3)]
intersection_table <- subset(mirwalk_pred, mirwalk_pred[,2] %in% deg$Gene.symbol)
#removing duplicate rows because some miRNA mRNA interactions are written multiple times (multiple binding sequences)
intersection_table <- distinct(intersection_table)

write.csv(intersection_table, "../files/intersection_table_miRNA_mRNA.csv", row.names = FALSE)

length(unique(intersection_table$target)) 
# Print 

miRNA_values <- c("hsa-miR-140-5p", "hsa-miR-532-5p", "hsa-miR-28-5p")
lncRNA_miRNA_filtered <- subset(lncmiR, miRNA %in% miRNA_values)
write.table(lncRNA_miRNA_filtered, "../files/lncRNA_miRNA_filtered.csv", row.names = FALSE)
colnames(lncRNA_miRNA_filtered)[1:2] <- c("source", "target")
colnames(intersection_table)[1:2] <- c("source", "target")
#adding column to see if up or down regulated for further visualization in cytoscape
deg_updated <- deg
deg_updated$regulation <- ifelse(deg_updated$logFC > 0, "up", "down")
dem_updated <- dem
dem_updated$regulation <- ifelse(dem_updated$logFC > 0, "up", "down")
#for loop to see if gene up or down regulated
for (i in 1:nrow(intersection_table)) {
  gene <- intersection_table[i, 2]
  match_row <- which(deg_updated$Gene.symbol == gene)
  if (length(match_row) > 0) {
    regulation <- deg_updated[match_row[1], "regulation"]
  } else {
    regulation <- "NA"
  }
  intersection_table[i, "regulation"] <- regulation
}
for (i in 1:nrow(lncRNA_miRNA_filtered)) {
  mirna <- lncRNA_miRNA_filtered[i, 2]
  match_row <- which(dem_updated$miRNA_ID == mirna)
  if (length(match_row) > 0) {
    regulation <- dem_updated[match_row[1], "regulation"]
  } else {
    regulation <- "NA"
  }
  lncRNA_miRNA_filtered[i, "regulation"] <- regulation
}
cytoscape <- rbind(intersection_table, lncRNA_miRNA_filtered)
write.table(cytoscape, "../files/cytoscape_to_import.csv", row.names = FALSE)


write.table(unique(intersection_table$target), "../files/deg_intersection.txt", row.names = FALSE, col.names = FALSE)

#getting hub genes from both PPI and ceRNA network, 
ppi_hub <- read.csv("../files/PPI_hub.csv", header = TRUE)#exported from cytoscape
cerna_hub <- read.csv("../files/ceRNA_hub.csv", header = TRUE)#exported from cytoscape
cerna_hub$name <- gsub('"', '', cerna_hub$name)
result <- merge(ppi_hub, cerna_hub, by.x = "display.name", by.y = "name")
result$display.name
#to use in venn diagram
write.csv(ppi_hub$display.name, "../files/ppi_hub_genes.txt", row.names = FALSE, col.names = FALSE)
write.csv(cerna_hub$name, "../files/cerna_hub_genes.txt", row.names = FALSE, col.names = FALSE)
```

## Discussion

Renal cell carcinoma (RCC) is a type of kidney cancer that is responsible for a significant number of cancer-related deaths worldwide (Hsieh et al., 2017). The dysregulation of various genes, including those involved in the regulation of metabolic pathways, membrane transport, and apoptosis, is known to contribute to the development and progression of RCC. In recent years, long non-coding RNAs (lncRNAs) have emerged as important regulators of gene expression and have been implicated in the pathogenesis of various cancers, including RCC (Slack & Chinnaiyan, 2019).

In a recent study, I investigated the ceRNA (competing endogenous RNA) network in RCC using bioinformatic analysis. Specifically, I focused on the interactions between lncRNAs, miRNAs, and mRNAs and how they contribute to the regulation of gene expression in RCC.

The study identified several lncRNAs that interacted with miRNAs and mRNAs in the RCC ceRNA network, including LINC00473 (down regulated), LINC00944 (up regulated), and LINC00982 (down regulated), as well as miRNAs hsa-miR-140-5p (up regulated), hsa-miR-28-5p (up regulated), and hsa-miR-532-5p (down regulated), which were found to interact with the genes RGMA (up regulated), BSND (down regulated), and THSD4 (down regulated).

This study is in-line with (Chen & Zheng, 2021) that provided clinicopathological and experimental evidence that lncRNA LINC00944 acts as an oncogene in RCC. Also according to (Li et al., 2021) LINC00473 have been found to play essential roles in the occurrence and development of multiple human cancers. Finally, LINC00982 has been found to regulate cell proliferation in gastric cancer (Fei et al., 2016), breast cancer (Chi et al., 2022), renal cancer (Zhang et al., 2019) and many others.

Same for the miRNAs obtained was found to be dysregulated and involved in many types of cancers including RCC and others, for example (Flamini et al., 2017) studied the involvement of miR-140-5p in the case of lung cancer invasion, (Xu et al., 2017) found that miR‑28‑5p upregulation promoted the progression of ovarian cancer cell cycle, proliferation, migration and invasion. Finally, (Zhai et al., 2018) revealed that miR-532-5p was one of the most heavily downregulated miRNAs (we find it also here). Overexpression of miR-532-5p inhibited RCC cell proliferation, while knockdown of miR-532-5p promoted cell proliferation. 

Considering the genes, The RGMA gene, also known as RGM domain family member A, is a gene that encodes the Repulsive Guidance Molecule A (RGMA) protein. RGMA plays a role in the development of the nervous system by providing guidance cues to growing axons and influencing their direction of growth. In addition to its role in axon guidance, RGMA has been implicated in other biological processes. It has been found to be involved in the regulation of cell migration, angiogenesis, and bone formation (Fujita & Yamashita, 2017). RGMA has been identified as a potential tumor suppressor in various types of cancer. Regarding its association with cancer, some studies have suggested that RGMA upregulation may have protumorigenic effects. For example, in colorectal cancer, RGMA upregulation has been associated with increased tumor growth, invasion, and metastasis. It has been proposed that the upregulated RGMA may disrupt the repulsive guidance cues and promote tumor cell migration and invasion. (Chan et al., 2008)

Regarding the BSND gene, the downregulation of BSND leads to a dysfunctional or absent barttin protein. This disrupts the normal function of chloride channels in the kidneys and inner ear, resulting in the characteristic symptoms of the syndrome, such as salt wasting, renal tubular defects, and hearing loss.\
Currently, there is no strong evidence linking BSND downregulation to the development or progression of renal cell carcinoma (RCC) or other types of cancer. However, it is worth noting that downregulation or altered expression of various genes, including those involved in ion transport and regulation, can be observed in cancer. These alterations can contribute to the dysregulation of cellular processes and potentially affect tumor development and progression.

THSD4, also known as Thrombospondin Type-1 Domain-Containing Protein 4, is a protein that is involved in various biological processes, including angiogenesis, cell adhesion, and extracellular matrix organization. (Tan et al., 2022) showed that the expression levels of THSD4 were downregulated in bladder cancer tissues compared with normal tissues.

These results provide insights into the regulatory mechanisms underlying RCC development and progression. The ceRNA network provides a novel perspective on the regulation of gene expression, particularly with the identification of lncRNAs as important players in the process. Additionally, the study highlights the potential for targeting these lncRNAs, miRNAs, and mRNAs as a therapeutic strategy for RCC treatment.

However, it is important to note that these findings are based on bioinformatic analysis, and further experimental validation is required to confirm the interactions between these RNAs and genes in RCC. Nonetheless, the results provide a valuable starting point for further investigation into the mechanisms underlying RCC pathogenesis and the development of new therapeutic strategies for this disease.

## References

Ala, U. (2020). Competing Endogenous RNAs, Non-Coding RNAs and Diseases: An Intertwined Story. *Cells*, *9*(7), 1574. <https://doi.org/10.3390/cells9071574>

Ashburner, M., Ball, C. G., Blake, J. B., Botstein, D., Butler, H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M. P., Hill, D. J., Issel-Tarver, L., Kasarskis, A., Lewis, S. R., Matese, J. C., Richardson, J. D., Ringwald, M., Rubin, G. J., & Sherlock, G. (2000). Gene Ontology: tool for the unification of biology. *Nature Genetics*, *25*(1), 25--29. <https://doi.org/10.1038/75556>

Chan, S. S. M., Griffith, O. L., Tai, I. T., & Jones, S. J. (2008). Meta-analysis of Colorectal Cancer Gene Expression Profiling Studies Identifies Consistently Reported Candidate Biomarkers. *Cancer Epidemiology, Biomarkers & Prevention*, *17*(3), 543--552. <https://doi.org/10.1158/1055-9965.epi-07-2615>

Chen, C., & Zheng, H. (2021). LncRNA LINC00944 Promotes Tumorigenesis but Suppresses Akt Phosphorylation in Renal Cell Carcinoma. *Frontiers in Molecular Biosciences*, *8*. <https://doi.org/10.3389/fmolb.2021.697962>

Chi, F., Qiu, F., Jin, X., Chen, L. Q., He, G., & Han, S. (2022). LINC00982 Inhibits the Proliferation, Migration, and Invasion of Breast Cancer Cells Through the miR-765/DPF3 Axis. *DNA And Cell Biology*, *41*(4), 424--436. <https://doi.org/10.1089/dna.2021.0866>

Dennis, G., Sherman, B., Hosack, D. A., Yang, J., Gao, W., Lane, H. C., & Lempicki, R. A. (2003). DAVID: Database for Annotation, Visualization, and Integrated Discovery. *Genome Biology*, *4*(9). <https://doi.org/10.1186/gb-2003-4-9-r60>

Fan, L., Huang, X., Chen, J. M., Zhang, K., Gu, Y., Sun, J., & Cui, S. (2020). Long Noncoding RNA MALAT1 Contributes to Sorafenib Resistance by Targeting miR-140-5p/Aurora-A Signaling in Hepatocellular Carcinoma. *Molecular Cancer Therapeutics*, *19*(5), 1197--1209. <https://doi.org/10.1158/1535-7163.mct-19-0203>

Fei, Z., Yu, X., Zhou, M., Su, H., Zheng, Z., & Xie, C. (2016). Upregulated expression of long non-coding RNA LINC00982 regulates cell proliferation and its clinical relevance in patients with gastric cancer. *Tumor Biology*, *37*(2), 1983--1993. <https://doi.org/10.1007/s13277-015-3979-9>

Flamini, V., Jiang, W. G., & Cui, Y. (2017). Therapeutic Role of MiR-140-5p for the Treatment of Non-small Cell Lung Cancer. *Anticancer Research*, *37*(8). <https://doi.org/10.21873/anticanres.11825>

Fujita, Y., & Yamashita, T. (2017). The roles of RGMa-neogenin signaling in inflammation and angiogenesis. *Inflammation and Regeneration*, *37*(1). <https://doi.org/10.1186/s41232-017-0037-6>

Haiyan, S., Luo, G., & Chen, Q. (2021). Long noncoding RNAs as tumorigenic factors and therapeutic targets for renal cell carcinoma. *Cancer Cell International*, *21*(1). <https://doi.org/10.1186/s12935-021-01805-2>

Hsieh, J. J., Purdue, M. P., Signoretti, S., Swanton, C., Albiges, L., Schmidinger, M., Heng, D. Y., Larkin, J., & Ficarra, V. (2017). Renal cell carcinoma. *Nature Reviews Disease Primers*, *3*(1). <https://doi.org/10.1038/nrdp.2017.9>

Kanehisa, M., Sato, Y., Furumichi, M., Morishima, K., & Tanabe, M. (2019). New approach for understanding genome variations in KEGG. *Nucleic Acids Research*, *47*(D1), D590--D595. <https://doi.org/10.1093/nar/gky962>

Karagkouni, D., Paraskevopoulou, M. D., Tastsoglou, S., Skoufos, G., Karavangeli, A., Pierros, V., Zacharopoulou, E., & Hatzigeorgiou, A. G. (2019). DIANA-LncBase v3: indexing experimentally supported miRNA targets on non-coding transcripts. *Nucleic Acids Research*. <https://doi.org/10.1093/nar/gkz1036>

Kung, J. T., Colognori, D., & Lee, J. T. (2013). Long Noncoding RNAs: Past, Present, and Future. *Genetics*, *193*(3), 651--669. <https://doi.org/10.1534/genetics.112.146704>

Li, L., Zhang, X., Liu, N., Chen, X., & Peng, C. (2021). LINC00473: A novel oncogenic long noncoding RNA in human cancers. *Journal of Cellular Physiology*, *236*(6), 4174--4183. <https://doi.org/10.1002/jcp.30176>

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, *43*(7), e47. <https://doi.org/10.1093/nar/gkv007>

Shen, J., Liang, C., Su, X., Wang, Q., Ke, Y., Fang, J., Zhang, D., & Duan, S. (2022). Dysfunction and ceRNA network of the tumor suppressor miR-637 in cancer development and prognosis. *Biomarker Research*, *10*(1). <https://doi.org/10.1186/s40364-022-00419-8>

Slack, F. J., & Chinnaiyan, A. M. (2019). The Role of Non-coding RNAs in Oncology. *Cell*, *179*(5), 1033--1055. <https://doi.org/10.1016/j.cell.2019.10.017>

Sticht, C., De La Torre, C., Parveen, A., & Gretz, N. (2018). miRWalk: An online resource for prediction of microRNA binding sites. *PLOS ONE*, *13*(10), e0206239. <https://doi.org/10.1371/journal.pone.0206239>

Tan, Z., Fu, S. H., Feng, R., Huang, Y., Li, N., Wang, H., & Wang, J. (2022). Identification of Potential Biomarkers for Progression and Prognosis of Bladder Cancer by Comprehensive Bioinformatics Analysis. *Journal of Oncology*, *2022*, 1--17. <https://doi.org/10.1155/2022/1802706>

Xu, J., Jiang, N., Shi, H., Zhao, S., Yao, S., & Shen, H. (2017). miR-28-5p promotes the development and progression of ovarian cancer through inhibition of N4BP1. *International Journal of Oncology*, *50*(4), 1383--1391. <https://doi.org/10.3892/ijo.2017.3915>

Zhai, W., Ma, J., Zhu, R., Xu, C., Zhang, J. Z., Chen, Y., Chen, Z., Gong, D., Zheng, J., Chen, C., Li, S., Li, B., Huang, Y., Xue, W., & Zheng, J. (2018). MiR-532-5p suppresses renal cancer cell proliferation by disrupting the ETS1-mediated positive feedback loop with the KRAS-NAP1L1/P-ERK axis. *British Journal of Cancer*, *119*(5), 591--604. <https://doi.org/10.1038/s41416-018-0196-5>

Zhang, C., Li, X. Y., Luo, Z. Z., Wu, T. W., & Hu, H. (2019). Upregulation of LINC00982 inhibits cell proliferation and promotes cell apoptosis by regulating the activity of PI3K/AKT signaling pathway in renal cancer. *Eur Rev Med Pharmacol Sci*, *23*(4), 1443-1450.
