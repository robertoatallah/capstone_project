# Capstone Project

## Microarray Data

### Materials & Methods

One dataset was chosen by setting the screening criteria for the species type as "Homo sapiens", from GEO database of National Center for Biotechnology Information (NCBI) ([https://www.ncbi.nlm.nih.gov/geo/),](https://www.ncbi.nlm.nih.gov/geo/),) with such keywords being searched as "**Renal Cell Carcinoma**", "miRNA", "mRNA". Then study type "Expression Profiling by Array" was selected. As a result, we included the dataset **GSE16441**, which is divided into two platforms. Platform GPL6480 (Agilent-014850 Whole Human Genome Microarray 4x44K G4112F), Samples GSM413237-GSM413270, contains the total RNA expression data from 17 RCC tumor samples and 17 corresponding non-tumor samples and was used to find differentially expressed lncRNA and mRNA. Platform GPL8659 (Agilent Human miRNA Microarray Rel12.0) Samples GSM413271-GSM413304, contains the miRNA microarray expression data from 17 RCC tumors and 17 corresponding non-tumor samples and was used to find differentially expressed microRNAs.

## DEGs Analysis

### Materials & Methods

EGs were analyzed by R package "Linear Models for Microarray Data (**limma**)" function for datasets and "**GEOquery**" R package was used to retrieve to GSE from GEO database. Let me note that the data was already LOWESS normalized and log2 transformed as indicated in the matrix file. For the GPL6480 platform (the total RNA one), **log Foldchange\>2 and adjust. p\<0.05** were regarded as threshold values for selecting DEGs and DELs. For the GPL8659 platform samples (the micro-RNA one), l**og Foldchange\>1 and adjust. p\<0.05** were regarded as threshold values for selecting DEGs and DElncRNAs. Statistical significance for the selection of this threshold was found, and those genes that were up- and down- regulated can also be selected for performing the subsequent analysis. R software was also used to draw volcano map of DELs, DEmiRs and DEmRNAs. This data would be used in the following ceRNA network construction and protein interaction network construction.

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
pData(gset)$data_processing[1] #Data has already been LOWESS normalized and log2 transformed, it is indicated in file.
# make proper column names to match toptable

fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000000000011111111111111111"
sml <- strsplit(gsms, split="")[[1]]
#Data has already been LOWESS normalized and log2 transformed, it is indicated in file.
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
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
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
#Data has already been LOWESS normalized and log2 transformed, it is indicated in file.
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
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
```

### Results

According to the cut-off criteria mentioned above, volcano plots were obtained. We finally found that in the samples of total RNA, 880 were differentially expressed, having 213 upregulated and 667 downregulated. Among these total RNA, we have 6 DElncRNAs, 5 downregulated and 1 upregulated. Also, in the miRNA sample we found 87 differentially expressed miRNAs, 55 being upregulated and 32 downregulated.

### Figures

![Figure 1: Volcano Plot for the DEmRNAs samples; Blue points are the ones selected](figures/figure1.png){width="342"}

![Figure 2: Box Plot for the log2FC expression levels in the samples for the DEmRNAs samples](figures/figure2.png){width="372"}

![Figure 3: Volcano Plot for the DEmiRs samples; Blue points are the ones selected](figures/figure3.png){width="337"}

![Figure 4: Box Plot for the log2FC expression levels in the samples for the DEmiRs samples](figures/figure4.png){width="371"}
