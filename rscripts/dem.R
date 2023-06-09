#   Differential expression analysis with limma
setwd("~/Documents/capstone/rscripts/")
#load packages
library(GEOquery)
library(limma)
library(dplyr)
# load series and platform data from GEO, this platform is for miRNA
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
# fit linear model
fit <- lmFit(gset, design) 
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
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

