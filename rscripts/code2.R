library(dplyr)
setwd("~/Documents/capstone/rscripts/")
#This is the code that I tried to find the genes that are targets to the three miRNA using 3 DB
#the search was done using the websites of the databases and the tables were downloaded and used here
pred_mirDB <- read.csv("../files/mirDB.csv")
pred_mirTarBase <- read.csv("../files/search_result-2.csv")
target1<- read.csv("../files/Copy of TargetScan8.0__miR-140-5p.predicted_targets.csv")
target2 <- read.csv("../files/Copy of TargetScan8.0__miR-28-5p_708-5p.predicted_targets.csv")
target3 <- read.csv("../files/Copy of TargetScan8.0__miR-532-5p.predicted_targets.csv")
target2 <- subset(target2, target2$Representative.miRNA == "hsa-miR-28-5p")

target1 <- target1[, -9]
target2 <- target2[, -9]
target3 <- target3[, -9]

pred_targetScan1 <- union(target1, target2)
pred_targetScan <- union(pred_targetScan1, target3)
#the one of target scan I performed it one by one so I had to do the union
#here renaming the columns to be the same so I could later join them
colnames(pred_mirDB)[colnames(pred_mirDB)=="Gene.Symbol"] <- "Target"
colnames(pred_mirDB)[colnames(pred_mirDB)=="miRNA.Name"] <- "miRNA"


colnames(pred_targetScan)[colnames(pred_targetScan)=="Target.gene"] <- "Target"
colnames(pred_targetScan)[colnames(pred_targetScan)=="Representative.miRNA"] <- "miRNA"
#taking only miRNA and target from tables because i am only interested in them
pred_mirDB <- pred_mirDB[, c("miRNA", "Target")]
pred_mirTarBase <- pred_mirTarBase[, c("miRNA", "Target")]
pred_targetScan <- pred_targetScan[, c("miRNA", "Target")]
#Combining all three tables together, without duplicates
combined_table <- rbind(pred_mirDB, pred_mirTarBase, pred_targetScan)
combined_table<-unique(combined_table)
#adding rows to take count if one row is found in a table or not
combined_table$mirDB <- integer(nrow(combined_table))
combined_table$targetScan <- integer(nrow(combined_table))
combined_table$mirTarBase <- integer(nrow(combined_table))
#conditions, 1 found, 0 not found
combined_table$mirDB <- ifelse(paste(combined_table$Target, combined_table$miRNA) %in% 
                                 paste(pred_mirDB$Target, pred_mirDB$miRNA), 1, 0)

combined_table$targetScan <- ifelse(paste(combined_table$Target, combined_table$miRNA) %in% 
                                      paste(pred_targetScan$Target, pred_targetScan$miRNA), 1, 0)

combined_table$mirTarBase <- ifelse(paste(combined_table$Target, combined_table$miRNA) %in% 
                                      paste(pred_mirTarBase$Target, pred_mirTarBase$miRNA), 1, 0)

#summing, then taking only ones with sum greater or equal than 2
combined_table$Sum <- rowSums(combined_table[,3:5])
subset_table <- subset(combined_table, Sum >= 2)
temp <- subset(head(subset_table, n=20))
#genes intersection predicted vs DEGs, found 17, not enough
subset_table_intersection <- subset(subset_table, subset_table[,2] %in% deg$Gene.symbol)
array<- unique(subset_table_intersection$Target)
write.table(unique(subset_table_intersection)$Target, "../files/genes_intersection.txt", row.names = FALSE, col.names = FALSE)

x <- subset_table_intersection
y <- lncRNA_miRNA_filtered[, 1:2]
merged <- merge(x, y, by.x = "miRNA", by.y = "target")
merged <- merged[, c("source", "miRNA","Target", "mirDB", "targetScan",  "mirTarBase", "Sum")]
colnames(merged)[1] <- "lncRNA"  
write.table(merged, "../files/cytoscape_to_import2.csv", row.names = FALSE)
#table to view network