setwd("~/Documents/capstone/rscripts/code/")
deg <- read.csv("../deg/DEG_table.csv")
dem <- read.csv("../dem/dem_table.csv")
#Reading file of lncRNA miRNA interactions then taking predicted miRNA ids without duplicates
lncmiR <- read.csv("lncRNA_miRNA.csv")
write.table(unique(lncmiR$miRNA), "miRNA_predicted.txt", row.names = FALSE, col.names = FALSE)
miRNA_predicted <- read.table("miRNA_predicted.txt")
#Finding intersections of miRNA in DEM file and predicted ones
miRNA_intersection <- intersect(dem$miRNA_ID, miRNA_predicted$V1)
write.table(miRNA_intersection, "miRNA_intersection.txt", row.names = FALSE, col.names=FALSE)
#We got 3 miRNAs found in both DEM and predicted.
#Now we want to find the genes that are targets to these three using 3 DB
pred_mirDB <- read.csv("mirDB.csv")
pred_mirTarBase <- read.csv("search_result-2.csv")
target1<- read.csv("Copy of TargetScan8.0__miR-140-5p.predicted_targets.csv")
target2 <- read.csv("Copy of TargetScan8.0__miR-28-5p_708-5p.predicted_targets.csv")
target3 <- read.csv("Copy of TargetScan8.0__miR-532-5p.predicted_targets.csv")
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

#genes intersection predicted vs DEGs
subset_table_intersection <- subset(subset_table, subset_table[,2] %in% deg$Gene.symbol)
array<- unique(subset_table_intersection$Target)
write.table(unique(subset_table_intersection)$Target, "genes_intersection.txt", row.names = FALSE, col.names = FALSE)


target_rows <- data.frame() # Create an empty data frame to store the target rows
for (target in array) {
  target_rows_subset <- subset(combined_table, Target == target)
  target_rows <- rbind(target_rows, target_rows_subset)
}

# Print 
write.table(target_rows, "target_rows.csv", row.names = FALSE)
miRNA_values <- c("hsa-miR-140-5p", "hsa-miR-532-5p", "hsa-miR-28-5p")
lncRNA_miRNA_filtered <- subset(lncmiR, miRNA %in% miRNA_values)
write.table(lncRNA_miRNA_filtered, "lncRNA_miRNA_filtered.csv", row.names = FALSE)
