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
#genes intersection predicted vs DEGs
subset_table_intersection <- subset(subset_table, subset_table[,2] %in% deg$Gene.symbol)
array<- unique(subset_table_intersection$Target)
write.table(unique(subset_table_intersection)$Target, "../files/genes_intersection.txt", row.names = FALSE, col.names = FALSE)


target_rows <- data.frame() # Create an empty data frame to store the target rows
for (target in array) {
  target_rows_subset <- subset(combined_table, Target == target)
  target_rows <- rbind(target_rows, target_rows_subset)
}
#MIRWalk Target of miRNA prediction
#mirnet prediction
mirnet_pred <- read.csv(file = "../files/miRWalk_miRNA_Targets-3.csv")
mirnet_pred <- mirnet_pred[, c(1,3)]
intersection_table <- subset(mirnet_pred, mirnet_pred[,2] %in% deg$Gene.symbol)
#removing duplicate rows
intersection_table <- distinct(intersection_table)

write.csv(intersection_table, "../files/intersection_table_miRNA_mRNA.csv", row.names = FALSE)

length(unique(intersection_table$genesymbol)) 
# Print 
write.table(target_rows, "../files/target_rows.csv", row.names = FALSE)
miRNA_values <- c("hsa-miR-140-5p", "hsa-miR-532-5p", "hsa-miR-28-5p")
lncRNA_miRNA_filtered <- subset(lncmiR, miRNA %in% miRNA_values)
write.table(lncRNA_miRNA_filtered, "../files/lncRNA_miRNA_filtered.csv", row.names = FALSE)
colnames(lncRNA_miRNA_filtered)[1:2] <- c("source", "target")
colnames(intersection_table)[1:2] <- c("source", "target")

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
intersection_table$regulation =="up" == TRUE

write.table(unique(intersection_table$target), "../files/deg_intersection.txt", row.names = FALSE, col.names = FALSE)


ppi_hub <- read.csv("../files/PPI_hub.csv", header = TRUE)
cerna_hub <- read.csv("../files/ceRNA_hub.csv", header = TRUE)
cerna_hub$name <- gsub('"', '', cerna_hub$name)
result <- merge(ppi_hub, cerna_hub, by.x = "display.name", by.y = "name")
result$display.name
write.csv(ppi_hub$display.name, "../files/ppi_hub_genes.txt", row.names = FALSE, col.names = FALSE)
write.csv(cerna_hub$name, "../files/cerna_hub_genes.txt", row.names = FALSE, col.names = FALSE)

ppi_hub_genes <- readLines("../files/ppi_hub_genes.txt")
cerna_hub_genes <- readLines("../files/cerna_hub_genes.txt")
library(VennDiagram)
venn.diagram(list(set1 = set1, set2 = set2), filename = "venn_diagram.png",
             col = c("cornflowerblue", "green"), fill = c("cornflowerblue", "green"),
             alpha = c(0.5, 0.5), label.col = c("white", "white"),
             cex = 1.5, fontfamily = "Helvetica", cat.cex = 1.8,
             margin = 0.05)
