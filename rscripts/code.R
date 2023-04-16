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


