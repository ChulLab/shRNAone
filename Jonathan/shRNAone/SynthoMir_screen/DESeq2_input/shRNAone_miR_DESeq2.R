# Load libraries
library(DESeq2)

# Step 1: Load count matrix
counts <- read.csv("shRNAone_miR_matrix.csv", row.names=1, check.names=FALSE)

# Verify data
print("Column types before conversion:")
print(sapply(counts, class))

# Round to integers
counts <- round(counts)

# Convert to matrix
counts <- as.matrix(counts)

# Step 2: Load metadata
colData <- read.csv("shRNAone_miR_metadata.csv", row.names=1, check.names=FALSE)
colData$condition <- factor(colData$condition, levels=c("day0", "day5", "DMSO", "Sorafenib"))
colData$replicate <- factor(colData$replicate)

# Verify column names match
if (!all(colnames(counts) == rownames(colData))) {
  stop("Column names in countData must match row names in colData.")
}

# Step 3: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=colData,
                              design=~replicate + condition)

# Step 4: Run DESeq2
dds <- DESeq(dds)

# Step 5: Extract results for Sorafenib vs. DMSO control
results <- results(dds, contrast=c("condition", "Sorafenib", "DMSO"))

# Step 6: Save results
write.csv(as.data.frame(results), "results_shRNAone_miR_DESeq2_output.csv")



# Step 7: Summarize results
print(summary(results))

#############################
##### DRAW PCA PLOT HERE#####
#############################

library(ggplot2)
library(locfit)

# Loading data
count_data <- read.csv("shRNAone_miR_matrix.csv", row.names = 1)
count_data <- as.matrix(round(count_data))
colData <- read.csv("shRNAone_miR_metadata.csv", row.names=1, check.names=FALSE)
colData$condition <- factor(colData$condition, levels=c("day0", "day5", "DMSO", "Sorafenib"))
colData$replicate <- factor(colData$replicate)

# Preparing DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = colData, 
                              design=~replicate + condition)

# Applying variance-stabilizing transformation with fitType='mean'
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = 'mean')

# Performing PCA
pca_data <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)

# Creating PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA Plot of Gene Expression Data",
       x = paste0("PC1: ", round(pca_data$PC1[1], 1), "% variance"),
       y = paste0("PC2: ", round(pca_data$PC2[1], 1), "% variance")) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")
