# Importing required libraries
library("DESeq2")
library("gplots")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("scatterplot3d")

# Importing the data from the count table
counts <- read.delim("path_to_your_file/abi_v3.txt", sep = "\t", header = TRUE, row.names = 1)
counts <- as.matrix(counts)

# Creating a data frame for the design
design <- data.frame(condition = factor(c("ABA", "ABA", "ABA", "Control", "Control", "Control", "ABA", "ABA", "ABA", "Control", "Control", "Control")))
rownames(design) <- colnames(counts)

# Creating DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~ condition)
dds$genotype <- factor(c("abi1", "abi1", "abi1", "abi1", "abi1", "abi1", "wt", "wt", "wt", "wt", "wt", "wt"))
dds$genotype <- relevel(dds$genotype, "wt")
dds$condition <- relevel(dds$condition, "Control")
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)

# Getting DESeq results
res1 <- results(dds, name = "genotypeabi1.conditionABA")
res1 <- res1[complete.cases(res1$padj), ]
res1 <- res1[order(res1$padj), ]
new_columns <- data.frame(GeneID = rownames(res1))
res1 <- cbind(new_columns, res1)

# Writing results to CSV
write.csv(res1, file = "abi1_interaction_effect.csv", quote = FALSE, row.names = FALSE)

# PCA plots
rld <- rlog(dds, blind = FALSE)
pdf("PCA_1&2_abi1_interaction_effect.pdf")
pcaData <- plotPCA(rld, intgroup = c("condition", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()

# pheatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[, c("condition", "genotype")])
pdf("pheatmap_abi1_interaction_effect.pdf")
pheatmap(assay(rld)[select, ], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df)
dev.off()

# Sample to sample distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
pdf("SampletoSampleDistance_abi1_interaction_effect.pdf")
rownames(sampleDistMatrix) <- paste(rld$condition, rld$genotype, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

# MA plot
pdf("MA_plot_abi1_interaction_effect.pdf")
plotMA(dds, main = "MA plot", ylim = c(-2, 2))
dev.off()

# Dispersion plot
pdf("Dispersionplot_abi1_interaction_effect.pdf")
plotDispEsts(dds)
dev.off()

# Scatterplot
rld <- rlog(dds, blind = FALSE)
rlogMat <- assay(rld)
PCA <- prcomp(t(rlogMat))
percentVar <- round(100 * PCA$sdev^2 / sum(PCA$sdev^2), 1)
dataGG <- data.frame(
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2],
  PC3 = PCA$x[, 3],
  PC4 = PCA$x[, 4],
  condition = colData(rld)$condition,
  genotype = colData(rld)$genotype
)
pdf("3D plot.pdf")
colors <- c("#bf0f0f", "#56B4E9")
colors <- colors[as.numeric(dataGG$condition)]
shapes <- c(15, 16, 17, 22)
shapes <- shapes[as.numeric(dataGG$genotype)]
s3d <- scatterplot3d(
  dataGG[, 1:3],
  pch = shapes,
  color = colors,
  main = "3D Plot",
  angle = 50
)
legend(
  "bottomright",
  legend = levels(dataGG$condition),
  pch = c(15, 16, 17, 22),
  col = c("#bf0f0f", "#56B4E9"),
  inset = 0.05,
  xpd = TRUE,
  horiz = TRUE
)
legend(
  "topleft",
  legend = levels(dataGG$genotype),
  pch = c(15, 16, 17, 22),
  inset = 0.05,
  xpd = TRUE,
  horiz = TRUE
)
dev.off()

# Alternative pheatmap
mat <- assay(rld)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[, c("condition")])
colnames(df) <- "condition"
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col = df)

# Plot counts
pdf("plotcounts.pdf")
res <- results(dds)
topGene <- rownames(res)[which.min(res$padj)]
data <- plotCounts(dds, gene = topGene, intgroup = c("condition", "genotype"), returnData = TRUE)
ggplot(data, aes(x = condition, y = count, color = genotype, group = genotype)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() +
  labs(title = topGene) + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Screeplot
pdf("Screeplot.pdf")
barplot(percentVar, cex.names = 1, xlab = paste("Principal component (PC), 1-", length(PCA$sdev)), ylab = "Proportion of variation (%)", main = "Scree plot", ylim = c(0, 100))
dev.off()

# Volcano plot
pdf("volcanoplot.pdf")
tab <- data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
head(tab)
par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2] ~ fold ~ change), ylab = expression(-log[10] ~ pvalue))
lfc <- 2
pval <- 0.05
signGenes <- (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
dev.off()
