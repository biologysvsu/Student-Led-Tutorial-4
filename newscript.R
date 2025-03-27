#Install and load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport", ask = FALSE, update = FALSE)
library(tximport)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# Point to Salmon quant files
files <- c("mock" = "quant_mock.sf",
           "covid" = "quant_covid.sf")

# Load tx2gene
tx2gene <- read.csv("tx2gene.csv", header=FALSE, col.names=c("TXNAME", "GENEID"))

# Import data
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# Optional: view matrix of expression
head(txi$abundance)


# Get top 50 most expressed transcripts
top_tx <- head(order(rowSums(txi$abundance), decreasing=TRUE), 50)
top_data <- log2(txi$abundance[top_tx, ] + 1)

##Create heatmaps for quantified transcripts (e.g., top 50):
pheatmap(top_data,
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         scale="row",
         fontsize=8,
         main="Top 50 Expressed Transcripts",
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50))
