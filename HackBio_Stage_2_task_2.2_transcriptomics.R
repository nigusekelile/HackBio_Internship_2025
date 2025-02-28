# Load necessary libraries
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Function to generate volcano plot and identify genes
analyze_rnaseq <- function(file_path) {
  # Read data
  data <- read.delim(file_path, sep = "   ", header = TRUE)
  
  # Add -log10(pvalue) column
  data$neg_log_pval <- -log10(data$pvalue)
  
  # Volcano plot
  volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = neg_log_pval)) +
    geom_point(aes(color = ifelse(abs(log2FoldChange) > 1 & pvalue < 0.01, "Significant", "Non-significant")), alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(p-value)", color = "Significance") +
    theme_minimal()
  
  # Identify up/downregulated genes
  upregulated <- data %>% filter(log2FoldChange > 1 & pvalue < 0.01) %>% arrange(desc(log2FoldChange))
  downregulated <- data %>% filter(log2FoldChange < -1 & pvalue < 0.01) %>% arrange(log2FoldChange))

# Get gene functions
get_gene_info <- function(genes) {
  gene_symbols <- genes$Gene
  gene_info <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "GENENAME", keytype = "SYMBOL")
  return(data.frame(Gene = gene_symbols, Description = unname(gene_info)))
}

top5_up <- head(upregulated, 5)
top5_down <- head(downregulated, 5)

top5_up_info <- get_gene_info(top5_up)
top5_down_info <- get_gene_info(top5_down)

return(list(
  volcano_plot = volcano_plot,
  upregulated = upregulated,
  downregulated = downregulated,
  top5_up_info = top5_up_info,
  top5_down_info = top5_down_info
))
}

# Usage
result <- analyze_rnaseq("results.txt")
print(result$volcano_plot)
