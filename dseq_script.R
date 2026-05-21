library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ── Загрузка данных ───────────────────────────────────────────────────────────
expr_raw <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/DEG_results.tsv",
                     sep='\t', row.names = 1)

colnames(expr_raw) <- gsub("-", ".", colnames(expr_raw))

sample_metadata <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/meta_responses.tsv", 
                            row.names = 1, sep='\t')
rownames(sample_metadata) <- gsub("-", ".", rownames(sample_metadata))
sample_metadata <- sample_metadata[sample_metadata$X0 %in% c('R', 'NR'), , drop = FALSE]

# ── Удаление выбросов ─────────────────────────────────────────────────────────
outliers <- c('Luc_65_S30_R1_001ReadsPerGene',
              'LuC_59_S1_R1_001ReadsPerGene', 
              'LuC_56_S19_R1_001ReadsPerGene',
              'Luc.1_S10_R1_001ReadsPerGene',
              'Luc.45_S5_R1_001ReadsPerGene',
              'LuC_72_S3_R1_001ReadsPerGene',
              'LuC_50_S5_R1_001ReadsPerGene',
              'Luc.92_S10_R1_001ReadsPerGene')
expr_raw        <- expr_raw[, !colnames(expr_raw) %in% outliers]
sample_metadata <- sample_metadata[!rownames(sample_metadata) %in% outliers, , drop = FALSE]

expr_raw        <- expr_raw[, intersect(colnames(expr_raw), rownames(sample_metadata))]
sample_metadata <- sample_metadata[colnames(expr_raw), , drop = FALSE]

# ── DESeq2 ────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = expr_raw,
                              colData   = sample_metadata,
                              design    = ~ X0)
dds$X0 <- relevel(dds$X0, ref = "NR")
dds    <- DESeq(dds)
res    <- results(dds)
summary(res)

# ── HGNC: Ensembl ID → gene symbol ───────────────────────────────────────────
hgnc_table <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/hgnc_complete_set.txt", 
                       row.names = 1, sep='\t')
symbol_map <- setNames(hgnc_table$symbol, hgnc_table$ensembl_gene_id)

rownames(res) <- sapply(rownames(res), function(id) {
  if (id %in% names(symbol_map) && !is.na(symbol_map[id]) && symbol_map[id] != "") {
    symbol_map[id]
  } else {
    id
  }
})

# ── Пороги ────────────────────────────────────────────────────────────────────
pvalue_threshold <- 1e-2
log2fc_threshold <- 2

# ── Таблица результатов с колонкой significance ───────────────────────────────
res_df <- as.data.frame(res) %>%
  mutate(
    gene_symbol  = rownames(.),
    significance = case_when(
      padj < pvalue_threshold & log2FoldChange >  log2fc_threshold ~ "Up in R",
      padj < pvalue_threshold & log2FoldChange < -log2fc_threshold ~ "Down in R",
      TRUE ~ "Not significant"
    )
  ) %>%
  select(gene_symbol, log2FoldChange, pvalue, padj, significance) %>%
  arrange(padj)

write.table(res_df,
            file      = "C:/Users/mefim/Downloads/rna_seq_diff_exp/dseq_table.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# ── Volcano plot ──────────────────────────────────────────────────────────────
palette <- c(
  "Up in R"       = "#E63946",
  "Down in R"     = "#457B9D",
  "Not significant" = "#CCCCCC"
)

significant_genes <- res_df %>%
  filter(significance != "Not significant") %>%
  arrange(padj) %>%
  slice_head(n = 20)   # подписываем топ-20, чтобы не было каши

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  
  # фоновые точки (not significant)
  geom_point(data = filter(res_df, significance == "Not significant"),
             alpha = 0.4, size = 1.2) +
  
  # значимые точки поверх
  geom_point(data = filter(res_df, significance != "Not significant"),
             alpha = 0.9, size = 1.8) +
  
  scale_color_manual(values = palette, name = NULL) +
  
  # пороговые линии
  geom_hline(yintercept = -log10(pvalue_threshold),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # подписи генов
  geom_text_repel(
    data          = significant_genes,
    aes(label     = gene_symbol),
    size          = 2.8,
    max.overlaps  = 15,
    segment.color = "grey50",
    segment.size  = 0.3,
    box.padding   = 0.4,
    force         = 2
  ) +
  
  # аннотации количества генов
  annotate("text", x = -9.5, y = max(-log10(res_df$padj), na.rm = TRUE) * 0.95,
           label = paste0("Down: ", sum(res_df$significance == "Down in R")),
           color = "#457B9D", size = 3.5, hjust = 0) +
  annotate("text", x =  9.5, y = max(-log10(res_df$padj), na.rm = TRUE) * 0.95,
           label = paste0("Up: ", sum(res_df$significance == "Up in R")),
           color = "#E63946", size = 3.5, hjust = 1) +
  
  xlim(-10, 10) +
  labs(
    title    = "Volcano plot: R vs NR",
    subtitle = paste0("padj < ", pvalue_threshold, "  |  |log2FC| > ", log2fc_threshold),
    x        = "Log2 Fold Change",
    y        = "-Log10 Adjusted p-value"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position   = c(0.88, 0.2),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size   = unit(0.4, "cm"),
    panel.grid.major  = element_line(color = "grey92"),
  )