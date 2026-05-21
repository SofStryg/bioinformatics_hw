library(limma)
library(ggplot2)
library(ggrepel)
library(Biobase)
library(dplyr)
library(patchwork)  # для объединения графиков

# ── Загрузка данных ───────────────────────────────────────────────────────────
exp <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/GSE63885/expression_for_limma2.csv", 
                header = TRUE, row.names = 'Gene.Symbol')
ann <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/GSE63885/annotation_for_limma.csv", 
                header = TRUE, row.names = 'X')
ann_1 <- read.csv("C:/Users/mefim/Downloads/rna_seq_diff_exp/GSE63885/limma_table.tsv", 
                  header = TRUE, sep = "\t")

# ── Фильтрация низкоэкспрессируемых генов (топ 50%) ──────────────────────────
medians <- apply(exp, 1, median)
exp     <- exp[medians > median(medians), ]

# ── ExpressionSet + дизайн ────────────────────────────────────────────────────
exp_set <- ExpressionSet(assayData  = as.matrix(exp),
                         phenoData  = AnnotatedDataFrame(ann))

slope  <- factor(ann$clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1,
                 levels = c("pNC", "pCR"))
design <- model.matrix(~ slope)
colnames(design) <- c("Intercept", "pCR")

# ── limma ─────────────────────────────────────────────────────────────────────
fit <- lmFit(exp_set, design)
fit <- eBayes(fit)
top <- topTable(fit, coef = "pCR", adjust = "BH", n = Inf)
top <- top[, c("logFC", "P.Value", "adj.P.Val")]
top$gene_symbol <- rownames(top)

# ── Пороги и столбцы significance ────────────────────────────────────────────
pvalue_threshold <- 0.05

for (lfc in c(1, 2, 3)) {
  # по adj.P.Val
  top[[paste0("sig_adjP_lfc", lfc)]] <- ifelse(
    top$adj.P.Val < pvalue_threshold & abs(top$logFC) > lfc, "+", "-"
  )
  # по P.Value
  top[[paste0("sig_pval_lfc", lfc)]] <- ifelse(
    top$P.Value < pvalue_threshold & abs(top$logFC) > lfc, "+", "-"
  )
}

# ── Подсчёт значимых генов ────────────────────────────────────────────────────
cat("════════════════════════════════════════\n")
cat("Количество значимых генов:\n\n")
for (lfc in c(1, 2, 3)) {
  n_adj  <- sum(top[[paste0("sig_adjP_lfc", lfc)]] == "+")
  n_pval <- sum(top[[paste0("sig_pval_lfc", lfc)]] == "+")
  cat(sprintf("  |logFC| > %d  |  adj.P.Val: %d  |  P.Value: %d\n", lfc, n_adj, n_pval))
}
cat("════════════════════════════════════════\n")

# ── Сохранение таблицы ────────────────────────────────────────────────────────
write.table(top,
            file      = "C:/Users/mefim/Downloads/rna_seq_diff_exp/GSE63885/limma_table.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# ── Функция volcano plot ──────────────────────────────────────────────────────
palette_volcano <- c("-" = "#CCCCCC", "+" = "#E63946")

plot_volcano <- function(df, lfc_threshold, sig_col, use_adjp = FALSE) {
  
  y_var   <- if (use_adjp) "adj.P.Val" else "P.Value"
  y_label <- if (use_adjp) "-Log10 Adjusted p-value" else "-Log10 p-value"
  p_line  <- -log10(pvalue_threshold)
  
  # разделяем на up/down/ns для цвета
  df <- df %>% mutate(
    direction = case_when(
      .data[[sig_col]] == "+" & logFC >  lfc_threshold ~ "Up in pCR",
      .data[[sig_col]] == "+" & logFC < -lfc_threshold ~ "Down in pCR",
      TRUE ~ "Not significant"
    )
  )
  
  pal <- c("Up in pCR"       = "#E63946",
           "Down in pCR"     = "#457B9D",
           "Not significant" = "#CCCCCC")
  
  sig_genes <- df %>%
    filter(direction != "Not significant") %>%
    arrange(.data[[y_var]]) %>%
    slice_head(n = 15)
  
  n_up   <- sum(df$direction == "Up in pCR")
  n_down <- sum(df$direction == "Down in pCR")
  y_max  <- max(-log10(df[[y_var]]), na.rm = TRUE)
  
  ggplot(df, aes(x = logFC, y = -log10(.data[[y_var]]), color = direction)) +
    
    geom_point(data = filter(df, direction == "Not significant"),
               alpha = 0.4, size = 1.2) +
    geom_point(data = filter(df, direction != "Not significant"),
               alpha = 0.9, size = 1.8) +
    
    scale_color_manual(values = pal, name = NULL) +
    
    geom_hline(yintercept = p_line,
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    
    geom_text_repel(
      data         = sig_genes,
      aes(label    = gene_symbol),
      size         = 2.5,
      max.overlaps = 15,
      segment.color = "grey50",
      segment.size  = 0.3,
      box.padding   = 0.4,
      force         = 2
    ) +
    
    annotate("text", x = -Inf, y = y_max * 0.97,
             label = paste0("Down: ", n_down),
             color = "#457B9D", size = 3.2, hjust = -0.1) +
    annotate("text", x =  Inf, y = y_max * 0.97,
             label = paste0("Up: ", n_up),
             color = "#E63946", size = 3.2, hjust = 1.1) +
    
    labs(
      title    = paste0("Volcano plot  |logFC| > ", lfc_threshold),
      subtitle = paste0(if (use_adjp) "adj." else "raw", 
                        " p-value < ", pvalue_threshold,
                        "  |  Up: ", n_up, "  Down: ", n_down),
      x = "Log Fold Change",
      y = y_label
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position   = c(0.88, 0.08),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key.size   = unit(0.4, "cm"),
      panel.grid.major  = element_line(color = "grey92")
    )
}

# ── Построение графиков (по raw p-value) ──────────────────────────────────────
p1 <- plot_volcano(top, 1, "sig_pval_lfc1")
p2 <- plot_volcano(top, 2, "sig_pval_lfc2")
p3 <- plot_volcano(top, 3, "sig_pval_lfc3")

# все три рядом
p1 + p2 + p3 + plot_layout(ncol = 3)
