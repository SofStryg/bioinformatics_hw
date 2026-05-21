genes <- c("BRCA1", "TP53", "EGFR")
expression <- c(12.5, 45.2, 30.1)
condition <- c("Control", "Treatment", "Treatment")

exp_data <- data.frame(genes, expression, condition)

print("Структура таблицы exp_data:")
str(exp_data)

png("r_task/expression_plot.png")
barplot(expression, names.arg = genes,
        xlab = "Genes", ylab = "Expression",
        col = "skyblue", main = "Gene Expression Levels")
dev.off()
