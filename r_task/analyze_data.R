data <- read.csv("r_task/sample_data.csv")

mean_score <- mean(data$Score)
cat("Среднее значение Score:", mean_score, "\n")

treatment_data <- subset(data, Group == "Treatment")
max_score_treatment <- max(treatment_data$Score)
cat("Максимальный Score в группе Treatment:", max_score_treatment, "\n")

png("r_task/score_boxplot.png")
boxplot(Score ~ Group, data = data,
        main = "Score Distribution by Group",
        xlab = "Group", ylab = "Score",
        col = c("lightgreen", "lightcoral"))
dev.off()
