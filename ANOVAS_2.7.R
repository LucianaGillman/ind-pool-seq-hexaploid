#This script is an example of Variance analysis and plot, you can substitute
# "CCC" and "Missing.Data" for the variables that you want to analyse
# Load required libraries
#library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)
library(car)
library(pgirmess)
library(ggpubr)

# Set working directory
setwd("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3")

# Load and preprocess data
load_data <- function(file_path) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  data <- data %>%
    mutate(
      Sample.Size = as.character(Sample.Size),
      MAF = as.character(MAF.pool),
      Accession = as.character(Accession),
      Missing.Data = as.character(Missing.Data),
      Profundidad = gsub("1M", "1.8Mr", gsub("2M", "3.0Mr", gsub("3M", "4.8Mr", Sequence.Depth))),
    )
  return(data)
}

result <- load_data("CCC_Rep.csv")

# Function to perform ANOVA and Tukey's HSD
perform_anova <- function(data, response, group_var, output_prefix) {
  levene_test <- leveneTest(as.formula(paste(response, "~", group_var)), data = data)
  write.csv(levene_test, paste0(output_prefix, "_levene.csv"))
  
  anova_model <- aov(as.formula(paste(response, "~", group_var)), data = data)
  anova_summary <- summary(anova_model)
  capture.output(anova_summary, file = paste0(output_prefix, "_anova.txt"))
  
  tukey <- TukeyHSD(anova_model)
  capture.output(tukey, file = paste0(output_prefix, "_tukey.txt"))
  
  cld <- multcompLetters4(anova_model, tukey)
  
  return(list(model = anova_model, tukey = tukey, cld = cld))
}

# Example: CCC vs Missing Data
ccc_results <- result %>%
  group_by(Missing.Data, Accession) %>%
  summarise(CCC = mean(CCC))

anova_results_ccc <- perform_anova(ccc_results, "CCC", "Missing.Data", "CCC_Missing.Data")

# Plot Function
plot_data <- function(data, x_var, y_var, cld, y_label, x_label, y_limits) {
  cld_table <- data %>%
    group_by(!!sym(x_var)) %>%
    summarise(
      mean = mean(!!sym(y_var)),
      quant = quantile(!!sym(y_var), probs = 0.75)
    ) %>%
    mutate(cld = cld$Letters)
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(data = cld_table, aes(x = !!sym(x_var), y = quant, label = cld), size = 4, vjust = -2) +
    geom_jitter(aes(color = Accession), width = 0.3) +
    stat_summary(fun.data = "mean_se", size = 0.0) +
    stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. - 0.25, yend = ..y..), size = 0.8) +
    stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. + 0.25, yend = ..y..), size = 0.8) +
    xlab(x_label) +
    ylab(y_label) +
    ylim(y_limits)
}

ccc_plot <- plot_data(ccc_results, "Missing.Data", "CCC", anova_results_ccc$cld$Missing.Data,
                      "CCC", "Missing data (%)", c(0.5, 1))
