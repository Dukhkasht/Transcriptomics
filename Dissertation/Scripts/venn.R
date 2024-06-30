library(ggvenn)

rm(list = ls())


data <- read.csv('maleinfertility.csv',sep = ',',header = TRUE)
data2 <- read.csv('diabetesgenes.csv',sep = ',',header = TRUE)

male_infertility_genes <- data$Maleinfertility
diabetes_genes <- data2$Diabetes



common_genes <- intersect(male_infertility_genes, diabetes_genes)

for (gene in common_genes) {
  cat(gene, "\n")
}

venn_list <- list(Male_Infertility = male_infertility_genes,
                  Diabetes = diabetes_genes)
ggvenn(venn_list)
