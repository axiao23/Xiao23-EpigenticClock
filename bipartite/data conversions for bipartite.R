#code for bipartite graphing
#code for TUBB2B vs all other genes 
library(broom)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
theme_set(theme_minimal())

methylation_data <- merged_data[-c(1:13)]
methylation_data <- methylation_data[-c(1:3)]
#scaling methylation data
methylation_data1 <- methylation_data %>% mutate_all(scale)
methylation_data1 <- cbind(methylation_data1, new_column = merged_data$submitter_id)
methylation_data1$submitter_id <- methylation_data1$new_column
gene_expression_data <- horvath_gene_expression[-c(1:8)]
#scaling gene expression data
gene_expression_data1 <- gene_expression_data %>% mutate_all(scale)
t_gene_expression_data <- t(gene_expression_data1)
t_gene_expression_data <- data.frame(t_gene_expression_data)
colnames(t_gene_expression_data) <- horvath_gene_expression$gene_name
t_gene_expression_data <- cbind(t_gene_expression_data, new_column = rownames(t_gene_expression_data))
t_gene_expression_data$submitter_id <- t_gene_expression_data$new_column
gene_expression_methylation <- merge(t_gene_expression_data, methylation_data1, by = "submitter_id",  all.x = FALSE, all.y = FALSE) 
new_gene_expression_methylation <- gene_expression_methylation[-c(7,8,18,20,21,22,25,29,30,31,32,33,34,39,40,46,49)]
new1_gene_expression_methylation <- new_gene_expression_methylation[-c(36 ,37 ,39 ,40 ,46 ,47 ,59 ,63 ,66 ,67 ,68 ,72 ,77 ,80 ,81 ,83 ,86 ,87, 88 ,89 ,90 ,91, 92 ,93 ,94, 95 ,96 ,97 ,98)]
new2_gene_expression_methylation <- new1_gene_expression_methylation[-c(73, 76, 77, 79, 81, 82, 83, 84, 85, 86, 80, 90, 93, 94, 96, 97, 98, 100, 101, 102, 103, 104, 105, 108, 109, 110, 122, 125, 126, 131, 133, 134, 135, 138, 139, 144, 146, 148)]
new3_gene_expression_methylation <- new2_gene_expression_methylation[-c(115, 118, 119, 122, 125, 127, 134, 135, 137, 139, 140, 141, 146, 147, 150, 152, 165, 170, 173, 174, 175, 178, 179, 187, 192, 200, 204, 205, 210, 213, 214, 218, 220, 223, 224, 226, 234, 236, 237, 238, 242, 253, 255, 256, 257, 262, 263, 265, 267, 271, 272, 276, 283, 288, 289, 291, 296, 299)]
new3_gene_expression_methylation <- new3_gene_expression_methylation[-c(78)]
gene_expression_methylation <- new3_gene_expression_methylation[-c(242, 244, 245)]
write.csv(gene_expression_methylation, "../data/TCGA_LUAD_EXPRESSION_METHYLATION_standardized_data.csv", row.names = F)
