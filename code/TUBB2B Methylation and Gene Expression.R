#merging and extracting TUBB2B methylation and expression data 
#get data - merged data, from Epigenetic-Clocks Audreydata, and Horvath_gene_expression 
library(tidyverse)
library(tibble)
TUBB2B_gene_expression <- horvath_gene_expression[190,]
TUBB2B_gene_expression_clean <- TUBB2B_gene_expression[-c(1:8)]
transpose_TUBB2B_expression <- t(TUBB2B_gene_expression_clean)
transpose_TUBB2B_expression <- data.frame(transpose_TUBB2B_expression)
colnames(transpose_TUBB2B_expression)[1]  <- "TUBB2B_expression"
transpose_TUBB2B_expression <- cbind(transpose_TUBB2B_expression, new_column = rownames(transpose_TUBB2B_expression))
colnames(transpose_TUBB2B_expression)[2] <- "submitter_id"
TUBB2B_methylation <- data.frame(merged_data$submitter_id, merged_data$cg22736354)
colnames(TUBB2B_methylation)[1]  <- "submitter_id"
colnames(TUBB2B_methylation)[2] <- "TUBB2B_methylation"
#merged TUBB2B Expression and Methylation 
TUBB2B_expression_methylation <- merge(transpose_TUBB2B_expression, TUBB2B_methylation, by = "submitter_id", all.x = FALSE, all.y = FALSE) 
#plotting TUBB2B values
install.packages("hexbin")
ggplot(TUBB2B_expression_methylation, aes(x = TUBB2B_methylation, y = TUBB2B_expression)) + 
  geom_point(alpha = 0.5, color = "blue") + labs(title = "TUBB2B gene expression vs methylation",
                                               x = "TUBB2B_methylation",
                                               y = "TUBBB_expression") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        strip.text = element_text(face = "italic"),
        text = element_text(size = 16))

ggplot(TUBB2B_expression_methylation, aes(x = TUBB2B_methylation, y = TUBB2B_expression)) + 
  geom_line( color = "blue") 

ggplot(TUBB2B_expression_methylation, aes(x = TUBB2B_methylation, y = TUBB2B_expression)) + 
  geom_hex() 

