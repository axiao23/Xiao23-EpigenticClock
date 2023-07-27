library(tidyverse)
library(broom)
#entire dataset linear regressions
tcga_dat = read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
dim(tcga_dat)

cg_indices = grep("cg",names(tcga_dat))
genes <- colnames(tcga_dat)[2:242]

# first, select just the methylation sites and TUBB2B expression
tcga_meth= tcga_dat %>% select(submitter_id,ZNF551, all_of(cg_indices))
dim(tcga_meth)

# now we are going to use a function called pivot_longer to convert the data so that
# there is one column that has all the cg probes in it

tcga_meth_long = tcga_meth %>% pivot_longer(cols = contains("cg"),
                                            names_to = "probeID", values_to="methylation")

#genes_loop <- for (colnam in colnames(tcga_meth_long[c(2:242)])) {
# print(colnam)
#}

model_results= tcga_meth_long %>%
  nest(data = -probeID) %>% # tell the code to work each probe
  mutate(
    fit = map(data, ~  lm(ZNF551 ~ methylation, data = .x)), # run the linear model
    tidied = map(fit, tidy) # make it look nice
  ) %>%
  unnest(tidied) 

methylation_coefs <- slice(model_results, seq(2, nrow(model_results), 2))
padj <- p.adjust(methylation_coefs$p.value, method = "BH")
methylation_coefs$padj <- padj

#regression matrix 
bipartite_matrix$ZNF551_estimate <-methylation_coefs$estimate
colnames(bipartite_matrix)[1]  <- "probeID"
colnames(bipartite_matrix)[2] <- "STRA6 estimate"