library(tidyverse)
library(broom)

tcga_dat = read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
dim(tcga_dat)

cg_indices = grep("cg",names(tcga_dat))

# first, select just the methylation sites and TUBB2B expression
tcga_meth= tcga_dat %>% select(submitter_id,TUBB2B, all_of(cg_indices))
dim(tcga_meth)

# now we are going to use a function called pivot_longer to convert the data so that
# there is one column that has all the cg probes in it

tcga_meth_long = tcga_meth %>% pivot_longer(cols = contains("cg"),
                                                names_to = "probeID", values_to="methylation")

model_results= tcga_meth_long %>%
  nest(data = -probeID) %>% # tell the code to work each probe
    mutate(
      fit = map(data, ~  lm(TUBB2B ~ methylation, data = .x)), # run the linear model
      tidied = map(fit, tidy) # make it look nice
    ) %>%
    unnest(tidied) 

# compare to one linear model just to show that we
# are getting the same results
lm(tcga_dat$TUBB2B ~ tcga_dat$cg02331561)

#graphing bipartite
library(igraph)
library(reshape2)
library(dplyr)

methylation_coefs <- slice(model_results, seq(2, nrow(model_results), 2))
padj <- p.adjust(methylation_coefs$p.value, method = "BH")
methylation_coefs$padj <- padj
df2 <- methylation_coefs[order(methylation_coefs$padj,decreasing=FALSE),]
methylation_sigvalues <- df2[1:555,]
linear_coefs <- data.frame(methylation_sigvalues$estimate)
g <- graph.incidence(linear_coefs, weighted = T)
V(g)$type <- bipartite_mapping(g)$type
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
E(g)$color <- "lightgray"
plot(g, layout = layout_as_bipartite)

#histogram 
padj <- p.adjust(methylation_coefs$p.value, method = "BH")
methylation_coefs$padj <- padj
df2 <- methylation_coefs[order(methylation_coefs$padj,decreasing=FALSE),]
methylation_sigvalues <- df2[1:555,]
ggplot(methylation_sigvalues, aes_string(x = "padj", size=5)) +
  geom_histogram()

#entire dataset linear regressions
tcga_dat = read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
dim(tcga_dat)

cg_indices = grep("cg",names(tcga_dat))
genes <- colnames(tcga_dat)[2:242]

# first, select just the methylation sites and TUBB2B expression
tcga_meth= tcga_dat %>% select(submitter_id,genes, all_of(cg_indices))
dim(tcga_meth)

# now we are going to use a function called pivot_longer to convert the data so that
# there is one column that has all the cg probes in it

tcga_meth_long = tcga_meth %>% pivot_longer(cols = contains("cg"),
                                            names_to = "probeID", values_to="methylation")

results_list = list()
counter = 1
for (i in colnames(tcga_meth_long[c(2:3)])){
  print(i)
  counter = counter + 1
  model1 = tcga_meth_long %>%
    select(submitter_id, probeID, i, methylation) %>%
  nest(data = -probeID) %>%
  mutate(
  fit = map(data, ~ lm(i ~ methylation, data = .x)),
    tidied = map(fit,tidy)
  )%>%
  unnest(tidied)
  results_list[[counter]] = model1%>% filter(term == "methylation")%>%
    select(-c(data,fit)) %>%
    mutate(geneName = i)
}
final_df = Reduce(rbind.data.frame, results_list)

#one regression padjust    
methylation_coefs <- slice(model_results, seq(2, nrow(model_results), 2))
padj <- p.adjust(methylation_coefs$p.value, method = "BH")
methylation_coefs$padj <- padj


