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
