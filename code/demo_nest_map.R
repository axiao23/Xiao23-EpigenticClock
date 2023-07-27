library(tidyverse)
library(broom)

tcga_dat = read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
dim(tcga_dat)

cg_indices = grep("cg",names(tcga_dat))
gene_indices = setdiff(2:ncol(tcga_dat),cg_indices)

output = matrix(nrow=length(gene_indices), ncol=length(cg_indices))
for(i in 1:length(gene_indices[1:3]))
{
  for(j in 1:length(cg_indices[1:3]))
  {
    print(paste(c("Gene:", names(tcga_dat)[gene_indices[i]], "cg:",names(tcga_dat)[cg_indices[j]])))
    # demo storing correlations
    output[i,j] = cor(tcga_dat[,gene_indices[i]], tcga_dat[,cg_indices[j]])
    
    # run the linear model with lm
    
    # store the coefficient for methylation in one matrix
    
    # store the p-value for methylation in another matrix
  }
}

## ALL THE CODE BELOW IS FOR NEST/MAP
## USE THIS LATER, ONCE YOU ARE COMFORTABLE WITH THE DOUBLE FOR LOOP 

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
