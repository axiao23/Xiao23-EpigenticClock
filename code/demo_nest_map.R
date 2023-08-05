library(tidyverse)
library(broom)

tcga_dat = read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
dim(tcga_dat)

#nested for loops here
cg_indices = grep("cg",names(tcga_dat))
length(cg_indices)
gene_indices = setdiff(2:ncol(tcga_dat),cg_indices)
length(gene_indices)
output = matrix(nrow=length(gene_indices), ncol=length(cg_indices))
#output2 = matrix(nrow=length(gene_indices), ncol=length(cg_indices))
#lmgenes <- data.frame(matrix(NA,    # Create empty data frame
                          #nrow = 2,
                          #ncol = 9))
#counter = 0
for(i in 1:length(gene_indices))
{
  for(j in 1:length(cg_indices))
  {
    counter = counter + 1
    #print(paste(c("Gene:", names(tcga_dat)[gene_indices[i]], "cg:",names(tcga_dat)[cg_indices[j]])))
    model1 = lm(tcga_dat[,gene_indices[i]] ~ tcga_dat[,cg_indices[j]], data = tcga_dat)
    #print(model1)
    #lmgenes[,counter] = coef(lm(tcga_dat[,gene_indices[i]] ~ tcga_dat[,cg_indices[j]], data = tcga_dat))
    #colnames(lmgenes) <- genes[1:3]
    #print(lmgenes)
    #beta = coef(lmgenes)
    #model1 = data.frame(beta)
    # demo storing correlations
    output[i,j] = summary(model1)$coef[2,1]
    #output2[i,j] = cor(tcga_dat[,gene_indices[i]], tcga_dat[,cg_indices[j]])
  
    # run the linear model with lm
    
    # store the coefficient for methylation in one matrix
    
    # store the p-value for methylation in another matrix
  }
}
output1 <- output[1:241,]
output2 <- data.frame(output1)
tcga_cgs <- tcga_meth[c(3:1798)]
colnames(output2) <- colnames(tcga_cgs)
rownames(output2)<- genes
write.csv(output2, "../data/regression_TCGA_LUAD.csv", row.names = F)


colnames(output) <-  
output$gene <- genes[1:241]
lmgenes1 <- lmgenes[2,]
lmgenes1[2,] <- lmgenes1[c(2,5,8)]
lmgenes1[3,] <- lmgenes1[c(3,6,9)]
lmgenes1 <- lmgenes1[c(1,4,7)]
colnames(lmgenes1) <- genes[1:3]
rownames(lmgenes1) <- tcga_meth_long$probeID[1:3]


# entire dataset nested loop
cg_indices = grep("cg",names(tcga_dat))
length(cg_indices)
gene_indices = setdiff(2:ncol(tcga_dat),cg_indices)
length(gene_indices)
#output = matrix(nrow=length(gene_indices), ncol=length(cg_indices))
#output2 = matrix(nrow=length(gene_indices), ncol=length(cg_indices))
lmgenes4 <- data.frame(matrix(NA,    # Create empty data frame
                             nrow = 2,
                             ncol = 243))
counter = 0
for(i in 1:length(gene_indices[201:243]))
{
  for(j in 1:length(cg_indices))
  {
    counter = counter + 1
    #print(paste(c("Gene:", names(tcga_dat)[gene_indices[i]], "cg:",names(tcga_dat)[cg_indices[j]])))
    #model1 = coef(lm(tcga_dat[,gene_indices[i]] ~ tcga_dat[,cg_indices[j]], data = tcga_dat))
    #print(model1)
    lmgenes4[,counter] = coef(lm(tcga_dat[,gene_indices[i]] ~ tcga_dat[,cg_indices[j]], data = tcga_dat))
    #colnames(lmgenes) <- genes[1:3]
    #print(lmgenes)
    #beta = coef(lmgenes)
    #model1 = data.frame(beta)
    # demo storing correlations
    #output[i,j] = coef(tcga_dat[,gene_indices[i]], tcga_dat[,cg_indices[j]])
    #output2[i,j] = cor(tcga_dat[,gene_indices[i]], tcga_dat[,cg_indices[j]])
    
    # run the linear model with lm
    
    # store the coefficient for methylation in one matrix
    
    # store the p-value for methylation in another matrix
  }
}

colnames(output) = 
lmgenes1 <- lmgenes[2,]
lmgenes1[2,] <- lmgenes1[c(2,5,8)]
lmgenes1[3,] <- lmgenes1[c(3,6,9)]
lmgenes1 <- lmgenes1[c(1,4,7)]
colnames(lmgenes1) <- genes[1:3]
rownames(lmgenes1) <- tcga_meth_long$probeID[1:3]

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
for (i in tcga_meth_long[c(2:3)]){
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

