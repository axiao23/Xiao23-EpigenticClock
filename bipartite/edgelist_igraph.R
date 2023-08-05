# clean nested loops and regressions for adjacency matrix
library(tidyverse)
library(broom)

tcga_dat= read.csv("../data/TCGA_LUAD_EXPRESSION_METHYLATION_standardized_data.csv")

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
    output[i,j] = summary(model1)$coef[2,3]
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
write.csv(output2, "../data/regression_t_stat_standardized_TCGA_LUAD.csv", row.names = F)


#creating an edgelist via igraph
library(igraph)
# creating a square matrix 
adj_matrix= read.csv("../data/regression_standardized_TCGA_LUAD.csv")
rownames(adj_matrix) = genes
df <- data.frame(matrix(ncol = 1796, nrow = 1555))
df1 <- data.frame(matrix(ncol = 241, nrow = 1796))
df2 <- data.frame(matrix(ncol = 2037, nrow = 241))
colnames(df) = colnames(adj_matrix)
rownames(df1) = rownames(adj_matrix)
colnames(df2) = colnames(adj_matrix2)
adj_matrix1 = rbind(adj_matrix,df)
adj_matrix2 = cbind(adj_matrix1,df1)
adj_matrix3 = rbind(adj_matrix2,df2)
adj_matrix3[is.na(adj_matrix3)] <- 0
adj_matrix3[1797:2037,242:2037] <- adj_matrix
adj_matrix3 = as.matrix(adj_matrix3)
igraph_regression<- graph_from_adjacency_matrix(adj_matrix3,
                            mode = "max",
                            weighted = TRUE,
                            diag = TRUE,
                            add.colnames = TRUE,
                            add.rownames = NA
)
dim(adj_matrix3)
#create an edgelist from igraph
edgelist_TCGA_LUAD <- as_edgelist(igraph_regression, names = FALSE)

