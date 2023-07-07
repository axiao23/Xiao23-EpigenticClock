#transpose and merging of the contributions and the corresponding genes
Horvath_genes_methylation <- merge(Horvath_outputs, coefHorvath, by = "CpGmarker", all.x = T, all.y = F)
View(Horvath_gene_methylation)
Horvath_gene_methylation_values <- Horvath_outputs[-c(1),-c(1,2)]


#gene enrichment analysis methylation 
library(fgsea)
library(tidyverse)
pathways_bp = gmtPathways("c5.go.bp.v7.5.1.symbols.gmt")
coefs = na.omit(Horvath_genes_methylation$CoefficientTraining.x)
names(coefs) = Horvath_genes_methylation$Symbol
gseaRes=fgsea(pathways_bp,coefs)

#plotting results(top)
topPathways <- gseaRes[head(order(pval), n=15)][order(NES), pathway]
plotGseaTable(pathways_bp[topPathways],coefs,gseaRes)


#plotting bottom
plotEnrichment(pathways_bp[["GOBP_METHYLATION"]],coefs)

