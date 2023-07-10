#transpose and merging of the contributions and the corresponding genes
Horvath_genes_methylation <- merge(Horvath_outputs, coefHorvath, by = "CpGmarker", all.x = T, all.y = F)
View(Horvath_gene_methylation)
Horvath_gene_methylation_values <- Horvath_outputs[-c(1),-c(1,2)]

#gene enrichment analysis methylation 
library(fgsea)
library(tidyverse)
library(ggplot2)
pathways_bp = gmtPathways("c5.go.bp.v7.5.1.symbols.gmt")
coefs = na.omit(Horvath_genes_methylation$CoefficientTraining.x)
names(coefs) = Horvath_genes_methylation$Symbol
gseaRes=fgsea(pathways_bp,coefs)
fgseaMultilevelRes <- fgseaMultilevel(pathways_bp, coefs, maxSize=500)

#visualization 


head(gseaRes[order(pval), ])
#top pathways
topPathwaysUp <- gseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- gseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways <- gseaRes[head(order(pval), n=15)][order(NES), pathway]
plotGseaTable(pathways_bp[topPathways],coefs,gseaRes, gseaParam=0.5)

#collapsed and main pathways
collapsedPathways <- collapsePathways(gseaRes[order(pval)][padj < 0.01], 
                                      pathways_bp, coefs)
mainPathways <- gseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways_bp[mainPathways], coefs, gseaRes, 
              gseaParam = 0.5)

#separate Enrichment plots
plotEnrichment(pathways_bp[["GOBP_METHYLATION"]],coefs)


#new code 







