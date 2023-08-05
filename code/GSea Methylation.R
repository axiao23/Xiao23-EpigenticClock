#transpose and merging of the contributions and the corresponding genes
Horvath_gene_methylation <- merge(Horvath_outputs, coefHorvath, by = "CpGmarker", all.x = T, all.y = F)
View(Horvath_gene_methylation)
Horvath_gene_methylation_values <- Horvath_outputs[-c(1),-c(1,2)]

#gene enrichment analysis methylation 
library(fgsea)
library(tidyverse)
library(ggplot2)
pathways_bp = gmtPathways("c5.go.bp.v7.5.1.symbols.gmt")
coefs = na.omit(Horvath_gene_methylation$CoefficientTraining.x)
names(coefs) = Horvath_gene_methylation$Symbol
gseaRes=fgsea(pathways_bp,coefs)
fgseaMultilevelRes <- fgseaMultilevel(pathways_bp, coefs, maxSize=500)
horvath_complete <- Horvath_gene_methylation[!is.na(Horvath_gene_methylation$CoefficientTraining.x),]
coefs2 = horvath_complete$CoefficientTraining.x
names(coefs2) = horvath_complete$Symbol
gseaRes1= fgsea(pathways_bp,coefs2)

B#visualization 
library(DOSE)
library(enrichplot)
setMethod("dotplot", signature(object = "enrichResult"),
          function(object, x = "GeneRatio", color = "p.adjust",
                   showCategory=10, size = NULL,
                   split = NULL, font.size=12, title = "",
                   label_format = 30, ...) {
            dotplot_internal(object, x, color, showCategory, size,
                             split, font.size, title, label_format, ...)
          })
df <- fortify(gseaRes[524:525,], showCategory = showCategory, split=split)
idx <- order(df[['pval']], decreasing = TRUE)
#df$pval <- factor(df$pval,
                         #levels=rev(unique(df$pval[idx])))
ggplot(df, aes(x=pathway, y=pval, size = size)) +
  geom_point() + coord_flip()
       
  

convert_data <- gseaRes[['pval']]
idx <- order(gseaRes[['pval']], decreasing = TRUE)
gseaRes$pval <- factor(gseaRes$pval,
                         levels=rev(unique(gseaRes$pval[idx])))
ggplot(gseaRes, aes_string(x = "pathway", y= "pval", size=5)) +
  geom_dotplot() 

barplot(edo, showCategory=20)
head(gseaRes[order(pval), ])
#top pathways
topPathwaysUp <- gseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- gseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways <- gseaRes[head(order(pval), n=15)][order(NES), pathway]
plotGseaTable(pathways_bp[topPathwaysDown],coefs,gseaRes, gseaParam=0.1)

# -------------------------------------------------------------------------


#collapsed and main pathways
collapsedPathways <- collapsePathways(gseaRes[order(pval)][padj < 0.01], 
                                      pathways_bp, coefs)
mainPathways <- gseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways_bp[mainPathways], coefs, gseaRes, 
              gseaParam = 0.5)

#separate Enrichment plots
plotEnrichment(pathways_bp[["GOBP_DNA_METHYLATION"]],coefs)

#visualization via bubbleplot
df <-gseaRes[order(gseaRes$pval),]
df <- df[1:30,]








