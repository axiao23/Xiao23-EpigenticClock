#condor and communities
library(netZooR)
install.packages("devtools")
library(devtools)
devtools::install_github("jplatig/condor")
library(condor)
elist = read.csv("../data/edgelist_TCGA_LUAD_data.csv")

condor.object <- createCondorObject(elist)
names(condor.object)
condor.object <- condorCluster(condor.object)
print(condor.object$red.memb)
print(condor.object$blue.memb)
condor.object <- condorQscore(condor.object)






#bipartite graph 
g <- graph.data.frame(elist, directed = F)
V(g)$type <- V(g)$name %in% elist[,2] #the second column of edges is TRUE type
E(g)$weight <- as.numeric(elist[,3])
V(g)$color <- V(g)$type
V(g)$color=gsub("FALSE","red",V(g)$color)
V(g)$color=gsub("TRUE","blue",V(g)$color)
plot(g, edge.color="gray30",edge.width=E(g)$weight, layout=layout_as_bipartite)

#CONDOR ANALYSIS
#enrich
q_methylation <- condor.object$qscores$red.qscore
core_stats <- condorCoreEnrich(test_nodes=c("cg23180365","cg00702638"),
                               q=q_methylation,perm=TRUE,plot.hist=TRUE)

q_genes <- condor.object$qscores$blue.qscore
core_stat_1 <- condorCoreEnrich(test_nodes=c("TUBB2B", "BAZ2A"),
                                q=q_genes,perm=TRUE,plot.hist=TRUE)


#heatmap
condorPlotHeatmap(condor.object)

#barplot
mod <- sort(condor.object$modularity, decreasing=T)
barplot(mod, horiz=T, main = 'Modularity',  xlab = 'Modularity',cex.names=0.4, xlim = c(0,0.15), col='orchid')

#communities
condorPlotCommunities(
  condor.object,
  color_list=c("darkgreen","darkorange", "darkviolet"),
  point.size = .01,
  xlab = "Methylation",
  ylab = "Gene"
)
