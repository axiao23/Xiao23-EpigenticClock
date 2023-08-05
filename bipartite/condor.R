#condor and communities
library(netZooR)
elist = read.csv("../data/edgelist_TCGA_LUAD_data.csv")
#gtoy = graph.edgelist(as.matrix(elist),directed=FALSE)
#plot(gtoy,vertex.label.dist=2)
condor.object <- createCondorObject(elist)
names(condor.object)
condor.object <- condorCluster(condor.object)
