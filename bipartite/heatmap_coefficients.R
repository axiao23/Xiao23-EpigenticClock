#graphing heatmap 
library(readr)
regression_TCGA_LUAD <- read_csv("regression_TCGA_LUAD.csv")
View(regression_TCGA_LUAD)
library("gplots")
mymat = as.matrix(regression_TCGA_LUAD[1:10,1:10])
heatmap.2(mymat,
          Rowv=T,
          Colv=T,
          trace="none",
          col = "redblue")
