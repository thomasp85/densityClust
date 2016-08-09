library(devtools)
load_all('~/Dropbox (Personal)/Projects/monocle-dev/')
setwd("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/densityClust")
load_all("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/densityClust")

irisDist <- dist(iris[,1:4])
irisClust <- densityClust(irisDist, gaussian=F)
plot(irisClust) # Inspect clustering attributes to define thresholds

irisClust <- findClusters(irisClust, rho=2, delta=2)
plotMDS(irisClust)
split(iris[,5], irisClust$clusters)

# library(proxy)
load("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/single-cell-3prime-paper/tsne_res")
max_components=2
tsne_data <- tsne_res$Y[, 1:max_components]
row.names(tsne_data) <- colnames(tsne_data)
tsne_data <- tsne_res$Y[, 1:max_components]
row.names(tsne_data) <- colnames(tsne_data)
dataDist <- dist(tsne_data[, ], upper = T, diag = T) #1:34000
#debug(densityClust)
dataClust <- densityClust::densityClust(dataDist, gaussian = F)

dataClust <- densityClust::findClusters(dataClust)

#the other approach: 
dataClust <- densityClust::findClusters(dataClust, rho = quantile(dataClust$rho, probs = 0.95), delta = quantile(dataClust$delta, probs = 0.95))

#
2147000000

adpclust(tsne_data)