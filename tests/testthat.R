# library(testthat)
# library(densityClust)
# library(ggplot2)
# # library(devtools)
# library(monocle)
# library(igraph)
# 
# load_all()
# 
# # test_check("densityClust")
# master_tsne <- readRDS('./tests/master_tsne.tsneresults.RDS')
# 
# dx <- FNN::get.knn(master_tsne$Y, k = 150)
# nn.index <- dx$nn.index
# nn.dist <- dx$nn.dist
# 
# rho <- apply(dx$nn.dist, 1, function(x) {
#   exp(-mean(x))
# })
# 
# qplot(master_tsne$Y[, 1], master_tsne$Y[, 2], color = rho)
# rho_order <- order(rho)
# 
# sample_ind <- sample(1:nrow(master_tsne$Y), 10000)
# delta_val <- smallest_dist_rho_order_coords(rho[rho_order[sample_ind]], master_tsne$Y[rho_order[sample_ind], ])
# 
# qplot(master_tsne$Y[rho_order[sample_ind], 1], master_tsne$Y[rho_order[sample_ind], 2], color = rho[rho_order[sample_ind]])
# qplot(master_tsne$Y[sample_ind, 1], master_tsne$Y[sample_ind, 2])
# 
# qplot(rho[rho_order[sample_ind]], delta_val, color = rho[rho_order[sample_ind]])
# 
# start_time <- Sys.time()
# master_tsneClust <- densityClust(master_tsne$Y[, ], gaussian=TRUE, k = 150)
# end_time <- Sys.time() 
# start_time - end_time 
# plot(master_tsneClust)
# 
# start_time <- Sys.time()
# master_tsneClust <- findClusters(master_tsneClust, rho=0.6, delta=5)
# end_time <- Sys.time() 
# start_time - end_time 
# 
# qplot(master_tsne$Y[, 1], master_tsne$Y[, 2], color = master_tsneClust$clusters)
# 
# start_time <- Sys.time()
# sample_master_tsneClust <- densityClust(master_tsne$Y, gaussian=TRUE, k = 150)
# end_time <- Sys.time() 
# start_time - end_time
# 
# test on lung dataset:

# qplot(lung@reducedDimA[1, ], lung@reducedDimA[2, ], color = pData(lung)$Time)
# dx <- FNN::get.knn(t(lung@reducedDimA), k = 15)
# nn.index <- dx$nn.index
# nn.dist <- dx$nn.dist
# 
# rho <- apply(dx$nn.dist, 1, function(x) {
#   exp(-mean(x))
# })
# 
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = rho)
# rho_order <- order(rho)
# 
# delta_val <- smallest_dist_rho_order_coords(rho[rho_order], master_tsne$Y[rho_order, ])
# qplot(rho[rho_order], delta_val)
# 
# irisDist <- dist(iris[,1:4])
# lungClust <- densityClust(t(lung@reducedDimA), gaussian=TRUE, k = 15)
# plot(lungClust) # Inspect clustering attributes to define thresholds
# 
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$delta > 3)
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$delta > 2.5)
# 
# lungClust <- findClusters(lungClust, rho=0.3, delta=10)
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$clusters)
# 
# plotMDS(irisClust)
# split(iris[,5], irisClust$clusters)
# 
# #
# filenames = c('Circle','two_moon','tree_300','Spiral_SUN','three_clusters','DistortedSShape')
# data_dir <- '~/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/toy/'
# X <- readMat(paste(data_dir, 'three_clusters', '.mat', sep =''))$X
# qplot(X[, 1], X[, 2])
# clusterClust <- densityClust(X, gaussian=TRUE)
# plot(clusterClust) # Inspect clustering attributes to define thresholds
# 
# qplot(X[, 1], X[, 2], color = clusterClust$delta > 6)
# 
# #artifical clusters: 
# X <- X[1:15, ]
# X[1:5, 1] <- 1 + c(-0.01, -0.02, 0, 0.01, 0.02) 
# X[1:5, 2] <- 1 + c(-0.01, -0.02, 0, 0.01, 0.02) 
# 
# X[6:10, 1] <- -1 + c(-0.01, -0.02, 0, 0.01, 0.02) 
# X[6:10, 2] <- 1 + c(-0.01, -0.02, 0, 0.01, 0.02)
# 
# X[11:15, 1] <- -1 + c(-0.01, -0.02, 0, 0.01, 0.02) 
# X[11:15, 2] <- -1 + c(-0.01, -0.02, 0, 0.01, 0.02)
# 
# clusterClust <- densityClust(X, gaussian=TRUE)
# plot(clusterClust) # Inspect clustering attributes to define thresholds
# 
# qplot(X[, 1], X[, 2], color = clusterClust$delta > 1.4)
# 
# 
# 
# irisClust <- densityClust(iris[,1:4], verbose=TRUE, k = 15)
# plot(irisClust)
# irisClust <- findClusters(irisClust, rho=0.7, delta=2)
# plotMDS(irisClust)
# 
# qplot(irisClust$distance[, 1], irisClust$distance[, 2], color = irisClust$clusters)
# #test why there is still NA values for the cluster assignment for certain cells: 
# NA_cells <- which(is.na(irisClust$clusters))
# 
# irisClust$nearest_higher_density_neighbor[NA_cells]
# 
# # lung <- load_lung()
# lung <- reduceDimension(lung, reduction_method='tSNE', num_dim = 3)
# lungClust <- densityClust(t(lung@reducedDimA), gaussian=TRUE, k = 10)
# plot(lungClust) # Inspect clustering attributes to define thresholds
# 
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$delta > 3)
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$delta > 2.5)
# 
# lungClust <- findClusters(lungClust, rho=0.3, delta=8)
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$clusters)
# plotMDS(lungClust)
# 
# lungClust <- densityClust(dist(t(lung@reducedDimA)), gaussian=TRUE, k = 7)
# plot(lungClust) # Inspect clustering attributes to define thresholds
# lungClust <- findClusters(lungClust, rho=0.3, delta=15)
# qplot(t(lung@reducedDimA)[, 1], t(lung@reducedDimA)[, 2], color = lungClust$clusters)
# plotMDS(lungClust)

# irisClust <- densityClust(iris[,1:4], verbose=TRUE, k = 5)
# plot(irisClust) # Inspect clustering attributes to define thresholds
# irisClust <- findClusters(irisClust, rho=0.7, delta=2)
# plotMDS(irisClust)




