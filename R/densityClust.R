#' Clustering by fast search and find of density peaks
#' 
#' This package implement the clustering algorithm described by Alex Rodriguez
#' and Alessandro Laio (2014). It provides the user with tools for generating 
#' the initial rho and delta values for each observation as well as using these
#' to assign observations to clusters. This is done in two passes so the user is
#' free to reassign observations to clusters using a new set of rho and delta 
#' thresholds, without needing to recalculate everything.
#' 
#' Plotting
#' 
#' Two types of plots are supported by this package, and both mimics the types of
#' plots used in the publication for the algorithm. The standard plot function
#' produces a decision plot, with optional colouring of cluster peaks if these 
#' are assigned. Furthermore \code{\link{plotMDS}} performs a multidimensional 
#' scaling of the distance matrix and plots this as a scatterplot. If clusters
#' are assigned observations are coloured according to their assignment.
#' 
#' Cluster detection
#' 
#' The two main functions for this package are \code{\link{densityClust}} and 
#' \code{\link{findClusters}}. The former takes a distance matrix and optionally
#' a distance cutoff and calculates rho and delta for each observation. The 
#' latter takes the output of \code{\link{densityClust}} and make cluster 
#' assignment for each observation based on a user defined rho and delta 
#' threshold. If the thresholds are not specified the user is able to supply 
#' them interactively by clicking on a decision plot.
#' 
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- densityClust(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{densityClust}}, \code{\link{findClusters}}, \code{\link{plotMDS}}
#' 
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @docType package
#' @name densityClust-package
#' 
#' @useDynLib densityClust
#' @importFrom Rcpp sourceCpp
NULL

#' Computes the local density of points in a distance matrix
#' 
#' This function takes a distance matrix and a distance cutoff and calculate the
#' local density for each point in the matrix. The computation can either be 
#' done using a simple summation of the points with the distance cutoff for each
#' observation, or by applying a gaussian kernel scaled by the distance cutoff 
#' (more robust for low-density data)
#' 
#' @param distance A distance matrix
#' 
#' @param dc A numeric value specifying the distance cutoff
#' 
#' @param gaussian Logical. Should a gaussian kernel be used to estimate the 
#' density (defaults to FALSE)
#' 
#' @return A vector of local density values, the index matching row and column 
#' indexes in the distance matrix
#' 
#' @noRd
#' 
localDensity <- function(distance, dc, gaussian=FALSE) {
   # These implementations are faster by virtue of being written in C++
   # They also avoid the need to convert `distance` to a matrix. 
   message("the length of the distance: ", length(distance))
  
   if(gaussian) {
      res <- gaussianLocalDensity(distance, attr(distance, "Size"), dc)
   } else {
      res <- nonGaussianLocalDensity(distance, attr(distance, "Size"), dc)
   }
   if(is.null(attr(distance, 'Labels'))) {
      names(res) <- NULL
   } else {
      names(res) <- attr(distance, 'Labels')
   }
   res
}
#' Calculate distance to closest observation of higher density
#' 
#' This function finds, for each observation, the minimum distance to an 
#' observation of higher local density.
#' 
#' @param distance A distance matrix
#' 
#' @param rho A vector of local density values as outputted by localDensity()
#' 
#' @return A vector of distances with index matching the index in rho
#' 
#' @noRd
#' 
distanceToPeak <- function(distance, rho) {
   # This implementation is faster by virtue of being written in C++.
   # It also avoids the need to convert `distance` to a matrix. 
   res <- distanceToPeakCpp(distance, rho);
   names(res) <- names(rho)
   res
}
#' Estimate the distance cutoff for a specified neighbor rate
#' 
#' This function calculates a distance cutoff value for a specific distance 
#' matrix that makes the average neighbor rate (number of points within the
#' distance cutoff value) fall between the provided range. The authors of the
#' algorithm suggests aiming for a neighbor rate between 1 and 2 percent, but 
#' also states that the algorithm is quite robust with regards to more extreme
#' cases.
#' 
#' @param distance A distance matrix
#' 
#' @param neighborRateLow The lower bound of the neighbor rate
#' 
#' @param neighborRateHigh The upper bound of the neighbor rate
#' 
#' @return A numeric value giving the estimated distance cutoff value
#' 
#' @examples
#' irisDist <- dist(iris[,1:4])
#' estimateDc(irisDist)
#' 
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @export
#' 
estimateDc <- function(distance, neighborRateLow=0.01, neighborRateHigh=0.02) {
   # This implementation uses binary search instead of linear search.
   
   size <- attr(distance, 'Size')
   # If size is greater than 448, there will be >100000 elements in the distance
   # object. Subsampling to 100000 elements will speed performance for very 
   # large dist objects while retaining good accuracy in estimating the cutoff
   if(size > 448) {
      distance <- distance[sample.int(length(distance), 100128)]
      size <- 448
   }
   
   low <- min(distance)
   high <- max(distance)
   dc <- 0
   while(TRUE) {
      dc <- (low + high) / 2
      # neighborRate = average of number of elements of comb per row that are 
      # less than dc minus 1 divided by size.
      # This implementation avoids converting `distance` to a matrix. The matrix is
      # symmetrical, so doubling the result from `distance` (half of the matrix) is
      # equivalent. The diagonal of the matrix will always be 0, so as long as dc 
      # is greater than 0, we add 1 for every element of the diagonal, which is
      # the same as size
      neighborRate <- (((sum(distance < dc) * 2 + (if (0 <= dc) size)) / size - 1)) / size 
      if(neighborRate >= neighborRateLow && neighborRate <= neighborRateHigh) break
      
      if(neighborRate < neighborRateLow) {
         low <- dc
      } else {
         high <- dc
      }
   }
   cat('Distance cutoff calculated to', dc, '\n')
   dc
}

#' Calculate clustering attributes based on the densityClust algorithm
#' 
#' This function takes a distance matrix and optionally a distance cutoff and 
#' calculates the values necessary for clustering based on the algorithm 
#' proposed by Alex Rodrigues and Alessandro Laio (see references). The actual 
#' assignment to clusters are done in a later step, based on user defined 
#' threshold values.
#' 
#' @details
#' The function calculates rho and delta for the observations in the provided 
#' distance matrix. If a distance cutoff is not provided this is first estimated
#' using \code{\link{estimateDc}} with default values.
#' 
#' The information kept in the densityCluster object is:
#' \describe{
#'   \item{rho}{A vector of local density values}
#'   \item{delta}{A vector of minimum distances to observations of higher density}
#'   \item{distance}{The initial distance matrix}
#'   \item{dc}{The distance cutoff used to calculate rho}
#'   \item{threshold}{A named vector specifying the threshold values for rho and delta used for cluster detection}
#'   \item{peaks}{A vector of indexes specifying the cluster center for each cluster}
#'   \item{clusters}{A vector of cluster affiliations for each observation. The clusters are referenced as indexes in the peaks vector}
#'   \item{halo}{A logical vector specifying for each observation if it is considered part of the halo}
#' }
#' Before running findClusters the threshold, peaks, clusters and halo data is 
#' NA.
#' 
#' @param distance A distance matrix
#' 
#' @param dc A distance cutoff for calculating the local density. If missing it
#' will be estimated with estimateDc(distance)
#' 
#' @param gaussian Logical. Should a gaussian kernel be used to estimate the 
#' density (defaults to FALSE)
#' 
#' @param verbose Logical. Should the running details be reported  
#'
#' @return A densityCluster object. See details for a description.
#' 
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- densityClust(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{estimateDc}}, \code{\link{findClusters}}
#' 
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @export
#' 
densityClust <- function(distance, dc, gaussian=FALSE, verbose = F, ...) {
  extra_arguments <- list(...)
  if(class(distance) %in% c('data.frame', 'matrix')) {
    dp_knn_args <- c(list(mat=distance, verbose = verbose),
                  extra_arguments[names(extra_arguments) %in% c("k", "use_dist")])
    #browser()
    res <- do.call(densityClust.knn, dp_knn_args)
    return(res)
  }
  else {
    if(missing(dc)) {
      if(verbose) {
          message('Calculating the distance cutoff')
      }
      dc <- estimateDc(distance)
    }
    if(verbose) {
      message('Calculating the local density for each sample based on distance cutoff')
    }
    rho <- localDensity(distance, dc, gaussian=gaussian)
    
    if(verbose) {
      message('Calculating the minimal distance of a sample to another sample with higher density')
    }
    delta <- distanceToPeak(distance, rho)
    
    if(verbose) {
      message('Returning result...')
    }
    res <- list(rho=rho, delta=delta, distance=distance, dc=dc, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA)
    class(res) <- 'densityCluster'
    res
  }
}
#' @export
#' @importFrom graphics plot points
#' @noRd
#' 
plot.densityCluster <- function(x, ...) {
    plot(x$rho, x$delta, main='Decision graph', xlab=expression(rho), ylab=expression(delta))
    if(!is.na(x$peaks[1])) {
        points(x$rho[x$peaks], x$delta[x$peaks], col=2:(1+length(x$peaks)), pch=19)
    }
}
#' Plot observations using multidimensional scaling and colour by cluster
#' 
#' This function produces an MDS scatterplot based on the distance matrix of the
#' densityCluster object, and, if clusters are defined, colours each observation
#' according to cluster affiliation. Observations belonging to a cluster core is
#' plotted with filled circles and observations belonging to the halo with
#' hollow circles.
#' 
#' @param x A densityCluster object as produced by \code{\link{densityClust}}
#' 
#' @param ... Additional parameters. Currently ignored
#' 
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- densityClust(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{densityClust}}
#' 
#' @export
#' 
plotMDS <- function (x, ...) {
    UseMethod("plotMDS", x)
}
#' @export
#' @importFrom stats cmdscale
#' @importFrom graphics plot points legend
#' @noRd
#' 
plotMDS.densityCluster <- function(x, ...) {
    mds <- cmdscale(x$distance)
    plot(mds[,1], mds[,2], xlab='', ylab='', main='MDS plot of observations')
    if(!is.na(x$peaks[1])) {
        for(i in 1:length(x$peaks)) {
            ind <- which(x$clusters == i)
            points(mds[ind, 1], mds[ind, 2], col=i+1, pch=ifelse(x$halo[ind], 1, 19))
        }
        legend('topright', legend=c('core', 'halo'), pch=c(19, 1), horiz=TRUE)
    }
}
#' @export
#' 
#' @noRd
#' 
print.densityCluster <- function(x, ...) {
    if(is.na(x$peaks[1])) {
        cat('A densityCluster object with no clusters defined\n\n')
        cat('Number of observations:', length(x$rho), '\n')
    } else {
        cat('A densityCluster object with', length(x$peaks), 'clusters defined\n\n')
        cat('Number of observations:', length(x$rho), '\n')
        cat('Observations in core:  ', sum(!x$halo), '\n\n')
        cat('Parameters:\n')
        cat('dc (distance cutoff)   rho threshold          delta threshold\n')
        cat(formatC(x$dc, width=-22), formatC(x$threshold[1], width=-22), x$threshold[2])
    }
}
#' Detect clusters in a densityCluster obejct
#' 
#' This function uses the supplied rho and delta thresholds to detect cluster
#' peaks and assign the rest of the observations to one of these clusters. 
#' Furthermore core/halo status is calculated. If either rho or delta threshold
#' is missing the user is presented with a decision plot where they are able to
#' click on the plot area to set the treshold. If either rho or delta is set, 
#' this takes presedence over the value found by clicking.
#' 
#' @param x A densityCluster object as produced by \code{\link{densityClust}}
#' 
#' @param ... Additional parameters passed on to \code{\link{findClusters.densityCluster}}
#' 
#' @return A densityCluster object with clusters assigned to all observations
#' 
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- densityClust(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{findClusters}}
#' 
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @export
#' 
findClusters <- function (x, ...) {
    UseMethod("findClusters", x)
}
#' @rdname findClusters
#' 
#' @param rho The threshold for local density when detecting cluster peaks
#' 
#' @param delta The threshold for minimum distance to higher density when detecting cluster peaks
#' 
#' @param plot Logical. Should a decision plot be shown after cluster detection
#' 
#' @param peaks A numeric vector indicates the index of density peaks used for clustering. This vector should be retrieved from the decision plot with caution. No checking involved.  
#'
#' @param verbose Logical. Should the running details be reported  
#'
#' @export
#' @importFrom graphics plot locator
findClusters.densityCluster <- function(x, rho, delta, plot=FALSE, peaks=NULL, verbose = F, ...) {
    if(class(x$distance) %in% c('data.frame', 'matrix')) {
      peak_ind <- which(x$rho > rho & x$delta > delta)
      x$peaks <- peak_ind

      knn_graph <- x$knn_graph
      shortest_dist <- shortest.paths(knn_graph, V(knn_graph), peak_ind)
      cluster <- apply(shortest_dist, 1, which.min)

      x$clusters <- factor(cluster)
      x$halo <- NULL #update this

      x$threshold['rho'] <- rho
      x$threshold['delta'] <- delta
    } 
    else {
      # Detect cluster peaks
      if(!is.null(peaks)) {
        if(verbose) {
          message('peaks are provided, clustering will be performed based on them')
        }
        x$peaks <- peaks
      } else {
        if(missing(rho) || missing(delta)) {
            x$peaks <- NA
            plot(x)
            cat('Click on plot to select thresholds\n')
            threshold <- locator(1)
            if(missing(rho)) rho <- threshold$x
            if(missing(delta)) delta <- threshold$y
            plot=TRUE
        }
        x$peaks <- which(x$rho > rho & x$delta > delta)
        x$threshold['rho'] <- rho
        x$threshold['delta'] <- delta
      }
      if(plot) {
          plot(x)
      }
      
      # Assign observations to clusters
      runOrder <- order(x$rho, decreasing = TRUE)
      cluster <- rep(NA, length(x$rho))
      if(verbose) {
        message('Assigning each sample to a cluster based on its nearest density peak')
      }
      for(i in runOrder) { 
        if((i %% round(length(runOrder) / 25)) == 0){
          if(verbose) {
            message(paste('the runOrder index is', i))
          }
        }

        if(i %in% x$peaks) {
          cluster[i] <- match(i, x$peaks)
        } else {
          higherDensity <- which(x$rho>x$rho[i])
          cluster[i] <- cluster[higherDensity[which.min(findDistValueByRowColInd(x$distance, attr(x$distance, 'Size'), i, higherDensity))]] 
        }
      }
      x$clusters <- cluster
      
      # Calculate core/halo status of observation
      border <- rep(0, length(x$peaks))
      if(verbose) {
        message('Identifying the border for each cluster')
      }
      for(i in 1:length(x$peaks)) {
          if(verbose) {
            message('the current index of the peak is ', i)
          }

          averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], '+')/2 
          index <- findDistValueByRowColInd(x$distance, attr(x$distance, 'Size'), which(cluster == i), which(cluster != i)) <= x$dc 
          if(any(index)) border[i] <- max(averageRho[index]) 
      }
      x$halo <- x$rho < border[cluster] 
    }
    x
}

#' Extract cluster membership from a densityCluster object
#' 
#' This function allows the user to extract the cluster membership of all the
#' observations in the given densityCluster object. The output can be formatted 
#' in two ways as described below. Halo observations can be chosen to be removed
#' from the output.
#' 
#' @details
#' Two formats for the output are available. Either a vector of integers 
#' denoting for each observation, which cluster the observation belongs to. If 
#' halo observations are removed, these are set to NA. The second format is a 
#' list with a vector for each group containing the index for the member 
#' observations in the group. If halo observations are removed their indexes are
#' omitted. The list format correspond to the following transform of the vector 
#' format \code{split(1:length(clusters), clusters)}, where \code{clusters} are 
#' the cluster information in vector format.
#' 
#' @param x The densityCluster object. \code{\link{findClusters}} must have
#' been performed prior to this call to avoid throwing an error.
#' 
#' @param ... Currently ignored
#' 
#' @return A vector or list with cluster memberships for the observations in the
#' initial distance matrix
#' 
#' @export
#' 
clusters <- function (x, ...) {
    UseMethod("clusters", x)
}
#' @rdname clusters
#' 
#' @param as.list Should the output be in the list format. Defaults to FALSE
#' 
#' @param halo.rm Logical. should halo observations be removed. Defaults to TRUE
#' 
#' @export
#' 
clusters.densityCluster <- function(x, as.list=FALSE, halo.rm=TRUE, ...) {
    if(!clustered(x)) stop('x must be clustered prior to cluster extraction')
    res <- x$clusters
    if(halo.rm) {
        res[x$halo] <- NA
    }
    if(as.list) {
        res <- split(1:length(res), res)
    }
    res
}

#' Check whether a densityCluster object have been clustered
#' 
#' This function checks whether \code{\link{findClusters}} has been performed on
#' the given object and returns a boolean depending on the outcome
#' 
#' @param x A densityCluster object
#' 
#' @return TRUE if \code{\link{findClusters}} have been performed, otherwise 
#' FALSE
#' 
#' @export
#' 
clustered <- function (x) {
    UseMethod("clustered", x)
}
#' @rdname clustered
#' 
#' @export
#' 
clustered.densityCluster <- function(x) {
    !any(is.na(x$peaks[1]), is.na(x$clusters[1]), is.na(x$halo[1]))
}

#' Extract labels
#' 
#' @noRd
#' 
#' @export
#' 
labels.densityCluster <- function(object, ...) {
    labels(object$distance)
}

#' Extract labels
#' 
#' @noRd
#' 
#' @export
#' 
densityClust.knn <- function(mat, k = 30, use_dist = T, verbose = F, ...) {
  if(verbose) {
    message('Calculating kNN graph')
  }
 
  dx <- FNN::get.knn(mat, k = k)

  nn.index <- dx$nn.index
  nn.dist <- dx$nn.dist
  N <- nrow(nn.index)

  edges <- reshape2::melt(t(nn.index)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
  edges_weight <- reshape2::melt(t(nn.dist)); #colnames(edges_weight) = c("B", "A", "C"); edges_weight = edges_weight[,c("A","B","C")]
  edges$B <- edges$C; 
  if(use_dist) 
    edges$C <- edges_weight$value
  else
    edges$C <- 1

  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

  Adj <- Matrix::sparseMatrix(i = c(edges$A, edges$B), j = c(edges$B, edges$A), x = c(edges$C, edges$C), dims = c(N, N))

  knn_graph <- graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = T)

  if(verbose) {
    message('Calculating the local density for each sample based on distance cutoff')
  }

  rho <- apply(dx$nn.dist, 1, function(x) {
      exp(-mean(x))
    })
  
  if(verbose) {
    message('Calculating the minimal distance of a sample to another sample with higher density')
  }
  
  rho_order <- order(rho)

  delta <- vector(mode = 'integer', length = N)
  for(i in 1:N) {
    dis <- igraph::shortest.paths(knn_graph, rho_order[i], rho_order[-c(1:i)])
    delta[rho_order[i]] <- min(dis)
  }

  delta[which(is.infinite(delta))] <- max(delta[which(!is.infinite(delta))]) #set the distance for the cell with highest density
 
  if(verbose) {
    message('Returning result...')
  }
  res <- list(rho=rho, delta=delta, distance=mat, dc=NULL, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA, knn_graph = knn_graph)
  class(res) <- 'densityCluster'
  res
}
