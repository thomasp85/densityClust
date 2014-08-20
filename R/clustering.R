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
#' @export
#' 
findClusters.densityCluster <- function(x, rho, delta, plot=FALSE, ...) {
  # Detect cluster peaks
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
  
  if(plot) {
    plot(x)
  }
  
  # Assign observations to clusters
  comb <- as.matrix(x$distance)
  runOrder <- order(x$rho, decreasing = TRUE)
  cluster <- rep(NA, length(x$rho))
  for(i in runOrder) {
    if(i %in% x$peaks) {
      cluster[i] <- match(i, x$peaks)
    } else {
      higherDensity <- which(x$rho>x$rho[i])
      cluster[i] <- cluster[higherDensity[which.min(comb[i, higherDensity])]]
    }
  }
  x$clusters <- cluster
  
  # Calculate core/halo status of observation
  border <- rep(0, length(x$peaks))
  for(i in 1:length(x$peaks)) {
    averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], function(X, Y){X})
    index <- comb[cluster == i, cluster != i] <= x$dc
    if(any(index)) border[i] <- max(averageRho[index])
  }
  x$halo <- x$rho < border[cluster]
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
