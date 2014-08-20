#' Calculate clusters in a densityCluster object
#' 
#' This function uses the supplied \eqn{\rho} and \eqn{\delta} thresholds to 
#' calculate cluster centers and assigns the rest of the points to one of these 
#' clusters, designating these points as belonging to either the core or the 
#' halo of the cluster. If either of the \eqn{\rho} or \eqn{\delta} thresholds 
#' are missing the user is presented with a decision plot where they are able to
#' click on the plot area to set the threshold. If either of the \eqn{\rho} or 
#' \eqn{\delta} thresholds is passed as an argument then this takes precedence
#' over the value set by clicking.
#' 
#' @param x a \code{densityCluster} object as produced by \code{\link{dclust}}.
#' @param ... additional parameters passed to pass on.
#'   
#' @return A \code{densityCluster} object with clusters assigned to all 
#'   observations.
#'   
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- dclust(irisDist, gaussian=TRUE)
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' @export
findClusters <- function (x, ...) {
  UseMethod("findClusters", x)
}


#' @rdname findClusters
#'   
#' @param rho the threshold for local density.
#' @param delta the threshold for minimum distance to higher density.
#' @param plot logical indicating whether a decision plot should be shown after
#'   cluster calculation.
#'   
#' @export
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
#' observations in the given \code{densityCluster} object. The output can be
#' formatted in two ways as described below. Halo observations can optionally be
#' removed from the output.
#' 
#' @details Two formats for the output are available, with the first being a vector of
#' integers denoting for each point which cluster the point belongs
#' to. If halo observations are removed, these are set to NA. The second format
#' is a list with a vector for each group containing the index for the member 
#' points in the group. If halo points are removed their indexes are
#' omitted. The list format corresponds to the transformation of the vector 
#' format \code{split(1:length(clusters), clusters)}, where \code{clusters} is 
#' the cluster information in vector format.
#' 
#' @param x the \code{densityCluster} object to inspect.
#' @param ... extra arguments to pass on.
#'   
#' @return A vector or list with cluster memberships for the observations in the
#'   initial distance matrix
#'   
#' @export
#' 
clusters <- function (x, ...) {
  UseMethod("clusters", x)
}


#' @rdname clusters
#' 
#' @param as.list logical indicating that the output should be in list format.
#' @param halo.rm logical indicating that halo points should be removed.
#' 
#' @export
clusters.densityCluster <- function(x, as.list=FALSE, halo.rm=TRUE, ...) {
  if(!is.clustered(x)) stop('Input data must be clustered prior to cluster extraction!')
  res <- x$clusters
  if(halo.rm) {
    res[x$halo] <- NA
  }
  if(as.list) {
    res <- split(1:length(res), res)
  }
  res
}


#' Check whether a densityCluster object has been clustered
#' 
#' @param x the \code{densityCluster} object to check.
#'
#' @return \code{TRUE} if the object has been clustered, otherwise \code{FALSE}.
#'   
#' @export
is.clustered <- function (x) {
  UseMethod("is.clustered", x)
}


#' @rdname is.clustered
#' @export
is.clustered.densityCluster <- function(x) {
  !any(is.na(x$peaks[1]), is.na(x$clusters[1]), is.na(x$halo[1]))
}


#' Extract labels
#' 
#' @noRd
#' 
#' @export
labels.densityCluster <- function(object, ...) {
  labels(object$distance)
}
