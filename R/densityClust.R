#' Computes the local density of points in a distance matrix
#' 
#' This function takes a distance matrix and a distance cutoff and calculates
#' the local density for each point in the matrix. The computation can either be
#' done using a simple summation of the points within the distance cutoff for
#' each point, or by applying a Gaussian kernel scaled by the distance cutoff 
#' (more robust for low-density data).
#' 
#' @param distance a distance matrix.
#' @param dc a numeric value specifying the distance cutoff.
#' @param gaussian logical indicating that a Gaussian kernel should be used to
#'   estimate the density.
#'   
#' @return A vector of local density values, one for each point.
#'   
#' @export
localDensity <- function(distance, dc, gaussian=FALSE) {
    comb <- as.matrix(distance)
    if(gaussian) {
        res <- apply(exp(-(comb/dc)^2), 1, sum)-1
    } else {
        res <- apply(comb < dc, 1, sum)-1
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
#' @param distance a distance matrix.
#' @param rho a vector of local density values as outputted by \code{\link{localDensity}}.
#' 
#' @return A vector of distances with index matching the index in \code{rho}.
#' 
#' @export 
distanceToPeak <- function(distance, rho) {
    comb <- as.matrix(distance)
    res <- sapply(1:length(rho), function(i) {
        peaks <- comb[rho>rho[i], i]
        if(length(peaks) == 0) {
            max(comb[,i])
        } else {
            min(peaks)
        }
    })
    names(res) <- names(rho)
    res
}


#' Estimate the distance cutoff for a specified neighbor rate
#' 
#' This function calculates a distance cutoff value for a specific distance 
#' matrix that makes the average neighbor rate (number of points within the 
#' distance cutoff value) fall between the provided range. The authors of the 
#' algorithm suggests aiming for a neighbor rate between 1 and 2 percent, but 
#' also state that the algorithm is quite robust with regards to more extreme 
#' cases.
#' 
#' @param distance a distance matrix.
#' @param neighborRateLow the lower bound of the neighbor rate.
#' @param neighborRateHigh the upper bound of the neighbor rate.
#'
#' @return A numeric value giving the estimated best distance cutoff value.
#'
#' @examples
#' irisDist <- dist(iris[,1:4])
#' estimateDC(irisDist)
#'   
#' @export
estimateDC <- function(distance, neighborRateLow=0.01, neighborRateHigh=0.02) {
    comb <- as.matrix(distance)
    size <- attr(distance, 'Size')
    dc <- min(distance)
    dcMod <- as.numeric(summary(distance)['Median']*0.01)
    while(TRUE) {
        neighborRate <- mean((apply(comb < dc, 1, sum)-1)/size)
        if(neighborRate > neighborRateLow && neighborRate < neighborRateHigh) break
        if(neighborRate > neighborRateHigh) {
            dc <- dc - dcMod
            dcMod <- dcMod/2
        }
        dc <- dc + dcMod
    }
    cat('Distance cutoff calculated as', dc,'\n')
    dc
}


#' Calculate clustering attributes
#' 
#' This function takes a distance matrix and optionally a distance cutoff and 
#' calculates the values necessary for clustering based on the algorithm 
#' proposed by Alex Rodrigues and Alessandro Laio (see \code{CITATION}). The 
#' actual assignment to clusters is done in a later step, based on user defined 
#' threshold values.
#' 
#' @details The function calculates \eqn{\rho} and \eqn{\delta} for the 
#'   observations in the provided distance matrix. If a distance cutoff is not 
#'   provided this is first estimated using \code{\link{estimateDC}} with 
#'   default values.
#'   
#'   The information kept in the densityCluster object is: \describe{ 
#'   \item{rho}{a vector of local density values.} \item{delta}{a vector of 
#'   minimum distances to observations of higher density.} \item{distance}{the 
#'   initial distance matrix.} \item{dc}{the distance cutoff used to calculate 
#'   rho.} \item{threshold}{a named vector specifying the threshold values for 
#'   rho and delta used for cluster detection.} \item{peaks}{a vector of indices
#'   specifying the cluster center for each cluster.} \item{clusters}{a vector 
#'   of cluster affiliations for each point, with the clusters being referenced as 
#'   indices in the peaks vector.} \item{halo}{a logical vector specifying for 
#'   each observation if it is considered part of the halo.} } Before running 
#'   \code{\link{findClusters}} the threshold, peaks, clusters and halo data is 
#'   \code{NA}.
#'   
#' @param distance a distance matrix.
#' @param dc a distance cutoff for calculating the local density. If missing it 
#'   will be estimated with \code{\link{estimateDC}(distance)}.
#' @param gaussian logical indicating that a Gaussian kernel should be used to
#'   estimate the density.
#'   
#' @return A densityCluster object. See Details for a description.
#'   
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- dclust(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{estimateDC}}, \code{\link{findClusters}}
#'   
#' @export
dclust <- function(distance, dc, gaussian=FALSE) {
    if(missing(dc)) {
        dc <- estimateDC(distance)
    }
    rho <- localDensity(distance, dc, gaussian=gaussian)
    delta <- distanceToPeak(distance, rho)
    res <- list(rho=rho, delta=delta, distance=distance, dc=dc, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA)
    class(res) <- 'densityCluster'
    res
}
