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
#' @return A densityCluster object. See details for a description.
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
#' @seealso \code{\link{estimateDc}}, \code{\link{findClusters}}
#' 
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @export
#' 
dclust <- function(distance, dc, gaussian=FALSE) {
    if(missing(dc)) {
        dc <- estimateDc(distance)
    }
    rho <- localDensity(distance, dc, gaussian=gaussian)
    delta <- distanceToPeak(distance, rho)
    res <- list(rho=rho, delta=delta, distance=distance, dc=dc, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA)
    class(res) <- 'densityCluster'
    res
}
