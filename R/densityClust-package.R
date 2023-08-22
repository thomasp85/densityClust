#' Clustering by fast search and find of density peaks
#'
#' This package implement the clustering algorithm described by Alex Rodriguez
#' and Alessandro Laio (2014). It provides the user with tools for generating
#' the initial rho and delta values for each observation as well as using these
#' to assign observations to clusters. This is done in two passes so the user is
#' free to reassign observations to clusters using a new set of rho and delta
#' thresholds, without needing to recalculate everything.
#'
#' @section Plotting:
#' Two types of plots are supported by this package, and both mimics the types of
#' plots used in the publication for the algorithm. The standard plot function
#' produces a decision plot, with optional colouring of cluster peaks if these
#' are assigned. Furthermore [plotMDS()] performs a multidimensional
#' scaling of the distance matrix and plots this as a scatterplot. If clusters
#' are assigned observations are coloured according to their assignment.
#'
#' @section Cluster detection:
#' The two main functions for this package are [densityClust()] and
#' [findClusters()]. The former takes a distance matrix and optionally
#' a distance cutoff and calculates rho and delta for each observation. The
#' latter takes the output of [densityClust()] and make cluster
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
#' @seealso [densityClust()], [findClusters()], [plotMDS()]
#' @aliases densityClust-package
#'
#' @references Rodriguez, A., & Laio, A. (2014). *Clustering by fast search and find of density peaks.* Science, **344**(6191), 1492-1496. doi:10.1126/science.1242072
#' 
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib densityClust, .registration = TRUE
## usethis namespace: end
NULL
