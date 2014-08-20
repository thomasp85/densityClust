#' Clustering by fast search and find of density peaks
#'
#' This package implements the clustering algorithm described by Alex Rodriguez 
#' and Alessandro Laio (2014, see references). It provides the user with tools 
#' for generating the \eqn{\rho} and \eqn{\delta} values for each observation, 
#' as used in the algorithm, as well as using these to assign observations to 
#' clusters. This is done in two passes so that the user is free to reassign 
#' observations to clusters using a new set of \eqn{\rho} and \eqn{\delta} 
#' thresholds without needing to recalculate everything.
#'
#' @section Approach: The algorithm is based on the assumptions that cluster 
#'   centers are surrounded by neighbors with lower local density and that they 
#'   are at a relatively large distance from any points with a higher local 
#'   density. For each data point \eqn{i}, a local density is calculated as 
#'   \deqn{ \rho_i = \sum_j \chi (d_{i,j} - d_c), }{ \rho[i] = \Sigma[j] \chi 
#'   (d[i,j] - d[c]), } with \eqn{ \chi(x) = 1 } if \eqn{ x < 0 } and \eqn{ 
#'   \chi(x) = 0 } otherwise, and where \eqn{d_{i,j}}{d[i,j]} is the distance 
#'   between points \eqn{i} and \eqn{j}, and \eqn{d_c}{d[c]} is a cutoff 
#'   distance. In addition, the minimum distance between the point and some 
#'   point of higher density is calculated as \deqn{ \delta_i = \min_{j:\rho_j >
#'   \rho_i}(d_{i,j}). }{ \delta[i] = min d[i,j], for j such that \rho[j] > 
#'   \rho[i]. }
#'
#'   Points of abnormally high \eqn{\rho} and \eqn{\delta} (as determined by 
#'   user inspection) are designated as cluster centers, with the remainder of 
#'   points each being assigned to the same cluster as that of its nearest 
#'   neighbour of higher density.
#'
#'   Points are designated as being part of the core of a cluster, or part of 
#'   its halo (i.e. the region in which points could be considered noise) by 
#'   considering a border region, defined as the set of points assigned to a 
#'   cluster within a distance \eqn{d_c}{d[c]} of points assigned to other 
#'   clusters. For each cluster, the point of highest density within the border 
#'   region is found and its density denoted by \eqn{\rho_b}{\rho[b]}. Points 
#'   within the cluster that have a density greater than \eqn{\rho_b}{\rho[b]} 
#'   are assigned to the core, with the remainder being assigned to the halo.
#'
#' @section Cluster detection: The two main functions for this package are 
#'   \code{\link{dclust}} and \code{\link{findClusters}}. The former takes a 
#'   distance matrix and optionally a distance cutoff and calculates \eqn{\rho} 
#'   and \eqn{\delta} for each observation. The latter takes the output of 
#'   \code{\link{dclust}} and makes cluster assignments for each observation 
#'   based on user defined \eqn{\rho} and \eqn{\delta} thresholds. If the 
#'   thresholds are not specified as arguments, the user is able to supply them 
#'   interactively by clicking on a decision plot of \eqn{\delta} against 
#'   \eqn{\rho}.
#'
#' @section Plotting: Two types of plots are supported by this package, and both
#'   mimic the types of plots used in the original paper. The standard plot 
#'   function produces a decision plot, with optional coloring of cluster 
#'   centers if these are assigned. Furthermore, \code{\link{plotMDS}} performs 
#'   a multidimensional scaling of the distance matrix and plots this as a 
#'   scatterplot. If clusters are assigned, points are colored according to 
#'   their assignment.
#'
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- dclust(irisDist, gaussian=TRUE)
#' # Display a decision plot for choosing rho and delta thresholds
#' plot(irisClust)
#' # We see two 'outliers' of rho and delta both greater than 2,
#' # so use these as thresholds for our clustering
#' irisClust <- findClusters(irisClust, rho=2, delta=2)
#' # Plot a multidimensional scaling of the distance matrix
#' plotMDS(irisClust)
#' # Look at the original labels for the iris dataset and compare to our clustering
#' split(iris[,5], irisClust$clusters)
#' 
#' @seealso \code{\link{dclust}}, \code{\link{findClusters}}, 
#'   \code{\link{plotMDS}}
#'
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and 
#'   find of density peaks. Science, 344 (6191), 1492-1496. 
#'   doi:10.1126/science.1242072
#'
#' @docType package
#' @name densityClust-package
#' @aliases densityClust
#'
NULL