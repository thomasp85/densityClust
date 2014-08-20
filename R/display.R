#' @export
#' 
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
#' @param x A densityCluster object as produced by \code{\link{dclust}}
#' 
#' @param ... Additional parameters. Currently ignored
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
#' @seealso \code{\link{dclust}}
#' 
#' @export
#' 
plotMDS <- function (x, ...) {
  UseMethod("plotMDS", x)
}
#' @export
#' 
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