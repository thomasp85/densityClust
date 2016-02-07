#' @name plotFuncs
#' @title Plotting Functions
#' @description Plotting functions for density cluster algorithm
#'
#' @param x A densityCluster object as produced by \code{\link{densityClust}}
#' @param dim.x,dim.y Vectors giving projections on first and second
#'   multi-dimensonal scaling (principal coordinates) axes. If not given,
#'   MDS is calculated and returned by the function.
#' @param n Number of observations to plot in Gamma graph.
#' @param col Vector of colors for clusters.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(iris)
#' data.dist <- dist(iris[, 1:4])
#' pca <- princomp(iris[, 1:4])
#'
#' op <- par(ask = TRUE)
#' dens.clust <- densityClust(data.dist)
#' plotDensClust(dens.clust)
#'
#' n.vec <- 2:3
#' clust.list <- lapply(n.vec, function(n) {
#'   new.clust <- findClusters(dens.clust, k = n)
#'   plotDensClust(new.clust, pca$scores)
#'   new.clust
#' })
#' names(clust.list) <- n.vec
#' par(op)
#'
#' @export
#'
plotDecisionGraph <- function(x, col = NULL) {
  plot(x$rho, x$delta, xlab = expression(rho), ylab = expression(delta))
  if(!is.na(x$peaks[1])) {
    if(is.null(col)) col <- 2:(1 + length(x$peaks))
    points(x$rho[x$peaks], x$delta[x$peaks], col = col, pch = 19)
    text(x$rho[x$peaks], x$delta[x$peaks], 1:length(x$peaks), adj = c(1, 1))
    rho <- x$threshold["rho"]
    delta <- x$threshold["delta"]
    segments(rho, delta, rho, par("usr")[4], col = "gray60")
    segments(rho, delta, par("usr")[2], delta, col = "gray60")
  }
}

#' @rdname plotFuncs
#' @export
#'
plotGammaGraph <- function(x, n = 20, col = NULL) {
  g <- sort(x$rho * x$delta, decreasing = TRUE)
  plot(g[1:n], type = "b", pch = 19, ylab = expression(gamma))
  if(!is.na(x$peaks[1])) {
    n.pks <- length(x$peaks)
    if(is.null(col)) col <- 2:(1 + n.pks)
    points(1:n.pks, g[1:n.pks], col = col, pch = 19)
    abline(v = n.pks + 0.5)
    text(1:n.pks, g[1:n.pks], 1:n.pks, adj = c(0, 0))
  }
}

#' @rdname plotFuncs
#' @export
#'
plotClustMDS <- function(x, dim.x = NULL, dim.y = NULL, col = NULL) {
  mds <- NULL
  if(is.null(dim.x) | is.null(dim.y)) {
    mds <- cmdscale(x$distance)
  } else if((length(dim.x) == length(dim.y)) & length(dim.x == nrow(x$distance))) {
    mds <- cbind(dim.x, dim.y)
  } else {
    warning("length of 'dim.x' or 'dim.y' is not same as number of observations in 'x'")
    return(NULL)
  }

  dim.x <- mds[, 1]
  dim.y <- mds[, 2]
  xlim <- range(dim.x)
  ylim <- range(dim.y)

  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  axis(1)
  axis(2)
  box()
  abline(h = 0, v = 0)
  pch <- 1
  pks <- x$peaks
  if(!all(is.na(pks))) {
    col <- if(is.null(col)) x$clusters + 1 else col[x$clusters]
    pch <- ifelse(x$halo, 1, 19)
  } else col <- 1
  points(dim.x, dim.y, col = col, pch = pch, cex = 0.7)
  if(!is.null(pks)) text(dim.x[pks], dim.y[pks], 1:length(pks), cex = 1.5)

  invisible(mds)
}


#' @rdname plotFuncs
#' @export
#'
plotDensClust <- function(x, dim.x = NULL, dim.y = NULL, n = 20, col = NULL) {
  op <- par(no.readonly = TRUE)
  par(mar = c(4, 5, 1, 1))
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = T), heights = c(1, 2))
  plotDecisionGraph(x = x, col = col)
  plotGammaGraph(x = x, n = n, col = col)
  mds <- plotClustMDS(x = x, dim.x = dim.x, dim.y = dim.y, col = col)
  layout(matrix(1))
  par(op)
  invisible(mds)
}