#' @name plotDensityClust
#' @title Plot densityCluster results
#' @description Generate a single panel of up to three diagnostic plots for a
#'   \code{densityClust} object.
#'
#' @param x A densityCluster object as produced by \code{\link{densityClust}}
#' @param type A character vector designating which figures to produce. Valid
#'   options include \code{"dg"} for a decision graph of \eqn{\delta} vs.
#'   \eqn{\rho}, \code{"gg"} for a gamma graph depicting the decrease of
#'   \eqn{\gamma} (= \eqn{\delta} * \eqn{\rho}) across samples, and \code{"mds"},
#'   for a Multi-Dimensional Scaling (MDS) plot of observations. Any combination
#'   of these three can be included in the vector, or to produce all plots,
#'   specify \code{type = "all"}.
#' @param n Number of observations to plot in the gamma graph.
#' @param mds A matrix of scores for observations from a Principal Components
#'   Analysis or MDS. If omitted, and a MDS plot has been requested, one will
#'   be calculated.
#' @param dim.x,dim.y The numbers of the dimensions to plot on the x and y
#'   axes of the MDS plot.
#' @param col Vector of colors for clusters.
#' @param alpha Value in \code{0:1} controlling transparency of points in the
#'   decision graph and MDS plot.
#'
#' @return A panel of the figures specified in \code{type} are produced.
#'   If designated, clusters are color-coded and labelled. If present in
#'   \code{x}, the rho and delta thresholds are designated in the
#'   decision graph by a set of solid black lines.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(iris)
#' data.dist <- dist(iris[, 1:4])
#' pca <- princomp(iris[, 1:4])
#'
#' # Run initial density clustering
#' dens.clust <- densityClust(data.dist)
#
#' op <- par(ask = TRUE)
#'
#' # Show the decision graph
#' plotDensityClust(dens.clust, type = "dg")
#'
#' # Show the decision graph and the gamma graph
#' plotDensityClust(dens.clust, type = c("dg", "gg"))
#'
#' # Cluster based on rho and delta
#' new.clust <- findClusters(dens.clust, rho = 4, delta = 2)
#'
#' # Show all graphs with clustering
#' plotDensityClust(new.clust, mds = pca$scores)
#'
#' par(op)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes_string geom_text geom_point geom_segment labs
#'   theme_bw theme scale_color_manual geom_line geom_label
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices rainbow
#' @export
#'
plotDensityClust <- function(x, type = "all", n = 20,
                             mds = NULL, dim.x = 1, dim.y = 2,
                             col = NULL, alpha = 0.8) {

  type <- tolower(type)
  if(any(pmatch(type, "all", nomatch = 0))) type <- c("dg", "gg", "mds")

  df <- data.frame(
    rho = x$rho, delta = x$delta, gamma = x$rho * x$delta,
    peaks = FALSE, cluster = factor(x$clusters), halo = x$halo
  )
  df$peaks[x$peaks] <- TRUE

  if(is.null(col)) {
    num.cols <- max(nlevels(df$cluster), 3)
    col <- if(num.cols <= 8) {
      brewer.pal(num.cols, "Set2")
    } else if(num.cols <= 12) {
      brewer.pal(num.cols, "Set3")
    } else rainbow(num.cols + 1)[1:num.cols]
  }

  plots <- list(dg = NULL, gg = NULL, mds = NULL)

  # Plot decision graph (dg)
  if(any(pmatch(type, "dg", nomatch = 0))) {
    plots$dg <- ggplot(df, aes_string(x = "rho", y = "delta"))
    if(!any(is.na(x$threshold))) {
      rho <- x$threshold["rho"]
      delta <- x$threshold["delta"]
      thresh.df <- data.frame(
        x = c(rho, rho),
        y = c(delta, delta),
        xend = c(rho, Inf),
        yend = c(Inf, delta)
      )
      plots$dg <- plots$dg +
        geom_segment(
          aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
          data = thresh.df, inherit.aes = F,
          lineend = "butt"
        )
    }
    if(any(df$peaks)) {
      plots$dg <- plots$dg +
        geom_label(
          aes_string(label = "cluster", color = "cluster"),
          data = df[df$peaks, ],
          fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col)
    }
    plots$dg <- plots$dg  +
      geom_point(
        data = df[!df$peaks, ],
        size = 3, color = "gray50", alpha = alpha
      ) +
      labs(x = expression(rho), y = expression(delta), color = "Cluster") +
      theme(legend.position = "none")
  }

  # Plot gamma graph (gg)
  if(any(pmatch(type, "gg", nomatch = 0))) {
    gg.df <- df[order(df$gamma, decreasing = TRUE), ]
    gg.df <- gg.df[1:n, , drop = FALSE]
    gg.df$Sample <- 1:nrow(gg.df)

    plots$gg <- ggplot(gg.df, aes_string(x = "Sample", y = "gamma")) + geom_line()
    if(any(gg.df$peaks)) {
      plots$gg <- plots$gg +
        geom_label(
          aes_string(label = "cluster", color = "cluster"),
          data = gg.df[gg.df$peaks, , drop = FALSE],
          fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col)
    }
    plots$gg <- plots$gg +
      geom_point(
        data = gg.df[!gg.df$peaks, , drop = FALSE],
        size = 3, color = "gray50"
      ) +
      labs(y = expression(gamma), color = "Cluster") +
      theme(legend.position = "none")
  }

  # Plot MDS (mds)
  if(any(pmatch(type, "mds", nomatch = 0))) {
    if(is.null(mds)) mds <- cmdscale(x$distance, k = max(dim.x, dim.y))
    df$x <- mds[, dim.x]
    df$y <- mds[, dim.y]

    plots$mds <- ggplot()
    plots$mds <- if(all(is.na(df$cluster))) {
      plots$mds +
        geom_point(
          aes_string(x = "x", y = "y"),
          data = df,
          size = 3, color = "gray50", alpha = alpha
        )
    } else {
      plots$mds +
        geom_point(
          aes_string(x = "x", y = "y", color = "cluster"),
          data = df[df$halo, , drop = FALSE],
          shape = 21, size = 3
        ) +
        geom_point(
          aes_string(x = "x", y = "y", color = "cluster"),
          data = df[!df$halo, , drop = FALSE],
          size = 3, alpha = alpha
        ) +
        geom_label_repel(
          aes_string(x = "x", y = "y", label = "cluster", color = "cluster"),
          data = df[df$peaks, , drop = FALSE],
          size = 6, fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col, na.value = "gray50")
    }
    plots$mds <- plots$mds +
      labs(x = paste("Dimension", dim.x), y = paste("Dimension", dim.y)) +
      theme(legend.position = "none")
  }

  has.plot <- !sapply(plots, is.null)
  switch(
    sum(has.plot),
    print(plots[[which(has.plot)]]),
    {
      plots <- plots[has.plot]
      if("mds" %in% names(plots)) plots$nrow <- 2 else plots$ncol <-2
      do.call(grid.arrange, plots)
    },
    {
      plots$layout_matrix <- matrix(c(1, 3, 2, 3), nrow = 2)
      do.call(grid.arrange, plots)
    }
  )
}