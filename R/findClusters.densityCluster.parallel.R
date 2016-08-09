#' @rdname findClusters
#' 
#' @param rho The threshold for local density when detecting cluster peaks
#' 
#' @param delta The threshold for minimum distance to higher density when detecting cluster peaks
#' 
#' @param cores Number of cores used in the mclapply function 
#'
#' @param plot Logical. Should a decision plot be shown after cluster detection
#' 
#' @export
#' @importFrom graphics plot locator

findClusters.densityCluster.ori <- function(x, rho, delta, plot=FALSE, ...) {
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
    #comb <- as.matrix(x$distance) ##huge matrix this step 
    runOrder <- order(x$rho, decreasing = TRUE)
    cluster <- rep(NA, length(x$rho))
    for(i in runOrder) { #can we parallel this part? 
      if((i %% 1000) == 0)
        message(paste('the runOrder index is ', i))

        if(i %in% x$peaks) {
            cluster[i] <- match(i, x$peaks)
        } else {
            higherDensity <- which(x$rho>x$rho[i])
            cluster[i] <- cluster[higherDensity[which.min(findDistValueByRowColInd(x$distance, i, higherDensity))]]#cluster[higherDensity[which.min(comb[i, higherDensity])]]
        }
    }
    x$clusters <- cluster
    
    # Calculate core/halo status of observation
    border <- rep(0, length(x$peaks))
    for(i in 1:length(x$peaks)) { #can we parallelize this part? 
        message(paste('the current index of the peak is ', i))

        averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], '+')/2 #this match the density of two cells as in the distance matrix 
        index <- findDistValueByRowColInd(x$distance, i, higherDensity) <= x$dc #comb[cluster == i, cluster != i] <= x$dc #distance matrix 
        if(any(index)) border[i] <- max(averageRho[index]) #calculate the matrix value 
    }
    x$halo <- x$rho < border[cluster] #should we do this for each cluster? 
    x
}

