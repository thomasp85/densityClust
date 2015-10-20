# This code generates expected results from the initial reference implementation 
# (i.e. https://github.com/thomasp85/densityClust/commit/b038fb30ea6f59d60a3a4b45eaa3ac9a504951f6)
# It is called by the testing code to compare the old results with the new.

set.seed(123)
dists <- list(
   dist(matrix(rnorm(1000), ncol = 4)),
   dist(matrix(rnorm(1000), ncol = 20)),
   dist(matrix(rnorm(10000), ncol = 40)),
   dist(matrix(sample(1:100000, 1000), ncol = 4)),
   dist(matrix(sample(1:100000, 1000), ncol = 20)),
   dist(matrix(sample(1:100000, 1000), ncol = 50))
)

referenceImplementation <- function(distance, dc, gaussian=FALSE) {
   if(missing(dc)) {
      dc <- reference_estimateDc(distance)
   }
   rho <- reference_localDensity(distance, dc, gaussian=gaussian)
   delta <- reference_distanceToPeak(distance, rho)
   res <- list(rho=rho, delta=delta, distance=distance, dc=dc, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA)
   class(res) <- 'densityCluster'
   res
}

reference_estimateDc <- function(distance, neighborRateLow=0.01, neighborRateHigh=0.02) {
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

reference_distanceToPeak <- function(distance, rho) {
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

reference_localDensity <- function(distance, dc, gaussian=FALSE) {
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

# Because the new implementation of estimateDc does not maintain equality with
# the previous implementation, calculate the cutoffs using the new 
# implementation. Then pass the calculated cutoffs into the reference 
# implementation and the new implementation to test that the rest of the 
# implementations are the same. 
dcs <- lapply(dists, estimateDc)

# Reference DCs for comparison 
referenceDcs <- lapply(dists, reference_estimateDc)

# non-Gaussian
densityClustReference <- Map(referenceImplementation, dists, dcs)

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# localDensityReference <- Map(reference_localDensity, dists, dcs)
# 
# distanceToPeakReference <- Map(reference_distanceToPeak, dists, localDensityReference)

# Gaussian
gaussianDensityClustReference <- Map(referenceImplementation, dists, dcs, TRUE) 

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# gaussianLocalDensityReference <- Map(f = function(x, y) reference_localDensity(x, y, gaussian = TRUE), dists, estimateDcReference)