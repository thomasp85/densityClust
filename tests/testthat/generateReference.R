# This code was run with the initial reference implementation (i.e. https://github.com/thomasp85/densityClust/commit/b038fb30ea6f59d60a3a4b45eaa3ac9a504951f6)
# It generates the .RData files with the expected results based on that 
# implementation, which are then tested against with new implementations. 

set.seed(123)
dists <- list(
   dist(matrix(rnorm(1000), ncol = 4)),
   dist(matrix(rnorm(1000), ncol = 20)),
   dist(matrix(rnorm(10000), ncol = 40)),
   dist(matrix(sample(1:100000, 1000), ncol = 4)),
   dist(matrix(sample(1:100000, 1000), ncol = 20)),
   dist(matrix(sample(1:100000, 1000), ncol = 50))
)

# non-Gaussian
densityClustReference <- lapply(dists, densityClust)
save(densityClustReference, file = "tests/testthat/testdata/densityClustReference.RData")

estimateDcReference <- lapply(dists, estimateDc)
save(estimateDcReference, file = "tests/testthat/testdata/estimateDcReference.RData")

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# localDensityReference <- Map(densityClust:::localDensity, dists, estimateDcReference)
# save(localDensityReference, file = "tests/testthat/testdata/localDensityReference.RData")
# 
# distanceToPeakReference <- Map(densityClust:::distanceToPeak, dists, localDensityReference)
# save(distanceToPeakReference, file = "tests/testthat/testdata/distanceToPeakReference.RData")

# Gaussian
gaussianDensityClustReference <- lapply(dists, FUN = function(x) densityClust(x, gaussian = TRUE))
save(gaussianDensityClustReference, file = "tests/testthat/testdata/gaussianDensityClustReference.RData")

gaussianLocalDensityReference <- Map(f = function(x, y) densityClust:::localDensity(x, y, gaussian = TRUE), dists, estimateDcReference)
save(gaussianLocalDensityReference, file = "tests/testthat/testdata/gaussianLocalDensityReference.RData")