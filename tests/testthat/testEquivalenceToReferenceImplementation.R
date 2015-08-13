set.seed(123)
dists <- list(
   dist(matrix(rnorm(1000), ncol = 4)),
   dist(matrix(rnorm(1000), ncol = 20)),
   dist(matrix(rnorm(10000), ncol = 40)),
   dist(matrix(sample(1:100000, 1000), ncol = 4)),
   dist(matrix(sample(1:100000, 1000), ncol = 20)),
   dist(matrix(sample(1:100000, 1000), ncol = 50))
)

context("Reference implementation")

load("testdata/densityClustReference.RData")
densityClustNewImp <- lapply(dists, densityClust)
test_that("Test equivalence to reference implementation of densityClust", {
   expect_identical(densityClustReference, densityClustNewImp)   
})

load("testdata/estimateDcReference.RData")
estimateDcNewImp <- lapply(dists, estimateDc)
test_that("Test equivalence to reference implementation of estimateDc", {
   expect_identical(estimateDcReference, estimateDcNewImp)   
})

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# load("testdata/localDensityReference.RData")
# localDensityNewImp <- Map(densityClust:::localDensity, dists, estimateDcNewImp)
# test_that("Test equivalence to reference implementation of localDensity", {
#    expect_identical(localDensityReference, localDensityNewImp)
# })
# 
# load("testdata/distanceToPeakReference.RData")
# distanceToPeakNewImp <- Map(densityClust:::distanceToPeak, dists, localDensityNewImp)
# test_that("Test equivalence to reference implementation of distanceToPeak", {
#    expect_identical(distanceToPeakReference, distanceToPeakNewImp)
# })

load("testdata/gaussianDensityClustReference.RData")
gaussianDensityClustNewImp <- lapply(dists, FUN = function(x) densityClust(x, gaussian = TRUE))
test_that("Test equivalence to reference implementation of gaussianDensityClust", {
   expect_identical(gaussianDensityClustReference, gaussianDensityClustNewImp)   
})

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# load("testdata/gaussianLocalDensityReference.RData")
# gaussianLocalDensityNewImp <- Map(f = function(x, y) densityClust:::localDensity(x, y, gaussian = TRUE), dists, estimateDcReference)
# test_that("Test equivalence to reference implementation of localDensity", {
#    expect_identical(gaussianLocalDensityReference, gaussianLocalDensityNewImp)
# })

