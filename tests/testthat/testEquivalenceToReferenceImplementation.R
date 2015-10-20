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

# get dcs and reference targets
source("generateReference.R")

dcComparison <- simplify2array(Map(function(x, y) abs(1 - x / y) <= 0.15, dcs, referenceDcs))
test_that("Reference DCs and new DCs are within 15% of each other", {
   expect_true(all(dcComparison))
})

densityClustNewImp <- lapply(dists, densityClust)
test_that("Test equivalence to reference implementation of densityClust", {
   expect_equal(densityClustReference, densityClustNewImp)   
})

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# localDensityNewImp <- Map(densityClust:::localDensity, dists, estimateDcNewImp)
# test_that("Test equivalence to reference implementation of localDensity", {
#    expect_equal(localDensityReference, localDensityNewImp)
# })
# 
# distanceToPeakNewImp <- Map(densityClust:::distanceToPeak, dists, localDensityNewImp)
# test_that("Test equivalence to reference implementation of distanc eToPeak", {
#    expect_equal(distanceToPeakReference, distanceToPeakNewImp)
# })

gaussianDensityClustNewImp <- lapply(dists, FUN = function(x) densityClust(x, gaussian = TRUE))
test_that("Test equivalence to reference implementation of gaussianDensityClust", {
   expect_equal(gaussianDensityClustReference, gaussianDensityClustNewImp)   
})

# convenient for debugging, but calling non-exported functions not allowed in CRAN
# gaussianLocalDensityNewImp <- Map(f = function(x, y) densityClust:::localDensity(x, y, gaussian = TRUE), dists, estimateDcReference)
# test_that("Test equivalence to reference implementation of localDensity", {
#    expect_equal(gaussianLocalDensityReference, gaussianLocalDensityNewImp)
# })