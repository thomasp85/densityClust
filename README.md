# Clustering by fast search and find of density peaks
This package implements the clustering algorithm described by Alex Rodriguez and Alessandro Laio (2014, see references). It provides the user with tools for generating the ρ and δ values for each observation, as used in the algorithm, as well as using these to assign observations to clusters. This is done in two passes so that the user is free to reassign observations to clusters using a new set of ρ and δ thresholds without needing to recalculate everything.


# Approach
The algorithm is based on the assumptions that cluster centers are surrounded by neighbors with lower local density and that they are at a relatively large distance from any points with a higher local density. For each data point i, a local density is calculated as

    ρ[i] = sum(χ(d[i,j] − d[c]), over j)

with ``χ(x) = 1`` if ``x < 0`` and ``χ(x) = 0`` otherwise, and where ``d[i,j]`` is the distance between points ``i`` and ``j``, and ``d[c]`` is a cutoff distance. In addition, the minimum distance between the point and some point of higher density is calculated as

    δ[i] = min(d[i,j]), such that ρ[j] > ρ[i]

Points of abnormally high ρ and δ (as determined by user inspection) are designated as cluster centers, with the remainder of points each being assigned to the same cluster as that of its nearest neighbour of higher density.

Points are designated as being part of the core of a cluster, or part of its halo (i.e. the region in which points could be considered noise) by considering a border region, defined as the set of points assigned to a cluster within a distance dc of points assigned to other clusters. For each cluster, the point of highest density within the border region is found and its density denoted by ρb. Points within the cluster that have a density greater than ρb are assigned to the core, with the remainder being assigned to the halo.


# Cluster detection
The two main functions for this package are ``dclust`` and ``findClusters``. The former takes a distance matrix and optionally a distance cutoff and calculates ρ and δ for each observation. The latter takes the output of ``dclust`` and makes cluster assignments for each observation based on user defined ρ and δ thresholds. If the thresholds are not specified as arguments, the user is able to supply them interactively by clicking on a decision plot of δ against ρ.


# Plotting
Two types of plots are supported by this package, and both mimic the types of plots used in the original paper. The standard plot function produces a decision plot, with optional coloring of cluster centers if these are assigned. Furthermore, plotMDS performs a multidimensional scaling of the distance matrix and plots this as a scatterplot. If clusters are assigned, points are colored according to their assignment.


## Usage
```R
irisDist <- dist(iris[,1:4])
irisClust <- dclust(irisDist, gaussian=TRUE)
# Display a decision plot for choosing rho and delta thresholds
plot(irisClust)
# We see two 'outliers' of rho and delta both greater than 2,
# so use these as thresholds for our clustering
irisClust <- findClusters(irisClust, rho=2, delta=2)
# Plot a multidimensional scaling of the distance matrix
plotMDS(irisClust)
# Look at the original labels for the iris dataset and compare to our clustering
split(iris[,5], irisClust$clusters)
```
Note that while the iris dataset contains information on three different species of iris, only two clusters are detected by the algorithm. This is because two of the species (versicolor and virginica) are not clearly separated by their data.


## References
Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344 (6191), 1492-1496. doi: [10.1126/science.1242072](http://dx.doi.org/10.1126/science.1242072)