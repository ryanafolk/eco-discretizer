<span style="font-variant:small-caps;">Eco-discretizer</span>
=========

Overview
---------
This script uses *k*-means clustering to discretize predicted niche occupancy profiles (PNOs) from niche models. The point is to convert continuous data to a discrete form that is useful for large numbers of comparative methods that use discrete data.

The script accepts as inputs a species list file and lists of extracted environmental data for points (format is with occurrences in rows and environmental data in columns, no headers). These must be in the same directory as the script. The output is a csv containing the species labels and a numeric character coding from 0 to *k* - 1. 

The distance matrices (raw and normalized) and k-means distortions are also saved in case they are useful for something. One idea is to plot the final averaged distance matrix in R, using a plot that makes sense like an MDS analysis. Classification surprises often happen when an individual is intermediate between an expected group and some neighbor.

The python libraries `numpy`, `scipy`, and `pandas` are required.

It is called like: 

```
./pno_discretization.py numberOfVariables numberOfCategories
```

where numberOfCategories is *k*, e.g.,

```
./pno_discretization.py 35 7
```

Approach
---------
Euclidean distances are calculated from all pairs of points between two species (100 randomly with replacement if the possible combinations are greater than this). The average of these is used to populated a species distance matrix. Finally, k-means clustering is applied to the matrix and a csv with the character coding (numerically from 0 to 1-k) is saved.

It is assumed that missing data is coded. as -9999 and points with missing data in any variable are discarded.
