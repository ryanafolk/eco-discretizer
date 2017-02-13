<span style="font-variant:small-caps;">Eco-discretizer</span>
=========

Overview
---------
This script uses k-means clustering to discretize predicted niche occupancy profiles (PNOs) from niche models. The point is to convert continuous data to a discrete form that is useful for large numbers of comparative methods that use discrete data.

The script accepts as inputs species list file and predicted niche occupancy profile files, in the exact same format as that used for `ambitus`. These must be in the same directory as the script. The output is a csv containing the species labels and a numeric character coding from 0 to k - 1. 

The distance matrices (raw and normalized) are also saved in case they are useful for something. One idea is to plot the final averaged distance matrix in R, using a plot that makes sense like an MDS analysis. Classification surprises often happen when an individual is intermediate between an expected group and some neighbor.

The python libraries `numpy`, `scipy`, and `pandas` are required.

It is called like: 

```
./pno_discretization.py numberOfVariables numberOfSamples numberOfCategories
```

E.g.,

```
./pno_discretization.py 12 10000 10
```

Approach
---------
PNOs are typically the way in which niche models are typically used on phylogenies. There is a problem in that in that the histograms output for different species do not have comparable bins; enforcing common bins across species is not necessarily desireable. Distance metrics for pdfs (probability density functions) can be applied to histograms but may have undesireable properties for non-overlapping pdfs, in that large numbers of ties are expected. Sampling a very large number of points (e.g., 10,000) from histograms yeilds comparable objects that are very close to the histograms they came from. 

From these pools of random points drawn from the histograms, the the cumulative Euclidean distance is calculated for every species pair in the tree, for each environmental predictor variable. Different scalings between variables would bias the result to favor clustering that follows variables with high magnitude and variance, so within each variable distance matrix, each is subtracted from the matrix mean and divided by the range to normalize. Then the individual distance matrices are averaged across variables. Finally, k-means clustering is applied and a csv with the character coding (numerically from 0 to 1-k) is saved.

The code takes about 2 minutes for 64 taxa and 12 environmental predictor variables on a Macbook Pro with a solid state drive.
