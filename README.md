<span style="font-variant:small-caps;">Eco-discretizer</span>
=========

Note 12/16/24: Make sure that the output is not missing any surprise species. Try removing the species the script mentions as having been dropped from the input species list to fix this (it may be the data filtering function is acting up.

Overview
---------
This script uses *k*-means clustering to discretize climatic datasets. The point is to convert continuous data to a discrete form that is useful for large numbers of comparative methods that use discrete data and also for exploratory analyses.

The script accepts as inputs a species list file `specieslist.csv` in the working directory (one species per line) and lists of extracted environmental data for points (format is with occurrences in rows and environmental data in columns, no headers). The naming format should be, e.g., `pno1_Species_a.csv` (numbering of variables starts at 1). These must be in the same directory as the script. The output is a csv containing the species labels and a numeric character coding from 0 to *k* - 1. The program will repeatedly test k values from 2 to the specified number.

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

Explanation of Files
---------
* pno_discretization.py&mdash;Python script for habitat coding.
* final_classification_k_7.csv&mdash;Result file for Saxifragales analysis.
* final_classification_k_7_biogeobearsformat.csv&mdash;Result file for Saxifragales analysis, in a format ready to use for BioGeoBEARS. For the BioGeoBEARS run script, see [https://github.com/ryanafolk/biogeographic_coder] and change file paths appropriately to point to habitat classifications.
* ultrametric_occur_matched_forcedultra.habitatclassificationmatched.tre&mdash;Tree used for BioGeoBEARS, with sampling matched to habitat classifications.

Possible errors
---------
`NameError: name 'process' is not defined` -- check for files that exist but are empty with `find . -size 0`.

 `FileNotFoundError: [Errno 2] No such file or directory:` -- check if the listed species only has some files. specieslist.csv should only contain entries that have data for all variables.


