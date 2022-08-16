#!/usr/bin/env python3
# Must be run on "pno" file from point extraction script
# Must use version WITH MISSING DATA
# Must be version of this WITHOUT sorting since we are calculating multidimensional distances and need to preserve point-order between variables
# Calculate a matrix of species distances for each variable
# The purpose of matrix normalization is to enforce equal weight among ecological variables

import csv # To process tab-delimited files
import numpy # Used for random sampling 
import sys # To process arguments
from scipy.spatial import distance # Distance function
from scipy.cluster.vq import kmeans, vq, whiten # Clustering
import pandas # Used to handle matrix of distances
import glob
import itertools # To get possible pairs
import random # Random list order

numpy.seterr(all='raise') # Exceptions instead of runtime warnings

try: # Parse command line arguments
	processtotal = int(sys.argv[1]) # Number of niche variables, = number of runs of BayesTraits per bootstrap
	clusternumber = int(sys.argv[2]) # Desired number of clusters
except: # Handle incorrect argument passing
	print("Error, two arguments must be provided: the number of variables and the desired number of clusters.")
	print("Example: python3 pno_discretization.py 12 10")
	sys.exit()

# Import list of species
specieslist = []
try:
	with open('specieslist.csv', 'r') as speciesfile:
		reader=csv.reader(speciesfile, delimiter='\t')
		for row in reader:
			specieslist.append(row[0])
except: 
	print("Error in parsing species list... Make sure there is a specieslist.csv in the working directory. If files were created on a Mac or Windows machine, verify that end-of-line characters are UNIX-compliant. Also ensure that you have the correct number of environmental variable files.")
	exit()

# Define dataframe containing distances
distanceDataFrame = pandas.DataFrame(data=numpy.zeros((len(specieslist), len(specieslist))), index=specieslist, columns=specieslist) # Empty dataframe

# Read input data in PNO format, return n-dimensional point data
def ReadItem(species):
	i = 1
	alldatalist = []
	while i <= processtotal:
		input = './pno{0}_{1}.csv'.format(i,species) # Iteratively open csv files for this species
		datalist = [] # List of environmental data points
		try:
			with open(input, 'r') as datafile:
				reader=csv.reader(datafile, delimiter=',')
				next(reader)
				for row in reader:
					datalist.append(float(row[1]))
		except: 
			print('PNO {0} file not found for the species {1}.'.format(process,species))
			print('Ensure that you have named these correctly; they will be treated as missing data.')
			return
		alldatalist.append(datalist) # Save the data column for this species and variable
		i+=1
	alldatalist = list(zip(*alldatalist)) # Zip returns a zip object in Python3; asterisk is to iterate over list elements for zipping	
	# Combine list of variable lists so that we have elements containing 
	# Alldatalist should contain a list of tuples of length processortotal; distance will be n-dimensional with n = processortotal
	return alldatalist

	
def ReadCleanedItem(species):
	input = './{0}_cleaned.csv'.format(species) # Iteratively open csv files for this species
	datalist = []
	try:
		with open(input, 'r') as datafile:
			reader=csv.reader(datafile, delimiter=',')
			for row in reader:
				row = list(map(float,row))
				datalist.append(row)
	except: 
		print('PNO {0} file not found for the species {1}.'.format(process,species))
		print('Ensure that you have named these correctly; they will be treated as missing data.')
		return
	return datalist

missing_species = [] # List to hold species for which the cleaning procedure removed everything
for species in specieslist:  # Clean a species and output n-dimensional point data
	print("Cleaning species {0}.".format(species))
	item = ReadItem(species) # Read in the species' data
	item_cleaned = []
	for tuple in item: # Clean out missing data values in ANY cell
		if -9999 in tuple:
			pass
		else:
			item_cleaned.append(tuple)
	if not item_cleaned: # Check if we removed all of the points
		missing_species.append(species)
		print("Lost everything for species {0}.".format(species))
		continue
	with open('./{0}_cleaned.csv'.format(species), "w+") as file:
		writer = csv.writer(file)
		writer.writerows(item_cleaned)	 
	missing_species = list(set(missing_species)) # Remove duplicate species -- needed for memory management

reduced_species = list(set(specieslist) - set(missing_species)) # Remove missing species
reduced_species.sort() # Alphabetize since we lose order by set operations
	
speciesDict = {} # Empty dictionary
for species in reduced_species: # Read in the cleaned species data as a matrix for easy access
	speciesDict[species] = ReadCleanedItem(species)

# Do the distance calculation, with various alternative actions for species with very large numbers of occurrences.
distanceDataFrame = pandas.DataFrame(data=numpy.zeros((len(reduced_species), len(reduced_species))), index=reduced_species, columns=reduced_species) # Empty dataframe to hold species pairwise distances 
for columnspecies in reduced_species:
	print("Distance matrix: on species {0}.".format(columnspecies))
	for rowspecies in reduced_species:
		rowsample = speciesDict[rowspecies]
		columnsample = speciesDict[columnspecies]
		speciesdistancelist = []
		if len(rowsample) * len(columnsample) <= 10000: # For moderate numbers of points
			allpairs = list(itertools.product(rowsample, columnsample)) # This yields the cartesian product (all possible pairs)
			if len(allpairs) > 100:
				random.shuffle(allpairs) # Randomize order before truncation
				del allpairs[100:] # Truncate list to 100 items
		else:  # For species that will produce very large cartesian products
			i = 1
			allpairs = []
			while i <= 100: # Randomly get a subset with replacement
				row_item = rowsample[random.randrange(0, len(rowsample))]
				column_item = columnsample[random.randrange(0, len(columnsample))]
				allpairs.append([row_item, column_item])
				i += 1
		for x in allpairs:
			speciesdistancelist.append(distance.euclidean(x[0], x[1])) # n-dimensional distance 
		speciesdistance = numpy.average(speciesdistancelist) # Average of pairwise distances
		distanceDataFrame.at[rowspecies, columnspecies] = speciesdistance # Populate dataframe with distance

print("Missing species are:")		
print(missing_species)

# Process and write the matrix
for species in reduced_species: 
	distanceDataFrame.at[species, species] = None # Clearing out the diagonal
distanceDataFrame.to_csv('./distancematrix.csv')

kmeans_list = []
for species in reduced_species: # Read in the cleaned species data as a matrix for easy access
	temp_list = []
	temp_list = speciesDict[species]
	samples_array = numpy.asarray(temp_list).T # Transpose the matrix
	samples = samples_array.tolist()
	median_list = []
	for i in samples:
		median_list.append(float(numpy.percentile(i, 50)))
	kmeans_list.append(median_list)
	
print("K-means for {0} species.".format(len(kmeans_list)))

with open("medians_data.csv", "w+") as file:
    writer = csv.writer(file)
    for item, species in zip(kmeans_list, reduced_species):
    	row = []
    	row.append(species)
    	for i in item:	
    		row.append(i)
    	writer.writerows([row])	

kmeans_array = numpy.asarray(kmeans_list)

kmeans_array = numpy.asmatrix(kmeans_array)

# RAW
# Perform the clustering
kmeans_array_normalized = whiten(kmeans_array)

distortion_list = []

for k in range(2, clusternumber + 1):
	print("K-means clustering: on k = {0}".format(k))
	
	centroids,distortion = kmeans(kmeans_array_normalized, k) # Won't work without ,_ since it returns a tuple

	classification,_ = vq(kmeans_array_normalized, centroids) # Won't work without ,_

	# Pair classification list with species list
	classificationFinal = []
	for x, y in zip (specieslist, classification):
		classificationFinal.append([x, y])
	
	# Save classification to csv
	with open("final_classification_k_{0}.csv".format(k), "w+") as file:
	    writer = csv.writer(file)
	    writer.writerows(classificationFinal)	
	
	distortion_list.append([k, distortion])

with open("distortions.csv".format(k), "w+") as file:
    writer = csv.writer(file)
    for row in distortion_list:
    	writer.writerows([row])	
