#!/usr/bin/env python3

import csv # To process tab-delimited files
import numpy # Used for random sampling 
import sys # To process arguments
from scipy.spatial import distance # Distance function
from scipy.cluster.vq import kmeans,vq # Clustering
import pandas # Used to handle matrix of distances
import glob

try: # Parse command line arguments
	processtotal = int(sys.argv[1]) # Number of niche variables, = number of runs of BayesTraits per bootstrap
	bootstraptotal = int(sys.argv[2]) # Number of times to sample random points in variable occupancy ranges
	clusternumber = int(sys.argv[3]) # Desired number of clusters
except: # Handle incorrect argument passing
	print("Error, three arguments must be provided: the number of variables, the number of samples for each variable range, and the desired number of clusters.")
	print("Example: python3 pno_discretization.py 12 10000 10")
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

# Define sampling procedure from distributions (identical to ambitus)
def SampleDistribution(process, bootstraptotal, species):
	input = './pno{0}_{1}.csv'.format(process,species)
	
	binlist = [] # List of bins created by phyloclim
	probabilitylist = [] # List of probabilities matched to bins
	
	try:
		with open(input, 'r') as datafile:
				reader=csv.reader(datafile, delimiter=',')
				next(reader)
				for row in reader:
					binlist.append(float(row[1]))
					probabilitylist.append(float(row[2]))
	except: 
		print('PNO {0} file not found for the species {1}.'.format(process,species))
		print('Ensure that you have named these correctly; they will be treated as missing data.')
		return
	
	# The following code handles rounding error of bin probabilities, which can be fairly high in phyloclim, and fixes probabilities by spreading the error evenly across bins
	sum = numpy.sum(probabilitylist)
	if sum > 1:
		correction = (sum - 1) / len(probabilitylist) # Calculate how much probabilities observed exceed one and calculate a correction value
		if correction > 1e-10: # Check if rounding error is outside of tolerances
			probabilitylistfixed = []
			for x in probabilitylist:
				y = x - correction
				probabilitylistfixed.append(y)
		else: 
			probabilitylistfixed = probabilitylist 
	elif sum < 1:
		correction = (1 - sum) / len(probabilitylist) # Calculate how much probabilities observed fall short of one and calculate a correction value
		if correction > 1e-10: # Check if rounding error is outside of tolerances
			probabilitylistfixed = []
			for x in probabilitylist:
				y = x + correction
				probabilitylistfixed.append(y)
		else:
			probabilitylistfixed = probabilitylist
	else:
		probabilitylistfixed = probabilitylist	
	
	for index, value in enumerate(probabilitylistfixed):
		if value < 0: 
			probabilitylistfixed[index] = 0
			probabilitylistfixed[probabilitylistfixed.index(max(probabilitylistfixed))] += value # This code is to deal with rounding error correction when it causes a negative value... add difference to highest bin
	
	sample = numpy.random.choice(binlist, size = bootstraptotal, replace=True, p = probabilitylistfixed) # This is the actual sampling procedure
	return sample

# Calculate a matrix of species distances for each variable
# The purpose of matrix normalization is to enforce equal weight among ecological variables
process = 1
while process <= processtotal:	
	distanceDataFrame = pandas.DataFrame(data=numpy.zeros((len(specieslist), len(specieslist))), index=specieslist, columns=specieslist) # Empty dataframe
	for columnspecies in specieslist:
		for rowspecies in specieslist:
			rowsample = SampleDistribution(process, bootstraptotal, rowspecies)
			columnsample = SampleDistribution(process, bootstraptotal, columnspecies)
			speciesdistance = distance.euclidean(rowsample,columnsample)
			distanceDataFrame.set_value(rowspecies, columnspecies, speciesdistance)
	for species in specieslist:
		distanceDataFrame.set_value(species, species, None)
	distanceDataFrame.to_csv('./distancematrix_raw_pno{0}.csv'.format(process))
	distanceDataFrameNormalized = (distanceDataFrame - distanceDataFrame.mean()) / (distanceDataFrame.max() - distanceDataFrame.min())
	distanceDataFrameNormalized.to_csv('./distancematrix_normalized_pno{0}.csv'.format(process))
	process += 1

# Define empty dataframe 
distanceAverageDataFrame = pandas.DataFrame(data=numpy.zeros((len(specieslist), len(specieslist))), index=specieslist, columns=specieslist) # Empty dataframe

# Get all output for next step
filelist = glob.glob('./distancematrix_normalized_pno*.csv')
framelist = []

# Open all normalized distance matrices, represent as a list
for file in filelist:
	try:
		distanceframe = pandas.DataFrame.from_csv(file)			
	except: 
		print("Error finding distance matrix {0}.".format(file))
		exit()
	framelist.append(distanceframe)
	distanceframe = []

# Get a new matrix that is the average of the normalized distance matrices	
templist = []
for columnspecies in specieslist:
	for rowspecies in specieslist:
		for frame in framelist:
			templist.append(frame.get_value(rowspecies, columnspecies))
		average = numpy.average(templist)
		distanceAverageDataFrame.set_value(rowspecies, columnspecies, average)
		templist = []

# Write final matrix to csv
distanceAverageDataFrame.to_csv('./distancematrix_normalized_averaged.csv')				

# Reopen final matrix (to avoid rerunning the entire script)
distanceAverageDataFrame = pandas.DataFrame.from_csv('./distancematrix_normalized_averaged.csv')

# Perform the clustering
fixedmatrix = distanceAverageDataFrame.fillna(0) # Scipy can't handle NaN for this
numpyMatrix = fixedmatrix.as_matrix()
centroids,_ = kmeans(numpyMatrix, clusternumber) # Won't work without ,_ since it returns a tuple
classification,_ = vq(numpyMatrix, centroids) # Won't work without ,_

# Pair classification list with species list
classificationFinal = []
for x, y in zip (specieslist, classification):
	classificationFinal.append([x, y])

# Save classification to csv
with open("final_classification.csv", "w+") as file:
    writer = csv.writer(file)
    writer.writerows(classificationFinal)	