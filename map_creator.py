############################################################
# INTESTAZIONE DA DECIDERE

# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	

# This library has been developed to generate a single .fits image containing gamma ray sources
############################################################


# IMPORTS
# All the libraries can be installed by pip except for:
#	ctools: download and install from
#	gammalib: download and install from
#	cscripts: download and install from 

import os 
import ctools
import cscripts
import gammalib
import sys
import argparse
import string
import random
import numpy
import math
from astropy.io import fits

# GLOBAL VARIABLES
min_background = 1e-3
max_background = 1
background_multiplicator = 1
sources_spread = 3

############################################################
# This library aims to create a .fits file representing a map of a limited region of sky possibly containing
#  gamma ray sources.
# The library is organized as follows:
#		1)	function main() :				this function is a sort of wrapper, allowing the user to   
#										create a .fits file through some parameters
#		2)	function usage() :				this function prints an inline help by calling 
#										"python map_creator.py -h"										
#		3)	function source_generation(): 	function devoted to randomly generate astronomic 
#										coordinates and intensity of gamma ray sources. Also the 
#										number of sources per map is randomly generated.
#		4)	function source_generation_INAF() : this function is imilar to the previous one but this 
#										one meets some specific requirements - INAF req.
#		5) 	function background_generation() :	function devoted to randomly generate the total 
#										background in the sky-region
#		6)	function sourcefile_generator() : function devoted to generate the source file that will 
#										be used by the image generation code  
############################################################			 	


############################################################
# This function works as a wrapper for the data analysis procedure. All the functions that will be used 
# here can also be launched as  standalone in a python environment. 
############################################################
def main():

	# # call the usage function in order to print the inline help if needed
	usage()

############################################################
#	Section devoted to: 
#		i) coordinates detecting: 	read the coordinates of the center (telescope pointing zone) 
#								from the .xml file
#	       ii) background generation: 	random generation of the background intensity within 
#								background_min and background_max (global vars)
#	      iii) source generation:		random generation of the number of sources per map, 
#								coordinates (RA and DEC) and intensity of the source itself.
#	      iv) logfile creation:			writing the obtained data on designed logfile (see also "recall.py")
#
#	NB:	now the source generation is done by the source_generation method. Remember that there's
#		also the source_generation_INAF one. 
############################################################
	log = sys.argv[4]
	open(log, 'w').close()
	logfile = open(log, 'a')
	coordinates = find_coordinates(logfile)
	background_prefactor = 10
	source_array = source_generation(coordinates, int(sys.argv[2]), background_prefactor, logfile)
	sourcefile_generator(source_array, background_prefactor, sys.argv[3])

############################################################
# function devoted to generate the parameters of the main program and the related --help function. 
# All the parameters needed by the main function are described below
############################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to create random gamma \
								    sources in a given sky region')
	parser.add_argument('input file', metavar = 'filename', type=str, help='source file used by \
						the program' )
	parser.add_argument('max sources number', metavar = 'n_sources_max', type=int, \
						help='maximum number of sources per image' )
	parser.add_argument('source file name', metavar = 'sourcefile', type=str, \
						help='generated source file name' )
	parser.add_argument('name of the logfile', metavar = 'logfile', type=str, \
						help='name of the logfile where generated data will be saved' )
	args = parser.parse_args()


 ###########################################################
# function devoted to randomly generate astronomic coordinates and intensity of gamma ray sources. 
# Also the number of sources per map is randomly generated.
############################################################
def source_generation( position, n_sources_max, background_prefactor, logfile):
	
############################################################
#	number of sources generation (random generation)
#	definition af the array where the generated values will be stored. This array has as many elements 
#as the number of sources in the map. Moreover, each element has 3 fields containing the RA, DEC and 
#intensity values definition of the spatial and intesity ranges within which the values will be generated
############################################################
	n_sources = random.randint(0,n_sources_max)	
	n_sources = 1
	source_array = numpy.zeros((n_sources, 3))	
	ra_range = [position[0] - 2.8, position[0] + 2.8]
	dec_range = [position[1] - 2, position[1] + 2]
	int_range = [background_prefactor*background_multiplicator, background_prefactor*background_multiplicator*sources_spread]
	
	# logfile writing
	logfile.write(str(n_sources) + '\n')
	
############################################################
#	i)	right ascension generation
#	ii)	declination generation
#	iii)	intensity generation
#	iv)	logfile writing
############################################################
	for i in range(0,n_sources):
		source_array[i][0] = round(random.uniform(ra_range[0], ra_range[1]),2)
		source_array[i][1] = round(random.uniform(dec_range[0], dec_range[1]),2)
		source_array[i][2] = round(random.uniform(int_range[0], int_range[1]),2)

	# coordinates ordering by RA
	for i in range(0,n_sources):
		for j in range(0,n_sources - 1 - i):
			if ( (source_array[j][1] < source_array[j+1][1] ) or (  (source_array[j][1] == \
				source_array[j+1][1] ) and (source_array[j][0] < source_array[j+1][0]) ) ):
				temp = source_array[j][0]
				source_array[j][0] = source_array[j+1][0]
				source_array[j+1][0] = temp
				temp = source_array[j][1]
				source_array[j][1] = source_array[j+1][1]
				source_array[j+1][1] = temp
				temp = source_array[j][2]
				source_array[j][2] = source_array[j+1][2]
				source_array[j+1][2] = temp

	# logfile writing
	for i in range(0,n_sources):
		logfile.write('\t\t\t\t\t\t\t\t\t\t' + str(source_array[i][0]) + '\t' + \
					str(source_array[i][1]) + '\t\t' + str(source_array[i][2]) + '\n')

	return source_array


############################################################
# function devoted to randomly generate astronomic coordinates and intensity of gamma ray sources 
# This functions follows some specific requirements - INAF REQUIREMENTS
############################################################
def source_generation_INAF(position, background_prefactor, sigma, logfile):

############################################################
#	number of sources is set to 1
#	definition af the array where the generated values will be stored. This array has as many elements 
#	as the number of sources (1 in this case) in the map. Moreover, each element has 3 fields containing 
#	the RA, DEC and intensity values definition of the spatial and intesity ranges within which the 
#	values will be generated
############################################################
	n_sources = random.randint(0,1)
	source_array = numpy.zeros((n_sources, 3))
	ra_range = [position[0] - 0.5, position[0] + 0.5]
	dec_range = [position[1] - 0.5, position[1] + 0.5]

	# logfile writing
	logfile.write(str(n_sources) + '\n')
	
############################################################
#	i)	right ascension generation
#	ii)	declination generation
#	iii)	intensity generation
#	iv)	logfile writing
############################################################
	for i in range(0,n_sources):
		source_array[i][0] = round(random.uniform(ra_range[0], ra_range[1]),2)
		source_array[i][1] = round(random.uniform(dec_range[0], dec_range[1]),2)
		source_array[i][2] = (sigma**2)*(1 + math.sqrt(sigma**2 + 4*background_prefactor))/2

	for i in range(0,n_sources):
		for j in range(0,n_sources - 1 - i):
			if ( (source_array[j][1] < source_array[j+1][1] ) or (  (source_array[j][1] == \
				source_array[j+1][1] ) and (source_array[j][0] < source_array[j+1][0]) ) ):
				temp = source_array[j][0]
				source_array[j][0] = source_array[j+1][0]
				source_array[j+1][0] = temp
				temp = source_array[j][1]
				source_array[j][1] = source_array[j+1][1]
				source_array[j+1][1] = temp
				temp = source_array[j][2]
				source_array[j][2] = source_array[j+1][2]
				source_array[j+1][2] = temp

	# logfile writing
	for i in range(0,n_sources):
		logfile.write('\t\t\t\t\t\t\t\t\t\t' + str(source_array[i][0]) + '\t' + \
					str(source_array[i][1]) + '\t\t' + str(source_array[i][2]) + '\n')

	return source_array	


############################################################
# function devoted to randomly generate the total background in the sky-region
############################################################
def background_generation(min_val, max_val, logfile):
	prefactor = round(random.uniform(min_val, max_val),4)
	logfile.write(str(prefactor) + '\t\t')	
	return prefactor

############################################################
# function devoted to generate the source file that will be used by the image generation code  
# ( see ExecuteCtools.py )
############################################################
def sourcefile_generator(source_array, background_prefactor, filename):

	# number of sources detection (from input source_array)
	size = len(source_array)

	# file writing procedure
	source = open(filename, 'w')
	source.write("<?xml version=\"1.0\" standalone=\"no\"?>\n")
	source.write("<source_library title=\"source library\">\n")
	
	# source data writing
	for i in range(0,size):
		source.write("  <source name=\"" + "Source_" + str(i+1) + "\" type=\"PointSource\">\n")
		source.write("    <spectrum type=\"PowerLaw\">\n")
		source.write("       <parameter name=\"Prefactor\"   \
					scale=\"1e-17\" value=\"2\" min=\"1e-07\" max=\"1000.0\" free=\"1\"/>\n")
		source.write("       <parameter name=\"Index\"       \
					scale=\"-1\"    value=\"2.48\" min=\"0.0\"   max=\"+5.0\"   free=\"1\"/>\n")
		source.write("       <parameter name=\"PivotEnergy\" scale=\"1e6\"   \
					value=\"0.3\"  min=\"0.01\"  max=\"1000.0\" free=\"0\"/>\n")
		source.write("    </spectrum>\n")
		source.write("    <spatialModel type=\"PointSource\">\n")
		source.write("      <parameter name=\"RA\"  scale=\"1.0\" value=\"" + \
str(source_array[i][0]) + "\" min=\"-360\" max=\"360\" free=\"0\"/>\n")
		source.write("      <parameter name=\"DEC\" scale=\"1.0\" value=\"" + \
str(source_array[i][1]) + "\" min=\"-90\"  max=\"90\"  free=\"0\"/>\n")
		source.write("    </spatialModel>\n")
		source.write("  </source>\n")

	# background data writing
	source.write("  <source name=\"CTABackgroundModel\" type=\"CTAIrfBackground\" \
				instrument=\"CTA\">\n")
	source.write("    <spectrum type=\"PowerLaw\">\n")
	source.write("      <parameter name=\"Prefactor\"   scale=\"1.0\"  value=\"1\"  \
				min=\"1e-3\" max=\"10\"   free=\"1\"/>\n")
	source.write("      <parameter name=\"Index\"       scale=\"1.0\"  value=\"0.0\"  \
				min=\"-5.0\" max=\"+5.0\"   free=\"1\"/>\n")
	source.write("      <parameter name=\"PivotEnergy\" scale=\"1e6\"  value=\"1.0\"  \
				min=\"0.01\" max=\"1000.0\" free=\"0\"/>\n")
	source.write("    </spectrum>\n")
	source.write("  </source>\n")
	source.write("</source_library>")

	# close the file
	source.close()


############################################################
# This function simply takes a text file name as input and stores in an array all the 
# words present in that file
############################################################
if __name__ == "__main__":
    main()


