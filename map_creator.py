########################################## MAP_CREATOR.PY ###########################################
# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	

# 1.7903
#################################################################################################

###########################################  IMPORTS ###############################################
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

############################################ GLOBAL VARIABLES ########################################
min_background = 1e-3
max_background = 1
background_multiplicator = 1
sources_spread = 3

######################################## FUNCTION  MAIN ############################################# 
def main():
	# Calling of the function usage()
	usage()

	'''
	Section devoted to: 
		i) coordinates detecting: 	read the coordinates of the center (telescope pointing zone) from the .xml file
	       ii) background generation: 	random generation of the background intensity within background_min and background_max (global vars)
	      iii) source generation:		random generation of the number of sources per map, coordinates (RA and DEC) and intensity of the source itself.
	      iv) logfile creation:			writing of the obtained data on the designed logfile (see also "recall.py")
	'''
	log = sys.argv[4]
	open(log, 'w').close()
	logfile = open(log, 'a')
	coordinates = find_coordinates(logfile)
	#background_prefactor = background_generation(min_background, max_background, logfile)
	background_prefactor = 10
	source_array = source_generation(coordinates, int(sys.argv[2]), background_prefactor, logfile)
	sourcefile_generator(source_array, background_prefactor, sys.argv[3])


######################################## FUNCTION USAGE #############################################
# function devoted to generate the parameters of the main program and the related --help function (call ">> python map_creator.py -h" ). 
#################################################################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to create random gamma sources in a given sky region')
	parser.add_argument('input file', metavar = 'filename', type=str, help='source file used by the program' )
	parser.add_argument('max sources number', metavar = 'n_sources_max', type=int, help='maximum number of sources per image' )
	parser.add_argument('source file name', metavar = 'sourcefile', type=str, help='generated source file name' )
	parser.add_argument('name of the logfile', metavar = 'logfile', type=str, help='name of the logfile where generated data will be saved' )
	args = parser.parse_args()

################################### FUNCTION FIND_COORDINATES #########################################
# function devoted to obtain the coordinates of the center of the sky-region pointed by the telescope  
#################################################################################################
def find_coordinates(logfile):
	
	# opening of the input file where center coordinates are stored (sys.argv[1])
	filename = sys.argv[1]
	obs_crab = open(filename, 'r')
	obs_crab_string = obs_crab.read()

	# writing logfile layout
	logfile.write('COORDINATES CENTER:\t')
	
	# RA detection in the opened logfile through the key word "ra=" + var ra_str init
	word = 'ra="'
	index = obs_crab_string.index(word)
	ra_str = ""

	for i in range (4,7):
		ra_str = ra_str + obs_crab_string[index + i]

	# integer casting of the returned value + logfile writing
	RA = int(ra_str)
	logfile.write('\t\t\tRA = ' + str(RA) + '\t')

	# DEC detection in the opened logfile through the key word "dec=" + var dec_str init
	word = 'dec="'
	index = obs_crab_string.index(word)	
	dec_str = ""

	for i in range (5,7):
		dec_str = dec_str + obs_crab_string[index + i]
	
	# integer casting of the returned value + var dec_str init
	DEC = int(dec_str)
	logfile.write('DEC = ' + str(DEC) + '\n\n')

	# telescope center coordinates assignement
	position = [RA, DEC]
	obs_crab.close()

	return position


###################################### FUNCTION SOURCE_GENERATION #####################################
# function devoted to randomly generate astronomic coordinates and intensity of gamma ray sources. Also the number of sources per map is randomly generated.
#################################################################################################
def source_generation( position, n_sources_max, background_prefactor, logfile):
	
	'''
	number of sources generation (random generation)
	definition af the array where the generated values will be stored. This array has as many elements as the number of sources in the map. Moreover, each element has 3 fields containing the RA, DEC and intensity values
	definition of the spatial and intesity ranges within which the values will be generated
	'''
	n_sources = random.randint(0,n_sources_max)	
	n_sources = 1
	source_array = numpy.zeros((n_sources, 3))	
	ra_range = [position[0] - 2.8, position[0] + 2.8]
	dec_range = [position[1] - 2, position[1] + 2]
	int_range = [background_prefactor*background_multiplicator, background_prefactor*background_multiplicator*sources_spread]
	
	# logfile writing
	logfile.write(str(n_sources) + '\n')
	
	'''
	i)	right ascension generation
	ii)	declination generation
	iii)	intensity generation
	iv)	logfile writing
	'''
	for i in range(0,n_sources):
		source_array[i][0] = round(random.uniform(ra_range[0], ra_range[1]),2)
		source_array[i][1] = round(random.uniform(dec_range[0], dec_range[1]),2)
		source_array[i][2] = round(random.uniform(int_range[0], int_range[1]),2)

	# coordinates ordering by RA
	for i in range(0,n_sources):
		for j in range(0,n_sources - 1 - i):
			if ( (source_array[j][1] < source_array[j+1][1] ) or (  (source_array[j][1] == source_array[j+1][1] ) and (source_array[j][0] < source_array[j+1][0]) ) ):
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
		logfile.write('\t\t\t\t\t\t\t\t\t\t' + str(source_array[i][0]) + '\t' + str(source_array[i][1]) + '\t\t' + str(source_array[i][2]) + '\n')

	return source_array


################################### FUNCTION SOURCE_GENERATION ##################################
# function devoted to randomly generate astronomic coordinates and intensity of gamma ray sources - INAF REQUIREMENTS
###########################################################################################
def source_generation_INAF(position, background_prefactor, sigma, logfile):

	'''
	number of sources is set to 1
	definition af the array where the generated values will be stored. This array has as many elements as the number of sources (1 in this case) in the map. Moreover, each element has 3 fields containing the RA, DEC and intensity values
	definition of the spatial and intesity ranges within which the values will be generated
	'''
	n_sources = random.randint(0,1)
	source_array = numpy.zeros((n_sources, 3))
	ra_range = [position[0] - 2, position[0] + 2]
	dec_range = [position[1] - 1.5, position[1] + 1.5]

	# logfile writing
	logfile.write(str(n_sources) + '\n')
	
	'''
	i)	right ascension generation
	ii)	declination generation
	iii)	intensity generation
	iv)	logfile writing
	'''
	for i in range(0,n_sources):
		source_array[i][0] = round(random.uniform(ra_range[0], ra_range[1]),2)
		source_array[i][1] = round(random.uniform(dec_range[0], dec_range[1]),2)
		source_array[i][2] = sigma*math.sqrt(background_prefactor)

	for i in range(0,n_sources):
		for j in range(0,n_sources - 1 - i):
			if ( (source_array[j][1] < source_array[j+1][1] ) or (  (source_array[j][1] == source_array[j+1][1] ) and (source_array[j][0] < source_array[j+1][0]) ) ):
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
		logfile.write('\t\t\t\t\t\t\t\t\t\t' + str(source_array[i][0]) + '\t' + str(source_array[i][1]) + '\t\t' + str(source_array[i][2]) + '\n')

	return source_array	

################################## FUNCTION BACKGROUND GENERATION ################################
# function devoted to randomly generate the total background in the sky-region
###########################################################################################
def background_generation(min_val, max_val, logfile):
	prefactor = round(random.uniform(min_val, max_val),4)
	logfile.write(str(prefactor) + '\t\t')	
	return prefactor

#################################### FUNCTION SOURCE_GENERATION #################################
# function devoted to generate the source file that will be used by the image generation code  ( see ExecuteCtools.py )
###########################################################################################
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
		#source.write("       <parameter name=\"Prefactor\"   scale=\"1e-17\" value=\"" + str(source_array[i][2]) + "\"  min=\"" + str(background_prefactor*background_multiplicator) + "\" max=\"1000.0\" free=\"1\"/>\n")
		source.write("       <parameter name=\"Prefactor\"   scale=\"1e-17\" value=\"2\" min=\"1e-07\" max=\"1000.0\" free=\"1\"/>\n")
		source.write("       <parameter name=\"Index\"       scale=\"-1\"    value=\"2.48\" min=\"0.0\"   max=\"+5.0\"   free=\"1\"/>\n")
		source.write("       <parameter name=\"PivotEnergy\" scale=\"1e6\"   value=\"0.3\"  min=\"0.01\"  max=\"1000.0\" free=\"0\"/>\n")
		source.write("    </spectrum>\n")
		source.write("    <spatialModel type=\"PointSource\">\n")
		source.write("      <parameter name=\"RA\"  scale=\"1.0\" value=\"" + str(source_array[i][0]) + "\" min=\"-360\" max=\"360\" free=\"0\"/>\n")
		source.write("      <parameter name=\"DEC\" scale=\"1.0\" value=\"" + str(source_array[i][1]) + "\" min=\"-90\"  max=\"90\"  free=\"0\"/>\n")
		source.write("    </spatialModel>\n")
		source.write("  </source>\n")

	# background data writing
	source.write("  <source name=\"CTABackgroundModel\" type=\"CTAIrfBackground\" instrument=\"CTA\">\n")
	source.write("    <spectrum type=\"PowerLaw\">\n")
	#source.write("      <parameter name=\"Prefactor\"   scale=\"1.0\"  value=\"" + str(background_prefactor) + "\"  min=\"1e-3\" max=\"10\"   free=\"1\"/>\n")
	source.write("      <parameter name=\"Prefactor\"   scale=\"1.0\"  value=\"1\"  min=\"1e-3\" max=\"10\"   free=\"1\"/>\n")
	source.write("      <parameter name=\"Index\"       scale=\"1.0\"  value=\"0.0\"  min=\"-5.0\" max=\"+5.0\"   free=\"1\"/>\n")
	source.write("      <parameter name=\"PivotEnergy\" scale=\"1e6\"  value=\"1.0\"  min=\"0.01\" max=\"1000.0\" free=\"0\"/>\n")
	source.write("    </spectrum>\n")
	source.write("  </source>\n")
	source.write("</source_library>")

	source.close()

if __name__ == "__main__":
    main()


