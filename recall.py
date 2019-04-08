############################################################
# INTESTAZIONE DA DECIDERE

# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	

# This library has been developed to generate a vast dataset of  .fits images
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
import pathlib
from astropy.io import fits


# SELF BUILT LIBRARIES
import map_creator

############################################################
# This library aims to create a vast dataset of .fits images containing gamma ray sources. Each map has
# at most n_sources placed randomly in the image area. 
# The library is organized as follows:
#		1)	function main() :				this function is a sort of iterative wrapper of all the
#										methods in map_creator.py, allowing the user to   
#										generate a dataset through some parameters
#		2)	function usage() :				this function prints an inline help by calling 
#										"python recall_inaf.py -h"										
############################################################


############################################################
# This function works as a wrapper for the data analysis procedure. All the functions that will be used 
# here can also be launched as  standalone in a python environment. 
############################################################
def main():

	# Calling of the function usage()
	usage()

	# Open a logfile where all the data referring to the generated sources will be stored. The file  
	# will be named "Map.log". Before writing on the file the programm will erase any previous data.
	
	logfile = "Map.log"
	open(logfile, 'w').close()
	log = open(logfile, 'a')
	
	
	# The program will now check if the folders where the data will be stored in are actually already 
	# existing and/or empty. If this is the case previous data will be erased, otherwise the folders 
	# will be created.	
	pwd = os.getcwd()
	gen_maps = pwd + '/Generated_Maps'
	gen_fits = pwd + '/Generated_Fits'
	gen_events = pwd + '/Generated_Events'
	gen_png = pwd + '/Generated_Png'	
	
	pathlib.Path(gen_maps).mkdir(parents=True, exist_ok=True)
	if not os.listdir(gen_maps) == []:
		os.system('rm -r Generated_Maps/*')

	pathlib.Path(gen_fits).mkdir(parents=True, exist_ok=True)
	if not os.listdir(gen_fits) == []:
		os.system('rm -r Generated_Fits/*')
	
	pathlib.Path(gen_events).mkdir(parents=True, exist_ok=True)
	if not os.listdir(gen_events) == []:
		os.system('rm -r Generated_Events/*')
	
	pathlib.Path(gen_png).mkdir(parents=True, exist_ok=True)
	if not os.listdir(gen_png) == []:
		os.system('rm -r Generated_Png/*')

	# file layout
	log.write('GENERATED COORDINATES AND INTENSITIES\n\n')
	
	# By calling "find_coordinates" function the program will retrieve the coordinates of the sky 
	# region the telescope is pointing at.
	n_maps = int(sys.argv[4])
	coordinates = map_creator.find_coordinates(log)

	log.write('TOTAL NUMBER OF MAPS:\t' + str(n_maps) + '\n\n')	
	log.write('FILE\t\t' + 'BACKGROUND\t' + 'N_SOURCES\t' + 'RA\t\t' + 'DEC\t\t' + \
			'INTENSITY\n\n')


############################################################
	# DATA GENERATION PROCEDURE
	# Each map will be namend ina progressive way and will contain a randomly generated set of 
	# sources, whose number is limited by the var "n_sources"  passed as the third parameter to the 
	# program. The background level is also generated. After that, the program "ExecuteCTools.py" will 
	# be invoked and an image file "cube.fits" for each map will be created and then stored in a 
	# dedicated folder. The files that will be stored in the those directorie are:
	# 1) Map_i.xml			(source file)
	# 2) Map_i.fits			(fits image - open with ds9)
	# 3) Events_i.fits		(original data structure - stored in a .fits object)
	# 4) SelectedEvents_i.fits	(data selected to create the final image - .fits object)
	# 5) skyi_image.png		(png version of the final image)	############################################################

	for i in range (0, n_maps): 
		filename = "Map_" + str(i+1) + '.xml'
		mapname = "Map_" + str(i+1) + ".fits"
		eventname = "Events_" + str(i+1) + ".fits"
		seleventname = "SelectedEvents_" + str(i+1) + ".fits"
		logname = "ExecuteCTools_" + str(i+1) + ".log"
		skyname = "sky" + str(i+1) + "_image.png"
		log.write('\n' + filename + '\t\t')
		background_prefactor = map_creator.background_generation(map_creator.min_background, map_creator.max_background, log)
		source_array = map_creator.source_generation(coordinates, int(sys.argv[3]), 
background_prefactor, log)
		map_creator.sourcefile_generator(source_array, background_prefactor,  filename)
		os.system("mv " + filename + " Generated_Maps/" + filename)
		print(mapname + ' : Generating\n')
		os.system("python ExecuteCTools.py -simmodel Generated_Maps/" + filename + \
				  " -observation " + str(sys.argv[1]) + " -runconf " + str(sys.argv[2]) + " &> " + logname )
		os.system("cp cube.fits Generated_Fits/" + mapname)
		os.system("mv " + logname + " Generated_Events/" + logname)
		os.system("cp sky1.png Generated_Png/" + skyname)


############################################################
		# this part is optional due to the big amount of memory occupied by events.fits files. Use only
		# if you have plenty of free space in your HDD/SSD (nearly 300 MB per file)
		#
		# CODE:
		#		
		# os.system("cp events.fits Generated_Events/" + eventname)
		# os.system("cp selected_events.fits Generated_Events/" + seleventname)	############################################################
		
	
	# directory management
	os.system("mv Map.log Generated_Maps/Map.log")
		
	os.system("rm cube.fits")
	os.system("rm events.fits")
	os.system("rm selected_events.fits")
	os.system("rm sky1.png")
	
	log.close()	

############################################################
# function devoted to generate the parameters of the main program and the related --help function. 
# All the parameters needed by the main function are described below
############################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to create random \
								     gamma sources in a given sky region')
	parser.add_argument('observation input file', metavar = 'filename', type=str, \
						help='source file used by the program' )
	parser.add_argument('configuration input file', metavar = 'fileconf', type=str, \
						help='configuration file used by the program' )
	parser.add_argument('max sources number', metavar = 'n_sources_max', type=int, \
						help='maximum number of sources per image' )
	parser.add_argument('number of generated maps', metavar = 'n_maps', type=int, \
						help='number of generated maps' )
	args = parser.parse_args()


if __name__ == "__main__":
    main()
