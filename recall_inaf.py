########################################## RECALL.PY ################################################
# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	
#################################################################################################
import map_creator
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

import time
from tqdm import tqdm



###################################### FUNCTION MAIN ###############################################
# The main function
#################################################################################################
def main():
	#################################
	# Calling of the function usage()
	#################################
	usage()

	###############################################################################################
	# Open a logfile where all the data referring to the generated sources will be stored. The file  
	# will be named "Map.log". Before writing on the file the programm will erase any previous data.
	###############################################################################################
	logfile = "Map.log"
	open(logfile, 'w').close()
	log = open(logfile, 'a')
	
	###############################################################################################
	# The program will now check if the folders where the data will be stored in are actually already 
	# existing and/or empty. If this is the case previous data will be erased, otherwise the folders 
	# will be created.
	###############################################################################################
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

	#############
	# file layout
	#############
	log.write('GENERATED COORDINATES AND INTENSITIES\n\n')

	#################################################################################################
	# By calling "find_coordinates" function the program will retrieve the coordinates of the sky 
	# region the telescope is pointing at.
	#################################################################################################
	n_maps = int(sys.argv[4])
	sigma = float(sys.argv[3])

	coordinates = numpy.zeros(2)
	coordinates[0] = float(sys.argv[1])
	coordinates[1] = float(sys.argv[2])
	log.write('TOTAL NUMBER OF MAPS:\t' + str(n_maps) + '\n\n')
	
	log.write('FILE\t\t' + 'BACKGROUND\t' + 'N_SOURCES\t' + 'RA\t\t' + 'DEC\t\t' + 'INTENSITY\n\n')

	#################################################################################################
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
	# 5) skyi_image.png		(png version of the final image)
	#################################################################################################
	for i in tqdm(range(0, n_maps)): 
		filename = "Map_" + str(i+1) + '.xml'
		mapname = "Map_" + str(i+1) + ".fits"
		eventname = "Events_" + str(i+1) + ".fits"
		seleventname = "SelectedEvents_" + str(i+1) + ".fits"
		ctobssim_logname = 'ctobssim_' + str(i+1) + '.log'
		logname = "ExecuteCTools_" + str(i+1) + ".log"
		skyname = "sky" + str(i+1) + "_image.png"
		log.write('\n' + filename + '\t\t')
		background_prefactor = map_creator.background_generation(map_creator.min_background, sigma**2, log)
		source_array = map_creator.source_generation_INAF(coordinates, background_prefactor, sigma, log)
		map_creator.sourcefile_generator(source_array, background_prefactor,  filename)
		os.system("mv " + filename + " Generated_Maps/" + filename)
		#print(mapname + ' : Generating\n')
		#os.system("python ExecuteCTools.py -simmodel Generated_Maps/" + filename + " -observation " + str(sys.argv[1]) + " -runconf " + str(sys.argv[2]) + ' -seed ' + str(i) + " &> " + logname )
	
		os.system( 'ctobssim inmodel=Generated_Maps/' + filename + ' ra=' + str(coordinates[0]) + ' dec=' + str(coordinates[1]) +  ' rad=10 tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:15:00 emin=0.1 emax=100.0 caldb=prod2 irf=South_0.5h outevents=events.fits seed=' + str(i) )
		os.system( 'ctbin inobs=events.fits coordsys=CEL proj=CAR xref=' + str(coordinates[0]) + ' yref=' + str(coordinates[1])  + ' binsz=0.02 nxpix=200 nypix=200 ebinalg=LOG emin=0.1 emax=100 enumbins=1 outcube=cube.fits'  )

		os.system("cp cube.fits Generated_Fits/" + mapname)
		#os.system("mv " + logname + " Generated_Events/" + logname)
		#os.system("cp sky1.png Generated_Png/" + skyname)

		#################################################################################################
		#this part is optional due to the big amount of memory occupied by events.fits files. Use only
		# if you have plenty of free space in your HDD/SSD (nearly 300 MB per file)
		#
		# CODE:
		#
		os.system("cp events.fits Generated_Events/" + eventname)
		os.system("cp ctobssim.log Generated_Events/" + ctobssim_logname )
		#os.system("cp selected_events.fits Generated_Events/" + seleventname)
		# 
		#################################################################################################
		
	
	os.system("mv Map.log Generated_Maps/Map.log")
	
	
	os.system("rm cube.fits")
	os.system("rm events.fits")
	os.system("rm ctobssim.log")
	os.system("rm ctbin.log")
	#os.system("rm selected_events.fits")
	#os.system("rm sky1.png")
	
	log.close()	


###################################### FUNCTION USAGE ############################################
# function devoted to generate the parameters of the main program and the related --help function. 
##################################################################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to create random gamma sources in a given sky region')
	parser.add_argument('RA source coordinate', metavar = 'RA', type=str, help='RA source coordinate - celestial coordinates' )
	parser.add_argument('DEC source coordinate', metavar = 'DEC', type=str, help='RA source coordinate - celestial coordinates' )
	parser.add_argument('source significativity', metavar = 'sigma', type=float, help='ratio between source intensity and the square root of the background intensity' )
	parser.add_argument('number of generated maps', metavar = 'n_maps', type=int, help='number of generated maps' )
	args = parser.parse_args()


if __name__ == "__main__":
    main()
