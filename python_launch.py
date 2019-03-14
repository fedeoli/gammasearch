##############################################################################
# INTESTAZIONE DA DECIDERE

# This library has been developed to perform statistical analysis on big amounts of source images
##############################################################################


# IMPORTS
# All the libraries can be installed by pip except for:
#	ctools: download and install from
#	gammalib: download and install from
#	cscripts: download and install from 

import matplotlib.pyplot as plt
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
import re
import scipy
import itertools as it
from astropy.io import fits
from subprocess import call
from astropy import wcs


# SELF BUILT LIBRARIES
import map_creator
import order
import data_analysis

def single_analysis():

	# parameters initialisation from argv (see usage() documentation)
	fits_file = sys.argv[1]
	dest_logfile = sys.argv[2]
	center_type = sys.argv[3]
	sigma_spa = sys.argv[4]
	accept_level = sys.argv[5]
	radius = sys.argv[6]
	binary_treshold = sys.argv[7]	
	intensity_treshold = sys.argv[8]	
	baricenter_distance = sys.argv[9]
	show = sys.argv[10]


	est_coord_array = []

	if os.path.isfile('logfile_map'):
		os.system('rm logfile_map' )

	flag = False #first row reading flag
	call(["./final_3.bin", fits_file, center_type, "logfile_map", show, "n", sigma_spa, accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance])
		
	with open('logfile_map', 'r') as f:
		for line in f:
			split_line = line.split()
			if ( not flag ) :
				flag = True
				n_sources = int(split_line[0])
			else :
				temp_dec = float(split_line[0])
				temp_ra = float(split_line[1])
				w  = wcs.WCS( fits_file)
				world = w.wcs_pix2world(temp_ra,200 - temp_dec, 1)
				est_coord_array.append(world[1])
				est_coord_array.append(world[0])
		f.close()
	
		if os.path.isfile(dest_logfile):
			os.system('rm ' + str(dest_logfile))		
	
		dest = open(dest_logfile, 'w+')
		dest.write("N_SOURCES = " + str(n_sources) + " \n\n")
		if (n_sources != 0):
			dest.write("SOURCE\t\t\tRA\t\t\t\tDEC\n")
			for i in range(n_sources):
				dest.write("SOURCE " + str(i) + "\t\t" + str(est_coord_array[2*i + 1]) + "\t\t" + str(est_coord_array[2*i]) + "\n")

		dest.close()

if __name__ == "__main__":
    single_analysis()





