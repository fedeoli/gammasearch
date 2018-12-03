##############################################################################
# INTESTAZIONE DA DECIDERE

# This library has been developed to perform statistical analysis on big amounts of source images
##############################################################################

# Imports

import data_analysis
import order
import map_creator
import os
import sys
import argparse

##############################################################################
# This function allows to analyze a general folder full of .fits file, writing the results on a logfile "est_map.log"
##############################################################################
def general_analysis():

	usage()

	# parameters initialisation from argv (see usage() documentation)
	fits_dir = sys.argv[1]
	dest_logfile = sys.argv[2]
	center_type = sys.argv[3]
	sigma_spa = sys.argv[4]
	accept_level = sys.argv[5]
	radius = sys.argv[6]
	binary_treshold = sys.argv[7]	
	intensity_treshold = sys.argv[8]	
	baricenter_distance = sys.argv[9]

	tmp = os.listdir(fits_dir)
	n_maps = len(tmp)

	[est_n_sources_array, est_coord_array] = data_analysis.convert_data(n_maps, center_type, sigma_spa, accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)

	os.remove('logfile_map')

	data_analysis.write_data(est_n_sources_array, est_coord_array, fits_dir)


##############################################################################
# function devoted to generate the parameters of the main program and the related --help function. All the parameters needed by the 
# main function are described below
##############################################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to analyze .fits maps and write the results in a logfile')
	parser.add_argument('Fits folder', metavar = 'fits_dir', type=str, help='Directory where .fits file are saved' )
	parser.add_argument('final result log file name', metavar = 'final_logfile', type=str, help='algorithm final result log file name' )
	parser.add_argument('type of centers computation', metavar = 'center type', type=str, help='type of centers computation. Even multiple input: (e.g.)  bim = baricenter + intensity + mean' )
	parser.add_argument('spatial gaussian sigma', metavar = 'sigma_spa', type=float, help='variance of the smoothing gaussian' )
	parser.add_argument('percentage of accepted area', metavar = 'accept_level', type=float, help='percentage to reach in circle detection area check' )
	parser.add_argument('mask circle radius', metavar = 'radius', type=int, help='template matching radius' )
	parser.add_argument('binarization treshold (255 valued)', metavar = 'binary_treshold', type=int, help='treshold used for image binarization' )
	parser.add_argument('intensity treshold ', metavar = 'intensity_treshold', type=float, help='background normalization treshold for intensity information handling' )
	parser.add_argument('blob centers minimum distance', metavar = 'baricenter_distance', type=int, help='minimum distance for two blob centers to be considered as separated' )
	args = parser.parse_args()

if __name__ == "__main__":
    general_analysis()

	
	
