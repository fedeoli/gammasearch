############################################################
# INTESTAZIONE DA DECIDERE

# This library has been developed to perform statistical analysis on big amounts of source images
############################################################


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
import math
from astropy.io import fits
from subprocess import call
from astropy import wcs
from tqdm import tqdm


# SELF BUILT LIBRARIES
import map_creator
import order

############################################################
# This library aims to process big amounts of source images in order to perform some statistical 
# analysis on the number of estimated sources and their estimated positions, expressed in RA and DEC. 
# The library is organized as follows:
#		1)	function main() :				this function is a sort of wrapper, allowing the user to   
#										analyze a bunch of images from scratch, through some 
#										parameters
#		2)	function usage() :				this function prints an inline help by calling 
#										"python data_analysis.py -h"										
#		3)	function words(): 				this function returns an array containing all the words in
#										a given file
#		4)	function find_n_maps() :		this function reads the number of maps from a  given 
#										logfile (see below)
#		5) 	function true_coordinates() :	this functin stores in an array all the coordinates written 
#										in a given logfile (see below)
#		6)	function convert_data() :		this function performs the analysis of the source images 
#										using "analyze_map.bin"  algorithm 
#										(see related documentation)
#		7)	function switcher() :			this function calls "convert_data()" with different 
#										parameters, according to argv[] values
#		8)	function difference() :			this function evaluates the error of the algorithm
#		9)	function statistics() :			this function computes the main statistics from the array 
#										returned by "difference()"
#		10)	function plot_result() :		this function plots the previously computed values
#		11)	function write_data() :		this function writes the estimated coordinates on a logfile
############################################################			 	


############################################################
# This function works as a wrapper for the data analysis procedure. All the functions that will be used 
# here can also be launched as  standalone in a python environment. 
############################################################
def main():

	# call the usage function in order to print the inline help if needed
	usage()

	# parameters initialisation from argv (see usage() documentation)
	logfile = sys.argv[1]
	fits_dir = sys.argv[2]
	dest_logfile = sys.argv[3]
	center_type = sys.argv[4]
	sigma_spa = sys.argv[5]
	accept_level = sys.argv[6]
	radius = sys.argv[7]
	binary_treshold = sys.argv[8]	
	intensity_treshold = sys.argv[9]	
	baricenter_distance = sys.argv[10]
	
	# screen flush
	os.system('clear')
	print('analysing data...')

	# open the logfile (created by recall.py)
	log = open(logfile, 'r')

	# retrieve the number of maps to be analyzed
	n_maps = find_n_maps(logfile)

	# retrieve the real sources coordinates from the given logfile
	[n_sources_array, coord_array] = true_coordinates(logfile, n_maps)

	# analyze all the maps depending on the given parameters (see switcher function)
	[est_n_sources_array_b, est_n_sources_array_i, est_n_sources_array_m, est_coord_array_b, 
est_coord_array_i, est_coord_array_m] = switcher(n_maps, center_type, sigma_spa, accept_level, 
radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)

	# set the number of real sources 
	true_n_sources = int(len(coord_array)/2)

	# error evaluation through function difference and depending on the input parameters
	if center_type.find('b') != -1:
		[diff_array_b, error_b] = difference(n_maps, n_sources_array, est_n_sources_array_b, 
coord_array, est_coord_array_b)
		[mean_dec_b, mean_ra_b, var_dec_b, var_ra_b] = statistics(diff_array_b)
	if center_type.find('i') != -1:
		[diff_array_i, error_i] = difference(n_maps, n_sources_array, est_n_sources_array_i, 
coord_array, est_coord_array_i)
		[mean_dec_i, mean_ra_i, var_dec_i, var_ra_i] = statistics(diff_array_i) 
	if center_type.find('m') != -1:
		[diff_array_m,error_m] = difference(n_maps, n_sources_array, est_n_sources_array_m, 
coord_array, est_coord_array_m) 
		[mean_dec_m, mean_ra_m, var_dec_m, var_ra_m] = statistics(diff_array_m)
	

	# final results presentation
	print('\n\nSTATISTICS EVALUATION\n\n')
	print('N_maps = ' + str(n_maps))
	print('\nN. TRUE SOURCES:	' + str(true_n_sources) + '\n')
	if center_type.find('b') != -1:
		print(str(mean_dec_b) + '\t' + str(mean_ra_b) + '\t' + str(var_dec_b) + '\t' + str(var_ra_b) \
+ '\n')
		est_n_sources_b = int(len(est_coord_array_b)/2)
		print('N. ESTIMATED SOURCES (baricenter):	' + str(est_n_sources_b))
		print('ERROR:\t' + str(error_b))
		write_data(est_n_sources_array_b, est_coord_array_b, fits_dir)
	if center_type.find('i') != -1:
		print(str(mean_dec_i) + '\t' + str(mean_ra_i) + '\t' + str(var_dec_i) + '\t' + str(var_ra_i) \
+ '\n')
		est_n_sources_i = int(len(est_coord_array_i)/2)
		print('N. ESTIMATED SOURCES (intensity):	' + str(est_n_sources_i))
		print('ERROR:\t' + str(error_i))
		write_data(est_n_sources_array_i, est_coord_array_i, fits_dir)
	if center_type.find('m') != -1:
		print(str(mean_dec_m) + '\t' + str(mean_ra_m) + '\t' + str(var_dec_m) + '\t' + str(var_ra_m) \
+ '\n') 
		est_n_sources_m = int(len(est_coord_array_m)/2)
		print('N. ESTIMATED SOURCES (mean):	' + str(est_n_sources_m))
		print('ERROR:\t' + str(error_m))
		write_data(est_n_sources_array_m, est_coord_array_m, fits_dir)


############################################################
# function devoted to generate the parameters of the main program and the related --help function. 
# All the parameters needed by the main function are described below
############################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to analyze generated maps and \
compare the obtained data to the theoretical expected results')
	parser.add_argument('theoretical map coordinates', metavar = 'logfile', type=str, \
help='theoretical map log (previously generated by recall.py)' )
	parser.add_argument('Generated Fits folder', metavar = 'fits_dir', type=str, help='Directory \
where generated fits file are saved' )
	parser.add_argument('final result log file name', metavar = 'final_logfile', type=str, \
help='algorithm final result log file name' )
	parser.add_argument('type of centers computation', metavar = 'center type', type=str, \
help='type of centers computation. Even multiple input: (e.g.)  bim = baricenter + intensity + mean' )
	parser.add_argument('spatial gaussian sigma', metavar = 'sigma_spa', type=float, \
help='variance of the smoothing gaussian' )
	parser.add_argument('percentage of accepted area', metavar = 'accept_level', type=float, \
help='percentage to reach in circle detection area check' )
	parser.add_argument('mask circle radius', metavar = 'radius', type=int, \
help='template matching radius' )
	parser.add_argument('binarization treshold (255 valued)', metavar = 'binary_treshold', \
type=int, help='treshold used for image binarization' )
	parser.add_argument('intensity treshold ', metavar = 'intensity_treshold', type=float, \
help='background normalization treshold for intensity information handling' )
	parser.add_argument('blob centers minimum distance', metavar = 'baricenter_distance', \
type=int, help='minimum distance for two blob centers to be considered as separated' )
	args = parser.parse_args()


############################################################
# This function simply takes a text file name as input and stores in an array all the 
# words present in that file
############################################################
def words(fileobj):
    for line in fileobj:
        for word in line.split():
            yield word

############################################################
# This function takes as input the logfile created by "recall.py" where all the maps coordinates are 
# stored. These coordinates will be used as benchmark to test the estimated ones
############################################################
def find_n_maps(logfile):

	# flag initialisation
	i = False;

	# open and split the logfile. Then look for the number of maps
	with open(logfile, 'r') as f:
		wordgen = words(f)
		for word in wordgen:
			if (i == True):
				n_maps = int(word)
				break
			if (word == 'MAPS:'):
				i = True

	return n_maps

############################################################
# This function takes in input the logfile generated by recall.py and the total number of maps. 
# It returns an array containing all the real coordinates in DEC and RA
############################################################
def true_coordinates(logfile, n_maps):

	# vars initialisation
	# array with the sources coordinates
	coord_array = [] 
	
	# array with the number of sources per Map
	n_sources_array = numpy.zeros((n_maps), dtype = 'int')	

	# flag initialisation
	j = True

	# coordinates reading and storage
	for i in range(0,n_maps):
		name = 'Map_'+ str(i+1) + '.xml'
		with open(logfile, 'r') as f:
			for line in f:
				split_line = line.split()
				if ( split_line != [] ):
					if ( str(split_line[0]) == name and j==True):
						n_sources_array[i] = int(split_line[2])
						if (n_sources_array[i] == 0):
							j = True
						else:
							j= False
					elif (j == False):
						coord_array.append(float(split_line[1]))
						coord_array.append(float(split_line[0]))
				else :
					j = True

	return n_sources_array, coord_array

############################################################
# This function performs the sources detection and the coordinates estimation by using the 
# analyze_map.bin algorithm (see documentation)
# It returns both the number of sources per map and the estimated coordinates in two separated arrays
############################################################
def convert_data(n_maps, center_type, sigma_spa, accept_level, radius, binary_treshold, \
intensity_treshold, baricenter_distance, fits_dir):

	# var initialisation
	est_n_sources_array = numpy.zeros(n_maps, dtype = 'int')
	est_coord_array = []

	name = order.natural_list(fits_dir)

	# map files analysis: the coordinates are evaluated in pixels and then converted in RA and DEC 
	# through the wcs:pix2world() function
	for i in tqdm(range(0,n_maps), ascii=True, desc="Map analysis:"):

		#first row reading flag
		flag = False 
		
		# C algorithm call
		call(["./final_3.bin", fits_dir + "/"+ str(name[i]), center_type, "logfile_map", "n", "n", \
sigma_spa, accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance])
		
		# result conversion from pixels to celestial coordinates
		with open('logfile_map', 'r') as f:
			for line in f:
				split_line = line.split()
				if ( flag == False ) :
					flag = True
					est_n_sources_array[i] = int(split_line[0]) 			
				else :
					temp_dec = float(split_line[0])
					temp_ra = float(split_line[1])
					w  = wcs.WCS( fits_dir + '/'+name[i])
					world = w.wcs_pix2world(temp_ra,200 - temp_dec, 1)
					est_coord_array.append(world[1])
					est_coord_array.append(world[0])
		f.close()
					
	return est_n_sources_array, est_coord_array

############################################################
# this function calls the analyze_map.bin algorithm (in the convert_data() function) changing the 
# parameters according to those entered by the user. It returns the estimated sources and 
# n_sources arrays, containg respectively the estimated coordinates and the number of source 
# per each map
############################################################
def switcher(n_maps, center_type,sigma_spa, accept_level, radius, binary_treshold, \
intensity_treshold, baricenter_distance, fits_dir) :

	# var initialisation
	est_coord_array_b = []
	est_coord_array_i = []
	est_coord_array_m = []
	est_n_sources_array_b = []
	est_n_sources_array_i = []
	est_n_sources_array_m = []


	# switch case driven by the entered parametrs and calling convert_data() 
	if (center_type == 'b') :
		[est_n_sources_array_b, est_coord_array_b] = convert_data(n_maps, 'b', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif (center_type == 'i') :
		[est_n_sources_array_i, est_coord_array_i] = convert_data(n_maps, 'i', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance , fits_dir )
	elif (center_type == 'm') :
		[est_n_sources_array_m, est_coord_array_m] = convert_data(n_maps, 'm', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif ( (center_type == 'bi') or (center_type == 'ib') ) :
		[est_n_sources_array_b, est_coord_array_b] = convert_data(n_maps, 'b', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir )
		[est_n_sources_array_i, est_coord_array_i] = convert_data(n_maps, 'i', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif ( (center_type == 'bm') or (center_type == 'mb') ) :
		[est_n_sources_array_b, est_coord_array_b] = convert_data(n_maps, 'b', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
		[est_n_sources_array_m, est_coord_array_m] = convert_data(n_maps, 'm', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif ( (center_type == 'mi') or (center_type == 'im') ) :
		[est_n_sources_array_i, est_coord_array_i] = convert_data(n_maps, 'i', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
		[est_n_sources_array_m, est_coord_array_m] = convert_data(n_maps, 'm', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif ( (center_type == 'bi') or (center_type == 'ib') ) :
		[est_n_sources_array_b, est_coord_array_b] = convert_data(n_maps, 'b', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
		[est_n_sources_array_i, est_coord_array_i] = convert_data(n_maps, 'i', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	elif ( (center_type == 'bim') or (center_type == 'bmi') or (center_type == 'ibm') or \
(center_type == 'imb') or (center_type == 'mbi') or (center_type == 'mib') ) :
		[est_n_sources_array_b, est_coord_array_b] = convert_data(n_maps, 'b', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
		[est_n_sources_array_i, est_coord_array_i] = convert_data(n_maps, 'i', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
		[est_n_sources_array_m, est_coord_array_m] = convert_data(n_maps, 'm', sigma_spa, \
accept_level, radius, binary_treshold, intensity_treshold, baricenter_distance, fits_dir)
	else :
		print('invalid center type\n')
		return -1

	return est_n_sources_array_b, est_n_sources_array_i, est_n_sources_array_m, \
est_coord_array_b, est_coord_array_i, est_coord_array_m


############################################################
# This functions takes as input the real and the estimated coordinates array and computes some 
# statistical analysis. More precisely it computes the mean and the variance of the error and the 
# number of wrongly detected sources. The wrong map number is stored in est_error.log
############################################################
def difference(n_maps, n_sources_array, est_n_sources_array, coord_array, est_coord_array) :

	# var initialisation
	difference = []
	init = 0
	est_init = 0
	error = 0	
	l = len(est_coord_array)
	errorlog = "est_error.log"

	# open the errorlog file
	open(errorlog, 'w').close()
	log = open(errorlog, 'a')

############################################################
	#  The estimation error is computed and stored in the difference array. The error on the estimated
	#  number of sources can be of 3 types:
	#	1) the estimated number is less than the real one
	#	2) the estimated number is greater than the real one
	#	3) the estimated number is the same but the estimated coordinates are far from the real ones
	# In both these cases not all the estimated coordinates are wrong. Therefore the estimated 
	# coordinates should be linked to the right source among the real coordinates array. In doing so 
	# the estimated coordinates are not wasted. 
	# This task is performed by the "fail_handle" function (see below). After this error habdling the 
	# results are stored in two arrays, difference and error, containing respectively the coordinates 
	# error and the number of missed sources per map.	############################################################

	for i in range(0,n_maps):
		
		if ( (n_sources_array[i] == 0) and (est_n_sources_array[i] > 0) ):
			error = error + est_n_sources_array[i]
			est_init = est_init + est_n_sources_array[i]
			log.write('Map ' + str(i+1) + '\n')
		elif ( (n_sources_array[i] == 1)  and ( (est_n_sources_array[i] == 0)  ) ):
			error = error + 1
			init = init + n_sources_array[i]
			log.write('Map ' + str(i+1) + '\n')
		elif ( (n_sources_array[i] == 1)  and ( (est_n_sources_array[i] == 1)  ) ):
			ra_diff = coord_array[2*init] - est_coord_array[2*est_init]
			dec_diff = coord_array[2*init+1] - est_coord_array[2*est_init+1]
			dist = math.sqrt(ra_diff**2 + dec_diff**2)
			if ( dist < 0.5 ):
				difference.append(abs(float(ra_diff)))
				difference.append(abs(float(dec_diff)))
			else:
				error = error + 2
				log.write('Map ' + str(i+1) + '\n')
			est_init = est_init + est_n_sources_array[i]
			init = init + n_sources_array[i]		
		elif ( (n_sources_array[i] == 1)  and ( (est_n_sources_array[i] > 1)  ) ):
			ra_diff = coord_array[2*init] - est_coord_array[2*(est_init)]
			dec_diff = coord_array[2*init+1] - est_coord_array[2*(est_init) + 1]
			pos = math.sqrt(ra_diff**2 + dec_diff**2)
			ra_ind = 0
			dec_ind = 0
			for j in range(1, est_n_sources_array[i]):
				ra_diff = coord_array[2*init] - est_coord_array[2*(est_init + j)]
				dec_diff = coord_array[2*init+1] - est_coord_array[2*(est_init + j) + 1]
				temp = math.sqrt(ra_diff**2 + dec_diff**2)
				if (temp < pos):
					pos = temp
					ra_ind = j
					dec_ind = j
			if (pos < 0.5):
				difference.append(abs(float(coord_array[2*init] - est_coord_array[2*(est_init + \
ra_ind)]))) 
				difference.append(abs(float(coord_array[2*init+1] - est_coord_array[2*(est_init + \
dec_ind)+1])))
				error = error + est_n_sources_array[i] - 1
			else:	
				error = error + est_n_sources_array[i] + 1
			log.write('Map ' + str(i+1) + '\n')
			init = init + 1
			est_init = est_init + est_n_sources_array[i]
			

	log.close()
	return difference, error
				

############################################################
# This function computes and returns the main statistics (mean and variance) on the array generated 
# by the "difference" funtion. 
############################################################
def statistics(difference) :
	
	# var initialisation
	mean_dec = 0
	mean_ra = 0
	var_dec = 0
	var_ra = 0

	length = int(len(difference)/2)

	#print('length = ' + str(length))

	# mean computation on both RA and DEC coordinates
	for i in range(0,length) :
		mean_dec = mean_dec + difference[2*i]
		mean_ra = mean_ra + difference[2*i + 1]
	
	if (length != 0):
		mean_dec = mean_dec/length
		mean_ra = mean_ra/length
	else:
		mean_dec = 0
		mean_ra = 0

	# variance computation on both RA and DEC coordinates
	for i in range(0,length) :
		var_dec = var_dec + (difference[2*i] - mean_dec)**2
		var_ra = var_ra + (difference[2*i+1] - mean_ra)**2
	
	if (length!=0):
		var_dec = var_dec/length
		var_ra = var_ra/length
	else:
		mean_dec = 0
		mean_ra = 0

	return mean_dec, mean_ra, var_dec, var_ra

############################################################
# This function handles the different kinds of error that can occur in the estimation if the number 
# of sources in a map. 
# The estimation error is computed and stored in the difference array. The error on the estimated 
# number of sources can be of 2 types:
#	1) the estimated number is less than the real one
#	2) the estimated number is greater than the real one
#	3) the estimated number is the same but the estimated coordinates are far from the real ones
# The function takes as input the real number of sources and the estimated one (in 2 different array) 
# and for each map generates all the possible combinations of estimeted/real indexes (through a 
# newton binomial). For each generated pair the error between the estimated coordinates and the real 
# ones is computed and compared to the error computed for the other pairs. The lowest error 
# therefore identifies the correct matching. The procedure is performed for both the error situations 
# previously described. The array of correct matchings is returned.
############################################################
def fail_handle(true_array, est_array):

	# var initilisation
	len_true = int(len(true_array)/2)
	len_est = int(len(est_array)/2)

	# handle the case in which the number of estimated sources il lower than the real one
	if len_true > len_est:
		
		# var initialisation
		diff_min_ra = 0
		diff_min_dec = 0
		pos = numpy.zeros(len_est)
		
		# computing the first coordinates error value
		for i in range(len_est):
			diff_min_ra = diff_min_ra + (true_array[2*i] - est_array[2*i])**2
			diff_min_dec = diff_min_dec + (true_array[2*i+1] - est_array[2*i+1])**2 
			pos[i] = i
		diff_min = diff_min_ra + diff_min_dec
		
		# creating the permutations and checking the error
		for data in it.permutations(range(len_true), len_est):
			tmp_ra = 0
			tmp_dec = 0
			for i in range(len_est):
				tmp_ra = tmp_ra + (true_array[2*data[i]] - est_array[2*i])**2
				tmp_dec = tmp_dec + (true_array[2*data[i]+1] - est_array[2*i+1])**2 
			tmp = tmp_ra+tmp_dec
			if diff_min > tmp:
				diff_min = tmp
				pos = data

	# handle the case in which the number of estimated sources il greater than the real one
	else:

		# var initialisation
		diff_min_ra = 0
		diff_min_dec = 0
		pos = numpy.zeros(len_true)

		# computing the first coordinates error value
		for i in range(len_true):
			diff_min_ra = diff_min_ra + (true_array[2*i] - est_array[2*i])**2
			diff_min_dec = diff_min_dec + (true_array[2*i+1] - est_array[2*i+1])**2 
			pos[i] = i
		diff_min = diff_min_ra + diff_min_dec
		
		# creating the permutations and checking the error
		for data in it.permutations(range(len_est), len_true):
			tmp_ra = 0
			tmp_dec = 0
			for i in range(len_true):
				tmp_ra = tmp_ra + (true_array[2*i] - est_array[2*data[i]])**2
				tmp_dec = tmp_dec + (true_array[2*i+1] - est_array[2*data[i]+1])** 2
			tmp = tmp_ra+tmp_dec
			if diff_min > tmp:
				diff_min = tmp
				pos = data
		
	return pos

############################################################
# This function orders all the statistical results obtained from previous functions and plots 
# them for both RA and DEC values
############################################################
def plot_result(diff_array):

	# var initialisation
	length = int(len(diff_array)/2)
	est_dec = numpy.zeros(length)
	est_ra = numpy.zeros(length)

	# vars for the gaussian representation 
	range_min = 0
	range_max = 0.5
	step = 0.0001
	n_bin = int(abs(range_max-range_min)/step)

	# splitting the statistical values in 2 arrays, one for RA and the other for DEC
	for i in range(0,length) :
		est_dec[i] = diff_array[2*i]
		est_ra[i] = diff_array[2*i + 1]

	# gaussian histogram creation
	hist_dec, junk_dec = numpy.histogram(est_dec, n_bin,(range_min,range_max))
	hist_ra, junk_ra = numpy.histogram(est_ra, n_bin,(range_min,range_max))
	x_array = numpy.arange(range_min,range_max,step)

	# plot if values are present
	if (diff_array != []) :
		plt.figure(1)
		plt.plot(x_array, hist_dec,'ro-')
		plt.title("Statistics DEC")
		plt.show()
		plt.figure(2)
		plt.plot(x_array, hist_ra, 'ro-')
		plt.title("Statistics RA")
		plt.show()
		plt.waitforbuttonpress()

############################################################
# This function writes all the estimated coordinates values in a logfile named "est_map.log"
############################################################
def write_data(est_n_sources_array, est_coord_array, fits_dir):
	
	# var initialisation
	logfile = "est_map.log"
	open(logfile, 'w').close()
	log = open(logfile, 'a')
	name = os.listdir(fits_dir)
	name = order.natural_list(fits_dir)
	n_maps = len(est_n_sources_array)

	# logfile header writing 
	log.write('ESTIMATED COORDINATES \n\n')
	log.write('TOTAL NUMBER OF ANALYZED MAPS:\t' + str(n_maps) + '\n\n')
	log.write('MAP\t\t\t\t\t\t\t\t\t\t\t' + 'N_SOURCES\t\t\t' + 'RA\t\t\t\t\t\t\t\t' + \
'DEC\n\n')

	# logfile data writing
	pos = 0
	for i in range(n_maps):
		log.write(name[i] + '\t\t\t\t\t\t\t\t\t\t\t' + str(est_n_sources_array[i]) + '\n')
		#log.write('Map_' + str(i) + '\t\t\t\t\t\t\t\t\t\t\t' + str(est_n_sources_array[i]) + '\n')
		for j in range(est_n_sources_array[i]):
			log.write('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t' + str(est_coord_array[2*(pos+j) + 1]) + \
'\t\t\t\t\t\t' + str(est_coord_array[2*(pos+j) ]) + '\n')
		
		pos = pos + est_n_sources_array[i]
		log.write('\n\n')

	log.close()

############################################################
 # if the file is run from terminal exeute function main() by default
############################################################
if __name__ == "__main__":
    main()

