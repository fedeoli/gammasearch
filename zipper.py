########################################## RECALL.PY ################################################
# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	
#################################################################################################
import map_creator
import data_analysis
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
import shutil
import math
from astropy.io import fits
from GammaPipeCommon.SkyImage import SkyImage
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import order

import time
from tqdm import tqdm


def map_foldering():
	coordinates = numpy.zeros(2,dtype=int)	
	
	xml_dir = sys.argv[1]
	#fits_dir = sys.argv[2]
	events_dir = sys.argv[2]
	coordinates[0] = float(sys.argv[3])
	coordinates[1] = float(sys.argv[4])
	
	dest_logfile = 'sigma_recap.log'
	dest = open(dest_logfile, 'w')

	if os.path.isdir('maps_temp'):
		shutil.rmtree('maps_temp')	

	os.system('cp -r ' + xml_dir + ' maps_temp')
	
	if os.path.isdir('zipdir'):
		shutil.rmtree('zipdir')	
	
	os.system('mkdir zipdir')
	os.system('mkdir zipdir/sigma3')
	os.system('mkdir zipdir/sigma4')
	os.system('mkdir zipdir/sigma5')
	os.system('mkdir zipdir/background')
	os.system('mkdir zipdir/total')

	os.system('mv maps_temp/Map.log zipdir/')

	name_array, n_files = order.directory_order_xml('maps_temp')
	n_sources_array, true_array = data_analysis.true_coordinates('zipdir/Map.log', n_files)
	
	coord_pos = 0

	
	for i in tqdm(range(n_files)):
		#os.system( 'ctobssim inmodel=maps_temp/' + name_array[i] + ' ra=' + str(coordinates[0]) + ' dec=' + str(coordinates[1]) +  ' rad=10 tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:15:00 emin=0.1 emax=100.0 caldb=prod2 irf=South_0.5h outevents=events.fits seed=' + str(i) )

		event_name = 'Events_' + str(i+1) + '.fits'
		ctobs_log_name = 'ctobssim_' + str(i+1) + '.log'
		cube_name = 'Map_' + str(i+1) + '.fits'

		os.system('cp ' + events_dir + '/'  + event_name + ' ./events.fits')
		os.system('cp ' + events_dir + '/'  + ctobs_log_name + ' ./ctobssim.log')
				
		if (n_sources_array[i] != 0):
			os.system( 'ctselect inobs=events.fits ra=' + str(true_array[2*coord_pos + 1]) + ' dec=' + str(true_array[2*coord_pos]) +  ' rad=0.2 tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:15:00 emin=0.1 emax=100.0 outobs=selected_events.fits')	
			coord_pos = coord_pos + 1
			S, B = find_counts('ctobssim.log', 'ctselect.log')
			sigma = float(S)/math.sqrt(float(B))
			dest.write('SIGMA ' + name_array[i] + ' :\t' + str(sigma) + '\n')
		else:
			os.system( 'ctselect inobs=events.fits ra=' + str(coordinates[0]) + ' dec=' + str(coordinates[1]) +  ' rad=0.2 tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:15:00 emin=0.1 emax=100.0 outobs=selected_events.fits')


		os.system( 'ctbin inobs=events.fits coordsys=CEL proj=CAR xref=' + str(coordinates[0]) + ' yref=' + str(coordinates[1])  + ' binsz=0.02 nxpix=200 nypix=200 ebinalg=LOG emin=0.1 emax=100 enumbins=1 outcube=cube.fits' )
		#os.system( 'ctskymap inobs=events.fits caldb=prod2 irf=South_0.5h outmap=sky.png emin=0.1 emax=100 nxpix=200 nypix=200 binsz=0.02 enumbins=1 outcube=cube.fits coordsys=CEL proj=CAR xref=' + str(coordinates[0]) + ' yref=' + str(coordinates[1]) )

		
		select_event_name = 'Selected_events_' + str(i+1) + '.fits'
		sky_name = 'Sky_' + str(i+1) + '.png'
		dir_name = 'Map_' + str(i+1)
		ctsel_log_name = 'ctselect_' + str(i+1) + '.log'

		# total folder construction
		os.system('cp cube.fits zipdir/total/' + cube_name)
		
		# partial folder construction
		if (n_sources_array[i] != 0):
			if (sigma < 4):
				os.system('cp cube.fits zipdir/sigma3/' + cube_name)
			elif (sigma < 5):
				os.system('cp cube.fits zipdir/sigma4/' + cube_name)
			else:
				os.system('cp cube.fits zipdir/sigma5/' + cube_name)
		else:
			os.system('cp cube.fits zipdir/background/' + cube_name)

		
		#os.system('mkdir ./zipdir/' + dir_name)
		#os.system('mv events.fits zipdir/' + dir_name + '/' + event_name )
		#os.system('mv selected_events.fits zipdir/' + dir_name + '/' + select_event_name )
		#os.system('mv cube.fits zipdir/' + dir_name + '/' + cube_name )
		#os.system('mv ctobssim.log zipdir/' + dir_name + '/' + ctobs_log_name )
		#os.system('mv ctselect.log zipdir/' + dir_name + '/' + ctsel_log_name )
		#os.system('mv sky1.png zipdir/' + dir_name + ' ' + sky_name )

	print('creating sigma3 log...')
	write_log('zipdir/sigma3', n_files ,'sigmalog', 'zipdir/Map.log' )
	name_correction('zipdir/sigma3', 'sigmalog')
	
	os.system('mv sigmalog zipdir/Map_sigma3.log')

	print('creating sigma4 log...')
	write_log('zipdir/sigma4', n_files ,'sigmalog', 'zipdir/Map.log' )
	name_correction('zipdir/sigma4', 'sigmalog')
	
	os.system('mv sigmalog zipdir/Map_sigma4.log')

	print('creating sigma5 log...')
	write_log('zipdir/sigma5', n_files ,'sigmalog', 'zipdir/Map.log' )
	name_correction('zipdir/sigma5', 'sigmalog')
	
	os.system('mv sigmalog zipdir/Map_sigma5.log')

	print('creating background log...')
	write_log('zipdir/background', n_files ,'sigmalog', 'zipdir/Map.log' )
	name_correction('zipdir/background', 'sigmalog')

	os.system('mv sigmalog zipdir/Background_sigma.log')


	if os.path.isdir('maps_temp'):
		shutil.rmtree('maps_temp')

	dest.close()

	if os.path.isfile('dataset.zip'):
		os.system('rm dataset.zip')
		
	print('zipping files...')
	os.system('zip -r dataset.zip zipdir &> trash.log')
	os.system('rm trash.log')

	return 0



def find_counts(logfile_ctobs, logfile_ctsel):
	# open and split the logfile. Then look for the number of maps
	line_count = 0
	temp_count = 0
	with open(logfile_ctobs, 'r') as f:
		for line in f:
			line_count = line_count + 1
	f.close()
	with open(logfile_ctobs, 'r') as f:
		for line in f:
			temp_count = temp_count+1
			if (temp_count == line_count -11):
				line_array = line.split()
				MC_sources = line_array[5]
	f.close()

	temp_count = 0
	line_count = 0
	with open(logfile_ctsel, 'r') as f:
		for line in f:
			line_count = line_count + 1
	f.close()
	with open(logfile_ctsel, 'r') as f:
		for line in f:
			temp_count = temp_count+1
			if (temp_count == line_count -8):
				line_array = line.split()
				obs_events = line_array[6]
	f.close()
	return MC_sources, obs_events

def write_log_old(fits_dir, true_array, n_sources_array, logfile):

	if os.path.isdir('log_maps_temp'):
		shutil.rmtree('log_maps_temp')	

	os.system('cp -r ' + fits_dir + ' log_maps_temp')

	#name_array, n_maps = order.directory_order('log_maps_temp')
	name_array = os.listdir('log_maps_temp')
	n_maps = len(name_array)

	log = open(logfile, 'w')
	log.write('GENERATED COORDINATES AND INTENSITIES\n\n')
	log.write('TOTAL NUMBER OF MAPS:\t' + str(n_maps) + '\n\n')
	log.write('FILE\t\t' +  'BACKGROUND\t' + 'N_SOURCES\t\t' + 'RA\t\t' + 'DEC\n\n')

	for i in range(0,n_maps):
		tmp = name_array[i]
		mapname = tmp.replace('fits', 'xml')

		tmp = tmp.replace('Map_', '')
		tmp = tmp.replace('.fits', '')

		map_index = int(tmp)
		n_sources = n_sources_array[map_index - 1 ]

		log.write(mapname + '\t\t' + str(n_sources) + '\t\t' + '\n\n')

		for i in range(0,n_sources):
			log.write('\t\t\t\t\t\t\t' + str(true_array[2*(map_index-1+i) + 1]) + '\t' + str(true_array[2*(map_index-1+i)]) + '\n\n')  

	log.close()
	shutil.rmtree('log_maps_temp')	

def write_log(fits_dir, n_maps_total, logfile, maplog):

	'''if os.path.isdir(fits_dir+'/../background_ord'):
		shutil.rmtree(fits_dir+'/../background_ord')	

	os.system('mkdir '+fits_dir+'/../background_ord')'''

	#name_array, n_maps = order.directory_order('log_maps_temp')
	name_array = os.listdir(fits_dir)
	n_maps = len(name_array)

	log = open(logfile, 'w')
	log.write('GENERATED COORDINATES AND INTENSITIES\n\n')
	log.write('TOTAL NUMBER OF MAPS:\t' + str(n_maps) + '\n\n')
	log.write('FILE\t\t' +  'BACKGROUND\t' + 'N_SOURCES\t\t' + 'RA\t\t' + 'DEC\t\t' + 'INTENSITY\n\n')

	found = False
	for i in range(0,n_maps):
		name = 'Map_' + str(i+1) + '.xml'
		curr_name = name_array[i]
		curr_name = curr_name.replace('.fits', '.xml')
		with open(maplog, 'r') as f:
			for line in f:
				split_line = line.split()
				if ( split_line != [] ):
					if ( (str(split_line[0]) == curr_name)  and (found == False) ):
						log.write(line)
						found = True
						n_temp=int(split_line[2])
					elif ((found == True) and (n_temp>0)):	
						log.write('\t' + line + '\n')
						found = False
					elif ((found == True) and (n_temp == 0)):
						found = False

def name_correction(fitsdir, maplog):
	
	nmaps = (len(os.listdir(fitsdir)))

	i = 1
	with open(maplog, 'r') as f:
		filedata = f.read()
	with open(maplog, 'r') as f:
		for line in f:
			split_line = line.split()
			if ( (split_line != []) and (split_line[0][0] == 'M') ):
				curr_name = split_line[0]
				filename = curr_name.replace( 'xml', 'fits')
				filename = fitsdir + '/' + filename
				tempname = fitsdir + '/Temp_' + str(i)
				tempname_log = 'tmp_' + str(i) + '.xml'
				os.rename(filename, tempname)
				filedata = filedata.replace(curr_name, tempname_log)
				i = i+1
	
	with open(maplog, 'w') as f:
		f.write(filedata)
	
	i = 1
	with open(maplog, 'r') as f:
		filedata = f.read()

	with open(maplog, 'r') as f:
		for line in f:
			split_line = line.split()
			if ( (split_line != []) and (split_line[0][0] == 't') ):
				curr_name = split_line[0]
				mapname = 'Map_' + str(i) + '.xml'
				filedata = filedata.replace(curr_name, mapname)
				i = i+1
	
	with open(maplog, 'w') as f:
		f.write(filedata)

	for i in range(nmaps):
		mapname = fitsdir +  '/Map_' + str(i+1) + '.fits'
		tempname = fitsdir + '/Temp_' + str(i+1)
		os.rename(tempname, mapname)
		
				
				
	


if __name__ == "__main__":
    map_foldering()


