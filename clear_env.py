############################################################
# INTESTAZIONE DA DECIDERE

# This library has been developed to clear the working directory from all junk/previous analysis files
############################################################


# IMPORTS
import os
import shutil
import argparse

############################################################
# This library aims to clear the working directory from all junk/previous analysis files.
# The library is organized as follows:
#		1)	function clear() :	This function checks if some files are present and if this is the
#							case it deletes them.
#		2)	function usage() :	this function prints an inline help by calling "python clear_env.py -h"
############################################################

############################################################
# This function simply checks if some files are present and if this is the case it deletes them. These 
# files are created by other libraries during the source generation and analysis and they can be deleted
# because their name is hard coded. 
############################################################
def clear():

	# call the usage function in order to print the inline help if needed
	usage()

	# delete log files 
	if os.path.isfile('est_error.log'):
		os.remove('est_error.log')
	if os.path.isfile('est_map.log'):
		os.remove('est_map.log')
	if os.path.isfile('ctbin.log'):
		os.remove('ctbin.log')
	if os.path.isfile('ctobssim.log'):
		os.remove('ctobssim.log')
	if os.path.isfile('ctselect.log'):
		os.remove('ctselect.log')
	if os.path.isfile('sigma_recap.log'):
		os.remove('sigma_recap.log')
	if os.path.isfile('single.log'):
		os.remove('single.log')
	if os.path.isfile('final.log'):
		os.remove('final.log')
	if os.path.isfile('general.log'):
		os.remove('general.log')
	if os.path.isfile('cube.fits'):
		os.remove('cube.fits')
	if os.path.isfile('events.fits'):
		os.remove('events.fits')
	if os.path.isfile('selected_events.fits'):
		os.remove('selected_events.fits')
	if os.path.isfile('logfile_map'):
		os.remove('logfile_map')
	
	# dlete directories 
	if os.path.isdir('__pycache__'):
		shutil.rmtree('__pycache__')
	if os.path.isdir('ANALYSIS3'):
		shutil.rmtree('ANALYSIS3')
	if os.path.isdir('Generated_Events'):
		shutil.rmtree('Generated_Events')
	if os.path.isdir('Generated_Fits'):
		shutil.rmtree('Generated_Fits')
	if os.path.isdir('Generated_Maps'):
		shutil.rmtree('Generated_Maps')
	if os.path.isdir('Generated_Png'):	
		shutil.rmtree('Generated_Png')


############################################################
# function devoted to generate the parameters of the main program and the related --help function. 
# All the parameters needed by the main function are described below
############################################################
def usage():
	parser = argparse.ArgumentParser(description='program devoted to clear the working directory \
from all junk/previous analysis files.\nJust run in a console "python clear_env.py"')
	args = parser.parse_args()


############################################################
 # if the file is run from terminal exeute function clear() by default
############################################################
if __name__ == "__main__":
    clear()     
