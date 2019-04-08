############################################################
# INTESTAZIONE DA DECIDERE

# Created by:	Elena Ruggiano, Federico Oliva
# Date: 	23-03-2018	

# This library has been developed to order the directory contents and rename them from 1 to n
############################################################


# IMPORTS
import os 
import re

############################################################
# This library aims to order the contents of a directory and to rename them from 1 to n. This has been 
# done only for the .fits files. 
# The library is organized as follows:
#		1)	function directory_order() :	this function sorts numerically the .fits images in a folder
#		2)	function natural_keys():		This function describes the natural sorting of a #										list				
#		3)	function atof():				This function converts a text into a float
#		4)	function natural_list():		This function uses the natural sorting on a directory 
#		5)	function directory_order_xml():This function sorts numerically the .xml files in a directory					
############################################################

############################################################
# This function sorts numerically the .fits images in a directory
############################################################
def directory_order(namedir):

	name_array = os.listdir(namedir)
	n_files = len(name_array)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]
		mapname = namedir + '/' + "Tmp_Name_" + str(i+1) + '.fits'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]	
		mapname = namedir + '/' + "Map_" + str(i+1) + '.fits'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	n_files = len(name_array)

	return name_array, n_files

############################################################
# This function describes the natural sorting of a list (numerically). To be used as parameter for 
# the list.sort() method
############################################################
def natural_keys(text):
	return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

############################################################
# This function converts a text into a float
############################################################
def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

############################################################
# This function uses the natural sorting on a directory 
############################################################
def natural_list(namedir):
	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	return name_array
	

############################################################
# This function sorts numerically the .xml files in a directory
############################################################
def directory_order_xml(namedir):

	name_array = os.listdir(namedir)
	n_files = len(name_array)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]
		mapname = namedir + '/' + "Tmp_Name_" + str(i+1) + '.xml'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]	
		mapname = namedir + '/' + "Map_" + str(i+1) + '.xml'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	n_files = len(name_array)

	return name_array, n_files
	
