import os 
import re

##############################################################################
# This function sorts numerically the .fits images in a directory
##############################################################################
def directory_order(namedir):

	name_array = os.listdir(namedir)
	#name_array.sort(key=natural_keys)
	n_files = len(name_array)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]
		mapname = namedir + '/' + "Tmp_Name_" + str(i+1) + '.fits'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	#name_array.sort(key=natural_keys)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]	
		mapname = namedir + '/' + "Map_" + str(i+1) + '.fits'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	n_files = len(name_array)

	return name_array, n_files

##############################################################################
# This function describes the natural sorting of a list (numerically). To be used as parameter for the list.sort() method
##############################################################################
def natural_keys(text):
	return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

##############################################################################
# This function converts a text into a float
##############################################################################
def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_list(namedir):
	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	return name_array
	

##############################################################################
# This function sorts numerically the .fits images in a directory
##############################################################################
def directory_order_xml(namedir):

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	n_files = len(name_array)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]
		mapname = namedir + '/' + "Tmp_Name_" + str(i+1) + '.xml'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)

	for i in range(n_files):
		sourcename = namedir + '/' + name_array[i]	
		mapname = namedir + '/' + "Map_" + str(i+1) + '.xml'
		os.rename(sourcename, mapname)

	name_array = os.listdir(namedir)
	name_array.sort(key=natural_keys)
	n_files = len(name_array)

	return name_array, n_files
	
