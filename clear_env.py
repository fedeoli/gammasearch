import os
import shutil

def clear():

	if os.path.isfile('est_error.log'):
		os.remove('est_error.log')
	if os.path.isfile('est_map.log'):
		os.remove('est_map.log')
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

if __name__ == "__main__":
    clear()     
