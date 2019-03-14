import os
import shutil

def clear():

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
