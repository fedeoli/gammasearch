# ==========================================================================
# Helper class for ctools integration into the Science Alert Generation system
#
# Copyright (C) 2018 Andrea Bulgarelli, Nicolo' Parmiggiani
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================

import os
import shutil
import gammalib
import ctools
import cscripts
import CToolsGammaPipe
from PipeConfiguration import PostAnalysisCopyFilesConfiguration
from PipeConfiguration import CToolsRunConfiguration


class PostAnalysisCopyFiles:
	def __init__(self, sessionconf):
		self.sessionconf = sessionconf
		return
		
	def execute(self):
		print('copy files')
		if self.sessionconf.WebImage == 1:
			#copy sky1.png
			os.system('sudo mkdir -p ' + self.sessionconf.WebImageDir)
			for entry in os.scandir(self.sessionconf.resdir):
				#print(entry.name)
				if entry.name.endswith('_sky1.png'):
					#print(entry.name)		
					shutil.copy(self.sessionconf.resdir + '/' + entry.name, self.sessionconf.WebImageDir + '/sky1.png')
				if entry.name.endswith('_cube.fits'):
					#print(entry.name)		
					shutil.copy(self.sessionconf.resdir + '/' + entry.name, self.sessionconf.WebImageDir + '/cube.fits')
			
		if self.sessionconf.TimeLine == 1:
			os.system('sudo mkdir -p ' + self.sessionconf.TimeLineDir)
			for entry in os.scandir(self.sessionconf.resdir):
				#print(entry.name)
				if  entry.name.endswith('_sky1.png'):
					shutil.copy(self.sessionconf.resdir + '/' + entry.name, self.sessionconf.TimeLineDir)
		
		return
		
		
		
def executePostAnalysis(filename):
	# Setup postanalysis
	if filename:
		sessionconf = PostAnalysisCopyFilesConfiguration(filename)
	
	if sessionconf.postanalysis == "copyfiles":
		copyfiles = PostAnalysisCopyFiles(sessionconf)
		copyfiles.execute()

def pipeline_binned():
	print('Run binned pipeline')
	"""
	Run binned pipeline
	"""
	# Set usage string
	usage = 'ExecuteCTools.py [-observation obsfilename] [-simmodel simmodelfilename] [-anamodel analysismodelfilename] [-confpipe configuration pipe][-seed seed]'

	# Set default options
	options = [{'option': '-observation', 'value': ''}, {'option': '-simmodel', 'value': ''}, {'option': '-anamodel', 'value': ''}, {'option': '-runconf', 'value': ''}, {'option': '-eventfilename', 'value': ''}, {'option': '-seed', 'value': '0'}, {'option': '-postanalysis', 'value': ''},]

	# Get arguments and options from command line arguments
	args, options = cscripts.ioutils.get_args_options(options, usage)

	# Extract script parameters from options
	obsfilename = options[0]['value']
	simfilename = options[1]['value']
	analysisfilename = options[2]['value']
	runconffilename = options[3]['value']
	eventfilename = options[4]['value']
	in_seed = int(options[5]['value'])
	postanalysis = options[6]['value']

	print('obsfilename: ' + obsfilename)
	print('simfilename: ' + simfilename)
	print('analysisfilename: ' + analysisfilename)
	print('runconffilename: ' + runconffilename)
	print('eventfilename: ' + eventfilename)
	print('postanalysis: ' + postanalysis)

	if not postanalysis:
		gp = CToolsGammaPipe.CToolsGammaPipe()
	
		gp.init(obsfilename, simfilename, analysisfilename, runconffilename, eventfilename)

		# Run analysis pipeline
		gp.run_pipeline(seed=in_seed)
		
		
		
	if postanalysis:
		
		executePostAnalysis(postanalysis)
		
	# Return
	return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run binned in-memory pipeline
    pipeline_binned()
