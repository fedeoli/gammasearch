# ==========================================================================
# Main class for ctools integration into the Science Alert Generation system
#
# Copyright (C) 2018 Andrea Bulgarelli, Nicolò Parmiggiani
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
from lxml import etree
import obsutils
from GammaPipeCommon.Configuration import ObservationConfiguration
from PipeConfiguration import CToolsRunConfiguration
from CTAGammaPipeCommon.create_fits import write_fits
from GammaPipeCommon.utility import Utility
from GammaPipeCommon.SkyImage import SkyImage
from conf import *

#import CTA3GHextractor_wrapper

class CToolsGammaPipe:

	def __init__(self):
		return

	def init(self, obsfilename, simfilename, analysisfilename, runconffilename, eventfilename):
		self.obsfilename = obsfilename
		self.simfilename = simfilename
		self.analysisfilename = analysisfilename
		self.runconffilename = runconffilename
		self.eventfilename = eventfilename

		# Setup run
		if self.runconffilename:
			self.runconf = CToolsRunConfiguration(self.runconffilename)


		#print all attribute of a class
		#attrs = vars(self.runconf)
		#print(attrs)

		# Setup observations
		if self.obsfilename:
			self.obs = self.open_observation(self.obsfilename)

		if self.runconf.skyframeref == 'fk5':
			print('pointing ra  ' + str(self.obsconf.point_ra) + ' dec ' + str(self.obsconf.point_dec) + ' frame ' + str(self.obsconf.point_frame))
			print('point roi ra ' + str(self.obsconf.roi_ra) + ' dec ' + str(self.obsconf.roi_dec) + ' frame ' + str(self.obsconf.roi_frame))
			print('run   roi ra ' + str(self.runconf.roi_ra) + ' dec ' + str(self.runconf.roi_dec) + ' frame ' + str(self.runconf.roi_frame))
		if self.runconf.skyframeref == 'galactic':
			print('pointing l  ' + str(self.obsconf.point_l) + ' b ' + str(self.obsconf.point_b) + ' frame ' + str(self.obsconf.point_frame))
			print('point roi l ' + str(self.obsconf.roi_l) + ' b ' + str(self.obsconf.roi_b) + ' frame ' + str(self.obsconf.roi_frame))
			print('run   roi l ' + str(self.runconf.roi_l) + ' b ' + str(self.runconf.roi_b) + ' frame ' + str(self.runconf.roi_frame))

		return

	def open_observation(self, obsfilename):

		print('Open observation')

		obs = gammalib.GObservations()

		self.obsconf = ObservationConfiguration(obsfilename, self.runconf.timesys, self.runconf.timeunit, self.runconf.skyframeref, self.runconf.skyframeunitref)
		print(self.obsconf.caldb)

		pntdir = gammalib.GSkyDir()

		in_pnttype = self.obsconf.point_frame
		print(in_pnttype)
		if in_pnttype == 'fk5' :
			pntdir.radec_deg(self.obsconf.point_ra, self.obsconf.point_dec)

		#if in_pnttype == 'equatorial' :
			#pntdir.radec_deg(self.obs_ra, self.obs_dec)

		#if in_pnttype == 'galactic' :
		#	pntdir.radec_deg(self.in_l, self.in_b)


		#pntdir.radec_deg(self.obsconf.obs_point_ra, self.obsconf.obs_point_dec)

		tstart = self.obsconf.tstart - self.runconf.timeref
		if self.runconf.timeref_timesys == 'mjd':
			tstart = tstart * 86400.

		print("TSTART " + str(tstart))

		obs1 = obsutils.set_obs(pntdir, tstart, self.obsconf.duration, 1.0, \
			self.obsconf.emin, self.obsconf.emax, self.obsconf.roi_fov, \
			self.obsconf.irf, self.obsconf.caldb, self.obsconf.id)

		print(obs1)

		obs.append(obs1)

		#print(obs1)
		return obs

	#rules
	#1) if simfilename perform simulation based on observation configuration. Do not select events based on run configuration
	#2) if eventfilename is present, read the file. Select events based on run configuration
	#3) if simfilename and eventfilename are not present, query the DB based on run and observation configuration

	#4) now that you have the event list in memory
	#4.A) make cts map
	#4.B) hypothesis generation (with spotfinders OR with the analysisfilename)
	#4.C) if you have an hypothesis, perform MLE

	def run_pipeline(self, debug=False, seed=0):
		"""
		Test unbinned pipeline with FITS file saving
		"""
		# Set script parameters
		events_name          = 'events.fits'
		cubefile_name 	     = ''
		selected_events_name = 'selected_events.fits'
		result_name          = 'results.xml'


		if self.simfilename and self.eventfilename:
			print('error')
			exit(10)

		############ Simulate events

		if self.simfilename:

			# Setup simulation model
			if self.simfilename:
				self.obs.models(gammalib.GModels(self.simfilename))

			if self.runconf.WorkInMemory == 0:
				print('# Generate simulated event list on disk')
				# Simulate events on disk
				sim = ctools.ctobssim(self.obs)
				sim['outevents'] = events_name
				sim['debug'] = debug
				sim['seed']    = seed
				sim['tmin'] = '2020-01-01T00:00:00'
				sim['tmax'] = '2020-01-01T00:15:00'

				print(sim)
				sim.execute()
				self.eventfilename = events_name

			if self.runconf.WorkInMemory == 1:
				print('# Generate simulated event list in memory')
				# Simulate events on memory
				sim = ctools.ctobssim(self.obs)
				sim['debug'] = debug
				sim['seed']    = seed
				sim.run()


		#if self.runfilename is not present
		#   import into DB if any
		#	exit(0)

		############ Get events from DB

		if not self.simfilename and not self.eventfilename:
			#load from DB
			print('# Load event list from DB')

			#write selected_events_name
			if self.runconf.timeref_timesys == 'tt':
				tstart_tt = self.runconf.tmin
				tstop_tt = self.runconf.tmax
			if self.runconf.timeref_timesys == 'mjd':
				tstart_tt = Utility.convert_mjd_to_tt(self.runconf.tmin)
				tstop_tt = Utility.convert_mjd_to_tt(self.runconf.tmax)

			datarepository_dictionary = get_data_repository()
			datarepositoryid = datarepository_dictionary['datarepositoryid']
			#tstart_tt = self.runconf.tmin
			#tstop_tt = self.runconf.tmax
			observationid = self.obsconf.id
			print("tstart"+str(tstart_tt))
			print("tstop"+str(tstart_tt))
			print(self.runconf.timeref_timesys)

			conf_dictionary = get_path_conf()

			path_base_fits = conf_dictionary['path_base_fits']
			tref_mjd = self.runconf.timeref

			if self.runconf.roi_frame == 'fk5':
				obs_ra = self.runconf.roi_ra
				obs_dec = self.runconf.roi_dec
			else:
				exit(10)

			emin = self.runconf.emin
			emax = self.runconf.emax
			fov = self.obsconf.roi_fov
			instrumentname = self.obsconf.instrument
			events_name = write_fits(tstart_tt,tstop_tt,observationid,datarepositoryid,path_base_fits, tref_mjd, obs_ra, obs_dec, emin, emax, fov, instrumentname)
			print("fitsname"+events_name)
			if self.runconf.WorkInMemory == 2:
				print('# Load event list from disk')
				#Load events from fits file on memory
				sim = ctools.ctobssim(self.obs)
				events = sim.obs()[0].events()
				for event in events:
					print(event)

				#events.load(events_name)
				events.load(selected_events_name)


		############ Select events

		if self.runconf.WorkInMemory == 1:
			print('# Select event list in memory')
			select = ctools.ctselect(sim.obs())

		#using file
		if self.runconf.WorkInMemory == 0:
			print('# Select event list from disk')
			select = ctools.ctselect()
			select['inobs']  = events_name
			select['outobs'] = selected_events_name

		select['ra']     = self.runconf.roi_ra
		select['dec']    = self.runconf.roi_dec
		select['rad']    = self.obsconf.roi_fov
		select['tmin']   = 'MJD ' + str(self.runconf.tmin)
		#select['tmin']   = 'INDEF'
		#select['tmin'] = 0.0
		#select['tmin'] = 'MJD 55197.0007660185'
		#select['tmax']   = self.runconf.tmax
		#select['tmax']   = 'INDEF'
		select['tmax']   = 'MJD ' + str(self.runconf.tmax)
		#select['tmax']   = 55197.0007660185 + self.runconf.tmax / 86400.
		select['emin']   = self.runconf.emin
		select['emax']   = self.runconf.emax

		if self.runconf.WorkInMemory == 1:
			select.run()

		if self.runconf.WorkInMemory == 0:
			select.execute()


		#print(self.runconf.roi_ra)
		#print(self.runconf.roi_dec)
		#print(self.obsconf.obs_fov)
		#print(self.runconf.tmin)
		#print(self.runconf.tmax)
		#print(self.runconf.emin)
		#print(self.runconf.emax)
		#print(select.obs()[0])
		print("select --")
		print(select)
		print("select --")

		if self.runconf.WorkInMemory == 2:
			print('# Load event list from disk')
			#Load events from fits file on memory
			sim = ctools.ctobssim(self.obs)
			events = sim.obs()[0].events()
			for event in events:
				print(event)

			#events.load(events_name)
			events.load(selected_events_name)

		print('Event list generated ----------------------')
		if self.runconf.WorkInMemory == 0:
			print(self.obs)
			print("obs --")
			print(self.obs[0])
			print("obs 0 --")
			localobs = self.obs

		if self.runconf.WorkInMemory == 1:
			print(select.obs())
			print(select.obs()[0])
			localobs = select.obs()

		for run in localobs:

			print('run ---'+selected_events_name)

			print(run)

			# Create container with a single observation
			container = gammalib.GObservations()
			container.append(run)


			if self.runconf.MakeCtsMap == 1:
				cubefile_name 	     = 'cube.fits'
				#event file in memory or read from fits file on memory
				if self.runconf.WorkInMemory == 1:
					bin = ctools.ctbin(container)

				if self.runconf.WorkInMemory == 0:
					#event file on disk
					bin = ctools.ctbin()
					bin['inobs']    = selected_events_name

				#make binned map on disk
				bin['outcube']  = cubefile_name

				#common configs
				bin['ebinalg']  = self.runconf.cts_ebinalg
				bin['emin']     = self.runconf.emin
				bin['emax']     = self.runconf.emax
				bin['enumbins'] = self.runconf.cts_enumbins
				bin['nxpix']    = self.runconf.cts_nxpix
				bin['nypix']    = self.runconf.cts_nypix
				bin['binsz']    = self.runconf.cts_binsz
				bin['coordsys'] = self.runconf.cts_coordsys
				bin['usepnt']   = self.runconf.cts_usepnt # Use pointing for map centre
				bin['proj']     = self.runconf.cts_proj

				#make binned map on disk
				if self.runconf.WorkInMemory == 0:
					bin.execute()
					# Set observation ID if make binned map on disk
					bin.obs()[0].id(cubefile_name)
					bin.obs()[0].eventfile(cubefile_name)
					print('track_1')
					print(os.getcwd())
					print(self.runconf.resdir)
					os.system('mkdir -p ' + os.getcwd() + '/'+ self.runconf.resdir)
					shutil.copy('./cube.fits', os.getcwd() + '/' + self.runconf.resdir + '/'+self.runconf.runprefix + '_cube.fits')
					
					if self.runconf.CtsMapToPng == 1:
						title = 'OBS ' + str(self.obsconf.id) + ' / MJD ' + str(self.runconf.tmin) + ' - ' + 'MJD ' + str(self.runconf.tmax)
						SkyImage.display(cubefile_name, "sky1.png", 3, title)
						#SkyImage.dispalywithds9_cts1(cubefile_name, "sky2", 10)
						#copy results
						print('track_2')
						os.system('mkdir -p ' + os.getcwd() + '/' + self.runconf.resdir)
						shutil.copy('./sky1.png', os.getcwd() + '/' + self.runconf.resdir + '/'+self.runconf.runprefix + '_sky1.png')

				#make binned map on memory
				if self.runconf.WorkInMemory == 1:
					bin.run()

				# Append result to observations
				#localobs.extend(bin.obs())

			#print(obs)
			#print(obs[0])
			print(str(self.obsconf.caldb))

			#hypothesis builders
			#3GHextractor
			# 3GH Extractor code
			if self.runconf.HypothesisGenerator3GH and cubefile_name:
				self.analysisfilename = CTA3GHextractor_wrapper.extract_source(cubefile_name)
				print(self.analysisfilename)
				#cv2.waitKey(0)
				print('HypothesisGeneratorEG3')

			#eseguire MLE
			print('analysisfilename: ' + self.analysisfilename)
			if self.analysisfilename:
				print('MLE')
				# Perform maximum likelihood fitting

				if self.runconf.WorkInMemory == 1:
					like = ctools.ctlike(select.obs())

				if self.runconf.WorkInMemory == 0:
					like = ctools.ctlike()
					like['inobs']    = selected_events_name

				like['inmodel']  = self.analysisfilename
				like['outmodel']  = result_name
				like['caldb'] = str(self.obsconf.caldb)
				like['irf']      = self.obsconf.irf
				like['statistic'] = 'DEFAULT'
				like.execute()
				logL = like.opt().value()
				print(logL)
				# insert logL into results.xml
				tree = etree.parse(result_name)
				contentnav = tree.find(".//source[@type='PointSource']")
				#contentdiv = contentnav.getparent()
				contentnav.set("ts",str(logL))
				#print(etree.tostring(tree,encoding="unicode",pretty_print=True))
				f = open(result_name, 'w')
				f.write(etree.tostring(tree,encoding="unicode",pretty_print=True))
				f.close()
				
				print('track_3')
				os.system('sudo mkdir -p ' + os.getcwd() + '/'+ self.runconf.resdir)
				shutil.copy('./results.xml', os.getcwd() + '/' + self.runconf.resdir + '/'+self.runconf.runprefix + '_results.xml')

			if self.runconf.HypothesisGenerator3GH:
				CTA3GHextractor_wrapper.print_graphs(self.simfilename, result_name, self.analysisfilename)
				#cv2.destroyAllWindows()
				print('HypothesisGeneratorEG3')

		# Return
		return
