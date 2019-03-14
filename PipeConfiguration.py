# ==========================================================================
# Utility functions to load configuration objects
#
# Copyright (C) 2018 Andrea Bulgarelli

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

from GammaPipeCommon.read_xml import read_xml
from GammaPipeCommon.utility import *
from GammaPipeCommon.Configuration import RunConfiguration
from GammaPipeCommon.Configuration import PostAnalysisConfiguration

class CToolsRunConfiguration(RunConfiguration):
	def __init__(self, filename):
		RunConfiguration.__init__(self,filename)
		
		#pipe configuration
		self.WorkInMemory = int(self.info_dict['run']['WorkInMemory']['value'])

		self.MakeCtsMap = int(self.info_dict['run']['MakeCtsMap']['value'])
		self.CtsMapToPng = int(self.info_dict['run']['CtsMapToPng']['value'])

		self.HypothesisGenerator3GH = int(self.info_dict['run']['HypothesisGenerator3GH']['value'])

		#counts map parameters
		self.cts_ebinalg = self.info_dict['run']['CountsMap']['ebinalg']
		self.cts_proj = self.info_dict['run']['CountsMap']['proj']
		self.cts_coordsys = self.info_dict['run']['CountsMap']['coordsys']
		self.cts_enumbins = int(self.info_dict['run']['CountsMap']['enumbins'])
		self.cts_nxpix = int(self.info_dict['run']['CountsMap']['nxpix'])
		self.cts_nypix = int(self.info_dict['run']['CountsMap']['nypix'])
		self.cts_binsz = float(self.info_dict['run']['CountsMap']['binsz'])
		self.cts_usepnt = bool(self.info_dict['run']['CountsMap']['usepnt'])
		
class PostAnalysisCopyFilesConfiguration(PostAnalysisConfiguration):
	def __init__(self, sessionconffilename):
		PostAnalysisConfiguration.__init__(self,sessionconffilename)
		self.resdir = self.info_dict['session']['DirectoryList']['results']
		self.WebImage = int(self.info_dict['session']['WebImage']['value'])
		self.WebImageDir = self.info_dict['session']['WebImage']['dir']
		self.TimeLine = int(self.info_dict['session']['TimeLine']['value'])
		self.TimeLineDir = self.info_dict['session']['TimeLine']['dir']