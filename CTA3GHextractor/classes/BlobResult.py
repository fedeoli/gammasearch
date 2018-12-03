from astropy.wcs import WCS
import numpy as np


class BlobResult:

	def __init__(self, fits_path, index, relative_path):
		"""
		Constructor
		"""
		self.path_map = fits_path
		self.index = index
		self.relative_path = relative_path
		self.bary = None
		self.diam = None
		self.radius = None
		self.radec = None
		self.mask = None
		return

	def set_bary(self, barycenter):
		self.bary = barycenter
		wcs = WCS(self.path_map)
		ra, dec = wcs.wcs_pix2world(self.bary[0], self.bary[1], 0)
		self.radec = [ra, dec]
		return

	def set_diameter(self, size):
		self.diam = size
		self.radius = size/2
		return

	def set_mask(self, shape):
		msk = np.zeros(shape)
		ykp = int(self.bary[0])  # why?
		xkp = int(self.bary[1])
		for i in range(-int(self.radius), int(self.radius)):
			for j in range(-int(self.radius), int(self.radius)):
				msk[xkp + i, ykp + j] = 1
		self.mask = msk.astype(np.uint8)
		return

	def make_xml_blob(self):
		with open(self.relative_path+"data/blob_model.xml", 'r') as model_xml:
			# read model
			parametrized = model_xml.read()

			# replace params
			parametrized = parametrized.replace("SPOT_NAME", "SPOT_{0}".format(self.index))
			parametrized = parametrized.replace("FOUND_RA", str(np.round(self.radec[0], 4)))
			parametrized = parametrized.replace("FOUND_DEC", str(np.round(self.radec[1], 4)))

		return parametrized

	def print_values(self):
		print("Barycenter: ({0}, {1})".format(int(round(self.bary[0])), int(round(self.bary[1]))))
		print("Radius: {0}".format(np.round(self.radius, 4)))
		print("RA,Dec: ({0}, {1})".format(np.round(self.radec[0], 4), np.round(self.radec[1], 4)))
		return
