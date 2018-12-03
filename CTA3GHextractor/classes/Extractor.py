from filters import kernelize as k, utils, centroid_extraction as ce
from classes import TimeChecker as timer, BlobResult as blob
import cv2
from astropy.io import fits
import numpy as np
from io import StringIO
from bs4 import BeautifulSoup


class Extractor:

	def __init__(self, fits_path=None, relative_path="../", debug_prints=True, prints=True):
		"""
		Constructor
		"""
		self.fits_path = fits_path
		self.relative_path = relative_path
		self.default_config = relative_path+"data/default.conf"
		self.config_loaded = False

		# instanzia il lettore
		self.median_iter = None
		self.median_ksize = None

		self.gaussian_iter = None
		self.gaussian_ksize = None
		self.gaussian_sigma = None

		self.local_mode = None
		self.threshold_mode = None

		self.local_stretch_ksize = None
		self.local_stretch_step_size = None
		self.local_stretch_min_bins = None

		self.local_eq_ksize = None
		self.local_eq_clip_limit = None

		self.adaptive_filtering = None
		self.adaptive_block_size = None
		self.adaptive_const = None

		self.morph_type = None
		self.morph_shape = None
		self.morph_size = None

		self.blob_filter_area = None
		self.blob_min_area = None
		self.blob_max_area = None
		self.blob_filter_circularity = None
		self.blob_min_circularity = None

		self.debug_images = None
		self.debug_prints = debug_prints
		self.prints = prints

		self.normalizedimg = None

		return


	def get_data(self, src, stretch_ksize, stretch_stepsize, stretch_minbin):
		data = fits.getdata(src)
		self.normalizedimg = utils.normalize(data).astype(np.uint8)
		data = k.local_stretching2(data, stretch_ksize, stretch_stepsize, stretch_minbin, False)
		data = utils.normalize(data)
		return data.astype(np.uint8)

	def perform_extraction(self):
		"""
		Identifies the gamma-ray source coordinates (RA,Dec) from the input map and creates a xml file for ctlike
		:return: (ra,dec) coordinates of the source and the path of the output xml
		"""
		#load default parameters
		if not self.config_loaded:
			self.load_config(self.default_config)

		time = timer.TimeChecker()
		# Open fits map
		img = self.get_data(self.fits_path, self.local_stretch_ksize, self.local_stretch_step_size, self.local_stretch_min_bins)
		if self.prints:
			print("loaded map: {0}".format(self.fits_path))
		time.toggle_time("read", self.debug_prints)

		# Filter map
		smoothed = k.median_gaussian(img, self.median_iter, self.median_ksize, self.gaussian_iter, self.gaussian_ksize, self.gaussian_sigma)
		time.toggle_time("smoothing", self.debug_prints)

		# Contrast Enhancing
		#localled = self.local_transformation(self.local_mode)(smoothed)
		time.toggle_time("contrast enhancing", self.debug_prints)

		# Binary segmentation and binary morphology
		segmented = self.segmentation(self.threshold_mode)(smoothed)
		segmented = self.morphology(segmented)
		time.toggle_time("segmentation", self.debug_prints)

		# Blob detection
		keypoints = self.blob_detection(segmented)

		# # Overlap keypoints and original image
		im_with_keypoints = cv2.drawKeypoints(img, keypoints, np.array([]), (0, 255, 0), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)

		# Create output file
		index = 0
		buffer = StringIO()

		if self.prints:
			for keyPoint in keypoints:
				print('----------------------------------')
				current_blob = blob.BlobResult(self.fits_path, index, self.relative_path)
				index = index + 1

				current_blob.set_bary(keyPoint.pt)
				current_blob.set_diameter(keyPoint.size)
				current_blob.set_mask(img.shape)
				current_blob.print_values()

				buffer.write(current_blob.make_xml_blob()+'\n')
			print('----------------------------------')

		time.toggle_time("blob extraction", self.debug_prints)
		if self.debug_prints:
			time.total()

		if self.debug_images:
			utils.show(Results=im_with_keypoints, Original=self.normalizedimg, Smoothed=smoothed, Segmented=segmented)
		else:
			utils.show(Results=im_with_keypoints, Original=self.normalizedimg)

		if self.prints:
			print('Done!')
			print('=================================')

		return utils.create_xml(buffer.getvalue(),self.relative_path)

	def local_stretching(self, img):
		return k.local_stretching(img, self.local_stretch_ksize, self.local_stretch_step_size, self.local_stretch_min_bins, self.debug_images)

	def local_equalization(self, img):
		return k.local_equalization(img, self.local_eq_ksize, self.local_eq_clip_limit, self.debug_images)

	def adaptive_threshold(self, img):
		return cv2.adaptiveThreshold(img, 255, self.adaptive_filtering, cv2.THRESH_BINARY, self.adaptive_block_size, self.adaptive_const)

	def morphology(self, img):
		return cv2.morphologyEx(img, self.morph_type, cv2.getStructuringElement(self.morph_shape, (self.morph_size, self.morph_size)))

	def blob_detection(self, img):
		# Detection parameters
		params = cv2.SimpleBlobDetector_Params()

		# Thresholds
		params.minThreshold = 0
		params.maxThreshold = 256

		# Filter blob by area
		params.filterByArea = self.blob_filter_area
		params.minArea = self.blob_min_area
		params.maxArea = self.blob_max_area

		# Filter blob by circularity
		params.filterByCircularity = self.blob_filter_circularity
		params.minCircularity = self.blob_min_circularity

		# # Remove unnecessary filters
		params.filterByConvexity = False
		params.filterByInertia = False

		detector = cv2.SimpleBlobDetector_create(params)

		# # Detect blobs
		reverse_segmented = 255 - img
		return detector.detect(reverse_segmented)

	def load_config(self, filepath):
		conf_file = open(filepath, "r")
		config = BeautifulSoup(conf_file, "html.parser")
		conf_file.close()
		root = config.gammarayanalysis
		if root.filtering.median:
			self.median_iter = utils.convert_node_value(root.filtering.median.iterations)
			self.median_ksize = utils.convert_node_value(root.filtering.median.kernelsize)
		if root.filtering.gaussian:
			self.gaussian_iter = utils.convert_node_value(root.filtering.gaussian.iterations)
			self.gaussian_ksize = utils.convert_node_value(root.filtering.gaussian.kernelsize)
			self.gaussian_sigma = utils.convert_node_value(root.filtering.gaussian.sigma)

		self.local_mode = utils.convert_node_value(root.localtransformation.type, "string")
		if self.local_mode == "Stretching":
			self.local_stretch_ksize = utils.convert_node_value(root.localtransformation.kernelsize)
			self.local_stretch_step_size = utils.convert_node_value(root.localtransformation.stepsize)
			self.local_stretch_min_bins = utils.convert_node_value(root.localtransformation.minbins)
		elif self.local_mode == "Equalization":
			self.local_eq_ksize = utils.convert_node_value(root.localtransformation.kernelsize)
			self.local_eq_clip_limit = utils.convert_node_value(root.localtransformation.stepsize)

		self.threshold_mode = utils.convert_node_value(root.segmentation.type, "string")
		if self.threshold_mode == "Adaptive":
			self.adaptive_filtering = self.adaptive_filter(utils.convert_node_value(root.segmentation.filter, "string"))
			self.adaptive_block_size = utils.convert_node_value(root.segmentation.blocksize)
			self.adaptive_const = utils.convert_node_value(root.segmentation.constant)

		if root.binarymorphology:
			self.morph_type = self.morph_operator(utils.convert_node_value(root.binarymorphology.type, "string"))
			self.morph_shape = self.set_morph_shape(utils.convert_node_value(root.binarymorphology.shape, "string"))
			self.morph_size = utils.convert_node_value(root.binarymorphology.size)

		if root.blobdetector:
			self.blob_filter_area = utils.convert_node_value(root.blobdetector.filterarea, "bool")
			if root.blobdetector.minarea:
				self.blob_min_area = utils.convert_node_value(root.blobdetector.minarea)
			if root.blobdetector.maxarea:
				self.blob_max_area = utils.convert_node_value(root.blobdetector.maxarea)
			self.blob_filter_circularity = utils.convert_node_value(root.blobdetector.filtercircularity, "bool")
			self.blob_min_circularity = utils.convert_node_value(root.blobdetector.mincircularity, "float")

		if root.debugimages:
			self.debug_images = utils.convert_node_value(root.debugimages, "bool")

		self.config_loaded = True
		return

	def local_transformation(self, x):
		return {
			"Stretching": self.local_stretching,
			"Equalization": self.local_equalization,
		}.get(x)

	def segmentation(self, x):
		return {
			"Adaptive": self.adaptive_threshold,
		}.get(x)

	def adaptive_filter(self, x):
		return {
			"Mean": cv2.ADAPTIVE_THRESH_MEAN_C,
			"Gaussian": cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
		}.get(x)

	def set_morph_shape(self, x):
		return {
			"Ellipse": cv2.MORPH_ELLIPSE,
			"Cross": cv2.MORPH_CROSS,
		}.get(x)

	def morph_operator(self, x):
		return {
			"Opening": cv2.MORPH_OPEN,
			"Closing": cv2.MORPH_CLOSE,
		}.get(x)
