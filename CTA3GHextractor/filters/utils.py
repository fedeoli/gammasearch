from astropy.io import fits
import numpy as np
import cv2
from filters import kernelize as k
import tkinter as tk


#  screen info
root = tk.Tk()
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
#print("Screen resolution: {0}x{1}".format(screen_width, screen_height))

dx = int(screen_width / 4)
dy = int(screen_height / 2)
label_bar_height = 65


def normalize(src):
	src_max = np.max(src)
	src_min = np.min(src)
	src = np.multiply(src-src_min, 255/(src_max - src_min))
	return src


def show(wait=False, **kwargs):
	x = 0
	y = 0
	for key in kwargs:
		label = key
		img = kwargs[key]
		cv2.namedWindow(label, cv2.WINDOW_NORMAL)
		cv2.imshow(label, img)
		cv2.resizeWindow(label, dx, dy - label_bar_height)
		cv2.moveWindow(label, x, y)

		screen_end = int(x / (screen_width - dx)) > 0

		x = x + dx if not screen_end else 0
		y += dy if screen_end else 0

		if wait:
			cv2.waitKey(0)
			cv2.destroyAllWindows()


def create_xml(sources, relative_path):
	output_file_path = "detected.xml"
	with open(relative_path+"data/output_model.xml", 'r') as model_xml:
		with open(output_file_path, 'w') as output_file:

			# read model
			parametrized = model_xml.read()

			# replace params
			parametrized = parametrized.replace("FOUND_SOURCES", str(sources))

			# write replaced
			output_file.write(parametrized)

	return output_file_path


def plot3d(src):
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import pylab as pl
	from PIL import Image
	import numpy as np

	img = Image.fromarray(src)

	z = np.asarray(img)
	mydata = z[::1, ::1]
	fig = pl.figure(facecolor='w')
	ax1 = fig.add_subplot(1, 2, 1)
	im = ax1.imshow(mydata, interpolation='nearest', cmap=pl.cm.jet)
	ax1.set_title('2D')

	ax2 = fig.add_subplot(1, 2, 2, projection='3d')
	x, y = np.mgrid[:mydata.shape[0], :mydata.shape[1]]
	ax2.plot_surface(x, y, mydata, cmap=pl.cm.jet, rstride=1, cstride=1, linewidth=0., antialiased=False)
	ax2.set_title('3D')
	ax2.set_zlim3d(0, 255)
	pl.show()


def sliding_window(image, stepSize, windowSize):
	# slide a window across the image
	for y in range(0, image.shape[0] - int(windowSize[0]) + 1, stepSize):
		for x in range(0, image.shape[1] - int(windowSize[1]) + 1, stepSize):
			# yield the current window
			yield (x, y, image[y:y + windowSize[1], x:x + windowSize[0]])


def perror(error):
	print("Error value: {0}".format(error))


def convert_node_value(value, type='int'):
	return {
		'int': lambda x: int(x),
		'float': lambda x: float(x),
		'string': lambda x: x,
		'bool': str2bool,
	}[type](value.string)


def str2bool(v):
	return v.lower() in ("yes", "true", "t", "1")
