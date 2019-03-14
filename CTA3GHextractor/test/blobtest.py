from classes import Extractor
import glob
import cv2
import collections
import argparse

parser = argparse.ArgumentParser(description='Detect Sources')
parser.add_argument('--filepath', type=str, help='Path of the fits file')
parser.add_argument('--config', type=str, help='Path of config file (default = data/default.conf)')
args = parser.parse_args()

fitspath = args.filepath if args.filepath else "../img/test/*.fits"
config = args.config if args.config else "../data/default.conf"


fits_names = glob.glob(fitspath)
index = len(fits_names) - 1
selected_mode = -1
selected_param = -1
run = True
debug_images = False

ext = Extractor.Extractor(fits_names[0], debug_prints=False, prints=False)
ext.load_config(config)
ext.debug_images = debug_images
ext.perform_extraction()


keys = {
	'right_arrow': 83,
	'left_arrow': 81,
	'up_arrow': 82,
	'down_arrow': 84,
	'esc': 27,
	'enter': 13,
	'w': 119,
	'a': 97,
	's': 115,
	'd': 100,
	'threshold': 116,  		# t
	'filter': 102,  		# f
	'stretch': 108,  		# s
	#  'equalization': 101,  	# e
	'r': 114,
	'v': 118,
	'p': 112,
	'1': 49,
	'2': 50,
	'3': 51,
	'4': 52,
	'5': 53,
	'anti-parallels': 127
}

mode = {
	'init': -1,
	'none': 0,
	'threshold': 1,
	'filter': 2,
	'stretch': 3,
	# 'equalization': 4,
}


# current parameters
params = collections.OrderedDict({
	'threshold': collections.OrderedDict(),
	'filter': collections.OrderedDict(),
	# 'equalization': collections.OrderedDict(),
	'stretch': collections.OrderedDict()
})

t = params['threshold']
f = params['filter']
# e = params['equalization']
s = params['stretch']
#params['local mode'] = ext.local_mode
t['adaptive kernel size'] = ext.adaptive_block_size
t['adaptive constant'] = ext.adaptive_const
f['number of median iterations'] = ext.median_iter
f['median kernel size'] = ext.median_ksize
f['number of gaussian iterations'] = ext.gaussian_iter
f['gaussian kernel size'] = ext.gaussian_ksize
# e['equalization kernel size'] = ext.local_eq_ksize
# e['clip limit'] = ext.local_eq_clip_limit
s['stretch kernel size'] = ext.local_stretch_ksize
s['stretch step size'] = ext.local_stretch_step_size
s['stretch min bins'] = ext.local_stretch_min_bins

# if not e['equalization kernel size']:
# 	e['equalization kernel size'] = 15
# 	e['clip limit'] = 2.0
#
# if not s['stretch kernel size']:
# 	s['stretch kernel size'] = 21
# 	s['stretch step size'] = 5
# 	s['stretch min bins'] = 3


def init_ext():
	ext.adaptive_block_size = t['adaptive kernel size']
	ext.adaptive_const = t['adaptive constant']
	ext.median_iter = f['number of median iterations']
	ext.median_ksize = f['median kernel size']
	ext.gaussian_iter = f['number of gaussian iterations']
	ext.gaussian_ksize = f['gaussian kernel size']
	# ext.local_eq_ksize = e['equalization kernel size']
	# ext.local_eq_clip_limit = e['clip limit']
	ext.local_stretch_ksize = s['stretch kernel size']
	ext.local_stretch_step_size = s['stretch step size']
	ext.local_stretch_min_bins = s['stretch min bins']
	#ext.local_mode = params['local mode']


def print_mode():
	print("---------------------------------"
		"\nSelect mode: \n\n"
		"t:\t\tadaptive threshold\n"
		"f:\t\tgaussian & median filter\n"
		# "e:\t\tlocal equalization\n"
		"l:\t\tlocal stretch\n\n"
		"r:\t\tprint results\n"
		"v:\t\tprint current values\n"
		"p:\t\tshow intermediate steps\n"
		"\n'w'-'s' to change map\n"
		"\nesc:\tquit\n"
		"---------------------------------")
	return


def print_values(params):
	print('\n==============VALUES==============')
	for x in params:
		print('\n--', x.upper(), '--\n')
		if x != 'local mode':
			for y in params[x]:
					print(params[x][y], '\t', y)
		else:
			print(params[x])
	print('==================================\n')


# print("\n=================================\nPlease keep focus on the maps\n
#  canvases when changing the parameters\n==================================")

while True:

	if run:
		init_ext()
		ext.debug_images = debug_images
		ext.perform_extraction()
		ext.prints = False
		ext.debug_prints = False


	if selected_mode == mode['init']:
		print_mode()
		selected_mode = mode['none']

	key = cv2.waitKey(0)
	run = True

	if key == keys['esc']:
		break

	elif key == keys['enter']:
		print_mode()

	elif key == keys['v']:
		print_values(params)

	# MODE
	#elif key in [keys['threshold'], keys['filter'], keys['equalization'], keys['stretch']]:
	elif key in [keys['threshold'], keys['filter'], keys['stretch']]:
		mode_key = [k for k, v in keys.items() if v == key][0]
		print("\nSelected mode:", mode_key.upper(), " - Select the parameter to change:\n")
		run = True
		selected_mode = mode[mode_key]

		# if key == keys['equalization']:
		# 	params['local mode'] = "Equalization"
		# elif key == keys['stretch']:
		# 	params['local mode'] = "Stretching"

		ord_dict = collections.OrderedDict(params[mode_key])
		idx = 0
		for x in ord_dict:
			idx += 1
			print(idx, ': ({0})\t'.format(ord_dict[x]), x)

	elif key == keys['r']:
		print("\n==============RESULTS=============")
		ext.prints = True

	elif key == keys['p']:
		debug_images = not debug_images
		if not debug_images:
			cv2.destroyWindow('Smoothed')
			cv2.destroyWindow('Segmented')

	# FITS maps
	elif key == keys['w'] and index < len(fits_names) - 1:
		index += 1
		ext = Extractor.Extractor(fits_names[index], debug_prints=False, prints=False)
		ext.load_config(config)
	elif key == keys['s'] and index > 0:
		index -= 1
		ext = Extractor.Extractor(fits_names[index], debug_prints=False, prints=False)
		ext.load_config(config)

	# PARAM
	elif key in [keys['1'], keys['2'], keys['3'], keys['4'], keys['5']]:
		if selected_mode != mode['none']:
			selected_param = int([k for k, v in keys.items() if v == key][0])
			print("\nSelected param:", selected_param, " - Use 'a'-'d' to change it.")
		else:
			print("\nNo mode selected!")
	elif key in [keys['a'], keys['d']]:
		if selected_param and selected_mode != ['none']:
			if key == keys['d']:
				sign = 1
			elif key == keys['a']:
				sign = -1

			if selected_mode == mode['threshold']:
				if selected_param == 1:
					t['adaptive kernel size'] += 2*sign
					if t['adaptive kernel size'] < 3:
						t['adaptive kernel size'] = 3

				elif selected_param == 2:
					t['adaptive constant'] += sign
				else:
					print("\nNo parameter selected!")

			elif selected_mode == mode['filter']:
				if selected_param == 1:
					f['number of median iterations'] += sign
					if f['number of median iterations'] < 0:
						f['number of median iterations'] = 0
				elif selected_param == 2:
					f['median kernel size'] += 2*sign
					if f['median kernel size'] < 3:
						f['median kernel size'] = 3
				elif selected_param == 3:
					f['number of gaussian iterations'] += sign
					if f['number of gaussian iterations'] < 0:
						f['number of gaussian iterations'] = 0
				elif selected_param == 4:
					f['gaussian kernel size'] += 2*sign
					if f['gaussian kernel size'] < 3:
						f['gaussian kernel size'] = 3
				else:
					print("\nNo parameter selected!")

			# elif selected_mode == mode['equalization']:
			# 	if selected_param == 1:
			# 		e['equalization kernel size'] += 2*sign
			# 		if e['equalization kernel size'] < 3:
			# 			e['equalization kernel size'] = 3
			# 	elif selected_param == 2:
			# 		e['clip limit'] += sign
			# 	else:
			# 		print("\nNo parameter selected!")

			elif selected_mode == mode['stretch']:
				if selected_param == 1:
					s['stretch kernel size'] += 2*sign
					if s['stretch kernel size'] < 2:
						s['stretch kernel size'] = 2
				elif selected_param == 2:
					s['stretch step size'] += sign
					if s['stretch step size'] < 1:
						s['stretch step size'] = 1
				elif selected_param == 3:
					s['stretch min bins'] += sign
					if s['stretch min bins'] < 1:
						s['stretch min bins'] = 1
				else:
					print("\nNo parameter selected!")

			else:
				print("\nNo mode selected!")
				run = False
		else:
			print("\nNo parameter selected!")
			run = False

	else:
		run = False
		# print(key)

cv2.destroyAllWindows()


