from classes import Extractor
import cv2
import argparse

parser = argparse.ArgumentParser(description='Detect Sources')
parser.add_argument('filepath', type=str, help='Path of the fits file')
parser.add_argument('--config', type=str, help='Path of config file (default = data/default.conf)')

args = parser.parse_args()

ext = Extractor.Extractor(args.filepath, "./", debug_prints=True, prints=True)

if args.config:
	ext.load_config(args.config)

ext.perform_extraction()
cv2.waitKey()
