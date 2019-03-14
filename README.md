# gammasearch
This repo develops an algorithm that locates gamma ray sources coordinates from source images in .fits format. The main algorithm is written in C while the analysis part is in python

0. REQUIREMENTS: in order to correctly run and use this library you need to install:

	-) gammalib ( git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git )

	-) ctools ( git clone https://cta-gitlab.irap.omp.eu/ctools/ctools.git )

	-) opencv library

	-) CFITSIO library ( https://heasarc.gsfc.nasa.gov/fitsio/ )
	

1. PROVIDED FILES:
	-) data_analysis.py : This function aims to process big amounts of source images in order to perform some statistical analysis on the number of estimated sources and their estimated positions, expressed in RA and DEC.

	-) map_creator.py : function devoted to randomly generate astronomic coordinates and intensity of gamma ray sources. Also the number of sources per map is randomly generated

	-) recall.py : function devoted to randomly generate folder of .fits file. To be used if you're not in possess of analyzable data yet

	-) recall_inaf.py : function devoted to randomly generate folder of .fits file. To be used if you're not in possess of analyzable data yet (with specific settings asked for INAF project)
  
	-) order.py : This library sorts numerically the .fits images in a directory
  
	-) param_tuning.py: This function aims to tune the parameters of the "analyze_map.bin" algorithm: it just consists of a series of  nested loops within which all the possible combinations of parameters are tested by calling "data_analysis.py"
  
	-) general_analysis.py: This library allows to analyze a general folder full of .fits file, writing the results on a logfile "est_map.log"

	-) analyze_map_source.bin: first version of the sources analyze algorithm	

	-) analyze_map_smooth.bin: second version of the sources analyze algorithm

	-) CToolsGammaPipe.py + ExecuteCTools.py + obsutils.py + PipeConfiguration.py + PostAnalysisCopyResults.py : gammalib library files

	-) gamma_algorithm_source.cpp : source code of the first version of the algorithm
	
	-) gamma_algorithm_source.cpp : source code of the fsecond version of the algorithm

	-) functions.h + functions.cpp : C library for the gamma_algorithm

	-) dataset.zip : zip file containing the zipdir directory described below

2. PROVIDED DIRECTORIES:

	-) CTA3GHextractor + CTAGammaPipeCommon + GammaPipeCommon : directories needed by gammalib 

	-) zipdir : testing directory with 10 maps to be analyzed

3. COMPILER OPTIONS :

	-) The C files have been compiled with the following command: 

		> g++ -O0 -g3 -Wall -fmessage-length=0 -g -pthread -I/path/to/include/opencv -I/path/to/fitsio/include $(pkg-config --cflags --libs opencv) -lcfitsio -o test gamma_algorithm.cpp functions.cpp

	WATCHOUT!! if you're building under archlinux there's an issue with the vtk library whose filesystem organization has been modified. Check "https://www.reddit.com/r/archlinux/comments/ams21f/linking_errors_with_opencv/"


FOR (ALMOST) EACH FUNCTION AN INLINE HELP IS AVAILABLE (except for order.py and param_ which is meant to be used as a library).

	> (e.g.) python file.py -h

TESTING COMMANDS:

MAP CREATION AND ZIPPING (not necessary - zipdir already provided)
		> python recall_inaf.py 220 45 4 10
		> python zipper.py Generated_Events Generated_Maps

SINGLE ANALYSIS
		> python -W ignore python_launch.py zipdir/sigma3/Map_1.fits single.log b 2.3 0.8 3 23 4.5 66 y

ANALYSIS ON THE GENERATED VALUES (generated files: final.log + est_error.log + est_map.log)
		> python -W ignore data_analysis.py zipdir/Map.log zipdir/total final.log b 2.3 0.8 3 23 4.5 66

ANALYSIS ON A GENERAL .fits FOLDER (WITHOUT KNOWN VALUES - generated files: general.log)
		> python -W ignore general_analysis.py zipdir/total final.log b 2.3 0.8 3 23 4.5 66

 


