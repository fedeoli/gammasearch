/*
 * utilizzo: ./binary "mapname.fits" + "b/i/m" + "logfilename" + "y/n"
 */

#define MAIN_FILE
#include "functions.h"

/***************************************************************************************************************************
 * This algorithm has been developed in order to estimate gamma ray sources position from an initial provided image.
 *
 * Algorithm goal: starting from a provided image "source.fits" it proceeds estimating the sources position. The initial images are in .fits
 * format (standard format used in astronomy data analysis) and will be handled with the CFITSIO C library developed by NASA
 * (available at https://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html).
 *
 * Used libraries:
 * 		i)   CFITSIO library
 * 		ii)  OpenCV C functions (version 2.47)
 * 		iii) functions.h library (self written - see attached file)
 *
 * 	Program usage (LINUX environment):
 *
 * 	user$user>> ./analyze_map.bin image.fits center_type logfile.log visual save sigma_spa accept_level radius bin_thresh int_thresh bar_dist
 *
 * 	parameters constraints:
 * 	image.fits = the image to be processed. Must be in .fits format
 * 	center_type = char var. Accepted values are "b", "i", "m".
 * 	logfile.log = the file where the retrieved data will be written. String var.
 * 	visual = char var. Accepted values are "y" or "n".
 * 	save = char var. Accepted values are "y" or "n".
 * 	sigma_spa = float var.
 * 	accept_level = float var.
 * 	radius = int var.
 * 	bin_thresh = int var.
 * 	int_thresh = int var.
 * 	bar_dist = int var.
 *
 ***************************************************************************************************************************/

int main(int argc, char *argv[]){

/***************************************************************************************************************************
 * VARIABLE DEFINITION:
 *
 * IplImage vars: OpenCV structures containing the processed image data. One var of this kind has been used for each algorithm step, allowing the
 * visualisation of the whole process ( if visual == y ).
 *
 * fitsfile vars: these vars are devoted to the .fits file management. They are usually arguments of the CFITSIO library functions.
 *
 * CvSize vars: structures used to initialise image dimensions and general info.
 ***************************************************************************************************************************/

//	IplImage vars
	IplImage* source_image, * elab_image[10], * elab_image_color[5];

//	fitsfile vars
	fitsfile * fits_pointer, *prova, *out, *open, *fits_dest;

//	CvSize vars
	CvSize gamma_image, double_gamma_image;

//	other vars
	int status =0, bitpix =0, dim =0, hdu_type =0, hdu_num=0, hdu_nums=0, *L, n_sources, *baricenter_int,
		 *area, *values, radius, * max_int, binarization_treshold, center_distance;
	long size[2] = {1,1}, fpixel[2] ={1,1};
	float sigma_spa, accept_level, intensity_treshold;
	unsigned int histogram_2[256];
	double *baricenter_r, *ratio, * max_real, *mean ;
	char * image_filepath, * center_type, * logfile, * visual, *save;

/****************************************************************************************************************************
 * parameter defintion from user input - for further info see below (argv[]). As input parameters are string, some conversions are needed for
 * the following processing.
 ****************************************************************************************************************************/
	image_filepath = argv[1];
	center_type = argv[2];
	logfile = argv[3];
	visual = argv[4];
	save = argv[5];
	sigma_spa = atof(argv[6]);
	accept_level = atof(argv[7]);
	radius = atoll(argv[8]);
	binarization_treshold = atoll(argv[9]);
	intensity_treshold = atof(argv[10]);
	center_distance = atoll(argv[10]);

// check if @center_type user input option is valid or not. Accepted values are "b", "i", "m". The parameters meaning will be clarified later on.

	if ((*center_type != 'b') && (*center_type != 'i') && (*center_type != 'm') ){
		printf("incorrect @center_type option\n");
		return -1;
	}

	else{

//	fits file opening and parameter reading procedure. See functions.cpp for further info.
	open = fits_image_open(open,image_filepath, status, hdu_type, hdu_num, hdu_nums);
	fits_pointer = fits_image_read_param(open, status, bitpix, dim, size);

//	CvSize structure initialization.
	gamma_image.height = size[0];
	gamma_image.width = size[1];
	double_gamma_image.width = 2*gamma_image.width;
	double_gamma_image.height = 2*gamma_image.height;

//	creation of a IplImage with the same size of the fits image, retrieved in the previous read_param function.
	source_image = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);

//	fits image data transposition in the created Iplimage structure.
	fits_pointer = fits_image_copy(fits_pointer, status, fpixel, size, source_image->imageData);

		
//	Image display if visual == "y"
	if(*visual == 'y')
		cvShowImage("Gamma Image", source_image);

//	Image Gaussian Smoothing: applied in order to reduce the noise caused by the background
	elab_image[0] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);
	gauss_smooth_fast(source_image, elab_image[0], sigma_spa);

//	Image display if visual == "y"
	if(*visual == 'y')
		cvShowImage("Smoothed Image", elab_image[0]);

// Image save if save == 'y'
		if(*save == 'y')
			out = fits_image_create(fits_pointer, fits_dest, "image_smooth.fits", elab_image[0]->imageData, bitpix, dim, size, status, gamma_image);

//	 Template matching image creation
		elab_image[1] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);

/*******************************************************************************************************************************
* Template matching: the following functions try to match a circle of radius = @radius within the source image. Firstly the grey-level histogram
*	  is created from the previously smoothed image and then the template matching procedure is run.
*******************************************************************************************************************************/
		buildHistogram(source_image, histogram_2, "smoothed image");
		values = (int*)malloc(gamma_image.height*gamma_image.width*sizeof(int));
		circle_detection(source_image, elab_image[1], values, binarization_treshold, accept_level, radius);
//		circle_detection(elab_image[0], elab_image[1], values, binarization_treshold, accept_level, radius);

//	Image display if visual == "y"
		if(*visual == 'y')
			cvShowImage("template matching", elab_image[1]);

//	Image save if save == 'y'
		if(*save == 'y')
			out = fits_image_create(fits_pointer, fits_dest, "image_cirlce.fits", elab_image[1]->imageData, bitpix, dim, size, status, gamma_image);

/******************************************************************************************************************************
* Labeling: this procedure labels in a different way all the blobs previously detected through the template matching. Final values are stored in L.
* For the labeling an 8-connectivity is used.
******************************************************************************************************************************/
		L = (int*)malloc((gamma_image.height * gamma_image.width)*sizeof(int));
		n_sources = labeling_8(elab_image[1], L);

/******************************************************************************************************************************
 * Area computation: in this procedure all the previously detected blobs areas are computed. Final values are stored in area. area[0] is the
 * background area.
 ******************************************************************************************************************************/
		area = (int*)malloc((n_sources+1)*sizeof(int));
		area_blobs(L, n_sources, area, gamma_image);

/******************************************************************************************************************************
 * Intensity ratio evaluation on smoothing: after the first blob detection some blobs could still be wrong background agglomerates which are
 * not as intense as real sources but still intense enough to be over-threshold. Therefore in the following procedure if a blob is not intense at
 * least a fixed times (@intensity_ratio) more than the background, it's killed.
 ******************************************************************************************************************************/
		ratio = (double*)malloc((n_sources+1) * sizeof(double));
		elab_image[2] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);
		intensity_ratio(elab_image[0],  elab_image[2] , L, area, intensity_treshold, n_sources, ratio);

//	Image display if visual == "y"
		if(*visual == 'y')
			cvShowImage("found sources", elab_image[2]);

//	Image save if save == 'y'
		if(*save == 'y')
			out = fits_image_create(fits_pointer, fits_dest, "image_ratio.fits", elab_image[2]->imageData, bitpix, dim, size, status, gamma_image);

//	Second labeling: after the wrong blob killing the labeling must be run again.
		n_sources = labeling_8(elab_image[2], L);

/*******************************************************************************************************************************
 * Second area computation: after the second labeling all the area must be computed again. Before proceeding the previous area vector must be
 * flushed from memory.
 *******************************************************************************************************************************/
		free(area);
		area = (int*)malloc((n_sources+1)*sizeof(int));
		area_blobs(L, n_sources, area, gamma_image);



/*******************************************************************************************************************************
 * Center detection: the following procedure computes the baricenter of the detected sources (blobs). Three different options are available and
 * defined by the user input parameter @center_type:
 * 		i)   baricenter: the center is assumed to be the baricenter of the blob.
 * 		ii)  intensity center: the center is assumed to be the pixel with the highest intensity within the blob.
 * 		ii)  mean: the center is assumed to be the mean of both the previous two methods.
 * 	The visualization of the centers is realized in colored images: red for @center_type = "b", blue for @center_type = "i", green for the
 * 	overlapping pixels of both the previous two (@center_type = "m").
 *
 * 	Lastly the separated blobs removing is executed: the labeling procedure could fail in detecting near pixels separated by a low intensity area
 * 	as part of the same blob. Therefore if two different blobs baricenters are distant less than a fixed threshold (@baricenter_dist) then the
 * 	smallest one is killed.
 *******************************************************************************************************************************/

		if (*center_type == 'b'){
//	Images creation
			elab_image[3] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);
			elab_image_color[0] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 3);
			elab_image_color[1] = cvCreateImage(double_gamma_image, IPL_DEPTH_8U, 3);


// Baricenter vector allocation (both for real and int values)
			baricenter_r = (double*) malloc(2 * n_sources * sizeof(double));
			baricenter_int = (int*) malloc(2 * n_sources * sizeof(int));

// Baricenter computation
			baricenter_real(L, n_sources, area, gamma_image,baricenter_r, baricenter_int);

//	Separated blobs removing procedure
			n_sources = remove_separated_blobs(elab_image[2], elab_image[3], L, baricenter_int, baricenter_r, n_sources, ratio, center_distance);

// Baricenters visualization
			visualize_centers(elab_image[3], elab_image_color[0], baricenter_int, n_sources, 'r');
			double_size_colour( elab_image_color[0],  elab_image_color[1]);

// The baricenter coordinates are sorted depending on the row pixel coord first and then on the col pixel coord
			order_coordinates(baricenter_r, n_sources);

//	The final baricenter coordinates are written in the logfile
			write_center_coordinates(baricenter_r, n_sources, logfile);

// Image display if visual == "y"
			if(*visual == 'y'){
				cvShowImage("removed blobs", elab_image[3]);
				cvShowImage("Colored both", elab_image_color[1]);
			}

// Image save if save == 'y'
			if(*save == 'y')
				out = fits_image_create(fits_pointer, fits_dest, "image_removed_blobs.fits", elab_image[3]->imageData, bitpix, dim, size, status, gamma_image);

// arrays are flushed from memory
			free(baricenter_r);
			free(baricenter_int);

// image release if visual == 'n'
			if(*visual== 'n'){
				cvReleaseImage(&elab_image[3]);
				cvReleaseImage(&elab_image_color[0]);
				cvReleaseImage(&elab_image_color[1]);
				}
		}

		else if (*center_type == 'i'){

//	Images creation
			elab_image[3] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);
			elab_image_color[0] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 3);
			elab_image_color[1] = cvCreateImage(double_gamma_image, IPL_DEPTH_8U, 3);

// Baricenter vector allocation (both for real and int values)
			max_real = (double*)malloc(2*n_sources*sizeof(double));
			max_int = (int*)malloc(2*n_sources*sizeof(int));

// Baricenter computation
			center_values(values, gamma_image, L, max_real, max_int, n_sources);

//	Separated blobs removing procedure
			n_sources = remove_separated_blobs(elab_image[2], elab_image[3], L, max_int, max_real, n_sources, ratio, center_distance);

// Baricenters visualization
			visualize_centers(elab_image[3], elab_image_color[0], max_int, n_sources, 'b');
			double_size_colour( elab_image_color[0],  elab_image_color[1]);

// The baricenter coordinates are sorted depending on the row pixel coord first and then on the col pixel coord
			order_coordinates(max_real, n_sources);

//	The final baricenter coordinates are written in the logfile
			write_center_coordinates(max_real, n_sources, logfile);

// Image display if visual == "y"
			if(*visual == 'y'){
				cvShowImage("removed blobs", elab_image[3]);
				cvShowImage("Colored both", elab_image_color[1]);
			}

// Image save if save == 'y'
			if(*save == 'y')
				out = fits_image_create(fits_pointer, fits_dest, "image_removed_blobs.fits", elab_image[3]->imageData, bitpix, dim, size, status, gamma_image);

// arrays are flushed from memory
			free(max_real);
			free(max_int);

// image release if visual == 'n'
			if(*visual == 'n'){
				cvReleaseImage(&elab_image[3]);
				cvReleaseImage(&elab_image_color[0]);
				cvReleaseImage(&elab_image_color[1]);
			}
		}

		else if (*center_type == 'm'){

//	Images creation
			elab_image[3] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 1);
			elab_image_color[0] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 3);
			elab_image_color[1] = cvCreateImage(gamma_image, IPL_DEPTH_8U, 3);
			elab_image_color[2] = cvCreateImage(double_gamma_image, IPL_DEPTH_8U, 3);


//	Baricenter vector allocation (both for real and int values) - baricenter mode
			baricenter_r = (double*) malloc(2 * n_sources * sizeof(double));
			baricenter_int = (int*) malloc(2 * n_sources * sizeof(int));

//	Baricenter computation
			baricenter_real(L, n_sources, area, gamma_image,baricenter_r, baricenter_int);

//	Separated blobs removing procedure
			n_sources = remove_separated_blobs(elab_image[2], elab_image[3], L, baricenter_int, baricenter_r, n_sources, ratio, center_distance);


//	The baricenter coordinates are sorted depending on the row pixel coord first and then on the col pixel coord
			order_coordinates(baricenter_r, n_sources);

// Baricenter vector allocation (both for real and int values) - Intensity mode
			max_real = (double*)malloc(2*n_sources*sizeof(double));
			max_int = (int*)malloc(2*n_sources*sizeof(int));

//	Baricenter computation
			center_values(values, gamma_image, L, max_real, max_int, n_sources);

//	Baricenters visualization
			visualize_both_centers(elab_image[3], elab_image_color[1], baricenter_int , max_int, n_sources);
			double_size_colour( elab_image_color[1],  elab_image_color[2]);

//	The baricenter coordinates are sorted depending on the row pixel coord first and then on the col pixel coord
			order_coordinates(max_real, n_sources);

// Baricenter vector allocation (both for real and int values) - Mean mode
			mean = (double*)malloc(2*n_sources*sizeof(double));

// Baricenter computation - mean mode
			mean_coordinates(baricenter_r, max_real, mean, n_sources);

//	The final baricenter coordinates are written in the logfile
			write_center_coordinates(mean, n_sources, logfile);

// Image display if visual == "y"
			if(*visual == 'y'){
				cvShowImage("removed blobs", elab_image[3]);
				cvShowImage("Colored both", elab_image_color[2]);
			}

// Image save if save == 'y'
			if(*save == 'y')
				out = fits_image_create(fits_pointer, fits_dest, "image_removed_blobs.fits", elab_image[3]->imageData, bitpix, dim, size, status, gamma_image);

//	arrays are flushed from memory
			free(max_real);
			free(max_int);
			free(baricenter_r);
			free(baricenter_int);
			free(mean);

// image release if visual == 'n'
			if(*visual == 'n'){
				cvReleaseImage(&source_image);
				cvReleaseImage(&elab_image[0]);
				cvReleaseImage(&elab_image[1]);
				cvReleaseImage(&elab_image[2]);
				cvReleaseImage(&elab_image[3]);
				cvReleaseImage(&elab_image_color[0]);
				cvReleaseImage(&elab_image_color[2]);
				cvReleaseImage(&elab_image_color[3]);
			}

		}

//	waiting for user input to close the opened images
		if(*visual == 'y')
			cvWaitKey();
		fits_close_file(fits_pointer, &status);

// final memory flush
		free(L);
		free(area);
		free(values);
		free(ratio);


		return 0;
}
}
