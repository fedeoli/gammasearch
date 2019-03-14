
#include <stdio.h>
#include "opencv/cv.h"
#include "opencv/highgui.h"
#include "fitsio.h"
#include "opencv2/core/types_c.h"
#include "opencv2/imgproc/imgproc_c.h"
#include <math.h>
#include <vtkInteractorObserver.h>


#ifndef SHAREFILE_INCLUDED
#define SHAREFILE_INCLUDED
#ifdef  MAIN_FILE
extern int glob_max;
#else
int glob_max;
#endif
#endif

//C interface function definition
fitsfile * fits_image_open (fitsfile * fits_pointer, char * image_filepath, int status, int hdu_type, int hdu_num, int hdu_nums);
fitsfile * fits_image_read_param(fitsfile * fits_pointer, int status, int bitpix, int dim, long * size);
fitsfile * fits_image_copy(fitsfile * fits_pointer, int status, long * fpixel, long * size, char * imagedata);
fitsfile * fits_image_copy_float(fitsfile * fits_pointer, int status, long * fpixel, long * size, char * imagedata);
fitsfile * fits_image_create(fitsfile * fits_source, fitsfile * fits_dest, char * filename, char * imagedata, int bitpix, int naxis, long *naxes, int status, CvSize size);
void buildHistogram(IplImage * input,unsigned int * histogram, char * window);
void visualizeHistogram(unsigned int * histogram, char * window);
double * pdf_create(unsigned int * histogram);
void imageEqualization(IplImage * input, IplImage * output, float * pdf);
void gauss_smooth_fast(IplImage * input, IplImage * output, float sigma);
void gauss_smooth(IplImage * input, IplImage * output, float sigma);
void bilateral_smooth(IplImage * input, IplImage * output, float sigma_spa, float sigma_int);
void bilateral_smooth_fast(IplImage * input, IplImage * output, float sigma_spa,  float sigma_int);
void gammaCorrection(IplImage * input, IplImage * output, double r);
void binarization_fixed_treshold(IplImage * input, IplImage * output, unsigned int treshold);
int binarization_otsu(IplImage * input, unsigned int * histogram);
int labeling(IplImage *input, int *L);
int labeling_8(IplImage *input, int *L);
void intensity_ratio(IplImage *input, IplImage *output, int *L, int *area, float treshold, int n_sources, double *ratio);
void different_intensity_sources(int n_sources, int *L, IplImage * output);
void area_blobs(int * L, int n_sources, int * area, CvSize size);
void baricenter_real(int * L, int n_sources, int * area, CvSize size, double * baricenter_real, int * baricenter_int);
void visualize_baricenter(IplImage * input, IplImage * output, int * baricenter_int, int n_sources);
void visualize_centers(IplImage *input, IplImage *output, int *baricenter_int, int n_sources, char colour);
void read_center_coordinate(char * filepath, float * coordinates);
void source_coordinate(float * center_coordinates, double * baricenter_real, float * source_coord, int n_sources);
void circle_detection(IplImage * input, IplImage * output, int * values, int treshold_int, float accept_level, int r);
void median_filter(IplImage *input, IplImage *output);
void median_cross_filter(IplImage *input, IplImage *output);
void mean_cross_filter(IplImage *input, IplImage *output);
void erase_detected (IplImage * source, IplImage * binary, IplImage * output);
void rescale(int *input, char *output, CvSize size);
void rescale_float(float *input, char *output, CvSize size);
void rescale_inv(char *input, int *output, CvSize size);
void rescale_circle(int *input, char* output, CvSize size);
void dilation(IplImage *input, IplImage *output);
void intensity_ratio(IplImage *input, IplImage *output, int *L, int *area, int treshold, int n_sources, double * ratio);
void visualize_both_centers(IplImage *input, IplImage *output, int *baricenter, int *max_int, int n_sources);
void center_values(int *values, CvSize size, int *L, double *max_real, int * max_int, int n_sources);
void double_size_colour(IplImage *input, IplImage *output);
int remove_separated_blobs(IplImage *input, IplImage *output, int *L, int *bar_int, double * bar_real, int n_sources, double *ratio, int center_distance) ;
void order_coordinates(double * array, int n_sources);
void mean_coordinates(double * bar, double * max, double * mena, int n_sources);
void write_center_coordinates(double * center_coordinates, int n_sources, char * filepath);
void change_frame(double * baricenter);
void back_divide(IplImage * input, IplImage* background, IplImage * output);


