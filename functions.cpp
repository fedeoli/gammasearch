#include "functions.h"

fitsfile * fits_image_open (fitsfile * fits_pointer, char * image_filepath, int status, int hdu_type, int hdu_num, int hdu_nums){

//	FITS file Image opening
	fits_open_file(&fits_pointer, image_filepath, READWRITE, &status);

//	 Move to the first header which contains the data
	fits_movabs_hdu(fits_pointer, 1, &hdu_type, &status);

//	 Check the position of the current header
	fits_get_hdu_num(fits_pointer, &hdu_num);
	fits_get_num_hdus(fits_pointer, &hdu_nums, &status);

	return fits_pointer;
}

fitsfile * fits_image_read_param(fitsfile * fits_pointer, int status, int bitpix, int dim, long * size) {
	int maxdim = 3;

//	FITS Image type
	fits_get_img_param(fits_pointer, maxdim, &bitpix, &dim, size, &status);
	return fits_pointer;
}

fitsfile * fits_image_copy(fitsfile * fits_pointer, int status, long * fpixel, long * size, char * imagedata) {

	int i = 0, j = 0, result = 0;
	int *pixels, *array;

//	FITS info retrieving
	pixels = (int *)malloc(size[1] * sizeof(int));
	array = (int *)malloc(size[0] * size[1] * sizeof(int));
	CvSize size_image = { size[0],size[1] };

	if (pixels == NULL) {
		printf("Memory allocation error\n");
		return(fits_pointer);
	}

	for (fpixel[1] = size[0]; fpixel[1] >= 1; fpixel[1]--) {
// 	extract row number fpixel[1] from fits pointer data and store it in array pixels
		result = fits_read_pix(fits_pointer, TINT, fpixel, size[0], NULL, pixels, NULL, &status);

//		jump out of loop on error
		if (result) {
			break;
		}

//		copy pixels content in array
		for (j = 0; j<size[1]; j++) {
			array[i*size[1] + j] = pixels[j];
		}

		i++;
	}

	rescale(array, imagedata, size_image);
//	rescale_circle(array, imagedata, size_image);

	free(pixels);
	free(array);

	return fits_pointer;
}

fitsfile * fits_image_copy_float(fitsfile * fits_pointer, int status, long * fpixel, long * size, char * imagedata) {

	int i = 0, j = 0, result = 0;
	float *pixels, *array;

//	FITS info retrieving
	pixels = (float *)malloc(size[1] * sizeof(float));
	array = (float *)malloc(size[0] * size[1] * sizeof(float));
	CvSize size_image = { size[0],size[1] };

	if (pixels == NULL) {
		printf("Memory allocation error\n");
		return(fits_pointer);
	}

	for (fpixel[1] = size[0]; fpixel[1] >= 1; fpixel[1]--) {
// 	extract row number fpixel[1] from fits pointer data and store it in array pixels
		result = fits_read_pix(fits_pointer, TFLOAT, fpixel, size[0], NULL, pixels, NULL, &status);

//		jump out of loop on error
		if (result) {
			break;
		}

//		copy pixels content in array
		for (j = 0; j<size[1]; j++) {
			array[i*size[1] + j] = pixels[j];
		}

		i++;
	}

	rescale_float(array, imagedata, size_image);

	free(pixels);
	free(array);

	return fits_pointer;
}

fitsfile * fits_image_create(fitsfile * fits_source, fitsfile * fits_dest, char * filename, char * imagedata, int bitpix, int naxis, long *naxes, int status, CvSize size){

	int maxdim = 3;
	long fpixel[2];
	int hdu_type;
	int max;

	int i,j, result;
	char *temp_image;
	char temp;

	int *pixels;
	//	 FITS info retrieving
	pixels = (int *)malloc(size.width * sizeof(int));
	temp_image = (char *)malloc(size.height*size.width * sizeof(char));

	i = 0;

	/*if (pixels == NULL) {
			printf("Memory allocation error\n");
			return(fits_source);
		}

		fpixel[0] = 1;
		fpixel[1] = 1;
		result = fits_read_pix(fits_source, TINT, fpixel, (size.height*size.width), NULL, temp_image, NULL, &status);
		if (result) {
			printf("error!\n");
		}*/

//	rescale_inv(imagedata,temp_image,size);

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			temp = imagedata[i * size.width + j];
			temp_image[(size.height - i -1) * size.width + j] = temp;
		}

	remove(filename);
	fits_create_file(&fits_dest, filename, &status);
	fits_copy_header(fits_source, fits_dest,  &status);
	fits_get_img_param(fits_source, maxdim, &bitpix, &naxis, naxes, &status);

	long nelements = size.height*size.width;
	fits_write_img(fits_dest, TBYTE, 1, nelements, temp_image, &status);

	fits_close_file(fits_dest, &status);
	free(temp_image);
	free(pixels);
	return fits_dest;
}


void buildHistogram(IplImage * input,unsigned int * histogram, char * window){
	int i=0, datavalue=0;

//	setting the histogram array to zero
	for (i = 0; i < 256; i++)
			histogram[i] = 0;

	for (int i = 0; i < input->height; i++)
			for (int j = 0; j < input->widthStep; j++)
				if (input->width != input->widthStep){
					if (((j + 1) % input->widthStep) != 0) {
						datavalue = ((uchar*) input->imageData)[i*input->widthStep + j];
						histogram[datavalue]++;
					}
				}
				else{
					datavalue = ((uchar*) input->imageData)[i*input->widthStep + j];
					histogram[datavalue]++;
				}
}

void visualizeHistogram(unsigned int * histogram, char * window){
	CvScalar val_struct;
	int outMargin = 20;
		int binw = 4;
		int binh = 3;
		int wout = (255 * binw) + (outMargin * 2);
		int hout = (100 * binh) + (outMargin * 2);
		CvSize size;
		size.width = wout;
		size.height = hout;
		IplImage * image_out = cvCreateImage(size, IPL_DEPTH_8U, 3);
		//draw output
		cvSet(image_out, cvScalar(255, 255, 255,0));
		//draw axes
		CvPoint pt1, pt2;
		int thickness = 1;
		pt1.x = outMargin;
		pt1.y = outMargin;
		pt2.x = outMargin;
		pt2.y = hout - outMargin;
		cvLine(image_out, pt1, pt2, cvScalar(0, 0, 0,0), thickness);
		pt1.x = wout - outMargin;
		pt1.y = hout - outMargin;
		cvLine(image_out, pt1, pt2, cvScalar(0, 0, 0,0), thickness);

		//draw histo
		thickness = CV_FILLED;
		unsigned int max = 0;
		for (int i = 0; i<255; i++) {
			if (histogram[i]>max)max = histogram[i];
		}
		for (int i = 0; i<255; i++) {
			if (histogram[i]>0) {
				pt1.x = outMargin + (i*binw);
				pt1.y = hout - outMargin;
				pt2.x = pt1.x + binw;
				pt2.y = pt1.y - ((histogram[i] * 100 / max)*binh);
				cvRectangle(image_out, pt1, pt2, cvScalar(255 - i, i / 3, i,0), thickness);
			}
		}

		//font init
		CvFont font;
		cvInitFont(&font, CV_FONT_HERSHEY_DUPLEX, 0.5f, 0.5f, 0, 2);
		pt1.x = outMargin;
		pt1.y = hout - (outMargin / 2);
		cvPutText(image_out, "0", pt1, &font, cvScalar(0, 0, 0,0));
		pt1.x = wout - (1.5*outMargin);
		cvPutText(image_out, "255", pt1, &font, cvScalar(0, 0, 0,0));

		//output visualization
		cvShowImage(window, image_out);
//		cvReleaseImage(&image_out);
//		cvDestroyWindow("Histogram");
		return;
}

double * pdf_create(unsigned int * histogram){

	static double pdf[256];
	int sum = 0, i=0;

	for (i = 0;i<256;i++)
		sum = sum + histogram[i];

	for (i = 0; i<256; i++){
		pdf[i]  =  ((long double)histogram[i]/sum);
	}

	return pdf;
}

void imageEqualization(IplImage * input, IplImage * output, float * pdf){

	int i=0;
	float cdf[256], temp=0;

//Cumulative distribution function creation
	for (i=0; i<256;i++ ){
		temp = temp + pdf[i];
		cdf[i] = (float) temp;
	}

//	Image equalization - double for
	for (int i = 0; i < input->height; i++)
				for (int j = 0; j < input->widthStep; j++){
					temp = (uchar) input->imageData[i*input->widthStep + j];
					output->imageData[i*output->widthStep + j] = (uchar) ((255)*cdf[(int) temp]);
				}

}

void gauss_smooth_fast(IplImage * input, IplImage * output, float sigma){

	int i, j, k, u, v, i_start, i_end, j_start, j_end, datavalue, G_smoothed;
	uchar temp;
	float G, n_f, * kernel_1D;

	IplImage * temp_image;
	CvSize size;

	size.height = input->height;
	size.width = input->width;

	temp_image = cvCreateImage(size, IPL_DEPTH_8U, 1);

	k = ceil(3*sigma);

//Kernel allocation
	kernel_1D = (float*) malloc((2*k+1)*(2*k+1)*sizeof(float));

//	Kernel assignment
	for (u = 0; u < (2 * k + 1); u++)
			kernel_1D[u] = exp(-(pow(u - k, 2)) / (2 * pow(sigma, 2))) / (sqrt(2 * M_PI) * sigma);

	//HORIZONTAL SMOOTHING
		for (i = 0; i < input->height ; i++)
			for (j = 0; j < input->width; j++) {
				G = 0;
				if (((j - k) >= 0) && ((j + k) <= input->width))

					//inside the image
					for (v = 0; v < (2 * k + 1); v++) {
						temp = ((uchar*)input->imageData)[i * input->width + (j + v - k)];
						datavalue = (int) temp;
						G = G + datavalue * kernel_1D[v];
					}
				else {

//					in the contours
//					left contour
					if ((j - k) < 0)
						j_start = -j + k;
					else
						j_start = 0;

//					right contour
					if ((j + k) >= input->width)
						j_end = input->width - j + k;
					else
						j_end = 2 * k + 1;

//					normalization factor
					n_f = 0;
					for (v = j_start; v < j_end; v++)
						n_f = n_f + kernel_1D[v];

					for (v = j_start; v < j_end; v++) {
						temp = ((uchar*) input->imageData)[i * input->width + (j + v - k)];
						datavalue = (int)temp;
						G = G + datavalue * kernel_1D[v] / n_f;
					}
				}
				G_smoothed = (uchar)G;
				temp_image->imageData[i*output->width + j] = (uchar) G_smoothed;
			}

//	VERTICAL SMOOTHING
			for (i = 0; i < input->height ; i++)
				for (j = 0; j < input->width ; j++) {
					G = 0;
					if (((i - k) >= 0) && ((i + k) <= input->height))
						//inside the image
						for (u = 0; u < (2 * k + 1); u++) {
							temp = ((uchar*)temp_image->imageData)[(i + u - k) * temp_image->width + j];
							datavalue = (int) temp;
							G = G + datavalue * kernel_1D[u];
						}
					else {
						if ((i - k) < 0)
							i_start = -i + k;
						else
							i_start = 0;

						if ((i + k) >= output->height)
							i_end = output->height  - i + k;
						else
							i_end = 2 * k + 1;

						n_f = 0; //normalization factor
						for (u = i_start; u < i_end; u++)
							n_f = n_f + kernel_1D[u];

						for (u = i_start; u < i_end; u++) {
							temp = ((uchar*) temp_image->imageData)[(i + u - k) * temp_image->width + j];
							datavalue = (int) temp;
							G = G + datavalue * kernel_1D[u] / n_f;
						}
					}
					G_smoothed = (uchar) G;
					output->imageData[i*output->width + j] = (uchar)G_smoothed;
				}

			free(kernel_1D);
			cvReleaseImage(&temp_image);

}

void gammaCorrection(IplImage * input, IplImage * output, double r){
	int i,j;

	for (i=0; i<input->height; i++)
		for (j=0;j<input->widthStep; j++){
			output->imageData[i*output->widthStep + j] = (uchar) pow(input->imageData[i*input->widthStep + j], r)*pow(255, 1-r);
		}
}

void bilateral_smooth(IplImage * input, IplImage * output, float sigma_spa, float sigma_int){

	int i, j, k, u, v, i_start, i_end, j_start, j_end, intensity_center, temp;
	float G, datavalue, n_f, G_val;

	k = ceil(3*sigma_spa);

	float* kernel_1D = (float*)malloc((2 * k + 1) * (2 *k + 1)*sizeof(float));

		for(u=0;u<(2*k+1);u++)
			for(v=0;v<(2*k+1);v++)
				kernel_1D[u*(2*k+1)+v]= exp(-(pow(u - k, 2) + pow(v - k, 2)) / (2 * pow(sigma_spa, 2))) / (2 * M_PI* pow(sigma_spa, 2));

	for (i = 0; i < input->height; i++)
			for (j = 0; j < input->widthStep; j++){

				G = 0;
				n_f = 0;
				temp = ((uchar*)input->imageData)[i * input->widthStep + j];
				intensity_center = (int)temp;

				if (((i - k) >= 0) && ((j - k) >= 0) && ((i + k) <= input->height) && ((j + k) <= input->widthStep)) {
					//inside the image

					for (u = 0; u < (2 * k + 1); u++)
						for (v = 0; v < (2 * k + 1); v++) {
							temp = ((uchar*)input->imageData)[(i + u - k)* input->widthStep + (j + v - k)];
							datavalue = (int)temp;
							n_f = n_f + kernel_1D[u * (2 * k + 1) + v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)* sigma_int);
						}

					for (u = 0; u < (2 * k + 1); u++)
						for (v = 0; v < (2 * k + 1); v++) {
							temp = ((uchar*)input->imageData)[(i + u - k)* input->widthStep + (j + v - k)];
							datavalue = (int)temp;
							G = G + datavalue * kernel_1D[u * (2 * k + 1) + v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)* sigma_int*n_f);
						}
				}
				else {
					//contours

					if ((i - k) < 0)
						i_start = -i + k;
					else
						i_start = 0;

					if ((i + k) >= input->height)
						i_end = k + input->height - i;
					else
						i_end = 2 * k + 1;

					if ((j - k) < 0)
						j_start = -j + k;
					else
						j_start = 0;

					if ((j + k) >= input->widthStep)
						j_end = k + input->widthStep - j;
					else
						j_end = 2 * k + 1;

					for (u = i_start; u < i_end; u++)
						for (v = j_start; v < j_end; v++) {
							temp = ((uchar*) input->imageData)[(i + u - k)* input->widthStep + (j + v - k)];
							datavalue = (int) temp;
							n_f = n_f + kernel_1D[u * (2 * k + 1) + v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)* sigma_int);
						}

					for (u = i_start; u < i_end; u++)
						for (v = j_start; v < j_end; v++) {
							temp = ((uchar*) input->imageData)[(i + u - k)*input->widthStep + (j + v - k)];
							datavalue = (int)temp;
							G = G + datavalue * kernel_1D[u * (2 * k + 1) + v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)* sigma_int*n_f) ;
						}
				}
					G_val = (int)G;
					output->imageData [i*output->widthStep+j]= (uchar)G_val;
				}

		free(kernel_1D);

}

void gauss_smooth(IplImage * input, IplImage * output, float sigma){

	int i, j, u, v,  datavalue;
	uchar temp;
	float k1,k2;

	CvSize size;

	size.height = input->height;
	size.width = input->width;


//	gaussian kernel evaluation
//	k = ceil(3*sigma);


//	kernel allocation
	float* kernel = (float*)malloc(input->height*input->width*sizeof(float));
	int* out = (int*)malloc(input->height*input->width*sizeof(int));
	k1 = (float)(input->height-1)/2;
	k2 = (float)(input->width-1)/2;
	int datavalue_final;


//	kernel assignement
		for(u=0;u<(input->height);u++)
			for(v=0;v<(input->width);v++)
				kernel[u*(input->width)+v]= exp(-(pow(u - k1, 2) + pow(v - k2, 2)) / (2 * pow(sigma, 2))) / (2 * M_PI* pow(sigma, 2));

//	gaussian smoothing through the image
	for (i = 0; i < input->height; i++)
			for (j = 0; j < input->width; j++) {
				temp = ((uchar*) input->imageData)[i* input->width + j];
				datavalue = (int) temp;
				datavalue_final = (int) (datavalue/kernel[i* input->width + j]);
				if (datavalue_final > 255)
//					printf("%d\n", datavalue_final);
				out[i*input->width + j] = datavalue_final;
			}

	rescale(out, output->imageData, size);

	free(kernel);
	free(out);

}

void bilateral_smooth_fast(IplImage * input, IplImage * output, float sigma_spa, float sigma_int){

	int i, j, k, u, v, i_start, i_end, j_start, j_end, intensity_center, temp, datavalue, G_val;
	float * kernel_1D, G, n_f;

	IplImage * temp_image;
	CvSize temp_size;

	temp_size.height = input->height;
	temp_size.width = input->width;

	temp_image = cvCreateImage(temp_size, IPL_DEPTH_8U, 1);

	k = ceil(3*sigma_spa);

	//Kernel allocation
		kernel_1D = (float*) malloc((2*k+1)*(2*k+1)*sizeof(float));

	//	Kernel computation
		for (u = 0; u < (2 * k + 1); u++)
				kernel_1D[u] = exp(-(pow(u - k, 2)) / (2 * pow(sigma_spa, 2))) / (sqrt(2 * M_PI) * sigma_spa);


	//Bilateral double 1D Gaussian smoothing (fast bilateral smoothing)

	//	HORIZONTAL SMOOTHING
	for (i = 0; i < input->height; i++)
		for (j = 0; j < input->widthStep; j++) {
			G = 0;
			n_f = 0;
			temp = ((uchar*) input->imageData)[i * input->widthStep + j];
			intensity_center = (int) temp;

			if (((j - k) >= 0) && ((j + k) <= input->widthStep)) {
				//inside the image
				for (v = 0; v < (2 * k + 1); v++) {
					temp = ((uchar*) input->imageData)[i * input->widthStep + (j + v - k)];
					datavalue = (int) temp;
					n_f = n_f + kernel_1D[v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int);
				}
				for (v = 0; v < (2 * k + 1); v++) {
					temp = ((uchar*)input->imageData)[i * input->widthStep + (j + v - k)];
					datavalue = (int) temp;
					G = G + datavalue * kernel_1D[v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int*n_f);
				}
			}
			else {

				//contours
				if ((j - k) < 0)
					j_start = -j + k;
				else
					j_start = 0;

				if ((j + k) >= input->widthStep)
					j_end = input->widthStep - j + k;
				else
					j_end = 2 * k + 1;

				for (v = j_start; v < j_end; v++) {
					temp = ((uchar*) input->imageData)[i * input->widthStep + (j + v - k)];
					datavalue = (int) temp;
					n_f = n_f + kernel_1D[v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int);
				}
				for (v = j_start; v < j_end; v++) {
					temp = ((uchar*) input->imageData)[i * input->widthStep + (j + v - k)];
					datavalue = (int) temp;
					G = G + datavalue * kernel_1D[v] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int*n_f);
				}
			}
			G_val = (int)G;
			temp_image->imageData[i*temp_image->widthStep + j] = (uchar)G_val;
		}
//	cvShowImage("temp", temp_image);


	//VERTICAL SMOOTHING
	for (i = 0; i < input->height; i++)
		for (j = 0; j < input->widthStep; j++) {
			G = 0;
			n_f = 0;
			temp = ((uchar*)temp_image->imageData)[i * temp_image->widthStep + j];
			intensity_center = (int) temp;

			if (((i - k) >= 0) && ((i + k) <= input->height)){
				//inside the image

				for (u = 0; u < (2 * k + 1); u++) {
					temp = ((uchar*) temp_image->imageData)[(i + u - k)* temp_image->widthStep + j];
					datavalue = (int) temp;
					n_f = n_f + kernel_1D[u] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int);
				}
				for (u = 0; u < (2 * k + 1); u++) {
					temp = ((uchar*) temp_image->imageData)[(i + u - k)* temp_image->widthStep + j];
					datavalue = (int) temp;
					G = G + datavalue * kernel_1D[u] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int*n_f);
				}
			}
			else {
				if ((i - k) < 0)
					i_start = -i + k;
				else
					i_start = 0;

				if ((i + k) >= output->height)
					i_end = output->height - i + k;
				else
					i_end = 2 * k + 1;

				n_f = 0; //normalization factor
				for (u = i_start; u < i_end; u++) {
					temp = ((uchar*)temp_image->imageData)[(i + u - k)* temp_image->widthStep + j];
					datavalue = (int) temp;
					n_f = n_f + kernel_1D[u] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int);
				}
				for (u = i_start; u < i_end; u++) {
					temp = ((uchar*) temp_image->imageData)[(i + u - k)* temp_image->widthStep + j];
					datavalue = (int) temp;
					G = G + datavalue * kernel_1D[u] * exp(-pow((datavalue - intensity_center), 2) / (2 * pow(sigma_int, 2))) / (sqrt(2 * M_PI)*sigma_int*n_f);
				}
			}
			G_val = (int)G;
			output->imageData[i*output->widthStep + j] = (uchar)G_val;
		}
//	cvShowImage("temp", temp_image);

	free(kernel_1D);
	cvReleaseImage(&temp_image);


}

void binarization_fixed_treshold(IplImage * input, IplImage * output, unsigned int treshold){

	int i,j;
	unsigned int datavalue;
	uchar temp;

	for(i=0; i< input->height;i++)
		for(j=0;j<input->widthStep;j++){
			temp = ((uchar*) input->imageData)[i*input->widthStep +j];
			datavalue = (int) temp;
			if(datavalue >= treshold)
				output->imageData[i*output->widthStep + j] = (uchar)255;
			else
				output->imageData[i*output->widthStep + j] = (uchar)0;
		}


}

int binarization_otsu(IplImage * input, unsigned int * histogram){

	int treshold = 0;
	double q1, mu_1, mu_2, mu, temp, var_b, temp_var_b;
	int i;
	double *  pmf;

	pmf = pdf_create(histogram);

	mu = 0;
	for (i = 0; i < 256; i++)
		mu = mu + pmf[i] * i;

	//	variables initialization
	q1 = pmf[0];
	mu_1 = 0;
	mu_2 = (mu - q1 * mu_1) / (1 - q1);
	var_b = q1 * (1 - q1)*pow((mu_1 - mu_2), 2);
	temp_var_b = var_b;


	for (i = 1; i < 256; i++) {
		if (1-q1 > 0.001){
		temp = q1;
		q1 = q1 + pmf[i];
		mu_1 = (temp * mu_1 + i * pmf[i]) / q1;
		mu_2 = (mu - q1 * mu_1) / (1 - q1);

		temp_var_b = q1*pow((mu_1 - mu), 2)/(1 - q1);
		}

		if (temp_var_b > var_b) {
			var_b = temp_var_b;
			treshold = i;
		}

	}

	return treshold;
}

int labeling(IplImage *input, int *L) {
	int i, j, k, datavalue, lp, lq, newlabel = 0, length_C, n, temp_label;
	uchar temp;
	int *C = (int*)malloc(input->height*input->widthStep * sizeof(int));

	//	Background Label
	C[0] = 0;

	//	First scan of the image
	for(i=0; i<input->height;i++)
		for (j = 0; j < input->widthStep; j++) {

			temp = ((uchar*)input->imageData)[i * input->widthStep + j];
			datavalue = (int)temp;

			if (datavalue == 0)
				//	Background
				L[i*input->widthStep + j] = 0;
			else {
				//	Foreground

				if (i == 0)
					lp = 0;
				else
					lp = L[(i - 1) * input->widthStep + j];

				if (j == 0)
					lq = 0;
				else
					lq = L[i * input->widthStep + (j - 1)];

				if ((lp == 0) && (lq == 0)) {
					//	both p and q are background pixels
					newlabel++;
					C[newlabel] = newlabel;
					L[i * input->widthStep + j] = newlabel;
				}
				else
					if(lp==0)
						//	p is a background pixel, q a foreground one
						L[i * input->widthStep + j] = lq;
					else
						if((lq==0)|(lp==lq))
							//	(q is a background pixel, p a foreground one) or (p and q has the same foreground label)
							L[i * input->widthStep + j] = lp;
						else {
							//	equivalence
							L[i * input->widthStep + j] = lp;
							for (k = 0; k <= newlabel; k++)
								if (C[k] == lq)
									C[k] = lp;
						}
			}
		}

	length_C = newlabel+1;

	//	n: number of different sources (the number of positive labels)
	n = 0;

	//	Remapping vector C with consecutive labels
	for(i=1; i < length_C; i++)
		if (C[i] >= (n+1)) {
			n++;
			temp_label = C[i];
			for (j = i; j < length_C; j++)
				if (C[j] == temp_label)
					C[j] = n;
		}

	//	Second scan of the image
	for (i = 0; i < input->height; i++)
		for (j = 0; j < input->widthStep; j++)
			if (L[i * input->widthStep + j] != 0)
				L[i * input->widthStep + j] = C[L[i * input->widthStep + j]];

	free(C);

	return n;
}

void different_intensity_sources(int n_sources, int *L, IplImage * output) {
	int i, j;
	for (i = 0; i < output->height; i++)
		for (j = 0; j < output->width; j++)
			if (n_sources){
				if (L[i*output->width + j] == 0)
					output->imageData[i*output->width + j] = (uchar)0;
				else
					output->imageData[i*output->width + j] = (uchar)(55 + (int)(200 * L[i*output->width + j] / n_sources));
			}
			else
				output->imageData[i*output->width + j] = (uchar)0;
}

void area_blobs(int * L, int n_sources, int * area, CvSize size){

	int i,j;

//	area array initialisation
	for(i=0;i<=n_sources;i++)
		area[i] =0;

	for (i=0;i<size.height; i++){
		for(j=0;j<size.width;j++){
				area[L[i*size.width + j]]++;
		}
	}
}

void baricenter_real(int * L, int n_sources, int * area, CvSize size, double * baricenter_real, int * baricenter_int){

	int i,j;

//	area array initialisation
	for(i=0;i<2*n_sources;i++)
		baricenter_real[i] =0;

	for (i=0;i<size.height; i++){
			for(j=0;j<size.width;j++){
				if (L[i*size.width + j] != 0){
					baricenter_real[2*(L[i*size.width + j] - 1)] = baricenter_real[2*(L[i*size.width + j] - 1)] + i;
					baricenter_real[2*(L[i*size.width + j] - 1)  + 1] = baricenter_real[2*(L[i*size.width + j] - 1)  + 1] + j;
				}
			}
		}

	for(i=0;i<n_sources;i++){
		baricenter_real[2*i] =baricenter_real[2*i]/area[i+1];
		baricenter_real[2*i + 1] = baricenter_real[2*i + 1]/area[i+1];
	}

//	printf("\n");
//	for(i=0;i<n_sources;i++)
//		printf("Baricenter blob %d = (%g,%g)\n", i+1, baricenter_real[2*i], baricenter_real[2*i+1]);

	for(i=0;i<n_sources;i++){
		baricenter_int[2*i] = round(baricenter_real[2*i]);
		baricenter_int[2*i + 1] = round(baricenter_real[2*i + 1]);
	}

//	printf("\n");

}

void visualize_baricenter(IplImage * input, IplImage * output, int * baricenter_int, int n_sources){

int i, i_b, j_b, j, k;

cvCvtColor(input, output, CV_GRAY2BGR);

for(i=0;i<n_sources;i++){
	i_b = baricenter_int[2*i];
	j_b = baricenter_int[2*i + 1];
	for (j=-1;j<2;j++){
		for(k=-1;k<2;k++){
			if( ((i_b + j < input->height ) && (i_b + j >= 0))  &&  ((j_b + k < input->width)  && (j_b + k >= 0))  ){
				output->imageData[3 * ((i_b + j)*output->width + (j_b + k)) ] = 0;
				output->imageData[3 * ((i_b + j)*output->width + (j_b + k)) + 1] = 0;
				output->imageData[3 * ((i_b + j)*output->width + (j_b + k)) + 2] = 255;
			}
		}
	}
}

}

void visualize_centers(IplImage *input, IplImage *output, int *baricenter_int, int n_sources, char colour) {
	int i, j, k, i_b, j_b;

	cvCvtColor(input, output, CV_GRAY2BGR);

	for (i = 0; i < n_sources; i++) {
		i_b = baricenter_int[2 * i];
		j_b = baricenter_int[2 * i + 1];

		for(j = -1; j <= 1; j++)
			for (k = -1; k <= 1; k++)
				if (((i_b + j) < input->height)&&((i_b + j) >= 0)&&((j_b + k) < input->width) && ((j_b + k) >= 0)){
					if (colour == 'r') {
						output->imageData[3 * ((i_b + j) * input->width + (j_b + k))] = 0;
						output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 1] = 0;
						output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 2] = 255;
					}
					else
						if (colour == 'b') {
							output->imageData[3 * ((i_b + j) * input->width + (j_b + k))] = 150;
							output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 1] = 0;
							output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 2] = 0;
						}
						else
							if (colour == 'g') {
								output->imageData[3 * ((i_b + j) * input->width + (j_b + k))] = 0;
								output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 1] = 150;
								output->imageData[3 * ((i_b + j) * input->width + (j_b + k)) + 2] = 0;
							}
				}
	}

}

void read_center_coordinate(char * filepath, float * coordinates){

	int cmp, i;
	char buff[50], ra[8], dec[7], temp;
	FILE * coord_file = fopen( filepath, "r");

//	RA initialisation
	for (i=0;i<8;i++)
		if (i == 3)
			ra[i] = '.';
		else
			ra[i] = 0;

//	 DEC initialisation
	for (i=0;i<7;i++)
		if (i == 2)
			dec[i] = '.';
		else
			dec[i] = 0;

//	printf("\n\n");
	do {
		fscanf(coord_file, "%s", buff);
		cmp = strncmp(buff, "ra=", 3);
	} while( cmp != 0 );


	i = 4;
	do {
		temp  = buff[i];
		if (temp != '"')
			ra[i-4] = temp;
		i++;
	} while( temp != '"' );
	coordinates[0] = atof(ra);


	fscanf(coord_file, "%s", buff);
	i = 5;
		do {
			temp  = buff[i];
			if (temp != '"')
				dec[i-5] = temp;
			i++;
		} while( temp != '"' );
		coordinates[1] = atof(dec);

	fclose(coord_file);

}

void source_coordinate(float * center_coordinates, double * baricenter_real, float * source_coord, int n_sources){

	int i;
	float i_bar, j_bar;


	for (i=0;i<n_sources;i++){
//		printf("%d\tok\n", i);
		i_bar = baricenter_real[2*i] - 100;
		j_bar = baricenter_real[2*i + 1] - 100;
		source_coord[2*i] = center_coordinates[0] - j_bar*(0.028);
		source_coord[2*i+1] = center_coordinates[1] - i_bar*(0.02);

//		printf("source[%d]:\tRA:\t%f\tDEC:\t%f\n", i+1, source_coord[2*i], source_coord[2*i+1]);
	}
//	printf("\n");

}

void circle_detection(IplImage * input, IplImage * output, int * values, int treshold_int, float accept_level, int r){
	int i, j, m, n, F, circ, intensity, A, accept_level_int, area_full=0;
	uchar temp;

	for (m = - r; m <= r; m++)
		for (n = -r; n <= r; n++) {
			circ = pow(m , 2) + pow(n , 2);
			if (circ <= r * r)
				area_full++;
		}

//	printf("Area = %d \n", area_full);
	accept_level_int = (int)(area_full * accept_level);
//	printf("Accept level = %d \n", accept_level_int);

	for (i = r; i < (input->height -r); i++)
		for (j = r; j < (input->width -r); j++) {
			F = 0;
			m = i-r;
			n = j-r;
			A = 0;

			do {
				circ = pow(m - i, 2) + pow(n - j, 2);
				if (circ <= pow(r, 2)) {
					A++;
					temp = ((uchar*) input->imageData)[m * input->width + n];
					intensity = (int) temp;
					if (intensity > treshold_int)
						F++;
				}
				if (n == (j + r)) {
					n = j - r;
					m++;
				} else
					n++;
			} while (((F + area_full - A) >= accept_level_int) && (m <= i + r) );




				if (F >= accept_level_int) {
					values[i * input->width + j] = F;
					output->imageData[i * input->width + j] = 255;
				}
				else{
					values[i * input->width + j] = 0;
					output->imageData[i * input->width + j] = 0;
				}
			}


	for (i = 0; i < r; i++)
			for (j = 0; j < input->width; j++) {
				F = 0;
				A = area_full;
				for (m = i - r; m <= i + r; m++)
					for (n = j - r; n <= j + r; n++) {
						circ = pow(m - i, 2) + pow(n - j, 2);
						if (circ <= pow(r,2))
							if ((m >= 0) && (n >= 0) && (m < input->height) && (n < input->width)) {
								temp = ((uchar*)input->imageData)[m * input->width + n];
								intensity = (int)temp;
								if (intensity > treshold_int)
									F++;
							}
							else
								A--;

						accept_level_int = (int)(A * accept_level);

						if (F >= accept_level_int) {
							values[i * input->width + j] = F;
							output->imageData[i * input->width + j] = 255;
						}
						else {
							values[i * input->width + j] = 0;
							output->imageData[i * input->width + j] = 0;
						}
					}
			}

	for ((i = input->height - r); i < input->height; i++)
				for (j = 0; j < input->width; j++) {
					F = 0;
					A = area_full;
					for (m = i - r; m <= i + r; m++)
						for (n = j - r; n <= j + r; n++) {
							circ = pow(m - i, 2) + pow(n - j, 2);
							if (circ <= pow(r,2))
								if ((m >= 0) && (n >= 0) && (m < input->height) && (n < input->width)) {
									temp = ((uchar*)input->imageData)[m * input->width + n];
									intensity = (int)temp;
									if (intensity > treshold_int)
										F++;
								}
								else
									A--;

							accept_level_int = (int)(A * accept_level);

							if (F >= accept_level_int) {
								values[i * input->width + j] = F;
								output->imageData[i * input->width + j] = 255;
							}
							else {
								values[i * input->width + j] = 0;
								output->imageData[i * input->width + j] = 0;
							}
						}
				}

	for (i = r; i < input->height-r; i++)
				for (j = 0; j < r; j++) {
					F = 0;
					A = area_full;
					for (m = i - r; m <= i + r; m++)
						for (n = j - r; n <= j + r; n++) {
							circ = pow(m - i, 2) + pow(n - j, 2);
							if (circ <= pow(r,2))
								if ((m >= 0) && (n >= 0) && (m < input->height) && (n < input->width)) {
									temp = ((uchar*)input->imageData)[m * input->width + n];
									intensity = (int)temp;
									if (intensity > treshold_int)
										F++;
								}
								else
									A--;

							accept_level_int = (int)(A * accept_level);

							if (F >= accept_level_int) {
								values[i * input->width + j] = F;
								output->imageData[i * input->width + j] = 255;
							}
							else {
								values[i * input->width + j] = 0;
								output->imageData[i * input->width + j] = 0;
							}
						}
				}

	for (i = r; i < input->height-r; i++)
					for (j = input->width -r ; j < input->width; j++) {
						F = 0;
						A = area_full;
						for (m = i - r; m <= i + r; m++)
							for (n = j - r; n <= j + r; n++) {
								circ = pow(m - i, 2) + pow(n - j, 2);
								if (circ <= pow(r,2))
									if ((m >= 0) && (n >= 0) && (m < input->height) && (n < input->width)) {
										temp = ((uchar*)input->imageData)[m * input->width + n];
										intensity = (int)temp;
										if (intensity > treshold_int)
											F++;
									}
									else
										A--;

								accept_level_int = (int)(A * accept_level);

								if (F >= accept_level_int) {
									values[i * input->width + j] = F;
									output->imageData[i * input->width + j] = 255;
								}
								else {
									values[i * input->width + j] = 0;
									output->imageData[i * input->width + j] = 0;
								}
							}
					}
}

void median_filter(IplImage *input, IplImage *output) {
	int i, j, m, n, k, datavalue, min;
	uchar temp, temp_central;
	int kernel[9], ord[9];

	for(i=0; i<input->height; i++)
		for (j = 0; j < input->width; j++) {
			k = 0;
			temp_central = ((uchar*)input->imageData)[ i* input->width + j];
			for(m = -1; m <= 1; m++)
				for (n = -1; n <= 1; n++) {
					if ( (i+m <0) || (i+m >= input->height)  || (j+n < 0) || (j+n >= input->width))
						datavalue = temp_central;
					else	{
						temp = ((uchar*)input->imageData)[(i + m)* input->width + (j + n)];
						datavalue = (int)temp;
					}
					kernel[k] = datavalue;
					k++;
				}

			for (m = 0; m < 8; m++) {
				min = kernel[m];
				for (n = m + 1; n < 9; n++)
					if (kernel[n] < min)
						min = kernel[n];
				ord[m] = min;
			}
			output->imageData[i*input->width + j] = (uchar)ord[7];
		}
}

void median_cross_filter(IplImage *input, IplImage *output) {
	int i, j, m, n, k, datavalue, min, central_pixel;
	uchar temp;
	int kernel[5], ord[5];

	for (i = 0; i<input->height; i++)
		for (j = 0; j < input->width; j++) {

			temp = ((uchar*)input->imageData)[i * input->width + j];
			central_pixel = (int)temp;

			if (i > 0) {
				temp = ((uchar*)input->imageData)[(i - 1)* input->width + j];
				datavalue = (int)temp;
				kernel[0] = datavalue;
			}
			else
				kernel[0] = central_pixel;

			if (j > 0) {
				temp = ((uchar*)input->imageData)[i * input->width + j - 1];
				datavalue = (int)temp;
				kernel[1] = datavalue;
			}
			else
				kernel[1] = central_pixel;

			kernel[2] = central_pixel;

			if (i < (input->height - 1)) {
				temp = ((uchar*)input->imageData)[(i + 1) * input->width + j];
				datavalue = (int)temp;
				kernel[3] = datavalue;
			}
			else
				kernel[3] = central_pixel;

			if (j < (input->width - 1)) {
				temp = ((uchar*)input->imageData)[i * input->width + j + 1];
				datavalue = (int)temp;
				kernel[4] = datavalue;
			}
			else
				kernel[4] = central_pixel;

			for (m = 0; m < 4; m++) {
				min = kernel[m];
				for (n = m + 1; n < 5; n++)
					if (kernel[n] < min)
						min = kernel[n];
				ord[m] = min;
			}

			output->imageData[i*input->width + j] = (uchar)ord[3];
		}
}

void mean_cross_filter(IplImage *input, IplImage *output) {
	int i, j, m, n, k, datavalue, min, central_pixel;
	uchar temp;
	int kernel[5], ord[5];
	float sum;

	for (i = 0; i<input->height; i++)
		for (j = 0; j < input->width; j++) {

			temp = ((uchar*)input->imageData)[i * input->width + j];
			central_pixel = (int)temp;

			if (i > 0) {
				temp = ((uchar*)input->imageData)[(i - 1)* input->width + j];
				datavalue = (int)temp;
				kernel[0] = datavalue;
			}
			else
				kernel[0] = central_pixel;

			if (j > 0) {
				temp = ((uchar*)input->imageData)[i * input->width + j - 1];
				datavalue = (int)temp;
				kernel[1] = datavalue;
			}
			else
				kernel[1] = central_pixel;

			kernel[2] = central_pixel;

			if (i < (input->height - 1)) {
				temp = ((uchar*)input->imageData)[(i + 1) * input->width + j];
				datavalue = (int)temp;
				kernel[3] = datavalue;
			}
			else
				kernel[3] = central_pixel;

			if (j < (input->width - 1)) {
				temp = ((uchar*)input->imageData)[i * input->width + j + 1];
				datavalue = (int)temp;
				kernel[4] = datavalue;
			}
			else
				kernel[4] = central_pixel;

			sum = 0;
			for (m = 0; m < 4; m++) {
				sum = sum + kernel[m];
			}
			sum = sum/5;

			output->imageData[i*input->width + j] = (uchar)sum;
		}
}

void erase_detected (IplImage * source, IplImage * binary, IplImage * output){
		int i, j, datavalue;
		uchar temp;
		CvSize size = { source->height , source->width };
		int *temp_image=(int*) malloc(source->height * source->width * sizeof(int));

		for(i = 0; i < source->height; i++)
			for (j = 0; j < source->width; j++) {
				temp = ((uchar*)binary->imageData)[i * binary->width + j ];
				datavalue = (int)temp;
				if (datavalue)
					temp_image[i * output->width + j] = 0;
				else {
					temp = ((uchar*)source->imageData)[i * source->width + j];
					datavalue = (int)temp;
					temp_image[i * output->width + j] = datavalue;
				}
			}

		rescale( temp_image, output->imageData, size);
		free(temp_image);
	}

void rescale(int *input, char *output, CvSize size) {
	int i, j;
	long int max, temp;

	max = 0;
	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			if (input[i * size.width + j]>max)
				max = input[i * size.width + j];
		}
//	printf("max input: %d\n", max);


	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			temp = (int)(input[i * size.width + j] * 255 / max);
			output[i * size.width + j] = (uchar)temp;
		}
}

void rescale_float(float *input, char *output, CvSize size) {
	int i, j;
	float max, temp;
	int temp_temp = 0;

	max = 0;
	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			if (input[i * size.width + j]>max)
				max = input[i * size.width + j];
		}


	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			temp = float(input[i * size.width + j] * 255/ max);
			temp_temp = round(temp);
			output[i * size.width + j] = (uchar)temp_temp;
		}
}

void rescale_inv(char *input, int *output, CvSize size) {
	int i, j, temp, max_dest = 0;
	char max_source = 0;

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			if (input[i * size.width + j]>max_source)
				max_source = input[i * size.width + j];
			if (output[i * size.width + j]>max_dest)
				max_dest = output[i * size.width + j];
		}

	printf("max imagedata: %d\n", max_source);
	printf("max fits before rescaling: %d\n", max_dest);

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			temp = (float)((uint8_t)input[i * size.width + j]*max_dest/(uint8_t)max_source);
			output[(size.height - i -1) * size.width + j] = temp;
		}

	for (i = 0; i < size.height; i++)
			for (j = 0; j < size.width; j++) {
				if (input[i * size.width + j]>max_source)
					max_source = input[i * size.width + j];
				if (output[i * size.width + j]>max_dest)
					max_dest = output[i * size.width + j];
			}

		printf("max imagedata: %d\n", max_source);
		printf("max fits after rescaling: %d\n", max_dest);

}

void dilation(IplImage *input, IplImage *output) {
	int i, j, m, n, datavalue;
	uchar temp;

	for (i = 0; i < input->height; i++)
		for (j = 0; j < input->width; j++)
			output->imageData[i * input->width + j] = input->imageData[i * input->width + j];

	for( i = 0; i < input->height; i++)
		for (j = 0; j < input->width; j++) {
			temp = ((uchar*)input->imageData)[i * input->width + j];
			datavalue = (int)temp;
			if (datavalue) {
				for(m=-4;m<=4;m++)
					for(n=-4;n<=4;n++)
						output->imageData[(i + m) * input->width + (j + n)] = 255;
			}

		}
}

void intensity_ratio(IplImage *input, IplImage *output, int *L, int *area, float treshold, int n_sources, double * ratio) {
	int i, j, k, datavalue;
	uchar temp;
	float *mu = (float*)malloc((n_sources+1)* sizeof(float));
	float compare;

	for (k = 0; k <= n_sources; k++)
		mu[k] = 0;

// mu = array whose elements consists of the summation of the blob's pixel intensities (from source image)
	for(i=0;i<input->height;i++)
		for (j = 0; j < input->width; j++) {
			temp = ((uchar*)input->imageData)[i * input->width + j];
			datavalue = (int)temp;
			mu[L[i*input->width + j]] = mu[L[i*input->width + j]] + datavalue;
		}



	// normalization wr2 blob area (in order to kill spot high-intensity background pixels)
	for (k = 0; k <= n_sources; k++)
		mu[k] = mu[k] / area[k];

	// normalization wr2 background average intensity
	for (k = 0; k <= n_sources; k++){
		ratio[k] = mu[k]/mu[0];
	}


	for (i = 0; i<input->height; i++)
		for (j = 0; j < input->width; j++) {
			compare = ratio[L[i*input->width + j]];
			if (compare > treshold)
				output->imageData[i*input->width + j] = (uchar)255;
			else
				output->imageData[i*input->width + j] = (uchar)0;
		}

	free(mu);
}

int labeling_8(IplImage *input, int *L) {
	int i, j, k, datavalue, lp, lq, lr, ls, newlabel = 0, length_C, n, temp_label;
	uchar temp;
	int *C = (int*)malloc(input->height*input->width * sizeof(int));

	//	Background Label
	C[0] = 0;

	//	First scan of the image
	for (i = 0; i<input->height; i++)
		for (j = 0; j < input->width; j++) {

			temp = ((uchar*)input->imageData)[i * input->width + j];
			datavalue = (int)temp;

			if (datavalue == 0)
				//	Background
				L[i*input->width + j] = 0;
			else {
				//	Foreground

				if (i == 0) {
					lp = 0;
					lq = 0;
					lr = 0;
				}
				else {
					lp = L[(i - 1) * input->width + (j - 1)];
					lq = L[(i - 1) * input->width + j];
					lr = L[(i - 1) * input->width + (j + 1)];
				}

				if (j == 0) {
					lp = 0;
					ls = 0;
				}
				else
					ls = L[i * input->width + (j - 1)];

				//	CLASS ASSIGNEMENT

				if ((lp == 0) && (lq == 0) && (lr == 0) && (ls == 0)) {
					//	p, q, r and s are background pixels
					newlabel++;
					C[newlabel] = newlabel;
					L[i * input->width + j] = newlabel;
				}
				else
					if ((lp == 0) && (lq == 0) && (lr == 0))
						//	s is the only foreground pixel
						L[i * input->width + j] = ls;
					else
						if ((lp == 0) && (lq == 0) && (ls == 0))
							//	r is the only foreground pixel
							L[i * input->width + j] = lr;
						else
							if ((lp == 0) && (lr == 0) && (ls == 0))
								//	q is the only foreground pixel
								L[i * input->width + j] = lq;
							else
								if ((lq == 0) && (lr == 0) && (ls == 0))
									//	p is the only foreground pixel
									L[i * input->width + j] = lp;
								else
									if ((lp == 0) && (lq == 0)) {
										L[i * input->width + j] = lr;
										if (lr != ls) {
											//	equivalence
											for (k = 0; k <= newlabel; k++)
												if (C[k] == ls)
													C[k] = lr;
										}
									}
									else
										if ((lp == 0) && (lr == 0)) {
											L[i * input->width + j] = lq;
											if (lq != ls) {
												//	equivalence
												for (k = 0; k <= newlabel; k++)
													if (C[k] == ls)
														C[k] = lq;
											}
										}
										else
											if ((lp == 0) && (ls == 0)) {
												L[i * input->width + j] = lq;
												if (lq != lr) {
													//	equivalence
													for (k = 0; k <= newlabel; k++)
														if (C[k] == lr)
															C[k] = lq;
												}
											}
											else
												if ((lq == 0) && (lr == 0)) {
													L[i * input->width + j] = lp;
													if (lp != ls) {
														//	equivalence
														for (k = 0; k <= newlabel; k++)
															if (C[k] == ls)
																C[k] = lp;
													}
												}
												else
													if ((lq == 0) && (ls == 0)) {
														L[i * input->width + j] = lp;
														if (lp != lr) {
															//	equivalence
															for (k = 0; k <= newlabel; k++)
																if (C[k] == lr)
																	C[k] = lp;
														}
													}
													else
														if ((lr == 0) && (ls == 0)) {
															L[i * input->width + j] = lp;
															if (lp != lq) {
																//	equivalence
																for (k = 0; k <= newlabel; k++)
																	if (C[k] == lq)
																		C[k] = lp;
															}
														}
														else
															if (lp == 0) {
																L[i * input->width + j] = lq;
																if (lq != lr) {
																	//	equivalence
																	for (k = 0; k <= newlabel; k++)
																		if (C[k] == lr)
																			C[k] = lq;
																}
																if (lq != ls) {
																	//	equivalence
																	for (k = 0; k <= newlabel; k++)
																		if (C[k] == ls)
																			C[k] = lq;
																}
															}
															else
																if (lq == 0) {
																	L[i * input->width + j] = lp;
																	if (lp != lr) {
																		//	equivalence
																		for (k = 0; k <= newlabel; k++)
																			if (C[k] == lr)
																				C[k] = lp;
																	}
																	if (lp != ls) {
																		//	equivalence
																		for (k = 0; k <= newlabel; k++)
																			if (C[k] == ls)
																				C[k] = lp;
																	}
																}
																else
																	if (lr == 0) {
																		L[i * input->width + j] = lp;
																		if (lp != lq) {
																			//	equivalence
																			for (k = 0; k <= newlabel; k++)
																				if (C[k] == lq)
																					C[k] = lp;
																		}
																		if (lp != ls) {
																			//	equivalence
																			for (k = 0; k <= newlabel; k++)
																				if (C[k] == ls)
																					C[k] = lp;
																		}
																	}
																	else
																		if (ls == 0) {
																			L[i * input->width + j] = lp;
																			if (lp != lq) {
																				//	equivalence
																				for (k = 0; k <= newlabel; k++)
																					if (C[k] == lq)
																						C[k] = lp;
																			}
																			if (lp != lr) {
																				//	equivalence
																				for (k = 0; k <= newlabel; k++)
																					if (C[k] == lr)
																						C[k] = lp;
																			}
																		}
																		else {
																			L[i * input->width + j] = lp;
																			if (lp != lq) {
																				//	equivalence
																				for (k = 0; k <= newlabel; k++)
																					if (C[k] == lq)
																						C[k] = lp;
																			}
																			if (lp != lr) {
																				//	equivalence
																				for (k = 0; k <= newlabel; k++)
																					if (C[k] == lr)
																						C[k] = lp;
																			}
																			if (lp != ls) {
																				//	equivalence
																				for (k = 0; k <= newlabel; k++)
																					if (C[k] == ls)
																						C[k] = lp;
																			}
																		}
			}
		}

	length_C = newlabel + 1;

	//	n: number of different sources (the number of positive labels)
	n = 0;

	//	Remapping vector C with consecutive labels
	for (i = 1; i < length_C; i++)
		if (C[i] >= (n + 1)) {
			n++;
			temp_label = C[i];
			for (j = i; j < length_C; j++)
				if (C[j] == temp_label)
					C[j] = n;
		}

	//	Second scan of the image
	for (i = 0; i < input->height; i++)
		for (j = 0; j < input->width; j++)
			if (L[i * input->width + j] != 0)
				L[i * input->width + j] = C[L[i * input->width + j]];

	free(C);

	return n;
}

void visualize_both_centers(IplImage *input, IplImage *output, int *baricenter, int *max_int, int n_sources) {
	int i, j, k, i_center, j_center, datavalue1, datavalue2, datavalue3;
	uchar temp1,temp2,temp3;

	cvCvtColor(input, output, CV_GRAY2BGR);


	for (i = 0; i < n_sources; i++) {
		i_center = baricenter[2 * i];
		j_center = baricenter[2 * i + 1];

		for (j = -1; j <= 1; j++)
			for (k = -1; k <= 1; k++)
				if (((i_center + j) < input->height) && ((i_center + j) >= 0) && ((j_center + k) < input->width) && ((j_center + k) >= 0)) {
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k))] = 0;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 1] = 0;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 2] = 255;
				}
	}

	for (i = 0; i < n_sources; i++) {
		i_center = max_int[2 * i];
		j_center = max_int[2 * i + 1];


		for (j = -1; j <= 1; j++)
			for (k = -1; k <= 1; k++)
				if (((i_center + j) < input->height) && ((i_center + j) >= 0) && ((j_center + k) < input->width) && ((j_center + k) >= 0)) {
					temp1 = ((uchar*)output->imageData)[3 * ((i_center + j) * input->width + (j_center + k))];
					datavalue1 = (int)temp1;
					temp2 = ((uchar*)output->imageData)[3 * ((i_center + j) * input->width + (j_center + k)) + 1];
					datavalue2 = (int)temp2;
					temp3 = ((uchar*)output->imageData)[3 * ((i_center + j) * input->width + (j_center + k)) + 2];
					datavalue3 = (int)temp3;
					if ((datavalue1 == 0)&& (datavalue2 == 0)&& (datavalue3 == 255)) {
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k))] = 0;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 1] = 150;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 2] = 0;
					}
					else {
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k))] = 150;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 1] = 0;
						output->imageData[3 * ((i_center + j) * input->width + (j_center + k)) + 2] = 0;
					}
				}
	}

}

void center_values(int *values, CvSize size, int *L, double *max_real, int * max_int,  int n_sources) {
	int i, j, m, n, k;
	int *max = (int*)malloc(n_sources * sizeof(int));
	int *n_max = (int*)malloc(n_sources * sizeof(int));

	for (i = 0; i < n_sources; i++) {
		max[i] = 0;
		n_max[i] = 0;
		max_real[2 * i] = 0;
		max_real[2 * i + 1] = 0;
	}

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++)
			if ((L[i*size.width + j] > 0) && (values[i*size.width + j] > max[L[i*size.width + j] - 1]))
				max[L[i*size.width + j] - 1] = values[i*size.width + j];

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++)
			if ((L[i*size.width + j] > 0) && (values[i*size.width + j] == max[L[i*size.width + j] - 1])) {
				max_real[2*(L[i*size.width + j] - 1)] = max_real[2 * (L[i*size.width + j] - 1)] + i;
				max_real[2 * (L[i*size.width + j] - 1) + 1] = max_real[2 * (L[i*size.width + j] - 1) + 1] + j;
				n_max[L[i*size.width + j] - 1]++;
			}

	for (i = 0; i < n_sources; i++) {
		max_real[2 * i] = max_real[2 * i] / n_max[i];
		max_real[2 * i + 1] = max_real[2 * i + 1] / n_max[i];
	}


	for(i=0;i<n_sources;i++){
		max_int[2*i] = round(max_real[2*i]);
		max_int[2*i + 1] = round(max_real[2*i + 1]);
	}

	free(max);
	free(n_max);

}

void double_size_colour(IplImage *input, IplImage *output) {
	int i, j, datavalue1, datavalue2, datavalue3;
	uchar temp1, temp2, temp3;

	for(i = 0; i < input->height; i++)
		for (j = 0; j < input->width; j++) {
			temp1 = ((uchar*)input->imageData)[3 * (i * input->width + j)];
			datavalue1 = (int)temp1;
			temp2 = ((uchar*)input->imageData)[3 * (i * input->width + j)+1];
			datavalue2 = (int)temp2;
			temp3 = ((uchar*)input->imageData)[3 * (i * input->width + j)+2];
			datavalue3 = (int)temp3;

			output->imageData[3 * (2 * i * output->width + 2 * j)] = datavalue1;
			output->imageData[3 * (2 * i * output->width + 2 * j) + 1] = datavalue2;
			output->imageData[3 * (2 * i * output->width + 2 * j) + 2] = datavalue3;

			output->imageData[3 * ((2 * i + 1) * output->width + 2 * j)] = datavalue1;
			output->imageData[3 * ((2 * i + 1) * output->width + 2 * j) + 1] = datavalue2;
			output->imageData[3 * ((2 * i + 1) * output->width + 2 * j) + 2] = datavalue3;

			output->imageData[3 * (2 * i * output->width + (2 * j + 1))] = datavalue1;
			output->imageData[3 * (2 * i * output->width + (2 * j + 1)) + 1] = datavalue2;
			output->imageData[3 * (2 * i * output->width + (2 * j + 1)) + 2] = datavalue3;

			output->imageData[3 * ((2 * i + 1) * output->width + (2 * j + 1))] = datavalue1;
			output->imageData[3 * ((2 * i + 1) * output->width + (2 * j + 1)) + 1] = datavalue2;
			output->imageData[3 * ((2 * i + 1) * output->width + (2 * j + 1)) + 2] = datavalue3;
		}
}

int remove_separated_blobs(IplImage *input, IplImage *output, int *L, int *bar_int, double * bar_real, int n_sources, double *ratio, int center_distance) {
	int i, j, m, n, temp;
	long int distance;

	cvCopy(input, output, 0);
	temp = n_sources;
//	printf("N_SOURCES_OLD: %d\n", temp);

	for(i=1; i<temp; i++)
		for (j = i + 1; j <= temp; j++) {

				distance = pow((bar_int[2 * (i - 1)] - bar_int[2 * (j - 1)]), 2) + pow((bar_int[2 * (i - 1)+1] - bar_int[2 * (j - 1)+1]), 2);

				if (distance <= center_distance) {
					if (ratio[i] < ratio[j]) {
						temp--;
						for (m = 0; m < input->height; m++)
							for (n = 0; n < input->width; n++)
								if (L[m*input->width + n] == i) {
									L[m*input->width + n] = 0;
									output->imageData[m*input->width + n] = (uchar)0;
								}
								else
									if (L[m*input->width + n] > i)
										L[m*input->width + n] = L[m*input->width + n] - 1;

						for (m=i-1; m<temp + 1; m++){
							bar_int[2*m] = bar_int[2*(m+1)];
							bar_int[2*m + 1] = bar_int[2*(m+1) + 1];
							bar_real[2*m] = bar_real[2*(m+1)];
							bar_real[2*m + 1] = bar_real[2*(m+1) + 1];
						}

						/*for (m = 2*(i - 1); m < 2*(n_sources-1); m++) {
							bar_int[m] = bar_int[(m + 2)];
							bar_int[m + 1] = bar_int[(m + 2) +1];
							bar_real[m] = bar_real[(m + 2)];
							bar_real[m + 1] = bar_real[(m + 2) +1];
						}*/
					}
					else {
							temp--;
							for (m = 0; m < input->height; m++)
								for (n = 0; n < input->width; n++)
									if (L[m*input->width + n] == j) {
										L[m*input->width + n] = 0;
										output->imageData[m*input->width + n] = (uchar)0;
									}
									else
										if (L[m*input->width + n] > j)
											L[m*input->width + n] = L[m*input->width + n] - 1;

							for (m=j-1; m<temp+1; m++){
								bar_int[2*m] = bar_int[2*(m+1)];
								bar_int[2*m + 1] = bar_int[2*(m+1) + 1];
								bar_real[2*m] = bar_real[2*(m+1)];
								bar_real[2*m + 1] = bar_real[2*(m+1) + 1];
							}

							/*for (m = 2*(j - 1); m < 2*(n_sources-1); m++) {
								bar_int[m] = bar_int[(m + 2)];
								bar_int[m + 1] = bar_int[(m + 2) +1];
								bar_real[m] = bar_real[(m + 2)];
								bar_real[m + 1] = bar_real[(m + 2) +1];
							}*/
					}
				}
			}

	return temp;
}

void order_coordinates(double * array, int n_sources){

	int i, j;
	double temp;

	for(i=0; i<2*n_sources; i=i+2)
		for(j=0;j<2*n_sources-i-2;j=j+2){
			if((array[j] > array[j+2]) || ((array[j] == array[j+2]) && (array[j+1] > array[j+3])) ){
				temp = array[j];
				array[j] = array[j+2];
				array[j+2] = temp;
				temp = array[j+1];
				array[j+1] = array[j+3];
				array[j+3] = temp;
			}
		}

//	for(i=0;i<n_sources;i++){
//		printf("%g\t%g\n", array[2*i], array[2*i+1]);
//	}

}

void mean_coordinates(double * bar, double * max, double * mean, int n_sources){
	int i;

	for(i=0;i<2*n_sources;i++)
		mean[i] = (bar[i] + max[i])/2;

//	for(i=0;i<n_sources;i++)
//			printf("%g\t%g\n", mean[2*i], mean[2*i+1]);

}

void write_center_coordinates(double * center_coordinates, int n_sources, char * filepath){
	FILE * coord_file = fopen( filepath, "w");
	int i;

	fprintf(coord_file, "%d\n", n_sources);
	for(i=0;i<n_sources;i++){
		fprintf(coord_file, "%7.5g\t%7.5g\n", center_coordinates[2*i], center_coordinates[2*i+1] );
	}
	fclose(coord_file);

}

void back_divide(IplImage * input, IplImage* background, IplImage * output){

	int i,j;
	int back, in_val, div_int;
	uchar back_char, in_val_char;
	float division, *array, in_val_float, back_float;

	CvSize size;

	size.height = input->height;
	size.width = input->width;

	array = (float*)malloc(input->height*input->width*sizeof(float));

	for(i=0; i<input->height; i++)
			for(j=0;j<input->width;j++){
				back_char = ((uchar*)background->imageData)[i * background->width + j];
				in_val_char = ((uchar*)input->imageData)[i * input->width + j];
				back = (int)back_char;
				back_float = (float)back_char;
				in_val = (int)in_val_char;
				in_val_float = (float)in_val_char;

				if ( (in_val >= 175) && (back <=80)  ) {
					in_val = 10;
					back = 1;
				}

				if (back == 0)
					back = 1;

				division = (float)in_val/(float)back;
//				div_int = round(division);
				array[i*output->width + j] = division;
			}

	rescale_float(array, output->imageData, size);
	free(array);

}

void rescale_circle(int *input, char* output, CvSize size){
	int i, j, circ;
	long int max, temp;

	max = 0;
	int r = 100;
	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			circ = pow(i-100 , 2) + pow(j-100 , 2);
			if ( (circ <= r*r) && (input[i * size.width + j]>max) )
				max = input[i * size.width + j];
		}

	for (i = 0; i < size.height; i++)
		for (j = 0; j < size.width; j++) {
			circ = pow(i-100 , 2) + pow(j-100 , 2);
			if (circ <= r*r){
				temp = (int)(input[i * size.width + j] * 255 / max);
				output[i * size.width + j] = (uchar)temp;
			}
			else
				output[i * size.width + j] = (uchar)0;
		}

}
