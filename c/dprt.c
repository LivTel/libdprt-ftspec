/* dprt.c
** Entry point for Data Pipeline Reduction Routines
** $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/dprt.c,v 0.12 2004-03-31 16:44:25 cjm Exp $
*/
/**
 * dprt.c is the entry point for the Data Reduction Pipeline (Real Time).
 * @author Chris Mottram, LJMU
 * @version $Revision: 0.12 $
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "dprt_jni_general.h"
#include "rjs_dprt.h"
#include "dprt.h"

/* ------------------------------------------------------- */
/* hash definitions */
/* ------------------------------------------------------- */
/**
 * This program only accepts FITS files with the bits per pixel of this value.
 */
#define FITS_GET_DATA_BITPIX		(16)
/**
 * This program only accepts FITS files with this number of axes.
 */
#define FITS_GET_DATA_NAXIS		(2)

/* ------------------------------------------------------- */
/* internal variables */
/* ------------------------------------------------------- */
/**
 * Revision Control System identifier.
 */
static char rcsid[] = "$Id: dprt.c,v 0.12 2004-03-31 16:44:25 cjm Exp $";

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */
static int Calibrate_Reduce_Fake(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts);
static int Expose_Reduce_Fake(char *input_filename,char **output_filename,double *seeing,double *counts,
	double *x_pix,double *y_pix,double *photometricity,double *sky_brightness,int *saturated);

/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
/**
 * This finction should be called when the library/DpRt is first initialised/loaded.
 * It allows the C layer to perform initial initialisation.
 * The function pointers to use a C routine to load the property from the config file are initialised.
 * Note these function pointers will be over-written by the functions in DpRtLibrary.c if this
 * initialise routine was called from the Java (JNI) layer.
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_General_Initialise
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_init
 */
int DpRt_Initialise(void)
{
	int retval,fake;


	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	if(!DpRt_JNI_Initialise())
		return FALSE;
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Initialise:Fake:%d\n",fake);
	if(fake == FALSE)
	{
		/* call real initialisation routine */
		fprintf(stdout,"Calling DpRt initialisation routine (dprt_init).\n");
		retval = dprt_init();
		fprintf(stdout,"DpRt initialisation routine (dprt_init) returned %d.\n",retval);
		if(retval == TRUE)
		{
			DpRt_JNI_Error_Number = dprt_err_int;
			strcpy(DpRt_JNI_Error_String,dprt_err_str);
			return FALSE;
		}
	}
	return TRUE;
}

/**
 * This finction should be called when the library/DpRt is about to be shutdown.
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_close_down
 */
int DpRt_Shutdown(void)
{
	int retval,fake;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Shutdown:Fake:%d\n",fake);
	if(fake == FALSE)
	{
		fprintf(stdout,"Calling DpRt shutdown routine (dprt_close_down).\n");
		retval = dprt_close_down();
		fprintf(stdout,"DpRt shutdown routine (dprt_lose_down) returned %d.\n",retval);
		if(retval != TRUE)
		{
			DpRt_JNI_Error_Number = dprt_err_int;
			strcpy(DpRt_JNI_Error_String,dprt_err_str);
			return FALSE;
		}
	}
	return TRUE;
}

/**
 * This routine does the real time data reduction pipeline on a calibration file. It is usually invoked from the
 * Java DpRtCalibrateReduce call in DpRtLibrary.java. If the DpRt_JNI_Get_Abort
 * routine returns TRUE during the execution of the pipeline the pipeline should abort it's
 * current operation and return FALSE.
 * @param input_filename The FITS filename to be processed.
 * @param output_filename The resultant filename should be put in this variable. This variable is the
 *       address of a pointer to a sequence of characters, hence it should be referenced using
 *       <code>(*output_filename)</code> in this routine.
 * @param meanCounts The address of a double to store the mean counts calculated by this routine.
 * @param peakCounts The address of a double to store the peak counts calculated by this routine.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ngat_dprt_ccs_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see #Calibrate_Reduce_Fake
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_process
 */
int DpRt_Calibrate_Reduce(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts)
{
	int fake,retval;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat,run_mode,full_reduction;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Fake:%d\n",fake);
	if(!DpRt_JNI_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Full Reduction Flag:%d\n",full_reduction);
	if(fake)
	{
		return Calibrate_Reduce_Fake(input_filename,output_filename,mean_counts,peak_counts);
	}
	else
	{
		if(full_reduction)
			run_mode = FULL_REDUCTION;
		else
			run_mode = QUICK_REDUCTION;
		fprintf(stdout,"DpRt_Calibrate_Reduce:Calling Calibration reduction routine (dprt_process(%d)).\n",
			run_mode);
		retval = dprt_process(input_filename,run_mode,output_filename,&l1mean,&l1seeing, 
			&l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
		fprintf(stdout,"DpRt_Calibrate_Reduce:Calibration reduction routine (dprt_process) returned %d.\n",
			retval);
		if(retval == TRUE) /* an error has occured */
		{
			DpRt_JNI_Error_Number = dprt_err_int;
			strcpy(DpRt_JNI_Error_String,dprt_err_str);
			(*output_filename) = NULL;
			(*mean_counts) = 0;
			(*peak_counts) = 0;
			return FALSE;
		}
		(*mean_counts) = (double)l1mean;
		(*peak_counts) = (double)l1counts;
	}
	return TRUE;
}

/**
 * This routine does the real time data reduction pipeline on an expose file. It is usually invoked from the
 * Java DpRtExposeReduce call in DpRtLibrary.java. If the <a href="#DpRt_Get_Abort">DpRt_Get_Abort</a>
 * routine returns TRUE during the execution of the pipeline the pipeline should abort it's
 * current operation and return FALSE.
 * @param input_filename The FITS filename to be processed.
 * @param output_filename The resultant filename should be put in this variable. This variable is the
 *       address of a pointer to a sequence of characters, hence it should be referenced using
 *       <code>(*output_filename)</code> in this routine.
 * @param seeing The address of a double to store the seeing calculated by this routine.
 * @param counts The address of a double to store the counts of th brightest pixel calculated by this
 *       routine.
 * @param x_pix The x pixel position of the brightest object in the field. Note this is an average pixel
 *       number that may not be a whole number of pixels.
 * @param y_pix The y pixel position of the brightest object in the field. Note this is an average pixel
 *       number that may not be a whole number of pixels.
 * @param photometricity In units of magnitudes of extinction. This is only filled in for standard field
 * 	reductions.
 * @param sky_brightness In units of magnitudes per arcsec&#178;. This is an estimate of sky brightness.
 * @param saturated This is a boolean, returning TRUE if the object is saturated.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ngat_dprt_ccs_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see #Expose_Reduce_Fake
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_process
 */
int DpRt_Expose_Reduce(char *input_filename,char **output_filename,double *seeing,double *counts,double *x_pix,
		       double *y_pix,double *photometricity,double *sky_brightness,int *saturated)
{
	int fake,retval;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat,run_mode,full_reduction;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Expose_Reduce:Fake:%d\n",fake);
	if(!DpRt_JNI_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Full Reduction Flag:%d\n",full_reduction);
	if(fake)
	{
		return Expose_Reduce_Fake(input_filename,output_filename,seeing,counts,x_pix,y_pix,
			photometricity,sky_brightness,saturated);
	}
	else
	{
		if(full_reduction)
			run_mode = FULL_REDUCTION;
		else
			run_mode = QUICK_REDUCTION;
		fprintf(stdout,"DpRt_Expose_Reduce:Calling Exposure reduction routine (dprt_process(%d)).\n",run_mode);
		retval = dprt_process(input_filename,run_mode,output_filename,&l1mean,&l1seeing, 
			&l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
		fprintf(stdout,"DpRt_Expose_Reduce:Exposure reduction routine (dprt_process) returned %d.\n",retval);
		if(retval == TRUE)
		{
			DpRt_JNI_Error_Number = dprt_err_int;
			strcpy(DpRt_JNI_Error_String,dprt_err_str);
			(*output_filename) = NULL;
			(*seeing) = 0.0;
			(*counts) = 0.0;
			(*x_pix) = 0.0;
			(*y_pix) = 0.0;
			(*photometricity) = 0.0;
			(*sky_brightness) = 0.0;
			(*saturated) = FALSE;
			return FALSE;
		}
		(*seeing) = (double)l1seeing;
		(*counts) = (double)l1counts;
		(*x_pix) = (double)l1xpix;
		(*y_pix) = (double)l1ypix;
		(*photometricity) = (double)l1photom;
		(*sky_brightness) = (double)l1skybright;
		(*saturated) = (int)l1sat;
	}
	return TRUE;
}

/**
 * This routine creates a master bias frame for each binning factor, created from biases in the specified
 * directory
 * @param directory_name A directory containing the  FITS filenames to be processed.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_process
 * @see ../../rjs/cdocs/rjs_dprt.html#MAKE_BIAS
 */
int DpRt_Make_Master_Bias(char *directory_name)
{
	int fake,retval,make_master_bias;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Bias:Fake:%d\n",fake);
	if(!DpRt_JNI_Get_Property_Boolean("dprt.make_master_bias",&make_master_bias))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias Flag:%d\n",make_master_bias);
	if(fake)
	{
		/* do nothing to fake this */
		return TRUE;
	}
	else
	{
		if(make_master_bias)
		{
			fprintf(stdout,"DpRt_Make_Master_Bias:Calling Make Master Bias routine (dprt_process).\n");
			retval = dprt_process(directory_name,MAKE_BIAS,NULL,&l1mean,&l1seeing, 
					      &l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
			fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias routine (dprt_process) returned %d.\n",
				retval);
			if(retval == TRUE)
			{
				DpRt_JNI_Error_Number = dprt_err_int;
				strcpy(DpRt_JNI_Error_String,dprt_err_str);
				return FALSE;
			}
		}
		else
		{
			fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias Flag was FALSE:"
				"Not making master bias.\n");
		}
	}
	return TRUE;
}

/**
 * This routine creates a master flat frame for each binning factor, created from flats in the specified
 * directory.
 * @param directory_name A directory containing the  FITS filenames to be processed.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see ../../rjs/cdocs/rjs_dprt.html#dprt_process
 * @see ../../rjs/cdocs/rjs_dprt.html#MAKE_FLAT
 */
int DpRt_Make_Master_Flat(char *directory_name)
{
	int fake,retval,make_master_flat;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
/* are we doing a fake reduction or a real one. */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Flat:Fake:%d\n",fake);
	if(!DpRt_JNI_Get_Property_Boolean("dprt.make_master_flat",&make_master_flat))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat Flag:%d\n",make_master_flat);
	if(fake)
	{
		/* do nothing to fake this */
		return TRUE;
	}
	else
	{
		if(make_master_flat)
		{
			fprintf(stdout,"DpRt_Make_Master_Flat:Calling Make Master Flat routine (dprt_process).\n");
			retval = dprt_process(directory_name,MAKE_FLAT,NULL,&l1mean,&l1seeing, 
					      &l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
			fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat routine (dprt_process) returned %d.\n",
				retval);
			if(retval == TRUE)
			{
				DpRt_JNI_Error_Number = dprt_err_int;
				strcpy(DpRt_JNI_Error_String,dprt_err_str);
				return FALSE;
			}
		}
		else
		{
			fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat Flag was FALSE:"
				"Not making master flat.\n");
		}
	}
	return TRUE;
}

/* ------------------------------------------------------- */
/* internal functions */
/* ------------------------------------------------------- */
/**
 * This routine does a fake real time data reduction pipeline on a calibration file. It is invoked from the
 * DpRt_Calibrate_Reduce routine.If the <a href="#DpRt_Get_Abort">DpRt_Get_Abort</a>
 * routine returns TRUE during the execution of the pipeline the pipeline should abort it's
 * current operation and return FALSE.
 * @param input_filename The FITS filename to be processed.
 * @param output_filename The resultant filename should be put in this variable. This variable is the
 *       address of a pointer to a sequence of characters, hence it should be referenced using
 *       <code>(*output_filename)</code> in this routine.
 * @param meanCounts The address of a double to store the mean counts calculated by this routine.
 * @param peakCounts The address of a double to store the peak counts calculated by this routine.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ngat_dprt_ccs_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Set_Abort
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Abort
 * @see #DpRt_Calibrate_Reduce
 */
static int Calibrate_Reduce_Fake(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts)
{
	fitsfile *fp = NULL;
	int retval=0,status=0,integer_value,naxis_one,naxis_two,i,j,value;
	unsigned short *data = NULL;
	double counts_count = 0.0;
	int max_count = 0;

/* set the error stuff to no error*/
	DpRt_JNI_Error_Number = 0;
	strcpy(DpRt_JNI_Error_String,"");
/* setup return values */
	(*mean_counts) = 0.0;
	(*peak_counts) = 0.0;
/* unset any previous aborts - ready to start processing */
	DpRt_JNI_Set_Abort(FALSE);
/* do processing  here */
	fprintf(stderr,"Calibrate_Reduce_Fake(%s).\n",input_filename);
/* open file */
	retval = fits_open_file(&fp,input_filename,READONLY,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 23;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Open failed.\n",input_filename);
		return FALSE;
	}
/* check bitpix */
	retval = fits_read_key(fp,TINT,"BITPIX",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 24;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get BITPIX.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_BITPIX)
	{
		DpRt_JNI_Error_Number = 25;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Wrong BITPIX value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* check naxis */
	retval = fits_read_key(fp,TINT,"NAXIS",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 26;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_NAXIS)
	{
		DpRt_JNI_Error_Number = 27;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Wrong NAXIS value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* get naxis1,naxis2 */
	retval = fits_read_key(fp,TINT,"NAXIS1",&naxis_one,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 28;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS1.\n",input_filename);
		return FALSE;
	}
	retval = fits_read_key(fp,TINT,"NAXIS2",&naxis_two,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 29;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS2.\n",input_filename);
		return FALSE;
	}
/* allocate data */
	data = (unsigned short *)malloc(naxis_one*naxis_two*sizeof(unsigned short));
	if(data == NULL)
	{
		DpRt_JNI_Error_Number = 30;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to allocate memory (%d,%d,%d).\n",
			input_filename,naxis_one,naxis_two,naxis_one*naxis_two*sizeof(short));
		return FALSE;
	}
/* read the data */
	retval = fits_read_img(fp,TUSHORT,1,naxis_one*naxis_two,NULL,data,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 31;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to read image(%d,%d).\n",
			input_filename,naxis_one,naxis_two);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* close file */
	retval = fits_close_file(fp,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 32;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Failed to close file.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* during processing regularily check the abort flag as below */
	if(DpRt_JNI_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*mean_counts) = 0.0;
		(*peak_counts) = 0.0;
		(*output_filename) = NULL;
		DpRt_JNI_Error_Number = 1;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* setup return values */
	(*mean_counts) = 0.0;
	(*peak_counts) = 0.0;
	counts_count = 0.0;
	max_count = 0;
	for(j=0;j<naxis_two;j++)
	{
		for(i=0;i<naxis_one;i++)
		{
			value = (int)(data[(naxis_one*j)+i]);
			counts_count += (double)value;
			if(value>max_count)
				max_count = value;
		}
		/* during processing regularily check the abort flag as below */
		if(DpRt_JNI_Get_Abort())
		{
			/* tidy up anything that needs tidying as a result of this routine here */
			(*mean_counts) = 0.0;
			(*peak_counts) = 0.0;
			(*output_filename) = NULL;
			DpRt_JNI_Error_Number = 45;
			sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
			if(data != NULL)
				free(data);
			return FALSE;
		}
	}
	if(data != NULL)
		free(data);
	(*mean_counts) = (float)(((double)counts_count)/((double)(naxis_one*naxis_two)));
	(*peak_counts) = (float)max_count;
/* setup filename - allocate space for string */
	(*output_filename) = (char*)malloc((strlen(input_filename)+1)*sizeof(char));
	/* if malloc fails it returns NULL - this is an error */
	if((*output_filename) == NULL)
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*mean_counts) = 0.0;
		(*peak_counts) = 0.0;
		(*output_filename) = NULL;
		DpRt_JNI_Error_Number = 2;
		sprintf(DpRt_JNI_Error_String,"Calibrate_Reduce_Fake(%s): Memory Allocation Error.\n",input_filename);
		return FALSE;
	}
/* set the filename to something more sensible here */
	strcpy((*output_filename),input_filename);
	return TRUE;
}

/**
 * This routine does the fake data reduction pipeline on an expose file. It is usually invoked from the
 * DpRt_Expose_Reduce routine. If the DpRt_Get_Abort routine returns TRUE during the execution of the pipeline 
 * the pipeline should abort it's current operation and return FALSE.
 * @param input_filename The FITS filename to be processed.
 * @param output_filename The resultant filename should be put in this variable. This variable is the
 *       address of a pointer to a sequence of characters, hence it should be referenced using
 *       <code>(*output_filename)</code> in this routine.
 * @param seeing The address of a double to store the seeing calculated by this routine.
 * @param counts The address of a double to store the counts of th brightest pixel calculated by this
 *       routine.
 * @param x_pix The x pixel position of the brightest object in the field. Note this is an average pixel
 *       number that may not be a whole number of pixels.
 * @param y_pix The y pixel position of the brightest object in the field. Note this is an average pixel
 *       number that may not be a whole number of pixels.
 * @param photometricity In units of magnitudes of extinction. This is only filled in for standard field
 * 	reductions.
 * @param sky_brightness In units of magnitudes per arcsec&#178;. This is an estimate of sky brightness.
 * @param saturated This is a boolean, returning TRUE if the object is saturated.
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see ngat_dprt_ccs_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Double
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Set_Abort
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Abort
 */
static int Expose_Reduce_Fake(char *input_filename,char **output_filename,double *seeing,double *counts,
	double *x_pix,double *y_pix,double *photometricity,double *sky_brightness,int *saturated)
{
	fitsfile *fp = NULL;
	int retval=0,status=0,integer_value,naxis_one,naxis_two,i,j,value;
	unsigned short *data = NULL;
	double telfocus,best_focus,fwhm_per_mm,atmospheric_seeing,atmospheric_variation,error;
	char *ch = NULL;

	/* set the error stuff to no error*/
	DpRt_JNI_Error_Number = 0;
	strcpy(DpRt_JNI_Error_String,"");

/* unset any previous aborts - ready to start processing */
	DpRt_JNI_Set_Abort(FALSE);
/* do processing  here */
	fprintf(stderr,"Expose_Reduce_Fake(%s).\n",input_filename);
/* setup return values */
	(*output_filename) = NULL;
	(*seeing) = 0.0;
	(*counts) = 0.0;
	(*x_pix) = 0.0;
	(*y_pix) = 0.0;
	(*photometricity) = 0.0;
	(*sky_brightness) = 0.0;
	(*saturated) = FALSE;
/* get parameters from config */
	if(!DpRt_JNI_Get_Property_Double("dprt.telfocus.best_focus",&best_focus))
		return FALSE;
	if(!DpRt_JNI_Get_Property_Double("dprt.telfocus.fwhm_per_mm",&fwhm_per_mm))
		return FALSE;
	if(!DpRt_JNI_Get_Property_Double("dprt.telfocus.atmospheric_seeing",&atmospheric_seeing))
		return FALSE;
	if(!DpRt_JNI_Get_Property_Double("dprt.telfocus.atmospheric_variation",&atmospheric_variation))
		return FALSE;
/* open file */
	retval = fits_open_file(&fp,input_filename,READONLY,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 33;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Open failed.\n",input_filename);
		return FALSE;
	}
/* check bitpix */
	retval = fits_read_key(fp,TINT,"BITPIX",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 34;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to get BITPIX.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_BITPIX)
	{
		DpRt_JNI_Error_Number = 35;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Wrong BITPIX value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* check naxis */
	retval = fits_read_key(fp,TINT,"NAXIS",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 36;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_NAXIS)
	{
		DpRt_JNI_Error_Number = 37;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Wrong NAXIS value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* get naxis1,naxis2 */
	retval = fits_read_key(fp,TINT,"NAXIS1",&naxis_one,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 38;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS1.\n",input_filename);
		return FALSE;
	}
	retval = fits_read_key(fp,TINT,"NAXIS2",&naxis_two,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 39;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS2.\n",input_filename);
		return FALSE;
	}
/* get telescope focus */
	retval = fits_read_key(fp,TDOUBLE,"TELFOCUS",&telfocus,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 40;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to get TELFOCUS.\n",input_filename);
		return FALSE;
	}
/* allocate data */
	data = (unsigned short *)malloc(naxis_one*naxis_two*sizeof(unsigned short));
	if(data == NULL)
	{
		DpRt_JNI_Error_Number = 41;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to allocate memory (%d,%d,%d).\n",
			input_filename,naxis_one,naxis_two,naxis_one*naxis_two*sizeof(short));
		return FALSE;
	}
/* read the data */
	retval = fits_read_img(fp,TUSHORT,1,naxis_one*naxis_two,NULL,data,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 42;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to read image(%d,%d).\n",
			input_filename,naxis_one,naxis_two);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* close file */
	retval = fits_close_file(fp,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_JNI_Error_Number = 43;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Failed to close file.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* during processing regularily check the abort flag as below */
	if(DpRt_JNI_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*output_filename) = NULL;
		DpRt_JNI_Error_Number = 44;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* get counts,x_pix,y_pix */
	for(j=0;j<naxis_two;j++)
	{
		for(i=0;i<naxis_one;i++)
		{
			value = (int)(data[(naxis_one*j)+i]);
			if(value> (*counts))
			{
				(*counts) = value;
				(*x_pix) = i;
				(*y_pix) = j;
			}
		}
	}
	if(data != NULL)
		free(data);
/* during processing regularily check the abort flag as below */
	if(DpRt_JNI_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*output_filename) = NULL;
		DpRt_JNI_Error_Number = 3;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
		return FALSE;
	}

	/* setup return values */
	ch = strstr(input_filename,"telFocus");
	if(ch != NULL)
	{
		error = (atmospheric_variation*((double)rand()))/((double)RAND_MAX);
		(*seeing) = (pow((telfocus-best_focus),2.0)*(fwhm_per_mm-atmospheric_seeing))+
				atmospheric_seeing+error;
		fprintf(stderr,"Expose_Reduce_Fake:telfocus %.2f:seeing set to %.2f.\n",telfocus,(*seeing));
	}
	else
	{
		(*seeing) = ((float)(rand()%50))/10.0;
	}
/* setup filename - allocate space for string */
	(*output_filename) = (char*)malloc((strlen(input_filename)+1)*sizeof(char));
/* if malloc fails it returns NULL - this is an error */
	if((*output_filename) == NULL)
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*output_filename) = NULL;
		DpRt_JNI_Error_Number = 4;
		sprintf(DpRt_JNI_Error_String,"Expose_Reduce_Fake(%s): Memory Allocation Error.\n",input_filename);
		return FALSE;
	}
/* set the filename to something more sensible here */
	strcpy((*output_filename),input_filename);
	return TRUE;
}

/*
** $Log: not supported by cvs2svn $
** Revision 0.11  2004/02/06 17:24:12  cjm
** Initialised DpRt_Error_String to blank at the start of each externally callable routine.
**
** Revision 0.10  2004/02/06 17:08:25  cjm
** Reset error number at start of routines.
**
** Revision 0.9  2004/01/15 16:04:13  cjm
** No change.
**
** Revision 0.8  2002/11/26 17:38:13  cjm
** Integrated librjs_dprt.
**
** Revision 0.7  2002/05/20 11:02:09  cjm
** Added photometricity, sky_brightness and saturated to DpRt_Expose_Reduce.
**
** Revision 0.6  2002/05/20 10:44:17  cjm
** Added property interaction routines.
**
** Revision 0.5  2001/05/17 10:04:35  cjm
** Added <math.h> include.
**
** Revision 0.4  2001/05/15 15:34:11  cjm
** telFocus simulation code added.
**
** Revision 0.3  1999/08/24 16:04:30  cjm
** Changed DpRt_Abort comment so CDoc could parse it correctly.
**
** Revision 0.2  1999/06/30 15:07:47  dev
** changes to return when error occurs in reduction code
** set return paramaters to 0/NULL.
**
** Revision 0.1  1999/06/24 11:06:42  dev
** initial revision
**
*/
