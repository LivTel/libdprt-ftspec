/* dprt.c
** Entry point for Data Pipeline Reduction Routines
** $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/dprt.c,v 0.13 2025-02-10 11:50:06 cjm Exp $
*/
/**
 * dprt.c is the entry point for the Data Reduction Pipeline (Real Time).
 * @author Chris Mottram, LJMU
 * @version $Revision: 0.13 $
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "dprt_jni_general.h"
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
static char rcsid[] = "$Id: dprt.c,v 0.13 2025-02-10 11:50:06 cjm Exp $";

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */

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
 */
int DpRt_Initialise(void)
{
	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	if(!DpRt_JNI_Initialise())
		return FALSE;
	/* call reduction library initialisation routine here. */
	return TRUE;
}

/**
 * This finction should be called when the library/DpRt is about to be shutdown.
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 */
int DpRt_Shutdown(void)
{
	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	/* call reduction library shutdown routine here. */
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
 * @see ngat_dprt_ftspec_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see #Calibrate_Reduce_Fake
 */
int DpRt_Calibrate_Reduce(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts)
{
	float l1mean,l1counts;
	int full_reduction;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	/* check parameters */
	if(input_filename == NULL)
	{
		DpRt_JNI_Error_Number = 100;
		strcpy(DpRt_JNI_Error_String,"DpRt_Calibrate_Reduce:input filename was NULL.");
		return FALSE;
	}
	if(output_filename == NULL)
	{
		DpRt_JNI_Error_Number = 101;
		strcpy(DpRt_JNI_Error_String,"DpRt_Calibrate_Reduce:output filename was NULL.");
		return FALSE;
	}
	if(mean_counts == NULL)
	{
		DpRt_JNI_Error_Number = 102;
		strcpy(DpRt_JNI_Error_String,"DpRt_Calibrate_Reduce:mean counts was NULL.");
		return FALSE;
	}
	if(peak_counts == NULL)
	{
		DpRt_JNI_Error_Number = 103;
		strcpy(DpRt_JNI_Error_String,"DpRt_Calibrate_Reduce:peak counts was NULL.");
		return FALSE;
	}
	/* do full reduction? */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Full Reduction Flag:%d\n",full_reduction);
	/*
	if(full_reduction)
		run_mode = FULL_REDUCTION;
	else
		run_mode = QUICK_REDUCTION;
	*/
	/* initialise return values */
	l1mean = 0.0f;
	l1counts= 0.0f;
	/* call reduction library routine here. */
	/* copy input filename to output - no reduction done */
	(*output_filename) = (char*)malloc((strlen(input_filename)+1)*sizeof(char));
	if((*output_filename) == NULL)
	{
		DpRt_JNI_Error_Number = 104;
		strcpy(DpRt_JNI_Error_String,"DpRt_Calibrate_Reduce:output filename was NULL.");
		return FALSE;
	}
	strcpy((*output_filename),input_filename);
	/* copy other return values */
	(*mean_counts) = (double)l1mean;
	(*peak_counts) = (double)l1counts;
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
 * @see ngat_dprt_ftspec_DpRtLibrary.html
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Get_Property_Boolean
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_Number
 * @see ../../jni_general/cdocs/dprt_jni_general.html#DpRt_JNI_Error_String
 * @see #Expose_Reduce_Fake
 */
int DpRt_Expose_Reduce(char *input_filename,char **output_filename,double *seeing,double *counts,double *x_pix,
		       double *y_pix,double *photometricity,double *sky_brightness,int *saturated)
{
	float l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat,full_reduction;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	/* check parameters */
	if(input_filename == NULL)
	{
		DpRt_JNI_Error_Number = 105;
		strcpy(DpRt_JNI_Error_String,"DpRt_Expose_Reduce:input filename was NULL.");
		return FALSE;
	}
	if(output_filename == NULL)
	{
		DpRt_JNI_Error_Number = 106;
		strcpy(DpRt_JNI_Error_String,"DpRt_Expose_Reduce:output filename was NULL.");
		return FALSE;
	}
	/* get whether to do full reduction */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Full Reduction Flag:%d\n",full_reduction);
	/* initialise return values */
	l1seeing = 0.0f;
	l1counts = 0.0f;
	l1xpix = 0.0f;
	l1ypix = 0.0f;
	l1photom = 0.0f;
	l1skybright = 0.0f;
	l1sat = 0;
	/* call reduction library routine here. */
	/* copy input filename to output - no reduction done */
	(*output_filename) = (char*)malloc((strlen(input_filename)+1)*sizeof(char));
	if((*output_filename) == NULL)
	{
		DpRt_JNI_Error_Number = 107;
		strcpy(DpRt_JNI_Error_String,"DpRt_Expose_Reduce:output filename was NULL.");
		return FALSE;
	}
	strcpy((*output_filename),input_filename);
	/* copy return values to function return values */
	(*seeing) = (double)l1seeing;
	(*counts) = (double)l1counts;
	(*x_pix) = (double)l1xpix;
	(*y_pix) = (double)l1ypix;
	(*photometricity) = (double)l1photom;
	(*sky_brightness) = (double)l1skybright;
	(*saturated) = (int)l1sat;
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
 */
int DpRt_Make_Master_Bias(char *directory_name)
{
	int make_master_bias;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	/* whether to do the make master bias or not */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.make_master_bias",&make_master_bias))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias Flag:%d\n",make_master_bias);
	if(make_master_bias)
	{
		fprintf(stdout,"DpRt_Make_Master_Bias:Calling Make Master Bias routine.\n");
		/* call make master bias routine here */
	}
	else
	{
		fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias Flag was FALSE:Not making master bias.\n");
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
 */
int DpRt_Make_Master_Flat(char *directory_name)
{
	int make_master_flat;

	DpRt_JNI_Error_Number = 0;
	DpRt_JNI_Error_String[0] = '\0';
	/* should we call make master flat or not */
	if(!DpRt_JNI_Get_Property_Boolean("dprt.make_master_flat",&make_master_flat))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat Flag:%d\n",make_master_flat);
	if(make_master_flat)
	{
		fprintf(stdout,"DpRt_Make_Master_Flat:Calling Make Master Flat routine.\n");
		/* call make master flat routine here */
	}
	else
	{
		fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat Flag was FALSE:Not making master flat.\n");
	}
	return TRUE;
}

/* ------------------------------------------------------- */
/* internal functions */
/* ------------------------------------------------------- */

/*
** $Log: not supported by cvs2svn $
*/
