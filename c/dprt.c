/* dprt.c
** Entry point for Data Pipeline Reduction Routines
** $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/dprt.c,v 0.8 2002-11-26 17:38:13 cjm Exp $
*/
/**
 * dprt.c is the entry point for the Data Reduction Pipeline (Real Time).
 * @author Chris Mottram, LJMU
 * @version $Revision: 0.8 $
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "dprt.h"
#include "rjs_dprt.h"

/* ------------------------------------------------------- */
/* hash definitions */
/* ------------------------------------------------------- */
/**
 * Property file name. This default value copied from the DpRtStatus.java source.
 */
#define PROPERTY_FILE_NAME "./dprt.properties"
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
static char rcsid[] = "$Id: dprt.c,v 0.8 2002-11-26 17:38:13 cjm Exp $";
/**
 * Internal Error Number - set this to a unique value for each location an error occurs.
 */
static int DpRt_Error_Number = 0;
/**
 * Internal Error String - set this to a descriptive string each place an error occurs.
 * Ensure the string is not longer than <a href="#DPRT_ERROR_STRING_LENGTH">DPRT_ERROR_STRING_LENGTH</a> long.
 * @see #DpRt_Error_Number
 * @see #DPRT_ERROR_STRING_LENGTH
 */
static char DpRt_Error_String[DPRT_ERROR_STRING_LENGTH] = "";

/**
 * Data type holding local data to dprt.c. This consists of the following:
 * <dl>
 * <dt>DpRt_Abort</dt><dd>Boolean value holding whether we should abort the current Data Reduction Pipeline. 
 * 	The reduction
 * 	routines should check this at regular intervals, and if it is TRUE
 * 	should abort processing and set the <a href="#DpRt_Error_Number">DpRt_Error_Number</a> and
 * 	<a href="#DpRt_Error_String">DpRt_Error_String</a> accordingly and return FALSE.</dd>
 * <dt>DpRt_Get_Property_Function_Pointer</dt><dd>Function pointer to actual routine that 
 * 	retrieves a keyword's string value from a Java property list Hashtable.</dd>
 * <dt>DpRt_Get_Property_Integer_Function_Pointer</dt><dd>Function pointer to actual routine that 
 * 	retrieves a keyword's integer value from a Java property list Hashtable.</dd>
 * <dt>DpRt_Get_Property_Double_Function_Pointer</dt><dd>Function pointer to actual routine that 
 * 	retrieves a keyword's double value from a Java property list Hashtable.</dd>
 * <dt>DpRt_Get_Property_Boolean_Function_Pointer</dt><dd>Function pointer to actual routine that retrieves 
 * 	a keyword's boolean value from a Java property list Hashtable.</dd>
 * <dt></dt><dd></dd>
 * </dl>
 * @see #DpRt_Calibrate_Reduce
 * @see #DpRt_Expose_Reduce
 * @see #DpRt_Error_Number
 * @see #DpRt_Error_String
 * @see #DpRt_Set_Abort
 * @see #DpRt_Get_Property
 * @see #DpRt_Get_Property_Integer
 * @see #DpRt_Get_Property_Double
 * @see #DpRt_Get_Property_Boolean
 */
struct DpRt_Struct
{
	volatile int DpRt_Abort;/* volatile as thread dependant */
	int (*DpRt_Get_Property_Function_Pointer)(char *keyword,char **value_string);
	int (*DpRt_Get_Property_Integer_Function_Pointer)(char *keyword,int *value);
	int (*DpRt_Get_Property_Double_Function_Pointer)(char *keyword,double *value);
	int (*DpRt_Get_Property_Boolean_Function_Pointer)(char *keyword,int *value);
};

/**
 * The single instance of struct DpRt_Struct, that holds local data to this source file.
 * Initialised to NULL/FALSE.
 * @see #DpRt_Struct
 */
static struct DpRt_Struct DpRt_Data = 
{
	FALSE,NULL,NULL,NULL,NULL
};

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */
static int Calibrate_Reduce_Fake(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts);
static int Expose_Reduce_Fake(char *input_filename,char **output_filename,double *seeing,double *counts,
	double *x_pix,double *y_pix,double *photometricity,double *sky_brightness,int *saturated);
static int DpRt_Get_Abort(void);
static int DpRt_Get_Property_From_C_File(char *keyword,char **value_string);
static int DpRt_Get_Property_Integer_From_C_File(char *keyword,int *value);
static int DpRt_Get_Property_Double_From_C_File(char *keyword,double *value);
static int DpRt_Get_Property_Boolean_From_C_File(char *keyword,int *value);

/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
/**
 * This finction should be called when the library/DpRt is first initialised/loaded.
 * It allows the C layer to perform initial initialisation.
 * The function pointers to use a C routine to load the property from the config file are initialised.
 * Note these function pointers will be over-written by the functions in DpRtLibrary.c if this
 * initialise routine was called from the Java (JNI) layer.
 * @see #DpRt_Set_Property_Function_Pointer
 * @see #DpRt_Set_Property_Integer_Function_Pointer
 * @see #DpRt_Set_Property_Double_Function_Pointer
 * @see #DpRt_Set_Property_Boolean_Function_Pointer
 * @see #DpRt_Get_Property_From_C_File
 * @see #DpRt_Get_Property_Integer_From_C_File
 * @see #DpRt_Get_Property_Double_From_C_File
 * @see #DpRt_Get_Property_Boolean_From_C_File
 * @see ../rjs/rjs_dprt.html#dprt_init
 */
int DpRt_Initialise(void)
{
	int retval,fake;

	if(DpRt_Data.DpRt_Get_Property_Function_Pointer == NULL)
		DpRt_Set_Property_Function_Pointer(DpRt_Get_Property_From_C_File);
	if(DpRt_Data.DpRt_Get_Property_Integer_Function_Pointer == NULL)
		DpRt_Set_Property_Integer_Function_Pointer(DpRt_Get_Property_Integer_From_C_File);
	if(DpRt_Data.DpRt_Get_Property_Double_Function_Pointer == NULL)
		DpRt_Set_Property_Double_Function_Pointer(DpRt_Get_Property_Double_From_C_File);
	if(DpRt_Data.DpRt_Get_Property_Boolean_Function_Pointer == NULL)
		DpRt_Set_Property_Boolean_Function_Pointer(DpRt_Get_Property_Boolean_From_C_File);
/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
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
			DpRt_Error_Number = dprt_err_int;
			strcpy(DpRt_Error_String,dprt_err_str);
			return FALSE;
		}
	}
	return TRUE;
}

/**
 * This finction should be called when the library/DpRt is about to be shutdown.
 * @see ../rjs/rjs_dprt.html#dprt_init
 */
int DpRt_Shutdown(void)
{
	int retval,fake;

/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Shutdown:Fake:%d\n",fake);
	if(fake == FALSE)
	{
		fprintf(stdout,"Calling DpRt shutdown routine (dprt_close_down).\n");
		retval = dprt_close_down();
		fprintf(stdout,"DpRt shutdown routine (dprt_lose_down) returned %d.\n",retval);
		if(retval != TRUE)
		{
			DpRt_Error_Number = dprt_err_int;
			strcpy(DpRt_Error_String,dprt_err_str);
			return FALSE;
		}
	}
	return TRUE;
}

/**
 * This routine does the real time data reduction pipeline on a calibration file. It is usually invoked from the
 * Java DpRtCalibrateReduce call in DpRtLibrary.java. If the <a href="#DpRt_Get_Abort">DpRt_Get_Abort</a>
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
 * @see DpRtLibrary.html
 * @see #DpRt_Get_Abort
 * @see #Calibrate_Reduce_Fake
 * @see ../rjs/rjs_dprt.html#dprt_process
 */
int DpRt_Calibrate_Reduce(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts)
{
	int fake,retval;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat,run_mode,full_reduction;

/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Calibrate_Reduce:Fake:%d\n",fake);
	if(!DpRt_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
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
		retval = dprt_process(input_filename,run_mode,&l1mean,&l1seeing, 
			&l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
		fprintf(stdout,"DpRt_Calibrate_Reduce:Calibration reduction routine (dprt_process) returned %d.\n",
			retval);
		if(retval == TRUE) /* an error has occured */
		{
			DpRt_Error_Number = dprt_err_int;
			strcpy(DpRt_Error_String,dprt_err_str);
			(*output_filename) = NULL;
			(*mean_counts) = 0;
			(*peak_counts) = 0;
			return FALSE;
		}
		(*output_filename) = NULL; /* diddly fix */
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
 * @see DpRtLibrary.html
 * @see #DpRt_Get_Abort
 * @see #Expose_Reduce_Fake
 * @see ../rjs/rjs_dprt.html#dprt_process
 */
int DpRt_Expose_Reduce(char *input_filename,char **output_filename,double *seeing,double *counts,double *x_pix,
		       double *y_pix,double *photometricity,double *sky_brightness,int *saturated)
{
	int fake,retval;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat,run_mode,full_reduction;

/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Expose_Reduce:Fake:%d\n",fake);
	if(!DpRt_Get_Property_Boolean("dprt.full_reduction",&full_reduction))
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
		retval = dprt_process(input_filename,run_mode,&l1mean,&l1seeing, 
			&l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
		fprintf(stdout,"DpRt_Expose_Reduce:Exposure reduction routine (dprt_process) returned %d.\n",retval);
		if(retval == TRUE)
		{
			DpRt_Error_Number = dprt_err_int;
			strcpy(DpRt_Error_String,dprt_err_str);
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
		(*output_filename) = NULL; /* diddly fix */
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
 * @see #DpRt_Get_Property_Boolean
 * @see ../rjs/rjs_dprt.html#dprt_process
 * @see ../rjs/rjs_dprt.html#MAKE_BIAS
 */
int DpRt_Make_Master_Bias(char *directory_name)
{
	int fake,retval,make_master_bias;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat;

/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Bias:Fake:%d\n",fake);
	if(!DpRt_Get_Property_Boolean("dprt.make_master_bias",&make_master_bias))
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
			retval = dprt_process(directory_name,MAKE_BIAS,&l1mean,&l1seeing, 
					      &l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
			fprintf(stdout,"DpRt_Make_Master_Bias:Make Master Bias routine (dprt_process) returned %d.\n",
				retval);
			if(retval == TRUE)
			{
				DpRt_Error_Number = dprt_err_int;
				strcpy(DpRt_Error_String,dprt_err_str);
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
 * @see #DpRt_Get_Property_Boolean
 * @see ../rjs/rjs_dprt.html#dprt_process
 * @see ../rjs/rjs_dprt.html#MAKE_FLAT
 */
int DpRt_Make_Master_Flat(char *directory_name)
{
	int fake,retval,make_master_flat;
	float l1mean,l1seeing,l1xpix,l1ypix,l1counts,l1photom,l1skybright;
	int l1sat;

/* are we doing a fake reduction or a real one. */
	if(!DpRt_Get_Property_Boolean("dprt.fake",&fake))
		return FALSE;
	fprintf(stdout,"DpRt_Make_Master_Flat:Fake:%d\n",fake);
	if(!DpRt_Get_Property_Boolean("dprt.make_master_flat",&make_master_flat))
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
			retval = dprt_process(directory_name,MAKE_FLAT,&l1mean,&l1seeing, 
					      &l1xpix,&l1ypix,&l1counts,&l1sat,&l1photom,&l1skybright);
			fprintf(stdout,"DpRt_Make_Master_Flat:Make Master Flat routine (dprt_process) returned %d.\n",
				retval);
			if(retval == TRUE)
			{
				DpRt_Error_Number = dprt_err_int;
				strcpy(DpRt_Error_String,dprt_err_str);
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

/**
 * A routine to set the value of the DpRt_Abort variable. This is used to
 * keep track of whether an abort event has occured.
 * @param value The value to set the abort value to, this should be TRUE of we want the data reduction
 * aborted, or FALSE if we don't or we are resetting the flag.
 * @see #DpRt_Data
 */
void DpRt_Set_Abort(int value)
{
	DpRt_Data.DpRt_Abort = value;
}

/**
 * This routine allows us to query the properties loaded into DpRt to get
 * the value associated with a keyword in the property Hashtable.
 * As libdprt can be called in two ways, from Java using JNI (the usual method) and
 * from a C test program (for testing), the mechanism for retrieving the value is different
 * in each case. A function pointer system is used.
 * @param keyword The keyword in the property file to look up.
 * @param value_string The address of a pointer to allocate and store the resulting value string in.
 * 	This pointer is dynamically allocated and must be freed using <b>free()</b>. 
 * @return The routine returns TRUE if it succeeds, FALSE if it fails.
 * @see #DpRt_Data
 */
int DpRt_Get_Property(char *keyword,char **value_string)
{
	if(keyword == NULL)
	{
		DpRt_Error_Number = 5;
		sprintf(DpRt_Error_String,"DpRt_Get_Property failed: Keyword was NULL.\n");
		return FALSE;
	}
	if(value_string == NULL)
	{
		DpRt_Error_Number = 6;
		sprintf(DpRt_Error_String,"DpRt_Get_Property failed: Value String Pointer was NULL.\n");
		return FALSE;
	}
	if(DpRt_Data.DpRt_Get_Property_Function_Pointer == NULL)
	{
		DpRt_Error_Number = 7;
		sprintf(DpRt_Error_String,"DpRt_Get_Property failed: Function Pointer was NULL.\n");
		return FALSE;
	}
	return DpRt_Data.DpRt_Get_Property_Function_Pointer(keyword,value_string);
}

/**
 * This routine allows us to query the properties loaded into DpRt to get
 * the integer value associated with a keyword in the property Hashtable.
 * As libdprt can be called in two ways, from Java using JNI (the usual method) and
 * from a C test program (for testing), the mechanism for retrieving the value is different
 * in each case. A function pointer system is used.
 * @param keyword The keyword in the property file to look up.
 * @param value The address of an integer to store the resulting value.
 * @return The routine returns TRUE if it succeeds, FALSE if it fails.
 * @see #DpRt_Data
 */
int DpRt_Get_Property_Integer(char *keyword,int *value)
{
	if(keyword == NULL)
	{
		DpRt_Error_Number = 8;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Integer failed: Keyword was NULL.\n");
		return FALSE;
	}
	if(value == NULL)
	{
		DpRt_Error_Number = 9;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Integer failed: Value Pointer was NULL.\n");
		return FALSE;
	}
	if(DpRt_Data.DpRt_Get_Property_Integer_Function_Pointer == NULL)
	{
		DpRt_Error_Number = 10;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Integer failed: Function Pointer was NULL.\n");
		return FALSE;
	}
	return DpRt_Data.DpRt_Get_Property_Integer_Function_Pointer(keyword,value);
}

/**
 * This routine allows us to query the properties loaded into DpRt to get
 * the double value associated with a keyword in the property Hashtable.
 * As libdprt can be called in two ways, from Java using JNI (the usual method) and
 * from a C test program (for testing), the mechanism for retrieving the value is different
 * in each case. A function pointer system is used.
 * @param keyword The keyword in the property file to look up.
 * @param value The address of an double to store the resulting value.
 * @return The routine returns TRUE if it succeeds, FALSE if it fails.
 * @see #DpRt_Data
 */
int DpRt_Get_Property_Double(char *keyword,double *value)
{
	if(keyword == NULL)
	{
		DpRt_Error_Number = 11;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Double failed: Keyword was NULL.\n");
		return FALSE;
	}
	if(value == NULL)
	{
		DpRt_Error_Number = 12;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Double failed: Value Pointer was NULL.\n");
		return FALSE;
	}
	if(DpRt_Data.DpRt_Get_Property_Double_Function_Pointer == NULL)
	{
		DpRt_Error_Number = 13;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Double failed: Function Pointer was NULL.\n");
		return FALSE;
	}
	return DpRt_Data.DpRt_Get_Property_Double_Function_Pointer(keyword,value);
}

/**
 * This routine allows us to query the properties loaded into DpRt to get
 * the boolean value associated with a keyword in the property Hashtable.
 * As libdprt can be called in two ways, from Java using JNI (the usual method) and
 * from a C test program (for testing), the mechanism for retrieving the value is different
 * in each case. A function pointer system is used.
 * @param keyword The keyword in the property file to look up.
 * @param value The address of an integer to store the resulting value, 1 is TRUE, 0 is FALSE.
 * @return The routine returns TRUE if it succeeds, FALSE if it fails.
 * @see #DpRt_Data
 */
int DpRt_Get_Property_Boolean(char *keyword,int *value)
{
	if(keyword == NULL)
	{
		DpRt_Error_Number = 14;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Boolean failed: Keyword was NULL.\n");
		return FALSE;
	}
	if(value == NULL)
	{
		DpRt_Error_Number = 15;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Boolean failed: Value Pointer was NULL.\n");
		return FALSE;
	}
	if(DpRt_Data.DpRt_Get_Property_Boolean_Function_Pointer == NULL)
	{
		DpRt_Error_Number = 16;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Boolean failed: Function Pointer was NULL.\n");
		return FALSE;
	}
	return DpRt_Data.DpRt_Get_Property_Boolean_Function_Pointer(keyword,value);
}

/**
 * Routine to set the function pointer that is called from <b>DpRt_Get_Property</b> .
 * @see #DpRt_Get_Property
 * @see #DpRt_Data
 */
void DpRt_Set_Property_Function_Pointer(int (*get_property_fp)(char *keyword,char **value_string))
{
	DpRt_Data.DpRt_Get_Property_Function_Pointer = get_property_fp;
}

/**
 * Routine to set the function pointer that is called from <b>DpRt_Get_Property_Integer</b> .
 * @see #DpRt_Get_Property_Integer
 * @see #DpRt_Data
 */
void DpRt_Set_Property_Integer_Function_Pointer(int (*get_property_integer_fp)(char *keyword,int *value))
{
	DpRt_Data.DpRt_Get_Property_Integer_Function_Pointer = get_property_integer_fp;
}

/**
 * Routine to set the function pointer that is called from <b>DpRt_Get_Property_Double</b> .
 * @see #DpRt_Get_Property_Double
 * @see #DpRt_Data
 */
void DpRt_Set_Property_Double_Function_Pointer(int (*get_property_double_fp)(char *keyword,double *value))
{
	DpRt_Data.DpRt_Get_Property_Double_Function_Pointer = get_property_double_fp;
}

/**
 * Routine to set the function pointer that is called from <b>DpRt_Get_Property_Boolean</b> .
 * @see #DpRt_Get_Property_Boolean
 * @see #DpRt_Data
 */
void DpRt_Set_Property_Boolean_Function_Pointer(int (*get_property_boolean_fp)(char *keyword,int *value))
{
	DpRt_Data.DpRt_Get_Property_Boolean_Function_Pointer = get_property_boolean_fp;
}

/**
 * A routine to return the current value of the error number. The error number is usually 0 for success,
 * and non-zero when an error occurs.
 * @return The current value of the error number.
 * @see #DpRt_Error_Number
 */
int DpRt_Get_Error_Number(void)
{
	return DpRt_Error_Number;
}

/**
 * A routine to get the current value of the error string. The current value is strcpyed into the 
 * supplied paramater. The error string is usually blank when no errors have occured, and non-blank when
 * an error has occured.
 * @param error_string A pointer to an area of memory reserved to accept the error string. This area of 
 * memory should be at least <a href="#DPRT_ERROR_STRING_LENGTH">DPRT_ERROR_STRING_LENGTH</a> in length.
 * @see #DpRt_Error_String
 */
void DpRt_Get_Error_String(char *error_string)
{
	strcpy(error_string,DpRt_Error_String);
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
 * @see DpRtLibrary.html
 * @see #DpRt_Get_Abort
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
	DpRt_Error_Number = 0;
	strcpy(DpRt_Error_String,"");
/* setup return values */
	(*mean_counts) = 0.0;
	(*peak_counts) = 0.0;
/* unset any previous aborts - ready to start processing */
	DpRt_Set_Abort(FALSE);
/* do processing  here */
	fprintf(stderr,"Calibrate_Reduce_Fake(%s).\n",input_filename);
/* open file */
	retval = fits_open_file(&fp,input_filename,READONLY,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 23;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Open failed.\n",input_filename);
		return FALSE;
	}
/* check bitpix */
	retval = fits_read_key(fp,TINT,"BITPIX",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 24;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get BITPIX.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_BITPIX)
	{
		DpRt_Error_Number = 25;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Wrong BITPIX value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* check naxis */
	retval = fits_read_key(fp,TINT,"NAXIS",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 26;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_NAXIS)
	{
		DpRt_Error_Number = 27;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Wrong NAXIS value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* get naxis1,naxis2 */
	retval = fits_read_key(fp,TINT,"NAXIS1",&naxis_one,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 28;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS1.\n",input_filename);
		return FALSE;
	}
	retval = fits_read_key(fp,TINT,"NAXIS2",&naxis_two,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 29;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to get NAXIS2.\n",input_filename);
		return FALSE;
	}
/* allocate data */
	data = (unsigned short *)malloc(naxis_one*naxis_two*sizeof(unsigned short));
	if(data == NULL)
	{
		DpRt_Error_Number = 30;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to allocate memory (%d,%d,%d).\n",
			input_filename,naxis_one,naxis_two,naxis_one*naxis_two*sizeof(short));
		return FALSE;
	}
/* read the data */
	retval = fits_read_img(fp,TUSHORT,1,naxis_one*naxis_two,NULL,data,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 31;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to read image(%d,%d).\n",
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
		DpRt_Error_Number = 32;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Failed to close file.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* during processing regularily check the abort flag as below */
	if(DpRt_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*mean_counts) = 0.0;
		(*peak_counts) = 0.0;
		(*output_filename) = NULL;
		DpRt_Error_Number = 1;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
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
		if(DpRt_Get_Abort())
		{
			/* tidy up anything that needs tidying as a result of this routine here */
			(*mean_counts) = 0.0;
			(*peak_counts) = 0.0;
			(*output_filename) = NULL;
			DpRt_Error_Number = 45;
			sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
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
		DpRt_Error_Number = 2;
		sprintf(DpRt_Error_String,"Calibrate_Reduce_Fake(%s): Memory Allocation Error.\n",input_filename);
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
 * @see DpRtLibrary.html
 * @see #DpRt_Get_Abort
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
	DpRt_Error_Number = 0;
	strcpy(DpRt_Error_String,"");

/* unset any previous aborts - ready to start processing */
	DpRt_Set_Abort(FALSE);
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
	if(!DpRt_Get_Property_Double("dprt.telfocus.best_focus",&best_focus))
		return FALSE;
	if(!DpRt_Get_Property_Double("dprt.telfocus.fwhm_per_mm",&fwhm_per_mm))
		return FALSE;
	if(!DpRt_Get_Property_Double("dprt.telfocus.atmospheric_seeing",&atmospheric_seeing))
		return FALSE;
	if(!DpRt_Get_Property_Double("dprt.telfocus.atmospheric_variation",&atmospheric_variation))
		return FALSE;
/* open file */
	retval = fits_open_file(&fp,input_filename,READONLY,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 33;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Open failed.\n",input_filename);
		return FALSE;
	}
/* check bitpix */
	retval = fits_read_key(fp,TINT,"BITPIX",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 34;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to get BITPIX.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_BITPIX)
	{
		DpRt_Error_Number = 35;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Wrong BITPIX value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* check naxis */
	retval = fits_read_key(fp,TINT,"NAXIS",&integer_value,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 36;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS.\n",input_filename);
		return FALSE;
	}
	if(integer_value != FITS_GET_DATA_NAXIS)
	{
		DpRt_Error_Number = 37;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Wrong NAXIS value(%d).\n",
			input_filename,integer_value);
		return FALSE;
	}
/* get naxis1,naxis2 */
	retval = fits_read_key(fp,TINT,"NAXIS1",&naxis_one,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 38;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS1.\n",input_filename);
		return FALSE;
	}
	retval = fits_read_key(fp,TINT,"NAXIS2",&naxis_two,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 39;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to get NAXIS2.\n",input_filename);
		return FALSE;
	}
/* get telescope focus */
	retval = fits_read_key(fp,TDOUBLE,"TELFOCUS",&telfocus,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 40;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to get TELFOCUS.\n",input_filename);
		return FALSE;
	}
/* allocate data */
	data = (unsigned short *)malloc(naxis_one*naxis_two*sizeof(unsigned short));
	if(data == NULL)
	{
		DpRt_Error_Number = 41;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to allocate memory (%d,%d,%d).\n",
			input_filename,naxis_one,naxis_two,naxis_one*naxis_two*sizeof(short));
		return FALSE;
	}
/* read the data */
	retval = fits_read_img(fp,TUSHORT,1,naxis_one*naxis_two,NULL,data,NULL,&status);
	if(retval)
	{
		fits_report_error(stderr,status);
		DpRt_Error_Number = 42;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to read image(%d,%d).\n",
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
		DpRt_Error_Number = 43;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Failed to close file.\n",input_filename);
		if(data != NULL)
			free(data);
		return FALSE;
	}
/* during processing regularily check the abort flag as below */
	if(DpRt_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*output_filename) = NULL;
		DpRt_Error_Number = 44;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
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
	if(DpRt_Get_Abort())
	{
		/* tidy up anything that needs tidying as a result of this routine here */
		(*output_filename) = NULL;
		DpRt_Error_Number = 3;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Operation Aborted.\n",input_filename);
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
		DpRt_Error_Number = 4;
		sprintf(DpRt_Error_String,"Expose_Reduce_Fake(%s): Memory Allocation Error.\n",input_filename);
		return FALSE;
	}
/* set the filename to something more sensible here */
	strcpy((*output_filename),input_filename);
	return TRUE;
}

/**
 * A routine to get the current value of the abort variable, to determine whether we should abort
 * processing the FITS file or not.
 * @return The current value of the DpRt_Data.DpRt_Abort variable, usually TRUE if we want to abort a reduction process
 * and FALSE if we don't.
 * @see #DpRt_Data
 */
static int DpRt_Get_Abort(void)
{
	return DpRt_Data.DpRt_Abort;
}

/**
 * Routine to get the value of the keyword from the property file.
 * This routine assumes keyword and value_string have been checked as being non-null.
 * @param keyword The keyword in the property file to look up.
 * @param value_string The address of a pointer to allocate and store the resulting value string in.
 * 	This pointer is dynamically allocated and must be freed using <b>free()</b>. 
 * @see #PROPERTY_FILE_NAME
 */
static int DpRt_Get_Property_From_C_File(char *keyword,char **value_string)
{
	char buff[256];
	FILE *fp = NULL;
	char *ch = NULL;
	int done;

	fp = fopen(PROPERTY_FILE_NAME,"r");
	if(fp == NULL)
	{
		DpRt_Error_Number = 17;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_From_C_File failed: File open (%s) failed.\n",
			PROPERTY_FILE_NAME);
		return FALSE;
	}
	done = FALSE;
	(*value_string) = NULL;
	while((done == FALSE)&&(fgets(buff,255,fp) != NULL))
	{
		if(strncmp(keyword,buff,strlen(keyword))==0)
		{
			ch = strchr(buff,'=');
			if(ch != NULL)
			{
				(*value_string) = (char*)malloc((strlen(ch)+1)*sizeof(char));
				if((*value_string) == NULL)
				{
					fclose(fp);
					DpRt_Error_Number = 18;
					sprintf(DpRt_Error_String,"DpRt_Get_Property_From_C_File failed: "
						"Memory Allocation error(%s,%s,%s,%d) failed.\n",
						PROPERTY_FILE_NAME,keyword,ch,(strlen(ch)+1));
					return FALSE;
				}
				strcpy((*value_string),ch+1);
			/* if the string terminates in a new-line, remove it */
				if((strlen((*value_string))>0) && ((*value_string)[strlen((*value_string))-1] == '\n'))
					(*value_string)[strlen((*value_string))-1] = '\0';
				done = TRUE;
			}
		}
	}
	fclose(fp);
	if(done == FALSE)
	{
		DpRt_Error_Number = 19;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_From_C_File failed:Failed to find keyword (%s,%s).\n",
			PROPERTY_FILE_NAME,keyword);
	}
	return done;
}

/**
 * Routine to get the integer value of the keyword from the property file.
 * This routine assumes keyword and value have been checked as being non-null.
 * DpRt_Get_Property_From_C_File is used to get the keyword's value, and sscanf used to convert it to
 * an integer.
 * @param keyword The keyword in the property file to look up.
 * @param value_string The address of an integer to store the resulting value string in.
 * @see #DpRt_Get_Property_From_C_File
 */
static int DpRt_Get_Property_Integer_From_C_File(char *keyword,int *value)
{
	char *value_string = NULL;
	int retval;

	if(!DpRt_Get_Property_From_C_File(keyword,&value_string))
		return FALSE;
	retval = sscanf(value_string,"%i",value);
	if(retval != 1)
	{
		DpRt_Error_Number = 20;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Integer_From_C_File failed:"
			"Failed to convert (%s,%s,%s).\n",PROPERTY_FILE_NAME,keyword,value_string);
		if(value_string != NULL)
			free(value_string);
		return FALSE;
	}
	if(value_string != NULL)
		free(value_string);
	return TRUE;
}

/**
 * Routine to get the double value of the keyword from the property file.
 * This routine assumes keyword and value have been checked as being non-null.
 * DpRt_Get_Property_From_C_File is used to get the keyword's value, and sscanf used to convert it to
 * a double.
 * @param keyword The keyword in the property file to look up.
 * @param value_string The address of an double to store the resulting value string in.
 * @see #DpRt_Get_Property_From_C_File
 */
static int DpRt_Get_Property_Double_From_C_File(char *keyword,double *value)
{
	char *value_string = NULL;
	int retval;

	if(!DpRt_Get_Property_From_C_File(keyword,&value_string))
		return FALSE;
	retval = sscanf(value_string,"%lf",value);
	if(retval != 1)
	{
		DpRt_Error_Number = 21;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Double_From_C_File failed:"
			"Failed to convert (%s,%s,%s).\n",PROPERTY_FILE_NAME,keyword,value_string);
		if(value_string != NULL)
			free(value_string);
		return FALSE;
	}
	if(value_string != NULL)
		free(value_string);
	return TRUE;
}

/**
 * Routine to get the boolean value of the keyword from the property file.
 * This routine assumes keyword and value have been checked as being non-null.
 * DpRt_Get_Property_From_C_File is used to get the keyword's value, and the string checked to see if it
 * contains <b>true</b> or <b>false</b>.
 * @param keyword The keyword in the property file to look up.
 * @param value_string The address of an integer to store the resulting value, either TRUE (1) or FALSE (0).
 * @see #DpRt_Get_Property_From_C_File
 */
static int DpRt_Get_Property_Boolean_From_C_File(char *keyword,int *value)
{
	char *value_string = NULL;

	if(!DpRt_Get_Property_From_C_File(keyword,&value_string))
		return FALSE;
	if((strcmp(value_string,"true")==0)||(strcmp(value_string,"TRUE")==0)||(strcmp(value_string,"True")==0))
		(*value) = TRUE;
	else if((strcmp(value_string,"false")==0)||(strcmp(value_string,"FALSE")==0)||
		(strcmp(value_string,"False")==0))
		(*value) = FALSE;
	else
	{
		DpRt_Error_Number = 22;
		sprintf(DpRt_Error_String,"DpRt_Get_Property_Boolean_From_C_File failed:"
			"Failed to convert (%s,%s,%s).\n",PROPERTY_FILE_NAME,keyword,value_string);
		if(value_string != NULL)
			free(value_string);
		return FALSE;
	}
	if(value_string != NULL)
		free(value_string);
	return TRUE;
}

/*
** $Log: not supported by cvs2svn $
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
