/* dprt.c
** Entry point for Data Pipeline Reduction Routines
** $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/dprt.c,v 0.3 1999-08-24 16:04:30 cjm Exp $
*/
/**
 * dprt.c is the entry point for the Data Reduction Pipeline (Real Time).
 * @author Lee Howells, LJMU
 * @version $Revision: 0.3 $
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dprt.h"

/* ------------------------------------------------------- */
/* internal variables */
/* ------------------------------------------------------- */
/**
 * Revision Control System identifier.
 */
static char rcsid[] = "$Id: dprt.c,v 0.3 1999-08-24 16:04:30 cjm Exp $";
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
 * Boolean value holding whether we should abort the current Data Reduction Pipeline. The reduction
 * routines should check this at regular intervals, and if it is TRUE
 * should abort processing and set the <a href="#DpRt_Error_Number">DpRt_Error_Number</a> and
 * <a href="#DpRt_Error_String">DpRt_Error_String</a> accordingly and return FALSE.
 * @see #DpRt_Calibrate_Reduce
 * @see #DpRt_Expose_Reduce
 * @see #DpRt_Error_Number
 * @see #DpRt_Error_String
 * @see #DpRt_Abort
 */
static volatile int DpRt_Abort = FALSE;/* diddly volatile as thread dependant */

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */
static int DpRt_Get_Abort(void);

/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
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
 */
int DpRt_Calibrate_Reduce(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts)
{
  /* set the error stuff to no error*/
  DpRt_Error_Number = 0;
  strcpy(DpRt_Error_String,"");

  /* unset any previous aborts - ready to start processing */
  DpRt_Set_Abort(FALSE);

  /* do processing  here */
  fprintf(stderr,"DpRt_Calibrate_Reduce(%s).\n",input_filename);

  /* during processing regularily check the abort flag as below */
  if(DpRt_Get_Abort())
  {
    /* tidy up anything that needs tidying as a result of this routine here */
    (*mean_counts) = 0.0;
    (*peak_counts) = 0.0;
    (*output_filename) = NULL;
    DpRt_Error_Number = 1;
    sprintf(DpRt_Error_String,"DpRt_Calibrate_Reduce(%s): Operation Aborted.\n",input_filename);
    return FALSE;
  }

  /* setup return values */
  (*mean_counts) = 1.0;
  (*peak_counts) = 2.0;

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
    sprintf(DpRt_Error_String,"DpRt_Reduce(%s): Memory Allocation Error.\n",input_filename);
    return FALSE;
  }
  /* set the filename to something more sensible here */
  strcpy((*output_filename),input_filename);

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
 * @return The routine should return whether it succeeded or not. TRUE should be returned if the routine
 *       succeeded and FALSE if they fail.
 * @see DpRtLibrary.html
 * @see #DpRt_Get_Abort
 */
int DpRt_Expose_Reduce(char *input_filename,char **output_filename,double *seeing,double *counts,double *x_pix,
		       double *y_pix)
{
  /* set the error stuff to no error*/
  DpRt_Error_Number = 0;
  strcpy(DpRt_Error_String,"");

  /* unset any previous aborts - ready to start processing */
  DpRt_Set_Abort(FALSE);

  /* do processing  here */
  fprintf(stderr,"DpRt_Expose_Reduce(%s).\n",input_filename);

  /* during processing regularily check the abort flag as below */
  if(DpRt_Get_Abort())
  {
    /* tidy up anything that needs tidying as a result of this routine here */
    (*output_filename) = NULL;
    (*seeing) = 0.0;
    (*counts) = 0.0;
    (*x_pix) = 0.0;
    (*y_pix) = 0.0;
    DpRt_Error_Number = 3;
    sprintf(DpRt_Error_String,"DpRt_Expose_Reduce(%s): Operation Aborted.\n",input_filename);
    return FALSE;
  }

  /* setup return values */
  (*seeing) = 1.0;
  (*counts) = 2.0;
  (*x_pix) = 3.0;
  (*y_pix) = 4.0;

  /* setup filename - allocate space for string */
  (*output_filename) = (char*)malloc((strlen(input_filename)+1)*sizeof(char));
  /* if malloc fails it returns NULL - this is an error */
  if((*output_filename) == NULL)
  {
    /* tidy up anything that needs tidying as a result of this routine here */
    (*output_filename) = NULL;
    (*seeing) = 0.0;
    (*counts) = 0.0;
    (*x_pix) = 0.0;
    (*y_pix) = 0.0;
    DpRt_Error_Number = 4;
    sprintf(DpRt_Error_String,"DpRt_Expose_Reduce(%s): Memory Allocation Error.\n",input_filename);
    return FALSE;
  }
  /* set the filename to something more sensible here */
  strcpy((*output_filename),input_filename);
  return TRUE;
}

/**
 * A routine to set the value of the <a href="#DpRt_Abort">DpRt_Abort</a> variable. This is used to
 * keep track of whether an abort event has occured.
 * @param value The value to set the abort value to, this should be TRUE of we want the data reduction
 * aborted, or FALSE if we don't or we are resetting the flag.
 * @see #DpRt_Abort
 */
void DpRt_Set_Abort(int value)
{
  DpRt_Abort = value;
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
 * A routine to get the current value of the abort variable, to determine whether we should abort
 * processing the FITS file or not.
 * @return The current value of the DpRt_Abort variable, usually TRUE if we want to abort a reduction process
 * and FALSE if we don't.
 * @see #DpRt_Abort
 */
static int DpRt_Get_Abort(void)
{
  return DpRt_Abort;
}

/*
** $Log: not supported by cvs2svn $
** Revision 0.2  1999/06/30 15:07:47  dev
** changes to return when error occurs in reduction code
** set return paramaters to 0/NULL.
**
** Revision 0.1  1999/06/24 11:06:42  dev
** initial revision
**
*/
