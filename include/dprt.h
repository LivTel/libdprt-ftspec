/* dprt.h
** $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/include/dprt.h,v 1.1 2005-03-11 11:56:53 cjm Exp $
*/
#ifndef DPRT_H
#define DPRT_H

/**
 * TRUE is the value usually returned from routines to indicate success.
 */
#ifndef TRUE
#define TRUE 1
#endif
/**
 * FALSE is the value usually returned from routines to indicate failure.
 */
#ifndef FALSE
#define FALSE 0
#endif

/* function declarations */
extern int DpRt_Initialise(void);
extern int DpRt_Shutdown(void);
extern int DpRt_Calibrate_Reduce(char *input_filename,char **output_filename,double *mean_counts,double *peak_counts);
extern int DpRt_Expose_Reduce(char *input_filename,char **output_filename,double *seeing,double *counts,double *x_pix,
		       double *y_pix,double *photometricity,double *sky_brightness,int *saturated);
extern int DpRt_Make_Master_Bias(char *directory_name);
extern int DpRt_Make_Master_Flat(char *directory_name);
#endif
