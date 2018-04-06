/*
 * ==========================  sfc =====================
 *   sfc shows 2-d field as a surface and shows the max/min, time & label.
 *   Updated to transpose X,Y for correct view when plotted
 *
 *   Note: this routine plots the entire field, including any ghost zones.
 *   All array contents must be set (including ghost zones, if any)
 *
 *   ATMS 502 / CSE 566, Spring 2016
 *
 *   Arguments:
 * 
 * 	q	input	real array	field to be displayed.
 * 	nx,ny	input	integers	dimensions of 'q', including ghost zones
 * 	time	input	real		time of integration
 * 	angh	input	real		horizontal viewing angle counter-
 * 					 clockwise from x-axis
 * 	angv	input	real		vertical viewing angle; >0=above plane
 * 					 of field average.
 * 	label	input	character*	character label
 *	name	input	character*	character label
 */

#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define IWTYPE 1
#define WKID   1

void sfc(q,nx,ny,time,angh,angv,label,name)
int nx,ny;
float q[][ny],time,angh,angv;
char *label,*name;
{

	int i,j;
	float qmin,qmax,*workarray;
	float *swapxy,*ptr;
        Gcolr_rep rgb1,rgb2;
	char mmlabel[80],tlabel[40];

	/*
	 * Find min,max
	 */

	qmin = q[0][0];
	qmax = q[0][0];
	for (j=0; j<ny; j++) {
	  for (i=0; i<nx; i++) {
	    if (q[i][j] < qmin) qmin=q[i][j];
	    if (q[i][j] > qmax) qmax=q[i][j];
	  }
  	}
	/*printf("Sfc: min %.2f, max %.2f\n",qmin,qmax);*/

        /*
         * Allocate and fill temporary array with 
         * dimensions reversed for NCAR routine
         */

        swapxy = (float *) malloc(nx*ny*sizeof(float));
	ptr = swapxy;

        for (i=0; i<nx; i++) {
          for (j=0; j<ny; j++) {
            *ptr++ = q[i][j];
          }
        }

	/*
	 * Create min/max and time labels
	 */

	sprintf(mmlabel,"MIN=%.3f  MAX=%.3f",qmin,qmax);
	sprintf( tlabel,"TIME = %.4f",time);

	/*
	 * Plot labels
	 */

	c_set(0.,1.,0.,1.,0.,1.,0.,1.,1);
	c_pcmequ(0.50,0.97,label  ,0.020,0.0, 0.0);
	c_pcmequ(0.95,0.02,mmlabel,0.015,0.0, 1.0);
	c_pcmequ(0.05,0.02, tlabel,0.015,0.0,-1.0);

	/*
	 * Additional labels
	 */

        c_pcmequ(0.02,0.99,name,0.01,90.,1.);
        c_pcmequ(0.98,0.06,"ATMS 502/CSE 566",0.01,0.,1.);
        c_pcmequ(0.02,0.06,"Spring 2016",0.01,0.,-1.);

	/*
	 * Set colors to default (commented out for now)
	 */

	/*
        rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
        rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
        gset_colr_rep(WKID,0,&rgb1);
        gset_colr_rep(WKID,1,&rgb2);
	*/

	/*
	 * Allocate scratch work array for ezsrfc
	 */

	workarray = (float *) malloc((2*nx*ny+nx+ny)*sizeof(float));

	/*
	 * Plot surface
	 */

	c_ezsrfc(swapxy,ny,nx,angh,angv,workarray);

	free(swapxy);
	free(workarray);
	return;
}
