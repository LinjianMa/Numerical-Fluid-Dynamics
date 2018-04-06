/*
 * ========================= plot1d ====================
 * Plot1d plots a one-dimensional function using NCAR graphics
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q	real array	array to plot (entire array passed here;
 *				  only i1 ... i2 is plotted)
 *	i1,i2	integers	indices bounding (physical portion of) array data
 *	nx	integer		size of array, not counting ghost points
 *	n	integer		time step counter, from main program
 *	qmax	real		max absolute value in s array - used in label
 *	overlay	integer		if nonzero, the plot will be an overlay of two fields
 *	qtrue	real array	[optional] overlay field, ignored if overlay=0
 *	title	char string	title (step, max info appended to this)
 *	name	char string	name to place on plot
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#define IWTYPE 1
#define WKID   1

void plot1d(q,i1,i2,nx,n,qmax,overlay,qtrue,title,name)
float q[],qmax,qtrue[];
int i1,i2,nx,n,overlay;
char *title,*name;
{
	int i,j;
	char label[200];
	float *both;

	sprintf(label,"%s: N=%d Max=%7.5f",title,n,qmax);
	printf("Plotting \"%s\"\n",label);

/*
 * ... If overlay=0, simply plot s starting from i1;
 * ... if overlay=1, fill the both[] array with both fields, and plot that.
 */
	c_agseti("frame.",2);	/* suppress automatic graphics frame advance for now */

	if (overlay == 0) {

	  c_ezy(&q[i1],nx,label);

	} else {

	  both = (float *) malloc(2*nx*sizeof(float));
	  for (i=i1,j=0; i<=i2; i++,j++) {
	    both[j]    = q[i];
	    both[nx+j] = qtrue[i];
	  }
	  c_agsetr("DASH/SELECTOR.",2.0);
	  c_agseti("DASH/LENGTH.",20);
	  c_agsetc(c_agdshn(1),"$$$$$$$$$$$$$$$$$$$$");
	  c_agsetc(c_agdshn(2),"$$''$$''$$''$$''$$''");
	  c_ezmy(both,nx,2,nx,label);
	  free(both);
	}
/*
 * ... Add additional labels
 */
	c_set(.05,.95,.05,.95,0.,1.,0.,1.,1);
        c_pcmequ(0.02,1.0,name,0.01,90.,1.);
        c_pcmequ(1.0,0.05,"ATMS 502/CSE 566",0.01,0.,1.);
        c_pcmequ(0.0,0.05,"Spring 2016",0.01,0.,-1.);

/*
 * ... Done
 */
	c_frame();
	c_agseti("frame.",1);
	return;
}
