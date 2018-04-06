/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm2:  Linear and nonlinear advection
 *  >>>>> Linjian Ma <<<<<
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#define IWTYPE 1
#define WKID   1

main()
{

/*
 * Definitions
 * NX is the number of physical grid points - Not including 'ghost points'
 * BC_WIDTH is the number of ghost points - on each side of the physical grid
 * I1 is the first C index to use - the left-most physical point
 * I2 is the last C index to use - the right-most physical point
 * NXDIM is the total array size, including grid points; use for declarations
 * MAXSTEP is the maximum number of time steps to take
 * "c" is [constant] flow speed (e.g. meters/sec) - used only in linear case
 * "dx" is grid spacing: distance (e.g. meters) between grid points in X, Y
 * "name" is a variable containing your name - used to label the plots.
 */

#define NX 121
#define NY 121

#define BC_WIDTH 1
#define I1 BC_WIDTH
#define I2 I1+NX-1
#define J1 BC_WIDTH
#define J2 J1+NY-1

#define NXDIM NX+2*BC_WIDTH
#define NYDIM NY+2*BC_WIDTH
#define x0 -0.5
#define y0 -0.5
#define xstart 0.0
#define ystart 0.3
#define MAXSTEP 1000

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
float dx   = 1.0/(NX-1);
float dy   = 1.0/(NY-1);
float dt   = M_PI/600;
int   nstep= 600;
int   nplot;
char *name = "Linjian Ma";

/* Arrays and other variables */

float q2[NXDIM][NYDIM],qtrue[NXDIM][NYDIM],qplot[NX][NY];
float u[NX+1][NY],v[NX][NY+1],uplot[NX][NY],vplot[NX][NY];
float q_av[MAXSTEP],qmin[MAXSTEP],qmax[MAXSTEP],courant,qtrue_av,qtrue_max,qtrue_min;
int i,j,n;
int method;
char label[200];
float cor_max,err_t,err_dissip,err_disper;
float cint    =  0.5;	
int colors    =  0;	
int pltzero   =  -1;	
float angh    =  45.0;
float angv    =  20.0;
/* Parameters and input .................................... */

	printf("Program #2       Numerical Fluid Dynamics, Spring 2016\n\n");

	printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",NX,BC_WIDTH,I1,I2,NXDIM);
	printf("NY=%d, BC_WIDTH=%d, J1=%d, J2=%d, NYDIM=%d\n",NY,BC_WIDTH,J1,J2,NYDIM);
	printf("dx=%8.5f, dy=%8.5f, dt=%8.5f, nstep=%d, nplot=%d\n",dx,dy,dt,nstep,nplot);

	printf(" Enter method, 1=lax 2=upstream \n");
	scanf("%d",&method);
	printf(" Enter the plotting interval \n");
	scanf("%d",&nplot);
/* Variables to reverse default black/white colors in NCAR Graphics */

	Gcolr_rep rgb1,rgb2;

/* Definitions for routines */

	void ic(),stats(),contr(),plot1d(),sfc(),bc(),integrate(),update(),advec1d(),advection(),error();

 /* Open the NCAR Graphics package and set colors.*/
	gopen_gks("stdout",0);
	gopen_ws(WKID, NULL, IWTYPE);
	gactivate_ws(WKID);

/* omit following four lines to invert black/white colors */
	rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
	rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
    gset_colr_rep(WKID,0,&rgb1);
    gset_colr_rep(WKID,1,&rgb2);

/* X-axis label*/
	c_anotat("step","S_value",0,0,0,0);

/* Set and plot the initial condition*/
    n = 0;
	ic(q2,u,v,dx,dy,I1,I2,J1,J2,NX,NY,x0,y0,xstart,ystart);
	stats(q2,I1,I2,J1,J2,q_av,qmax,qmin,n,NX,NY);

/*  plot the initial condition  */
	for ( j=J1; j<=J2; j++ ) 
    for ( i=I1; i<=I2; i++) 
	{
		qplot[i-I1][j-J1] = q2[i][j];
		uplot[i-I1][j-J1] = u[i-I1][j-J1];
		vplot[i-I1][j-J1] = v[i-I1][j-J1];
    }
	printf("Plotting contours.\n");
	sprintf(label,"Contour plot at n = %d",n);
	contr(uplot,NX,NY,0.1,0.0,label,colors,pltzero,-1,0,0,0,name); 
	contr(vplot,NX,NY,0.1,0.0,label,colors,pltzero,-1,0,0,0,name); 
	contr(qplot,NX,NY,cint,0.0,label,colors,pltzero,-1,0,0,0,name); 
	printf("Plotting surface.\n");
	sprintf(label,"Surface plot at n = %d",n);
	sfc(qplot,NX,NY,(dt*(float)n),angh,angv,label,name);

/* Copy the initial condition to the "qtrue" array. */
    for (i=I1; i<=I2; i++)
	for (j=J1; j<=J2; j++)
	{ 	
		qtrue[i][j]=q2[i][j];
	}
	qtrue_av=q_av[0];
	qtrue_max=qmax[0];
	qtrue_min=qmin[0];

/* Integrate */
 	for (n=1; n<=nstep; n++) {
/* Compute values at next step */
		advection(q2,u,v,dx,dy,dt,I1,I2,J1,J2,NXDIM,NYDIM,method);
/* Stats */
		stats(q2,I1,I2,J1,J2,q_av,qmax,qmin,n,NX,NY);
/* plot */
		if ( (n % nplot)==0)
		{
			for ( j=J1; j<=J2; j++ ) 
			for ( i=I1; i<=I2; i++) 
			{
				qplot[i-I1][j-J1] = q2[i][j];
			}
			printf("Plotting contours.\n");
			sprintf(label,"Contour plot at n = %d",n);
			contr(qplot,NX,NY,cint,(dt*(float)nstep),label,colors,pltzero,-1,0,0,0,name); 
			printf("Plotting surface.\n");
			sprintf(label,"Surface plot at n = %d",n);
			sfc(qplot,NX,NY,(dt*(float)nstep),angh,angv,label,name);
		}

	}	/* end of time loop n = 1,...,nstep */

/* calculate the error */
	error(q2,qtrue,I1,I2,J1,J2,NXDIM,NYDIM,nstep,q_av[nstep],qtrue_av,&cor_max,&err_t,&err_dissip,&err_disper);
	printf("cor = %8.5f\n",cor_max); 
	printf("total error = %8.5f\n",err_t); 
	printf("dissipation error = %8.5f\n",err_dissip); 
	printf("dispersion error = %8.5f\n",err_disper); 

	sprintf(label,"Smin vs. Time");
	c_ezy( qmin, nstep, label);
	sprintf(label,"Smax vs. Time");
	c_ezy( qmax, nstep, label);

/*  Close the graphics package and stop.*/
	gdeactivate_ws(WKID);
	gclose_ws(WKID);
	gclose_gks();

	exit;
}

