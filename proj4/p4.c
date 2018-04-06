/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm4:  Linear and nonlinear advection
 *  >>>>> Linjian Ma <<<<<
 */
#include <time.h> 
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

time_t time1,time2;
float elapsed;
time1 = time(0);

int ratio;        /* nest ratio */
float dx   = 1.0/(NX-1);
float dy   = 1.0/(NY-1);
int   nstep;
int   nplot;
int   nmove;
char *name = "Linjian Ma";
int nest_option;
int nestX1,nestX2,nestY1,nestY2,nestX1old,nestX2old,nestY1old,nestY2old;
nestX1 = -1;
nestX2 = 0;
nestY1 = 0;
nestY2 = 0;
float q4[NXDIM][NYDIM],q4_max[NXDIM][NYDIM],qtrue[NXDIM][NYDIM],qplot[NX][NY],q4nest[NXDIM][NXDIM];
float u[NX+1][NY],v[NX][NY+1],unest[NX+1][NY],vnest[NX][NY+1],uplot[NX][NY],vplot[NX][NY];
float q_av[MAXSTEP],qmin[MAXSTEP],qmax[MAXSTEP],courant,qtrue_av,qtrue_max,qtrue_min;
int trunc_xcenter,trunc_ycenter,trunc_x1,trunc_x2,trunc_y1,trunc_y2;
float trunc_emax;
int i,j,n;
int method = 1;
char label[200];
float cor_max,err_t,err_dissip,err_disper;
float cint    =  0.5;	
int colors    =  0;	
int pltzero   =  -1;	
float angh    =  45.0;
float angv    =  20.0;
float r_cone;
/* Parameters and input .................................... */

	printf("Program #4       Numerical Fluid Dynamics, Spring 2016\n\n");

	printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",NX,BC_WIDTH,I1,I2,NXDIM);
	printf("NY=%d, BC_WIDTH=%d, J1=%d, J2=%d, NYDIM=%d\n",NY,BC_WIDTH,J1,J2,NYDIM);
	printf("dx=%8.5f, dy=%8.5f\n",dx,dy);

	printf(" Enter the number of running time steps \n");
	scanf("%d",&nstep);
	printf(" Enter the plotting interval \n");
	scanf("%d",&nplot);
	printf(" Enter the nest movement interval \n");
	scanf("%d",&nmove);
	printf(" Feedback or not, 1=add feedback, 0=without feedback \n");
	scanf("%d",&nest_option);
	printf(" Input the nest ratio = coarse grid size/nest size \n");
	scanf("%d",&ratio);
	printf(" Input the initial cone radius \n");
	scanf("%f",&r_cone);

float dt       = M_PI/nstep;
float nestsize = (NX-1)/ratio;
float dx_nest   = 1.0/(NX-1)/ratio;
float dy_nest   = 1.0/(NY-1)/ratio;
float dt_nest   = dt/ratio;
/* Variables to reverse default black/white colors in NCAR Graphics */

	Gcolr_rep rgb1,rgb2;

/* Definitions for routines */

	void ic(),stats(),contr(),plot1d(),sfc(),bc(),integrate(),update(),\
		advec1d(),advection(),error(),dointerp(),nestwind();

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
	ic(q4,u,v,dx,dy,I1,I2,J1,J2,NX,NY,x0,y0,xstart,ystart,r_cone);
/* initial condition for q4_max */
	for ( j=J1; j<=J2; j++ ) 
    for ( i=I1; i<=I2; i++) 
	{
		q4_max[i][j] = q4[i][j];
    }
	stats(q4,I1,I2,J1,J2,q_av,qmax,qmin,n,NX,NY);

/*  plot the initial condition  */
	for ( j=J1; j<=J2; j++ ) 
    for ( i=I1; i<=I2; i++) 
	{
		qplot[i-I1][j-J1] = q4[i][j];
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
		qtrue[i][j]=q4[i][j];
	}
	qtrue_av=q_av[0];
	qtrue_max=qmax[0];
	qtrue_min=qmin[0];
/* Integrate */
 	for (n=1; n<=nstep; n++) {
			/* Move the nest */
			if ( n == 1)
			{
				/* Calculate truncation error */
				truncation(q4,u,v,dx,dy,I1,I2,J1,J2,NXDIM,NYDIM,BC_WIDTH,&trunc_emax,\
					&trunc_xcenter,&trunc_ycenter,&trunc_x1,&trunc_x2,&trunc_y1,&trunc_y2,dt);
				printf("Max trunc error = %8.5f, bounded by ( %d: %d, %d: %d); center = ( %d, %d)",\
					trunc_emax,trunc_x1,trunc_x2,trunc_y1,trunc_y2,trunc_xcenter,trunc_ycenter);
				nestX1 = trunc_xcenter - nestsize/2;
				nestX2 = nestX1 + nestsize;
				nestY1 = trunc_ycenter - nestsize/2;
				nestY2 = nestY1 + nestsize;
				if (nestX1 <= I1)
				{
					nestX1 = I1;
					nestX2 = nestX1 + nestsize;
				}
				if (nestX2 >= I2)
				{
					nestX2 = I2;
					nestX1 = nestX2 - nestsize;
				}
				if (nestY1 <= J1)
				{
					nestY1 = J1;
					nestY2 = nestY1 + nestsize;
				}
				if (nestY2 >= J2)
				{
					nestY2 = J2;
					nestY1 = nestY2 - nestsize;
				}
				/* move */
				dointerp(q4,q4nest,NX,NY,NXDIM,NYDIM,I1,I2,J1,J2,nestX1,nestX2,nestY1,nestY2,\
					n,ratio,1);
				nestX1old = nestX1;
				nestX2old = nestX2;
				nestY1old = nestY1;
				nestY2old = nestY2;
				/* plot */
				for ( j=J1; j<=J2; j++ ) 
				for ( i=I1; i<=I2; i++) 
				{
					qplot[i-I1][j-J1] = q4[i][j];
				}
				printf("Plotting contours.\n");
				sprintf(label,"Contour plot at n = %d",n);
				contr(qplot,NX,NY,cint,(dt*(float)n),label,colors,pltzero,nestX1,nestX2,nestY1,nestY2,name);
				/* nest plot */
				for ( j=J1; j<=J2; j++ ) 
				for ( i=I1; i<=I2; i++) 
				{
					qplot[i-I1][j-J1] = q4nest[i][j];
				}
				printf("Plotting contours.\n");
				sprintf(label,"nested contour plot at n = %d",n);
				contr(qplot,NX,NY,cint,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
			}
			if ( (n % nmove)==0)
			{
				/* Calculate truncation error */
				truncation(q4,u,v,dx,dy,I1,I2,J1,J2,NXDIM,NYDIM,BC_WIDTH,&trunc_emax,\
					&trunc_xcenter,&trunc_ycenter,&trunc_x1,&trunc_x2,&trunc_y1,&trunc_y2,dt);
				printf("Max trunc error = %8.5f, bounded by ( %d: %d, %d: %d); center = ( %d, %d)",\
					trunc_emax,trunc_x1,trunc_x2,trunc_y1,trunc_y2,trunc_xcenter,trunc_ycenter);
				nestX1 = trunc_xcenter - nestsize/2;
				nestX2 = nestX1 + nestsize;
				nestY1 = trunc_ycenter - nestsize/2;
				nestY2 = nestY1 + nestsize;
				if (nestX1 <= I1)
				{
					nestX1 = I1;
					nestX2 = nestX1 + nestsize;
				}
				if (nestX2 >= I2)
				{
					nestX2 = I2;
					nestX1 = nestX2 - nestsize;
				}
				if (nestY1 <= J1)
				{
					nestY1 = J1;
					nestY2 = nestY1 + nestsize;
				}
				if (nestY2 >= J2)
				{
					nestY2 = J2;
					nestY1 = nestY2 - nestsize;
				}
				/* move */
				dointerp(q4,q4nest,NX,NY,NXDIM,NYDIM,I1,I2,J1,J2,nestX1,nestX2,nestY1,nestY2,\
					n,ratio,-1,nestX1old,nestX2old,nestY1old,nestY2old);
				nestX1old = nestX1;
				nestX2old = nestX2;
				nestY1old = nestY1;
				nestY2old = nestY2;
			}
			/* Compute values at next step */
			nestwind(unest,vnest,NX,NY,nestX1,nestY1,dx,ratio,1);
			for (i=1;i<=ratio;i++)
			{
				bc(q4,I1,I2,J1,J2,NXDIM,NYDIM,BC_WIDTH,1,q4nest,NX,NY,nestX1,nestX2,nestY1,nestY2,n,ratio);
				advection(q4nest,unest,vnest,dx_nest,dy_nest,dt_nest,I1,I2,J1,J2,NXDIM,NYDIM,method,BC_WIDTH,0);
			}
		bc(q4,I1,I2,J1,J2,NXDIM,NYDIM,BC_WIDTH,0);
		advection(q4,u,v,dx,dy,dt,I1,I2,J1,J2,NXDIM,NYDIM,method,BC_WIDTH,1);
		if (nest_option == 1)
		{
			dointerp(q4,q4nest,NX,NY,NXDIM,NYDIM,I1,I2,J1,J2,nestX1,nestX2,nestY1,nestY2,\
				n,ratio,2);
		}
/* q4_max */
        for ( j=J1; j<=J2; j++ ) 
		for ( i=I1; i<=I2; i++) 
		{
			if ( q4_max[i][j] < q4[i][j] )
			{
				q4_max[i][j] = q4[i][j];
			}
		}
/* Stats */
		stats(q4,I1,I2,J1,J2,q_av,qmax,qmin,n,NX,NY);
/* plot */
		if ( (n % nplot)==0)
		{
			for ( j=J1; j<=J2; j++ ) 
			for ( i=I1; i<=I2; i++) 
			{
				qplot[i-I1][j-J1] = q4[i][j];
			}
			printf("Plotting contours.\n");
			sprintf(label,"Contour plot at n = %d",n);
			if (nest_option == 1)
			{
				contr(qplot,NX,NY,cint,(dt*(float)n),label,colors,pltzero,nestX1,nestX2,nestY1,nestY2,name); 
				printf("Plotting surface.\n");
				sprintf(label,"Surface plot at n = %d",n);
				sfc(qplot,NX,NY,(dt*(float)n),angh,angv,label,name);
				/* nested plot */
				for ( j=J1; j<=J2; j++ ) 
				for ( i=I1; i<=I2; i++) 
				{
					qplot[i-I1][j-J1] = q4nest[i][j];
				}
				printf("Plotting contours.\n");
				sprintf(label,"nested contour plot at n = %d",n);
				contr(qplot,NX,NY,cint,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
			}
			if (nest_option == 0)
			{
				contr(qplot,NX,NY,cint,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
				printf("Plotting surface.\n");
				sprintf(label,"Surface plot at n = %d",n);
				sfc(qplot,NX,NY,(dt*(float)n),angh,angv,label,name);
			}
		}

	}	/* end of time loop n = 1,...,nstep */

/* plot q4-max */
    for ( j=J1; j<=J2; j++ ) 
	for ( i=I1; i<=I2; i++) 
	{
		qplot[i-I1][j-J1] = q4_max[i][j];
	}
	printf("Plotting contours.\n");
	sprintf(label,"2D max_q field");
	contr(qplot,NX,NY,cint,(dt*(float)nstep),label,colors,pltzero,-1,0,0,0,name); 
	printf("Plotting surface.\n");
	sprintf(label,"2D max_q field");
	sfc(qplot,NX,NY,(dt*(float)nstep),angh,angv,label,name);

/* calculate the error */
	error(q4,qtrue,I1,I2,J1,J2,NXDIM,NYDIM,nstep,q_av[nstep],qtrue_av,&cor_max,&err_t,&err_dissip,&err_disper);
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

time2 = time(0);
elapsed = (float) (time2 - time1);
printf("wallclock time = %8.5f\n",elapsed); 

	exit;
}

