/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm1:  Linear and nonlinear advection
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

#define NX 75
#define BC_WIDTH 1
#define I1 BC_WIDTH
#define I2 I1+NX-1
#define NXDIM NX+2*BC_WIDTH
#define MAXSTEP 200

float c    = 1.0;
float dx   = 0.1;
char *name = "Linjian Ma";

/* Arrays and other variables */

	float q1[NXDIM],q2[NXDIM],qtrue[NXDIM];
	float qtrace[MAXSTEP],dt,courant,qmax;
	int i,j,n,nstep,nplot;
	char advection_type, reply[10];

/* Variables for run history */

	float history[MAXSTEP][NX];

/* Variables to reverse default black/white colors in NCAR Graphics */

	Gcolr_rep rgb1,rgb2;

/* Definitions for routines */

	void ic(),stats(),plot1d(),sfc(),bc(),integrate(),update();

/* Parameters and input .................................... */

	printf("Program #1       Numerical Fluid Dynamics, Spring 2016\n\n");

	printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",
		NX,BC_WIDTH,I1,I2,NXDIM);

	printf("Enter desired time step: ");
	scanf("%f",&dt);
	courant = c*dt/dx;

	printf("For time step dt = %6.3f, Courant number = %5.2f\n",  dt, courant);
	printf("Number of steps for a complete loop = %.0f\n",( dx*(float)NX / c / dt ));

	printf("Enter number of time steps to take: ");
	scanf("%d",&nstep);

	printf("Enter plot interval, in steps (1=each step): ");
	scanf("%d",&nplot);

	printf("Enter L for linear, N for nonlinear advection: ");
	scanf("%s",reply);
	advection_type = reply[0];
/*
 * Open the NCAR Graphics package and set colors.
 */
	gopen_gks("stdout",0);
	gopen_ws(WKID, NULL, IWTYPE);
	gactivate_ws(WKID);

	/* omit following four lines to invert black/white colors */
	rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
	rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
        gset_colr_rep(WKID,0,&rgb1);
        gset_colr_rep(WKID,1,&rgb2);
/*
 * X-axis label
 */
	c_anotat("I","Q",0,0,0,0);
/*
 * Set default Y plot bounds
 */
	c_agsetf("Y/MINIMUM.",-1.2);
	c_agsetf("Y/MAXIMUM.", 1.2);

/*
 * Set and plot the initial condition
 */
	ic(q1,dx,I1,I2,NX);
	stats(q1,I1,I2,NX,0,&qmax);
	plot1d(q1,I1,I2,NX,0,qmax,0,qtrue,"Initial condition",name);
/*
 * Copy the initial condition to the "qtrue" array.
 */
	for (i=I1; i<=I2; i++) qtrue[i]=q1[i];

/*
 * .. Integrate .....................................................
 */

	for (n=1; n<=nstep; n++) {

/*  . . . Set boundary conditions				*/
	  bc(q1,I1,I2,NX);

/*  . . . Compute values at next step				*/
	  advection(q1,q2,c,dx,dt,I1,I2,NX,advection_type);

/*  . . . Do array update at end of time step			*/
 	  update(q1,q2,I1,I2,NX);

/*  . . . Copy latest solution q2() to history array		*/
	  for (i=I1,j=0; i<=I2; i++,j++)
	    history[n-1][j] = q2[i];

/*  . . . Stats							*/
	  stats(q2,I1,I2,NX,n,&qmax);
	  qtrace[n-1] = qmax;

/*  . . . Plot fields when needed				*/
	  if (n == nstep) {	/* integration is done, showing final state */
	    if (advection_type == 'L') {
	      plot1d(q2,I1,I2,NX,n,qmax,1,qtrue,"Final solution",name);
	    } else {
	      plot1d(q2,I1,I2,NX,n,qmax,0,qtrue,"Final solution",name);
	    }
	  } else if (n%nplot == 0) {	/* intermediate plot before done */
	    plot1d(q2,I1,I2,NX,n,qmax,0,qtrue,"Solution",name);
	  }

/*  . . . Check if problem out of bounds			*/
	  if (qmax > 1.5) {
	    printf("Stop - solution blowing up at step %d\n",n);
	    plot1d(q2,I1,I2,NX,n,qmax,0,qtrue,"Solution blowing up",name);
	    nstep = n;
	    break;
	  }

	}	/* end of time loop n = 1,...,nstep */

/*
 * Run complete - do final plots
 */

/* . . Plot Qmax(t)						*/
	c_agsetf("Y/MINIMUM.",0.0);
	c_agsetf("Y/MAXIMUM.",1.5);
	c_anotat("N","Max value",0,0,0,0);
	plot1d(qtrace,0,nstep-1,nstep,0,qmax,0,qtrue,"Qmax vs. time",name);

/* . . Plot history surface s(x,t)				*/
	sfc(history,nstep,NX,(dt*(float)nstep),-90., 0.,
          "Time-Space evolution (view: -90,  0)",name);
	sfc(history,nstep,NX,(dt*(float)nstep),-75., 5.,
	  "Time-Space evolution (view: -75, +5)",name);
	sfc(history,nstep,NX,(dt*(float)nstep),-30.,20.,
	  "Time-Space evolution (view: -30,+20)",name);

/*
 *  Close the graphics package and stop.
 */

	gdeactivate_ws(WKID);
	gclose_ws(WKID);
	gclose_gks();

	exit;
}

