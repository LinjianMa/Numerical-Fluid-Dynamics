/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm5:  Linear and nonlinear advection
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

#define NX 401
#define NY 2
#define NZ 71

#define BC_WIDTH 2
#define I1 BC_WIDTH
#define I2 I1+NX-1
#define J1 BC_WIDTH
#define J2 J1+NY-1
#define K1 BC_WIDTH
#define K2 K1+NZ-1

#define NXDIM NX+2*BC_WIDTH
#define NYDIM NY+2*BC_WIDTH
#define NZDIM NZ+2*BC_WIDTH
#define x0  0
#define y0  0
#define z0  0
#define MAXSTEP 1000

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
float dx   = 50;
float dy   = 1.0/(NY-1);
float dz   = 50;
float dt   = 0.2;
float tstep;
int   nstep;
int   nplot;
char *name = "Linjian Ma";

/* Arrays and other variables */

float theta[NXDIM][NZDIM],theta_d2[NXDIM][NZDIM],theta_d1[NXDIM][NZDIM];
float p1[NXDIM][NZDIM],p2[NXDIM][NZDIM],p3[NXDIM][NZDIM];
float thetaplot[NX][NZ],pplot[NX][NZ];
float u1[I2+2][NZDIM],u2[I2+2][NZDIM],u3[I2+2][NZDIM];
float w1[NXDIM][K2+2],w2[NXDIM][K2+2],w3[NXDIM][K2+2];
float uplot[NX+1][NZ],wplot[NX][NZ+1];
float rho[K2+1];
/*int delta_theta[2],xstart[2],zstart[2],xradius[2],zradius[2];*/
int K_u,K_w,K_theta;
float q_av[MAXSTEP],qmin[MAXSTEP],qmax[MAXSTEP],courant,qtrue_av,qtrue_max,qtrue_min;
int i,j,k,n;
int method = 1;
char label[200];
float cor_max,err_t,err_dissip,err_disper;
float cint    =  2.0;	
int colors    =  0;	
int pltzero   =  -1;	
float angh    =  45.0;
float angv    =  20.0;
/* Parameters and input .................................... */

	printf("Program #5       Numerical Fluid Dynamics, Spring 2016\n\n");

	printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",NX,BC_WIDTH,I1,I2,NXDIM);
	printf("NY=%d, BC_WIDTH=%d, J1=%d, J2=%d, NYDIM=%d\n",NY,BC_WIDTH,J1,J2,NYDIM);
	printf("NZ=%d, BC_WIDTH=%d, K1=%d, K2=%d, NZDIM=%d\n",NZ,BC_WIDTH,K1,K2,NZDIM);
	printf("dx=%8.5f, dy=%8.5f,dx=%8.5f, dt=%8.5f, nstep=%d, nplot=%d\n",dx,dy,dz,dt,nstep,nplot);

	printf(" Enter the number of steps \n");
	scanf("%d",&nstep);  
	printf(" Enter the plotting interval \n");
	scanf("%d",&nplot);  

	printf(" Enter the diffusivity of u \n");
	scanf("%d",&K_u);
	printf(" Enter the diffusivity of w \n");
	scanf("%d",&K_w);
	printf(" Enter the diffusivity of theta \n");
	scanf("%d",&K_theta);
	/*K_u = 50;
	K_w = 50;
	K_theta = 40;*/
	/*delta_theta[0] = -15;
	delta_theta[1] = 0;
	xstart[0] = 10050;
	xstart[1] = 20050;
	zstart[0] = 2050;
	zstart[1] = 2050;
	xradius[0] = 4000;
	xradius[1] = 4000;
	zradius[0] = 1000;
	zradius[1] = 1000;*/
/* Variables to reverse default black/white colors in NCAR Graphics */

	Gcolr_rep rgb1,rgb2;

/* Definitions for routines */

	void ic(),stats(),contr(),plot1d(),sfc(),bc(),integrate(),update(),advec1d(),advection(),PGF(),error(),diffusion();

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
	ic(rho,theta_d1,u1,w1,p1,dx,dz,I1,I2,K1,K2,NXDIM,NZDIM,x0,z0,BC_WIDTH);
	for (i=I1;i<=I2+1;i++)  /*u2 = u1*/
	for (k=K1;k<=K2;k++)
	{
		u2[i][k] = u1[i][k];
	}
	for (i=I1;i<=I2;i++)   /*w2 = w1*/
	for (k=K1;k<=K2+1;k++)
	{
		w2[i][k] = w1[i][k];
	}
	for (i=I1;i<=I2;i++)   /*p2 = p1*/
	for (k=K1;k<=K2;k++)
	{
		p2[i][k] = p1[i][k];
	}

/*  plot the initial condition  */
	/* plot p */
	for ( k=K1; k<=K2; k++ ) 
	for ( i=I1; i<=I2; i++ ) 
	{
		pplot[i-I1][k-K1] = p1[i][k];
	}
	printf("Plotting contours.\n");
	sprintf(label,"ploting p,   n = %d",n);
	contr(pplot,NX,NZ,25.0,0.0,label,colors,pltzero,-1,0,0,0,name); 
	/* plot theta */
	for ( k=K1; k<=K2; k++) 
    for ( i=I1; i<=I2; i++) 
	{
		thetaplot[i-I1][k-K1] = theta_d1[i][k];
    }
	printf("Plotting contours.\n");
	sprintf(label,"ploting theta,   n = %d",n);
	contr(thetaplot,NX,NZ,1.0,0.0,label,colors,pltzero,-1,0,0,0,name); 
	/* plot u */
	for ( k=K1; k<=K2; k++ ) 
	for ( i=I1; i<=I2+1; i++ ) 
	{
		uplot[i-I1][k-K1] = u1[i][k];
	}
	printf("Plotting contours.\n");
	sprintf(label,"ploting u,   n = %d",n);
	contr(uplot,NX+1,NZ,2.0,0.0,label,colors,pltzero,-1,0,0,0,name); 
	/* plot w */
	for ( k=K1; k<=K2+1; k++ ) 
	for ( i=I1; i<=I2; i++ ) 
	{
		wplot[i-I1][k-K1] = w1[i][k];
	}
	printf("Plotting contours.\n");
	sprintf(label,"ploting w,   n = %d",n);
	contr(wplot,NX,NZ+1,2.0,0.0,label,colors,pltzero,-1,0,0,0,name); 

/* Integrate */
 	for (n=1; n<=nstep; n++) {
		if (n == 1)
		{
			tstep = dt;
		}
		else 
		{
			tstep = 2*dt;
		}
		bc(theta_d1,p1,u1,w1,I1,I2,K1,K2,NXDIM,NZDIM,BC_WIDTH);
		for (i=I1-1;i<=I2+1;i++)
		for (k=K1-1;k<=K2+1;k++)
		{
			p3[i][k] = p1[i][k];
			theta_d2[i][k] = theta_d1[i][k];
		}
		for (i=I1;i<=I2+1;i++)
		for (k=K1-1;k<=K2+1;k++)
		{
			u3[i][k] = u1[i][k];
		}
		for (i=I1-1;i<=I2+1;i++)
		for (k=K1;k<=K2+1;k++)
		{
			w3[i][k] = w1[i][k];
		}
		
/* advect */
		advection(theta_d2,u3,u2,w3,w2,dx,dz,dt,tstep,I1,I2,K1,K2,NXDIM,NZDIM,BC_WIDTH);
/* diffusion */
		diffusion(theta_d2,theta_d1,u3,u1,w3,w1,dx,dz,tstep,I1,I2,K1,K2,NXDIM,NZDIM,BC_WIDTH,K_u,K_w,K_theta);
/*  PGF  */
		PGF(theta_d1,p3,u3,w3,rho,dx,dz,tstep,I1,I2,K1,K2,NXDIM,NZDIM,BC_WIDTH);

/*Return the values back*/
		for (i=I1-1;i<=I2+1;i++)
		for (k=K1-1;k<=K2+1;k++)
		{
			p1[i][k] = p2[i][k];
			p2[i][k] = p3[i][k];
			theta_d1[i][k] = theta_d2[i][k];
		}
		for (i=I1;i<=I2+1;i++)
		for (k=K1-1;k<=K2+1;k++)
		{
			u1[i][k] = u2[i][k];
			u2[i][k] = u3[i][k];
		}
		for (i=I1-1;i<=I2+1;i++)
		for (k=K1;k<=K2+1;k++)
		{
			w1[i][k] = w2[i][k];
			w2[i][k] = w3[i][k];
		}
/* plot */
		if ( (n % nplot)==0)
		{
			/* plot p */
			for ( k=K1; k<=K2; k++ ) 
			for ( i=I1; i<=I2; i++ ) 
			{
				pplot[i-I1][k-K1] = p3[i][k];
			}
			printf("Plotting contours.\n");
			sprintf(label,"ploting p,   n = %d",n);
			contr(pplot,NX,NZ,25.0,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
			/* plot theta_d */
			for ( k=K1; k<=K2; k++ ) 
			for ( i=I1; i<=I2; i++ ) 
			{
				thetaplot[i-I1][k-K1] = theta_d2[i][k];
			}
			printf("Plotting contours.\n");
			sprintf(label,"ploting theta,   n = %d",n);
			contr(thetaplot,NX,NZ,1.0,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
			/* plot u */
			for ( k=K1; k<=K2; k++ ) 
			for ( i=I1; i<=I2+1; i++ ) 
			{
				uplot[i-I1][k-K1] = u3[i][k];
			}
			printf("Plotting contours.\n");
			sprintf(label,"ploting u,   n = %d",n);
			contr(uplot,NX+1,NZ,2.0,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
			/* plot w */
			for ( k=K1; k<=K2+1; k++ ) 
			for ( i=I1; i<=I2; i++ ) 
			{
				wplot[i-I1][k-K1] = w3[i][k];
			}
			printf("Plotting contours.\n");
			sprintf(label,"ploting w,   n = %d",n);
			contr(wplot,NX,NZ+1,2.0,(dt*(float)n),label,colors,pltzero,-1,0,0,0,name); 
		}

	}	/* end of time loop n = 1,...,nstep */

/*  Close the graphics package and stop.*/
	gdeactivate_ws(WKID);
	gclose_ws(WKID);
	gclose_gks();

	exit;
}

