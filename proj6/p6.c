/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm6
 *  >>>>> Linjian Ma <<<<<
 */
#include <time.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#define IWTYPE 1
#define WKID   1*/

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

#define NX 300
#define NY 300
#define NZ 75

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
float dy   = 50;
float dz   = 50;
float dt   = 0.1;
float tstep;
int   nstep;
int   nplot;
char *name = "Linjian Ma";

/* Variables for runtime timing */
time_t tstart,time1,time2,tlast[10],elapsed;
float estimate1,estimate2;
int min1,sec1,min2,sec2;

tstart = time(0);
/* Arrays and other variables */

float theta[NXDIM][NYDIM][NZDIM],theta_d2[NXDIM][NYDIM][NZDIM],theta_d1[NXDIM][NYDIM][NZDIM];
float p1[NXDIM][NYDIM][NZDIM],p2[NXDIM][NYDIM][NZDIM],p3[NXDIM][NYDIM][NZDIM];
float thetaplot[NX][NY][NZ],pplot[NX][NY][NZ];
float pplot_xz[NX][NZ],pplot_xy[NX][NY],pplot_yz[NY][NZ];
float u1[I2+2][NYDIM][NZDIM],u2[I2+2][NYDIM][NZDIM],u3[I2+2][NYDIM][NZDIM];
float v1[NXDIM][J2+2][NZDIM],v2[NXDIM][J2+2][NZDIM],v3[NXDIM][J2+2][NZDIM];
float w1[NXDIM][NYDIM][K2+2],w2[NXDIM][NYDIM][K2+2],w3[NXDIM][NYDIM][K2+2];
float uplot[NX][NY][NZ],vplot[NX][NY][NZ],wplot[NX][NY][NZ];
float rho[K2+2];

/*int delta_theta[2],xstart[2],zstart[2],xradius[2],zradius[2];*/

float K_u,K_v,K_w,K_theta;
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

	printf("Program #6       Numerical Fluid Dynamics, Spring 2016\n\n");

	printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",NX,BC_WIDTH,I1,I2,NXDIM);
	printf("NY=%d, BC_WIDTH=%d, J1=%d, J2=%d, NYDIM=%d\n",NY,BC_WIDTH,J1,J2,NYDIM);
	printf("NZ=%d, BC_WIDTH=%d, K1=%d, K2=%d, NZDIM=%d\n",NZ,BC_WIDTH,K1,K2,NZDIM);
	printf("dx=%8.5f, dy=%8.5f,dx=%8.5f, dt=%8.5f, nstep=%d, nplot=%d\n",dx,dy,dz,dt,nstep,nplot);

	printf(" Enter the number of steps \n");
	scanf("%d",&nstep);  
	printf(" Enter the plotting interval \n");
	scanf("%d",&nplot);  

	/*printf(" Enter the diffusivity of u \n");
	scanf("%d",&K_u);
	printf(" Enter the diffusivity of w \n");
	scanf("%d",&K_w);
	printf(" Enter the diffusivity of theta \n");
	scanf("%d",&K_theta);*/
	K_u = 40.0;
	K_v = 40.0;
	K_w = 40.0;
	K_theta = 5.0;

/* Definitions for routines */
	void ic(),stats(),plot1d(),bc(),integrate(),update(),advec1d(),advection(),PGF(),error(),diffusion(),putfield();

/* Set and plot the initial condition*/
    n = 0;
	ic(rho,theta_d1,u1,v1,w1,p1,dx,dy,dz,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,x0,y0,z0,BC_WIDTH);
	bc(theta_d1,p1,u1,v1,w1,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH);
	for (i=I1;i<=I2+1;i++)  /*u2 = u1*/
	for (j=J1-1;j<=J2+1;j++)
	for (k=K1-1;k<=K2+1;k++)
	{
		u2[i][j][k] = u1[i][j][k];
	}
	for (i=I1-1;i<=I2+1;i++)    /*v2 = v1*/
	for (j=J1-1;j<=J2+1;j++)
	for (k=K1-1;k<=K2+1;k++)
	{
		v2[i][j][k] = v1[i][j][k];
	}
	for (i=I1-1;i<=I2+1;i++)    /*w2 = w1*/
	for (j=J1-1;j<=J2+1;j++)
	for (k=K1;k<=K2+1;k++)
	{
		w2[i][j][k] = w1[i][j][k];
	}
	for (i=I1-1;i<=I2+1;i++)    /*p2 = p1*/
	for (j=J1-1;j<=J2+1;j++)
	for (k=K1-1;k<=K2+1;k++)
	{
		p2[i][j][k] = p1[i][j][k];
	}

/*  plot the initial condition  */
/* ... NOTE: pass field name in double quotes (string, not char expected) */
	for (i=I1;i<=I2;i++)   
	for (j=J1;j<=J2;j++)
	for (k=K1;k<=K2;k++)
	{
		thetaplot[i-I1][j-J1][k-K1] = theta_d1[i][j][k];
		pplot[i-I1][j-J1][k-K1] = p1[i][j][k];
		uplot[i-I1][j-J1][k-K1] = (u1[i][j][k]+u1[i+1][j][k])/2;
		vplot[i-I1][j-J1][k-K1] = (v1[i][j][k]+v1[i][j+1][k])/2;
		wplot[i-I1][j-J1][k-K1] = (w1[i][j][k]+w1[i][j][k+1])/2;
	}

	putfield("T",dt*(float)n,thetaplot,NX,NY,NZ);
	putfield("P",dt*(float)n,pplot,NX,NY,NZ);
	putfield("U",dt*(float)n,uplot,NX,NY,NZ);
	putfield("V",dt*(float)n,vplot,NX,NY,NZ);
	putfield("W",dt*(float)n,wplot,NX,NY,NZ);

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
		bc(theta_d1,p1,u1,v1,w1,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH);
		bc(theta_d1,p2,u2,v2,w2,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH);

	#pragma omp parallel for shared(p3,p1) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			p3[i][j][k] = p1[i][j][k];
		}
	#pragma omp parallel for shared(theta_d2,theta_d1) private(i,j,k)
		for (i=I1-2;i<=I2+2;i++)
		for (j=J1-2;j<=J2+2;j++)
		for (k=K1-2;k<=K2+2;k++)
		{
			theta_d2[i][j][k] = theta_d1[i][j][k];
		}
	#pragma omp parallel for shared(u3,u1) private(i,j,k)
		for (i=I1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			u3[i][j][k] = u1[i][j][k];
		}
	#pragma omp parallel for shared(v3,v1) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			v3[i][j][k] = v1[i][j][k];
		}
	#pragma omp parallel for shared(w3,w1) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1;k<=K2+1;k++)
		{
			w3[i][j][k] = w1[i][j][k];
		}
		
/* advect */
		advection(theta_d2,p3,u3,u2,v3,v2,w3,w2,dx,dy,dz,dt,tstep,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH);
/* diffusion */
		diffusion(theta_d2,theta_d1,u3,u1,v3,v1,w3,w1,dx,dy,dz,tstep,dt,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH,K_u,K_v,K_w,K_theta);
/*  PGF  */
		PGF(theta_d1,p3,u3,v3,w3,rho,dx,dy,dz,tstep,I1,I2,J1,J2,K1,K2,NXDIM,NYDIM,NZDIM,BC_WIDTH);

/*Return the values back*/
	#pragma omp parallel for shared(p1,p2,p3,theta_d1,theta_d2) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			p1[i][j][k] = p2[i][j][k];
			p2[i][j][k] = p3[i][j][k];
			theta_d1[i][j][k] = theta_d2[i][j][k];
		}
	#pragma omp parallel for shared(u1,u2,u3) private(i,j,k)
		for (i=I1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			u1[i][j][k] = u2[i][j][k];
			u2[i][j][k] = u3[i][j][k];
		}
	#pragma omp parallel for shared(v1,v2,v3) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1-1;k<=K2+1;k++)
		{
			v1[i][j][k] = v2[i][j][k];
			v2[i][j][k] = v3[i][j][k];
		}
	#pragma omp parallel for shared(w1,w2,w3) private(i,j,k)
		for (i=I1-1;i<=I2+1;i++)
		for (j=J1-1;j<=J2+1;j++)
		for (k=K1;k<=K2+1;k++)
		{
			w1[i][j][k] = w2[i][j][k];
			w2[i][j][k] = w3[i][j][k];
		}
/* plot */
		if ( (n % nplot)==0)
		{

		#pragma omp parallel for shared(thetaplot,pplot,uplot,vplot,wplot,theta_d2,p3,u3,v3,w3) private(i,j,k)
			for (i=I1;i<=I2;i++)   
			for (j=J1;j<=J2;j++)
			for (k=K1;k<=K2;k++)
			{
				thetaplot[i-I1][j-J1][k-K1] = theta_d2[i][j][k];
				pplot[i-I1][j-J1][k-K1] = p3[i][j][k];
				uplot[i-I1][j-J1][k-K1] = (u3[i][j][k]+u3[i+1][j][k])/2;
				vplot[i-I1][j-J1][k-K1] = (v3[i][j][k]+v3[i][j+1][k])/2;
				wplot[i-I1][j-J1][k-K1] = (w3[i][j][k]+w3[i][j][k+1])/2;
			}

			putfield("T",dt*(float)n,thetaplot,NX,NY,NZ);
			putfield("P",dt*(float)n,pplot,NX,NY,NZ);
			putfield("U",dt*(float)n,uplot,NX,NY,NZ);
			putfield("V",dt*(float)n,vplot,NX,NY,NZ);
			putfield("W",dt*(float)n,wplot,NX,NY,NZ);
		}
		/*
		* TIMING; report wallclock per-time-step based on
		* last 10 steps, and estimates to complete run.
		* first estimate is based on last 10 steps,
		* second on time since start of run
		*/
        /*time2 = time(0);
        if (n>10) {
			time1     = tlast[n%10];
            elapsed   = time2-time1;
            estimate1 = (float)(nstep-n) * ( (float)(time2-time1)/10.0 );
            min1 = (int) (estimate1 / 60.0);
            sec1 = (int) (estimate1 - 60.0*(float)min1 );
            estimate2 = (float)(nstep-n)  * ( (float)(time2-tstart) / (float)n );
            min2 = (int) (estimate2 / 60.0);
            sec2 = (int) (estimate2 - 60.0*(float)min2 );
            printf("Step %4d: time/step %5.2f sec, est. time remaining %03d:%02d / %03d:%02d\n",n,(float)elapsed/10.0,min1,sec1,min2,sec2);
        }
        tlast[n%10] = time2;*/
	}	/* end of time loop n = 1,...,nstep */

	exit;
}

