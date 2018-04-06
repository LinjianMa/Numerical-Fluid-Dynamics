/*
 * ============================ ic =====================
 * IC sets the initial condition
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	IC data. Set 1..nx here;
 *				  [0],[nx+1] = ghost zones
 *				  if 1 ghost point on each side
 *	dx	real		grid spacing
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ic(rho,theta_d,u,w,p,dx,dz,i1,i2,k1,k2,nx,nz,x0,z0,BC_WIDTH)
int i1,i2,k1,k2,nx,nz,BC_WIDTH;
float dx,dz,rho[],u[][nz],w[][k2+2],p[][nz],theta_d[][nz],x0,z0;
{
	#ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
	int i,k,m;
	float x[i2+1],z[k2+1],d[i2+1][k2+1];
	float theta_0 = 300;
	float g = 9.81;
	float cp = 1004;
	float Rd = 287;
	float P0 = 100000;
	float T,P,rm;
	int delta_theta[2],xstart[2],zstart[2],xradius[2],zradius[2];

	printf(" Enter the first temperature perturbation value \n");
	scanf("%d",&delta_theta[0]); 
	printf(" Enter the first temperature perturbation x coordinate \n");
	scanf("%d",&xstart[0]); 
	printf(" Enter the first temperature perturbation z coordinate \n");
	scanf("%d",&zstart[0]); 
	printf(" Enter the first temperature perturbation x radius \n");
	scanf("%d",&xradius[0]); 
	printf(" Enter the first temperature perturbation z radius \n");
	scanf("%d",&zradius[0]);   

	printf(" Enter the second temperature perturbation value \n");
	scanf("%d",&delta_theta[1]); 
	printf(" Enter the second temperature perturbation x coordinate \n");
	scanf("%d",&xstart[1]); 
	printf(" Enter the second temperature perturbation z coordinate \n");
	scanf("%d",&zstart[1]); 
	printf(" Enter the second temperature perturbation x radius \n");
	scanf("%d",&xradius[1]); 
	printf(" Enter the second temperature perturbation z radius \n");
	scanf("%d",&zradius[1]);
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
	for (i=i1;i<=i2+1;i++)  /*u*/
	for (k=k1;k<=k2;k++)
	{
		/*if (k <= (k2-k1)/2)
		{
			u[i][k] = 20;
		}
		if (k >= (k2-k1)/2+1)
		{
			u[i][k] = -20;
		}*/
		u[i][k] = 0;
		
	}
	for (i=i1;i<=i2;i++)   /*w*/
	for (k=k1;k<=k2+1;k++)
	{
		/*if (i <= (i2-i1)/2)
		{
			w[i][k] = 20;
		}
		if (i >= (i2-i1)/2+1)
		{
			w[i][k] = -20;
		}*/
		w[i][k] = 0;
	}
	for (i=i1;i<=i2;i++)   /*p*/
	for (k=k1;k<=k2;k++)
	{
		p[i][k] = 0;
	}

	for (i=i1;i<=i2;i++)  
	for (k=k1;k<=k2;k++)
	{
		x[i] = dx/2+dx*(i-i1);
		z[k] = dz/2+dz*(k-k1);
		theta_d[i][k] = 0;

		for (m=0;m<2;m++)
		{
			rm = sqrt(pow((x[i]-xstart[m])/xradius[m],2.0)+pow((z[k]-zstart[m])/zradius[m],2.0));
			if (rm <= 1.0)
			{
				theta_d[i][k] = theta_d[i][k] + delta_theta[m]/2.0*(cos(rm*M_PI)+1);
			}
		}
	}

	for (k=k1;k<=k2;k++ )  /*rho*/
	{
		T = 300.0-g/cp*z[k];
		P = P0*pow(T/theta_0,cp/Rd);
		rho[k] = P/Rd/T;
	}
	return;
}

