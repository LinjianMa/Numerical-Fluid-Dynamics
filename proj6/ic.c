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

void ic(rho,theta_d,u,v,w,p,dx,dy,dz,i1,i2,j1,j2,k1,k2,nx,ny,nz,x0,y0,z0,BC_WIDTH)
int i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH;
float dx,dy,dz,rho[],u[][ny][nz],v[][j2+2][nz],w[][ny][k2+2],p[][ny][nz],theta_d[][ny][nz],x0,y0,z0;
{
	#ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
	int i,j,k,m;
	float x[nx],y[ny],z[nz],d[nx][ny][nz];
	float theta_0 = 300;
	float g = 9.81;
	float cp = 1004;
	float Rd = 287;
	float P0 = 100000;
	float T,P,rm;
	int delta_theta[2],delta_v[2],xstart[2],ystart[2],zstart[2],xradius[2],yradius[2],zradius[2];
	float upertur = 2.0;
	srand(0.0); 

	/*printf(" Enter the first temperature perturbation value \n");
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
	scanf("%d",&zradius[1]);*/
	delta_theta[0] = -25;
	delta_theta[1] = -25;
	delta_v[0] = -40;
	delta_v[1] = 40;
	xstart[0] = 25;
	xstart[1] = 14975;
	ystart[0] = 7525;
	ystart[1] = 7525;
	zstart[0] = 1525;
	zstart[1] = 1525;
	xradius[0] = 3500;
	xradius[1] = 3500;
	yradius[0] = 999999;
	yradius[1] = 999999;
	zradius[0] = 1750;
	zradius[1] = 1750;
#pragma omp parallel for shared(u) private(i,j,k)
	for (i=i1;i<=i2+1;i++)  /*u*/
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		u[i][j][k] = 0;	
	}
	for (i=i1+2;i<=i2-1;i++)  /*u*/
	for (j=j1+1;j<=j2-1;j++)
	for (k=k1+1;k<=k2-1;k++)
	{
		u[i][j][k] = u[i][j][k] + upertur * ( rand() / (RAND_MAX + 1.0) ) - upertur/2.0;	
	}
#pragma omp parallel for shared(v) private(i,j,k)
	for (i=i1;i<=i2;i++)  /*v*/
	for (j=j1-1;j<=j2+1;j++)
	for (k=k1;k<=k2;k++)
	{
		v[i][j][k] = 0;	
	}
#pragma omp parallel for shared(w) private(i,j,k)
	for (i=i1;i<=i2;i++)    /*w*/
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2+1;k++)
	{
		w[i][j][k] = 0;
	}
#pragma omp parallel for shared(p) private(i,j,k)
	for (i=i1;i<=i2;i++)   /*p*/
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		p[i][j][k] = 0;
	}

	for (i=i1;i<=i2;i++)  
	for (j=j1;j<=j2;j++) 
	for (k=k1;k<=k2;k++)
	{
		x[i] = dx/2+dx*(i-i1);
		y[j] = dy/2+dy*(j-j1);
		z[k] = dz/2+dz*(k-k1);
		theta_d[i][j][k] = 0;

		for (m=0;m<2;m++)
		{
			rm = sqrt(pow((x[i]-xstart[m])/xradius[m],2.0)+pow((y[j]-ystart[m])/yradius[m],2.0)+pow((z[k]-zstart[m])/zradius[m],2.0));
			if (rm <= 1.0)
			{
				theta_d[i][j][k] = theta_d[i][j][k] + delta_theta[m]/2.0*(cos(rm*M_PI)+1);
			}
		}
	}

	for (i=i1;i<=i2;i++)  
	for (j=j1;j<=j2+1;j++) 
	for (k=k1;k<=k2;k++)
	{
		x[i] = dx/2+dx*(i-i1);
		y[j] = dy/2+dy*(j-j1);
		z[k] = dz/2+dz*(k-k1);

		for (m=0;m<2;m++)
		{
			rm = sqrt(pow((x[i]-xstart[m])/xradius[m],2.0)+pow((y[j]-ystart[m])/yradius[m],2.0)+pow((z[k]-zstart[m])/zradius[m],2.0));
			if (rm <= 1.0)
			{
				v[i][j][k] = v[i][j][k] + delta_v[m]/2.0*(cos(rm*M_PI)+1);
			}
		}
	}

	for (k=k1-1;k<=k2+1;k++ )  /*rho*/
	{
		z[k] = dz/2+dz*(k-k1);
		T = 300.0-g/cp*z[k];
		P = P0*pow(T/theta_0,cp/Rd);
		rho[k] = P/Rd/T;
	}
	return;
}

