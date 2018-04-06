/*
 * ======================= PGF ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
 *	c	real		true speed of wave
 *	dx	real		grid spacing
 *	dt	real		time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */
#include <stdio.h>
#include <stdlib.h>

void PGF(theta_d,p,u,w,rho,dx,dz,tstep,i1,i2,k1,k2,nx,nz,BC_WIDTH)
int i1,i2,k1,k2,nx,nz,BC_WIDTH;
float rho[],p[][nz],theta_d[][nz],u[][nz],w[][k2+2],dx,dz,tstep;
{
	int i,k;
	float p1d[nx],p1d_2[nz],u1d[i2+2],u1d_2[nz],w1d[nx],w1d_2[k2+2],courant;
	float theta_0 = 300.0;
	float g = 9.81;
	float Cs = 80;

	for (i=i1+1;i<=i2;i++)    /* u */
	for (k=k1;k<=k2;k++)
	{
		u[i][k] = u[i][k] - tstep/rho[k]*(p[i][k]-p[i-1][k])/dx;
	}

	for (i=i1;i<=i2;i++)    /* w */
	for (k=k1+1;k<=k2;k++)
	{
		w[i][k] = w[i][k] - tstep*2.0/(rho[k-1]+rho[k])*(p[i][k]-p[i][k-1])/dz;
		w[i][k] = w[i][k] + tstep*g*(theta_d[i][k]/theta_0+theta_d[i][k-1]/theta_0)/2;
	} 
	
	for (i=i1;i<=i2;i++)    /* p */
	for (k=k1;k<=k2;k++)
	{
		p[i][k] = p[i][k] - tstep*Cs*Cs*rho[k]*(u[i+1][k]-u[i][k])/dx;
		p[i][k] = p[i][k] - tstep*Cs*Cs*((rho[k]+rho[k+1])/2*w[i][k+1]-(rho[k-1]+rho[k])/2*w[i][k])/dz;
	} 

	return;
}

