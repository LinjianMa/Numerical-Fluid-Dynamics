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

void diffusion(theta_d2,theta_d1,u3,u1,w3,w1,dx,dz,tstep,i1,i2,k1,k2,nx,nz,BC_WIDTH,K_u,K_w,K_theta)
int i1,i2,k1,k2,nx,nz,BC_WIDTH;
float theta_d2[][nz],theta_d1[][nz],u3[][nz],w3[][k2+2],u1[][nz],w1[][k2+2],dx,dz,tstep;
int K_u,K_w,K_theta;
{
	int i,k;

	for (i=i1+1;i<=i2;i++)    /* u */
	for (k=k1;k<=k2;k++)
	{
		u3[i][k] = u3[i][k] + tstep*K_u*(u1[i+1][k]-2*u1[i][k]+u1[i-1][k])/dx/dx;
		u3[i][k] = u3[i][k] + tstep*K_u*(u1[i][k+1]-2*u1[i][k]+u1[i][k-1])/dz/dz;
	}

	for (i=i1;i<=i2;i++)    /* w */
	for (k=k1+1;k<=k2;k++)
	{
		w3[i][k] = w3[i][k] + tstep*K_w*(w1[i+1][k]-2*w1[i][k]+w1[i-1][k])/dx/dx;
		w3[i][k] = w3[i][k] + tstep*K_w*(w1[i][k+1]-2*w1[i][k]+w1[i][k-1])/dz/dz;
	} 
	
	for (i=i1;i<=i2;i++)    /* theta */
	for (k=k1;k<=k2;k++)
	{
		theta_d2[i][k] = theta_d2[i][k] + tstep*K_theta*(theta_d1[i+1][k]-2*theta_d1[i][k]+theta_d1[i-1][k])/dx/dx;
		theta_d2[i][k] = theta_d2[i][k] + tstep*K_theta*(theta_d1[i][k+1]-2*theta_d1[i][k]+theta_d1[i][k-1])/dz/dz;
	} 

	return;
}

