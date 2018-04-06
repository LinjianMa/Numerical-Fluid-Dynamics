/*
 * ======================= advection ====================
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

void advection(theta_d,u3,u2,w3,w2,dx,dz,dt,tstep,i1,i2,k1,k2,nx,nz,BC_WIDTH)
int i1,i2,k1,k2,nx,nz,BC_WIDTH;
float theta_d[][nz],u3[][nz],u2[][nz],w3[][k2+2],w2[][k2+2],dx,dz,dt,tstep;
{
	int i,k;
	float theta1d_x[nx],theta1d_x2[nx],theta1d_z[nz],theta1d_z2[nz];
	float u1d[nx],w1d[nz],courant;

	for (i=i1+1;i<=i2;i++)    
	for (k=k1;k<=k2;k++)
	{
		u3[i][k] = u3[i][k] - tstep*((u2[i-1][k]+u2[i][k])/2*(u2[i][k]-u2[i-1][k])/dx + (u2[i][k]+u2[i+1][k])/2*(u2[i+1][k]-u2[i][k])/dx)/2;
		u3[i][k] = u3[i][k] - tstep*((w2[i-1][k+1]+w2[i][k+1])/2*(u2[i][k+1]-u2[i][k])/dz + (w2[i-1][k]+w2[i][k])/2*(u2[i][k]-u2[i][k-1])/dz)/2;
	}

	for (i=i1;i<=i2;i++)    
	for (k=k1+1;k<=k2;k++)
	{
		w3[i][k] = w3[i][k] - tstep*((u2[i][k-1]+u2[i][k])/2*(w2[i][k]-w2[i-1][k])/dx + (u2[i+1][k-1]+u2[i+1][k])/2*(w2[i+1][k]-w2[i][k])/dx)/2;
		w3[i][k] = w3[i][k] - tstep*((w2[i][k+1]+w2[i][k])/2*(w2[i][k+1]-w2[i][k])/dz + (w2[i][k]+w2[i][k-1])/2*(w2[i][k]-w2[i][k-1])/dz)/2;
	} 
/* x array */
	for (k=k1;k<=k2;k++)                           
	{
		for (i=i1-BC_WIDTH;i<=i2+BC_WIDTH;i++)
		{
			theta1d_x[i]=theta_d[i][k];
		}
		for (i=i1;i<=i2+1;i++)
		{
			u1d[i]=u2[i][k];
		}

		advect1d(theta1d_x,theta1d_x2,u1d,dx,dt,i1,i2,nx);

		for (i=i1;i<=i2;i++)
		{
			theta_d[i][k]=theta1d_x2[i];
		}
	}
/* z array */
	for (i=i1;i<=i2;i++)                           
	{
		for (k=k1-BC_WIDTH;k<=k2+BC_WIDTH;k++)
		{
			theta1d_z[k]=theta_d[i][k];
		}
		for (k=k1;k<=k2+1;k++)
		{
			w1d[k]=w2[i][k];
		}

		advect1d(theta1d_z,theta1d_z2,w1d,dz,dt,k1,k2,nz);

		for (k=k1;k<=k2;k++)
		{
			theta_d[i][k]=theta1d_z2[k];
		}
	}
	return;
}

