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

void diffusion(theta_d2,theta_d1,u3,u1,v3,v1,w3,w1,dx,dy,dz,tstep,dt,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH,K_u,K_v,K_w,K_theta)
int i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH;
float theta_d2[][ny][nz],theta_d1[][ny][nz],u3[][ny][nz],w3[][ny][k2+2],v3[][j2+2][nz],u1[][ny][nz],w1[][ny][k2+2],v1[][j2+2][nz],dx,dy,dz,tstep,dt;
float K_u,K_v,K_w,K_theta;
{
	int i,j,k;

#pragma omp parallel for shared(u3,u1) private(i,j,k)
	for (i=i1+1;i<=i2;i++)
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		u3[i][j][k] = u3[i][j][k] + tstep*K_u*(u1[i+1][j][k]-2*u1[i][j][k]+u1[i-1][j][k])/dx/dx;
		u3[i][j][k] = u3[i][j][k] + tstep*K_u*(u1[i][j+1][k]-2*u1[i][j][k]+u1[i][j-1][k])/dy/dy;
		u3[i][j][k] = u3[i][j][k] + tstep*K_u*(u1[i][j][k+1]-2*u1[i][j][k]+u1[i][j][k-1])/dz/dz;
	}
#pragma omp parallel for shared(v3,v1) private(i,j,k)
	for (i=i1;i<=i2;i++)
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		v3[i][j][k] = v3[i][j][k] + tstep*K_v*(v1[i+1][j][k]-2*v1[i][j][k]+v1[i-1][j][k])/dx/dx;
		v3[i][j][k] = v3[i][j][k] + tstep*K_v*(v1[i][j+1][k]-2*v1[i][j][k]+v1[i][j-1][k])/dy/dy;
		v3[i][j][k] = v3[i][j][k] + tstep*K_v*(v1[i][j][k+1]-2*v1[i][j][k]+v1[i][j][k-1])/dz/dz;
	}
#pragma omp parallel for shared(w3,w1) private(i,j,k)
	for (i=i1;i<=i2;i++) 
	for (j=j1;j<=j2;j++)
	for (k=k1+1;k<=k2;k++)
	{
		w3[i][j][k] = w3[i][j][k] + tstep*K_w*(w1[i+1][j][k]-2*w1[i][j][k]+w1[i-1][j][k])/dx/dx;
		w3[i][j][k] = w3[i][j][k] + tstep*K_w*(w1[i][j+1][k]-2*w1[i][j][k]+w1[i][j-1][k])/dy/dy;
		w3[i][j][k] = w3[i][j][k] + tstep*K_w*(w1[i][j][k+1]-2*w1[i][j][k]+w1[i][j][k-1])/dz/dz;
	} 
#pragma omp parallel for shared(theta_d2,theta_d1) private(i,j,k)
	for (i=i1;i<=i2;i++)    /* theta */
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		theta_d2[i][j][k] = theta_d2[i][j][k] + dt*K_theta*(theta_d1[i+1][j][k]-2*theta_d1[i][j][k]+theta_d1[i-1][j][k])/dx/dx;
		theta_d2[i][j][k] = theta_d2[i][j][k] + dt*K_theta*(theta_d1[i][j+1][k]-2*theta_d1[i][j][k]+theta_d1[i][j-1][k])/dy/dy;
		theta_d2[i][j][k] = theta_d2[i][j][k] + dt*K_theta*(theta_d1[i][j][k+1]-2*theta_d1[i][j][k]+theta_d1[i][j][k-1])/dz/dz;
	} 

	return;
}

