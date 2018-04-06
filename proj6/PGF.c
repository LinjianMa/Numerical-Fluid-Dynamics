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

void PGF(theta_d,p,u,v,w,rho,dx,dy,dz,tstep,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH)
int i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH;
float rho[],p[][ny][nz],theta_d[][ny][nz],u[][ny][nz],v[][j2+2][nz],w[][ny][k2+2],dx,dy,dz,tstep;
{
	int i,j,k;
	float theta_0 = 300.0;
	float g = 9.81;
	float Cs = 60;

#pragma omp parallel for shared(u,p) private(i,j,k)
	for (i=i1+1;i<=i2;i++)    /* u */
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		u[i][j][k] = u[i][j][k] - tstep/rho[k]*(p[i][j][k]-p[i-1][j][k])/dx;
	}
#pragma omp parallel for shared(v,p) private(i,j,k)
	for (i=i1;i<=i2;i++)    /* v */
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		v[i][j][k] = v[i][j][k] - tstep/rho[k]*(p[i][j][k]-p[i][j-1][k])/dy;
	}
#pragma omp parallel for shared(w,p,theta_d) private(i,j,k)
	for (i=i1;i<=i2;i++)    /* w */
	for (j=j1;j<=j2;j++)
	for (k=k1+1;k<=k2;k++)
	{
		w[i][j][k] = w[i][j][k] - tstep*2.0/(rho[k-1]+rho[k])*(p[i][j][k]-p[i][j][k-1])/dz;
		w[i][j][k] = w[i][j][k] + tstep*g*(theta_d[i][j][k]/theta_0+theta_d[i][j][k-1]/theta_0)/2;
	} 
	bc(theta_d,p,u,v,w,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);
#pragma omp parallel for shared(p,u,v,w) private(i,j,k)	
	for (i=i1;i<=i2;i++)   
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		p[i][j][k] = p[i][j][k] - tstep*Cs*Cs*rho[k]*(u[i+1][j][k]-u[i][j][k])/dx;
		p[i][j][k] = p[i][j][k] - tstep*Cs*Cs*rho[k]*(v[i][j+1][k]-v[i][j][k])/dy;
		p[i][j][k] = p[i][j][k] - tstep*Cs*Cs*((rho[k]+rho[k+1])/2*w[i][j][k+1]-(rho[k-1]+rho[k])/2*w[i][j][k])/dz;
	}

	return;
}

