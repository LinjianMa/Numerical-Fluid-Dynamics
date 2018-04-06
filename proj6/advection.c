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

void advection(theta_d,p,u3,u2,v3,v2,w3,w2,dx,dy,dz,dt,tstep,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH)
int i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH;
float theta_d[][ny][nz],p[][ny][nz],u3[][ny][nz],u2[][ny][nz],v3[][j2+2][nz],v2[][j2+2][nz],w3[][ny][k2+2],w2[][ny][k2+2],dx,dy,dz,dt,tstep;
{
	int i,j,k;
	float theta1d_x[nx],theta1d_x2[nx],theta1d_y[ny],theta1d_y2[ny],theta1d_z[nz],theta1d_z2[nz];
	float u1d[nx],v1d[ny],w1d[nz],courant;

	#pragma omp parallel for shared(u3,u2,v2,w2) private(i,j,k)
	for (i=i1+1;i<=i2;i++)    
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		u3[i][j][k] = u3[i][j][k] - tstep*((u2[i-1][j][k]+u2[i][j][k])/2*(u2[i][j][k]-u2[i-1][j][k])/dx + (u2[i][j][k]+u2[i+1][j][k])/2*(u2[i+1][j][k]-u2[i][j][k])/dx)/2;
		u3[i][j][k] = u3[i][j][k] - tstep*((v2[i-1][j+1][k]+v2[i][j+1][k])/2*(u2[i][j+1][k]-u2[i][j][k])/dy + (v2[i-1][j][k]+v2[i][j][k])/2*(u2[i][j][k]-u2[i][j-1][k])/dy)/2;
		u3[i][j][k] = u3[i][j][k] - tstep*((w2[i-1][j][k+1]+w2[i][j][k+1])/2*(u2[i][j][k+1]-u2[i][j][k])/dz + (w2[i-1][j][k]+w2[i][j][k])/2*(u2[i][j][k]-u2[i][j][k-1])/dz)/2;
	}
	#pragma omp parallel for shared(v3,u2,v2,w2) private(i,j,k)
	for (i=i1;i<=i2;i++)  
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)
	{
		v3[i][j][k] = v3[i][j][k] - tstep*((u2[i][j-1][k]+u2[i][j][k])/2*(v2[i][j][k]-v2[i-1][j][k])/dx + (u2[i+1][j-1][k]+u2[i+1][j][k])/2*(v2[i+1][j][k]-v2[i][j][k])/dx)/2;
		v3[i][j][k] = v3[i][j][k] - tstep*((v2[i][j+1][k]+v2[i][j][k])/2*(v2[i][j+1][k]-v2[i][j][k])/dy + (v2[i][j][k]+v2[i][j-1][k])/2*(v2[i][j][k]-v2[i][j-1][k])/dy)/2;
		v3[i][j][k] = v3[i][j][k] - tstep*((w2[i][j-1][k]+w2[i][j][k])/2*(v2[i][j][k]-v2[i][j][k-1])/dz + (w2[i][j-1][k+1]+w2[i][j][k+1])/2*(v2[i][j][k+1]-v2[i][j][k])/dz)/2;
	} 
	#pragma omp parallel for shared(w3,u2,v2,w2) private(i,j,k)
	for (i=i1;i<=i2;i++)  
	for (j=j1;j<=j2;j++)
	for (k=k1+1;k<=k2;k++)
	{
		w3[i][j][k] = w3[i][j][k] - tstep*((u2[i][j][k-1]+u2[i][j][k])/2*(w2[i][j][k]-w2[i-1][j][k])/dx + (u2[i+1][j][k-1]+u2[i+1][j][k])/2*(w2[i+1][j][k]-w2[i][j][k])/dx)/2;
		w3[i][j][k] = w3[i][j][k] - tstep*((v2[i][j][k-1]+v2[i][j][k])/2*(w2[i][j][k]-w2[i][j-1][k])/dy + (v2[i][j+1][k-1]+v2[i][j+1][k])/2*(w2[i][j+1][k]-w2[i][j][k])/dy)/2;
		w3[i][j][k] = w3[i][j][k] - tstep*((w2[i][j][k+1]+w2[i][j][k])/2*(w2[i][j][k+1]-w2[i][j][k])/dz + (w2[i][j][k]+w2[i][j][k-1])/2*(w2[i][j][k]-w2[i][j][k-1])/dz)/2;
	} 
/* x array */
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)                           
	{
		for (i=i1-BC_WIDTH;i<=i2+BC_WIDTH;i++)
		{
			theta1d_x[i]=theta_d[i][j][k];
		}
		for (i=i1;i<=i2+1;i++)
		{
			u1d[i]=u2[i][j][k];
		}

		advect1d(theta1d_x,theta1d_x2,u1d,dx,dt/2,i1,i2,nx);

		for (i=i1;i<=i2;i++)
		{
			theta_d[i][j][k]=theta1d_x2[i];
		}
	}
	bc(theta_d,p,u3,v3,w3,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);
/* y array */
	for (i=i1;i<=i2;i++)
	for (k=k1;k<=k2;k++)                           
	{
		for (j=j1-BC_WIDTH;j<=j2+BC_WIDTH;j++)
		{
			theta1d_y[j]=theta_d[i][j][k];
		}
		for (j=j1;j<=j2+1;j++)
		{
			v1d[j]=v2[i][j][k];
		}

		advect1d(theta1d_y,theta1d_y2,v1d,dy,dt/2,j1,j2,ny);

		for (j=j1;j<=j2;j++)
		{
			theta_d[i][j][k]=theta1d_y2[j];
		}
	}
	bc(theta_d,p,u3,v3,w3,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);
/* z array */
	for (i=i1;i<=i2;i++)
	for (j=j1;j<=j2;j++)
	{
		for (k=k1-BC_WIDTH;k<=k2+BC_WIDTH;k++)
		{
			theta1d_z[k]=theta_d[i][j][k];
		}
		for (k=k1;k<=k2+1;k++)
		{
			w1d[k]=w2[i][j][k];
		}

		advect1d(theta1d_z,theta1d_z2,w1d,dz,dt,k1,k2,nz);

		for (k=k1;k<=k2;k++)
		{
			theta_d[i][j][k]=theta1d_z2[k];
		}
	}
	bc(theta_d,p,u3,v3,w3,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);
/* y array */
	for (i=i1;i<=i2;i++)
	for (k=k1;k<=k2;k++)                           
	{
		for (j=j1-BC_WIDTH;j<=j2+BC_WIDTH;j++)
		{
			theta1d_y[j]=theta_d[i][j][k];
		}
		for (j=j1;j<=j2+1;j++)
		{
			v1d[j]=v2[i][j][k];
		}

		advect1d(theta1d_y,theta1d_y2,v1d,dy,dt/2,j1,j2,ny);

		for (j=j1;j<=j2;j++)
		{
			theta_d[i][j][k]=theta1d_y2[j];
		}
	}
	bc(theta_d,p,u3,v3,w3,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);
/* x array */
	for (j=j1;j<=j2;j++)
	for (k=k1;k<=k2;k++)                           
	{
		for (i=i1-BC_WIDTH;i<=i2+BC_WIDTH;i++)
		{
			theta1d_x[i]=theta_d[i][j][k];
		}
		for (i=i1;i<=i2+1;i++)
		{
			u1d[i]=u2[i][j][k];
		}

		advect1d(theta1d_x,theta1d_x2,u1d,dx,dt/2,i1,i2,nx);

		for (i=i1;i<=i2;i++)
		{
			theta_d[i][j][k]=theta1d_x2[i];
		}
	}
	bc(theta_d,p,u3,v3,w3,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH);

	return;
}

