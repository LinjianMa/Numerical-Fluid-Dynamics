/*
 * ============================ bc =====================
 * BC sets the boundary conditions
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		main array size, not including
 *				extra 'ghost' zones/points
 */
#include <stdio.h>
#include <stdlib.h>

void bc(theta_d,p,u,w,i1,i2,k1,k2,nx,nz,BC_WIDTH)
int i1,i2,k1,k2,nx,nz,BC_WIDTH;
float theta_d[][nz],p[][nz],u[][nz],w[][k2+2];
{
	int i,k;
	for (i=i1;i<=i2;i++)   /*up and bottom bc */
	{
		u[i][k2+1] = u[i][k2];
		u[i][k1-1] = u[i][k1];

		theta_d[i][k2+1] = theta_d[i][k2];
		theta_d[i][k1-1] = theta_d[i][k1];
		theta_d[i][k2+2] = theta_d[i][k2];
		theta_d[i][k1-2] = theta_d[i][k1];

		p[i][k2+1] = p[i][k2];
		p[i][k1-1] = p[i][k1];

		w[i][k2+1] = 0;
		w[i][k1] = 0;
	}
	for (k=k1;k<=k2;k++)    /*left and right bc*/
	{
		u[i1][k] = -u[i1+1][k];
		u[i2+1][k] = -u[i2][k];

		theta_d[i1-1][k] = theta_d[i1+1][k];
		theta_d[i2+1][k] = theta_d[i2-1][k];
		theta_d[i1-2][k] = theta_d[i1+2][k];
		theta_d[i2+2][k] = theta_d[i2-2][k];

		p[i1-1][k] = p[i1+1][k];
		p[i2+1][k] = p[i2-1][k];

		w[i1-1][k] = w[i1+1][k];
		w[i2+1][k] = w[i2-1][k];
	}
	return;
}

