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
#include <math.h>

void truncation(q1,u,v,dx,dy,i1,i2,j1,j2,nx,ny,BC_WIDTH,trunc_emax,trunc_xcenter,trunc_ycenter,trunc_x1,trunc_x2,trunc_y1,trunc_y2,dt)
int i1,i2,j1,j2,nx,ny,BC_WIDTH,*trunc_xcenter,*trunc_ycenter,*trunc_x1,*trunc_x2,*trunc_y1,*trunc_y2;
float q1[][ny],u[][ny-2*BC_WIDTH],v[][ny-2*BC_WIDTH+1],dx,dy,dt,*trunc_emax;
{
	int i,j;
	float q1d[ny],trunc_e[nx][ny],trunc_e1d[ny],u1d[ny],v1d[ny],courant;
    *trunc_emax = 0;
	*trunc_x1 = nx;
	*trunc_x2 = 0;
	*trunc_y1 = ny;
	*trunc_y2 = 0;

/* x array */
	for (j=j1+2;j<=j2-2;j++)                           
	{
		for (i=i1;i<=i2;i++)
		{
			q1d[i]=q1[i][j];
		}
		for (i=i1;i<=i2+1;i++)
		{
			u1d[i-i1]=u[i-i1][j-j1];
		}

		truncation1d(q1d,trunc_e1d,u1d,dx,i1,i2,nx,dt);

		for (i=i1+2;i<=i2-2;i++)
		{
			trunc_e[i][j]=trunc_e1d[i];
		}
	}

/* y array */
	for (i=i1+2;i<=i2-2;i++)                           
	{
		for (j=j1;j<=j2;j++)
		{
			q1d[j]=q1[i][j];
		}
		for (j=j1;j<=j2+1;j++)
		{
			v1d[j-j1]=v[i-i1][j-j1];
		}

		truncation1d(q1d,trunc_e1d,v1d,dy,j1,j2,ny,dt);

		for (j=j1+2;j<=j2-2;j++)
		{
			if ( trunc_e1d[j] > trunc_e[i][j] )
			{
				trunc_e[i][j] = trunc_e1d[j];
			}
		}
	}
/* error max, position */
	for (i=i1+2;i<=i2-2;i++)                           
	for (j=j1+2;j<=j2-2;j++)
	{
		if (trunc_e[i][j] > *trunc_emax)
		{
			*trunc_emax = trunc_e[i][j];
		}
	}
/* region calculation*/
	for (i=i1+2;i<=i2-2;i++)                           
	for (j=j1+2;j<=j2-2;j++)
	{
		if (trunc_e[i][j] >= *trunc_emax/2)
		{
			if ((i-i1+1) < *trunc_x1)
			{	
				*trunc_x1 = i-i1+1;
			}
			if ((i-i1+1) > *trunc_x2)
			{	
				*trunc_x2 = i-i1+1;
			}
			if ((j-j1+1) < *trunc_y1)
			{	
				*trunc_y1 = j-j1+1;
			}
			if ((j-j1+1) > *trunc_y2)
			{	
				*trunc_y2 = j-j1+1;
			}
		}
	}
	*trunc_xcenter = (*trunc_x1 + *trunc_x2)/2;
	*trunc_ycenter = (*trunc_y1 + *trunc_y2)/2;
	return;
}



