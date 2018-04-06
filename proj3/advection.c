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

void advection(q1,u,v,dx,dy,dt,i1,i2,j1,j2,nx,ny,method,BC_WIDTH)
int i1,i2,j1,j2,nx,ny,BC_WIDTH;
float q1[][ny],u[][ny-2*BC_WIDTH],v[][ny-2*BC_WIDTH+1],dx,dy,dt;
{
	int i,j;
	float q1d[ny],q1d_2[ny],u1d[ny],v1d[ny],courant;

/* Set boundary conditions */
	bc(q1,i1,i2,j1,j2,nx,ny,BC_WIDTH);
/* x array */
	for (j=j1;j<=j2;j++)                           
	{
		for (i=i1-BC_WIDTH;i<=i2+BC_WIDTH;i++)
		{
			q1d[i]=q1[i][j];
		}
		for (i=i1;i<=i2+1;i++)
		{
			u1d[i-i1]=u[i-i1][j-j1];
		}

		advect1d(q1d,q1d_2,u1d,dx,dt,i1,i2,nx,method);

		for (i=i1;i<=i2;i++)
		{
			q1[i][j]=q1d_2[i];
		}
	}
/* Set boundary conditions */
	bc(q1,i1,i2,j1,j2,nx,ny,BC_WIDTH);
/* y array */
	for (i=i1;i<=i2;i++)                           
	{
		for (j=j1-BC_WIDTH;j<=j2+BC_WIDTH;j++)
		{
			q1d[j]=q1[i][j];
		}
		for (j=j1;j<=j2+1;j++)
		{
			v1d[j-j1]=v[i-i1][j-j1];
		}

		advect1d(q1d,q1d_2,v1d,dy,dt,j1,j2,ny,method);

		for (j=j1;j<=j2;j++)
		{
			q1[i][j]=q1d_2[j];
		}
	}
	return;
}

