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

void advect1d(q1d,q1d_2,u1d,dx,dt,i1,i2,nx,method,flag)
int i1,i2,nx;
float q1d[],q1d_2[],u1d[],dx,dt;
{
	int i;
	float courant;
	if (flag == 1)
	{
		for (i=i1;i<=i2;i++)
		{
			courant=(u1d[i-i1]+u1d[i-i1+1])/2*dt/dx;
			if (method == 1)
			{
				q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+\
					courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1]);
			}
		}
	}
	if (flag == 0)
	{
		for (i=i1+1;i<=i2-1;i++)
		{
			courant=(u1d[i-i1]+u1d[i-i1+1])/2*dt/dx;
			if (method == 1)
			{
				q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+\
					courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1]);
			}
		}
	}
	return;
}

