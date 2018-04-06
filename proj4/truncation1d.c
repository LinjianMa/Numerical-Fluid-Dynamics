
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

void truncation1d(q1d,trunc_e1d,u1d,dx,i1,i2,nx,dt)
int i1,i2,nx;
float q1d[],trunc_e1d[],u1d[],dx,dt;
{
	int i;
	float courant,q_xxx;
	for (i=i1+2;i<=i2-2;i++)
	{
        courant=(u1d[i-i1]+u1d[i-i1+1])/2*dt/dx;
        q_xxx = (q1d[i+2]-2*q1d[i+1]+2*q1d[i-1]-q1d[i-2])/(2*dx*dx*dx);
		trunc_e1d[i] = fabs(q_xxx*(courant*courant-1.0)*(u1d[i-i1]+u1d[i-i1+1])/2*dx*dx/6);
	}
	return;
}

