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

void advect1d(q1d,q1d_2,u1d,dx,dt,i1,i2,nx)
int i1,i2,nx;
float q1d[],q1d_2[],u1d[],dx,dt;
{
	int i;
	float courant,F2,F1,dq;
	/* piecewise linear */
	for (i=i1;i<=i2;i++)
	{        
		if (u1d[i] >= 0)
		{
			courant = u1d[i]*dt/dx;	
			dq = (q1d[i]-q1d[i-2])/2;
			F1 = courant*(q1d[i-1]+(1-courant)/2*dq);
		}
		if (u1d[i] < 0)
		{
			courant = -u1d[i]*dt/dx;	
			dq = (q1d[i+1]-q1d[i-1])/2;
			F1 = courant*(-q1d[i]+(1-courant)/2*dq);
		}
        
		if (u1d[i+1] >= 0)
		{
			courant = u1d[i+1]*dt/dx;
			dq = (q1d[i+1]-q1d[i-1])/2;
			F2 = courant*(q1d[i]+(1-courant)/2*dq);
		}
		if (u1d[i+1] < 0)
		{
			courant = -u1d[i+1]*dt/dx;
			dq = (q1d[i+2]-q1d[i])/2;
			F2 = courant*(-q1d[i+1]+(1-courant)/2*dq);
		}
		q1d_2[i]=q1d[i]-(F2-F1)+dt/dx*q1d[i]*(u1d[i+1]-u1d[i]);
	}
	/* lax wendroff */
	/*for (i=i1;i<=i2;i++)
	{
        courant=(u1d[i]+u1d[i+1])/2*dt/dx;
		q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1]);
	}*/
	return;
}

