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

void advect1d(q1d,q1d_2,u1d,dx,dt,i1,i2,nx,method)
int i1,i2,nx;
float q1d[],q1d_2[],u1d[],dx,dt;
{
	int i;
	float courant;
	for (i=i1;i<=i2;i++)
	{
        courant=(u1d[i-i1]+u1d[i-i1+1])/2*dt/dx;
		if (method == 1)
		{
			q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+\
				courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1]);
		}
		if (method == 2)
		{
			q1d_2[i]=q1d[i]+courant/60*(q1d[i-3]-9*q1d[i-2]+45*q1d[i-1]-45*q1d[i+1]+9*q1d[i+2]-q1d[i+3])+\
				pow(courant,2.0)/360*(2*q1d[i-3]-27*q1d[i-2]+270*q1d[i-1]-490*q1d[i]+270*q1d[i+1]-27*q1d[i+2]+2*q1d[i+3])+\
				pow(courant,3.0)/48*(-q1d[i-3]+8*q1d[i-2]-13*q1d[i-1]+13*q1d[i+1]-8*q1d[i+2]+q1d[i+3])+\
				pow(courant,4.0)/144*(-q1d[i-3]+12*q1d[i-2]-39*q1d[i-1]+56*q1d[i]-39*q1d[i+1]+12*q1d[i+2]-q1d[i+3])+\
				pow(courant,5.0)/240*(q1d[i-3]-4*q1d[i-2]+5*q1d[i-1]-5*q1d[i+1]+4*q1d[i+2]-q1d[i+3])+\
				pow(courant,6.0)/720*(q1d[i-3]-6*q1d[i-2]+15*q1d[i-1]-20*q1d[i]+15*q1d[i+1]-6*q1d[i+2]+q1d[i+3]);
		}
		if (method == 3)
		{
			if (courant >= 0)
			{
				q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+\
					courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1])-\
					(1+courant)/6*courant*(courant-1)*(q1d[i+1]-3*q1d[i]+3*q1d[i-1]-q1d[i-2]);
			}
			if (courant < 0)
			{
				q1d_2[i]=q1d[i]-courant/2*(q1d[i+1]-q1d[i-1])+\
					courant*courant/2*(q1d[i+1]-2*q1d[i]+q1d[i-1])-\
					(1+fabs(courant))/6*courant*(courant+1)*(q1d[i-1]-3*q1d[i]+3*q1d[i+1]-q1d[i+2]);
			}
		}
	}
	return;
}

