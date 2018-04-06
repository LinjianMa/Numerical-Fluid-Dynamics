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

void advection(q1,q2,c,dx,dt,i1,i2,nx,advection_type)
int i1,i2,nx;
char advection_type;
float q1[],q2[],c,dx,dt;
{
	int i;
	float courant;
    courant=c*dt/dx;
	if (advection_type == 'L') {
     printf(" >Put advection code here for linear advection.\n");
	  for (i=i1; i<=i2; i++)
	  q2[i]=q1[i]-courant/2*(q1[i+1]-q1[i-1])+courant*courant/2*(q1[i+1]-2*q1[i]+q1[i-1]);
	} 
	else if (advection_type == 'N') {
			  printf(" >Put advection code here for nonlinear advection.\n");
	  for (i=i1; i<=i2; i++)
	  q2[i]=q1[i]-q1[i]*dt/dx/2*(q1[i+1]-q1[i-1])+q1[i]*dt/dx*q1[i]*dt/dx/2*(q1[i+1]-2*q1[i]+q1[i-1]);
	} else {
	  printf("Advection: Error, unrecognized advection type '%c'\n",advection_type);
	  exit(1);
	}

	return;
}

