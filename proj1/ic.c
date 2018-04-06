/*
 * ============================ ic =====================
 * IC sets the initial condition
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	IC data. Set 1..nx here;
 *				  [0],[nx+1] = ghost zones
 *				  if 1 ghost point on each side
 *	dx	real		grid spacing
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *
 */

void ic(q1,dx,i1,i2,nx)
int i1,i2,nx;
float dx,q1[];
{
	int i,j;
	float x,pi,length;

	pi = 4.0*atan(1.0);
 	length = dx * (float) nx;

	for (i=i1,j=1; i<=i2; i++,j++) {
	  x = dx * (float)(j-1);
	  q1[i] = sin( 2.*pi/length*x );
	}

	return;
}

