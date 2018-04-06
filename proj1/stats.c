/*
 * ========================= stats =====================
 * Stats computes and prints out the max Q values
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q2	real array	Latest data. Check 1..nx;
 *				  [0],[nx+1] = ghost zones
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		size of data array
 *	n	integer		time step counter
 *	qmax	real		holds max absolute
 *				  value of q2
 */

void stats(q2,i1,i2,nx,n,qmax)
int i1,i2,nx,n;
float q2[],*qmax;
{
	int i;
	float qtmp;

	qtmp = fabs(q2[i1]);
	for (i=i1+1; i<=i2; i++) {
	  if (fabs(q2[i]) > qtmp) qtmp = fabs(q2[i]);
	}

	printf("Step %3d, Max = %8.5f\n",n,qtmp);
	*qmax = qtmp;

	return;
}

