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

void bc(q1,i1,i2,nx)
int i1,i2,nx;
float q1[];
{
	printf(" >Set boundary conditions here.\n");
    q1[i1-1]=q1[i2];
	q1[i2+1]=q1[i1];
	return;
}

