/*
 * ========================= update ====================
 * Update: replace old values with new ones
 * We are NOT copying ghost points here.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1,q2	real arrays	old, new data arrays
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		size of arrays
 *
 */

void update(q1,q2,i1,i2,nx)
int i1,i2,nx;
float q1[],q2[];
{
	int i;

	for (i=i1; i<=i2; i++)
	  q1[i] = q2[i];

	return;
}

