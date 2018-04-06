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

void bc(q1,i1,i2,j1,j2,nx,ny)
int i1,i2,j1,j2,nx,ny;
float q1[][ny];
{
	int i,j;
	/*printf(" >Set boundary conditions here.\n");*/
	for (i=i1;i<=i2;i++)
	{
		q1[i][j1-1]=q1[i][j1];
		q1[i][j2+1]=q1[i][j2];
	}
	for (j=j1;j<=j2;j++)
	{
		q1[i1-1][j]=q1[i1][j];
		q1[i2+1][j]=q1[i2][j];
	}
	return;
}

