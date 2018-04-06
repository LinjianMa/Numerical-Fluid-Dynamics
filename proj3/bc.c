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

void bc(q1,i1,i2,j1,j2,nx,ny,BC_WIDTH)
int i1,i2,j1,j2,nx,ny,BC_WIDTH;
float q1[][ny];
{
	int i,j;
	/*printf(" >Set boundary conditions here.\n");*/
	for (j=1;j<=BC_WIDTH;j++)
	for (i=i1;i<=i2;i++)
	{
		q1[i][j1-j]=q1[i][j1];
		q1[i][j2+j]=q1[i][j2];
	}
	for (i=1;i<=BC_WIDTH;i++)
	for (j=j1;j<=j2;j++)
	{
		q1[i1-i][j]=q1[i1][j];
		q1[i2+i][j]=q1[i2][j];
	}
	return;
}

