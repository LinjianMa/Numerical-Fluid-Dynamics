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
#include <stdio.h>
#include <stdlib.h>

void bc(q1,i1,i2,j1,j2,nxdim,nydim,BC_WIDTH,mode,q1nest,nx,ny,nestX1,nestX2,nestY1,nestY2,n,ratio)
int i1,i2,j1,j2,nxdim,nydim,BC_WIDTH,nx,ny,mode,nestX1,nestX2,nestY1,nestY2,n,ratio;
/*mode=1 nest,mode=0 no gradient*/
float q1[][nydim],q1nest[][nydim];
{
	int i,j;
	if (mode == 0)
	{
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
	}
	if (mode == 1)
	{
		dointerp(q1,q1nest,nx,ny,nxdim,nydim,i1,i2,j1,j2,nestX1,nestX2,nestY1,nestY2,\
				n,ratio,10);
	}

	return;
}

