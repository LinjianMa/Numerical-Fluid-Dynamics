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
#include <stdio.h>
#include <stdlib.h>

void ic(q,u,v,dx,dy,i1,i2,j1,j2,nx,ny,x0,y0,xstart,ystart,r)
int i1,i2,j1,j2,nx,ny;
float dx,dy,q[][ny+2],u[][ny],v[nx][ny+1],x0,y0,xstart,ystart,r;
{
	#ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
	int i,j;
	float x[nx],y[ny],d[nx][ny];
	for (i=0;i<i2 ;i++ )  /*x*/
	{
		x[i]=x0+dx*(i);
	}
	for (j=0;j<j2 ;j++ )  /*y*/
	{
		y[j]=y0+dy*(j);
	}

	for (i=0;i<=i2;i++ )  /*u*/
		for (j=0;j<j2;j++ )
	{
    u[i][j]=-2*y[j];
	}
	for (j=0;j<=j2 ;j++ )  /*v*/
		for (i=0;i<i2 ;i++ )  
	{
    v[i][j]=2*x[i];
	}

	for (i=i1;i<=i2 ;i++ )    /*q*/
	{
		for (j=j1;j<=j2 ;j++ )
	{
    d[i][j]=sqrt(pow(x[i-i1]-xstart,2.0)+pow(y[j-j1]-ystart,2.0));
	if (d[i][j]<r)  {q[i][j]=5*(1+cos(M_PI*d[i][j]/r));}
	else {q[i][j]=0;}/*;}*/
	}
	}
	return;
}

