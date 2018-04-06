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

void ic(q,u,v,dx,dy,i1,i2,j1,j2,nx,ny,x0,y0,BC_WIDTH)
int i1,i2,j1,j2,nx,ny,BC_WIDTH;
float dx,dy,q[][ny],u[][ny-2*BC_WIDTH],v[][ny-2*BC_WIDTH+1],x0,y0;
{
	#ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
	#define xstart 0.0
	#define ystart -0.05
	int i,j;
	float x[nx-2*BC_WIDTH+1],y[ny-2*BC_WIDTH+1],d[nx][ny];
    float r=0.165;
	for (i=0;i<=i2-i1+1;i++ )  /*x*/
	{
		x[i]=x0+dx*(i)-dx/2;
	}
	for (j=0;j<=j2-j1+1;j++ )  /*y*/
	{
		y[j]=y0+dy*(j)-dy/2;
	}

	for (i=0;i<=i2-i1+1;i++ )  /*u*/
		for (j=0;j<j2-j1+1;j++ )
	{
    u[i][j]=sin(5*M_PI*x[i])*sin(5*M_PI*(y[j]+dy/2));
	}
	for (j=0;j<=j2-j1+1;j++ )  /*v*/
		for (i=0;i<i2-i1+1;i++ )  
	{
    v[i][j]=cos(5*M_PI*(x[i]+dx/2))*cos(5*M_PI*y[j]);
	}

	for (i=i1;i<=i2;i++)    /*q*/
	{
		for (j=j1;j<=j2;j++)
	{
    /*d[i][j]=sqrt(pow(x[i-i1]+dx/2-xstart,2.0)+pow(y[j-j1]+dy/2-ystart,2.0));
	if (d[i][j]<r)  {q[i][j]=5*(1+cos(M_PI*d[i][j]/r));}*/
	if (((x[i-i1]+dx/2)>=-0.020 && (x[i-i1]+dx/2)<=0.020 && (y[j-j1]+dy/2)>=-0.25 && (y[j-j1]+dy/2)<=0.25) ||\
		((x[i-i1]+dx/2)>=-0.30 && (x[i-i1]+dx/2)<=0.30 && (y[j-j1]+dy/2)>=-0.25 && (y[j-j1]+dy/2)<=-0.22) ||\
		((x[i-i1]+dx/2)>=-0.30 && (x[i-i1]+dx/2)<=0.30 && (y[j-j1]+dy/2)>=0.22 && (y[j-j1]+dy/2)<=0.25))
	{
		q[i][j]=10.0;
	}
	else {q[i][j]=0;}
	}
	}
	return;
}

