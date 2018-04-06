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

void bc(theta_d,p,u,v,w,i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH)
int i1,i2,j1,j2,k1,k2,nx,ny,nz,BC_WIDTH;
float theta_d[][ny][nz],p[][ny][nz],u[][ny][nz],v[][j2+2][nz],w[][ny][k2+2];
{
	int i,j,k;

	/*** U:  X (Symmetry ... but asymmetry for U) Boundaries ***/
#pragma omp parallel for shared(u) private(j,k)
    for (j=i1; j<=j2; j++)  
    for (k=k1; k<=k2; k++) 
	{ 
        u[i1][j][k] = -u[i1+1][j][k];
        u[i2+1][j][k] = -u[i2][j][k];
	}
    /*** U:  Z (0-gradient) Boundaries ***/
#pragma omp parallel for shared(u) private(i,j)
    for (i=i1; i<=i2+1; i++)
    for (j=j1; j<=j2; j++) 
	{
        u[i][j][k1-1] = u[i][j][k1];
        u[i][j][k2+1] = u[i][j][k2];
    } 
	/*** U:  Y (Periodic) Boundaries ***/
#pragma omp parallel for shared(u) private(i,k)
    for (i=i1; i<=i2+1; i++)
    for (k=k1; k<=k2; k++) 
	{
        u[i][j1-1][k] = u[i][j2][k];
        u[i][j2+1][k] = u[i][j1][k];
    }  

    /*** W:  Z (Rigid upper/lower lid) Boundaries ***/
#pragma omp parallel for shared(w) private(i,j)
	for (i=i1; i<=i2; i++)
    for (j=j1; j<=j2; j++) 
	{
        w[i][j][k1]   = 0;
        w[i][j][k2+1] = 0;
    }
    /*** W:  X (Symmetry) Boundaries ***/
#pragma omp parallel for shared(w) private(j,k)
    for (j=j1; j<=j2; j++)
    for (k=k1; k<=k2+1; k++) 
	{
        w[i1-1][j][k] = w[i1+1][j][k];
        w[i2+1][j][k] = w[i2-1][j][k];					
    }
    /*** W:  Y (Periodic) Boundaries ***/
#pragma omp parallel for shared(w) private(i,k)
	for (i=i1; i<=i2; i++)
    for (k=k1; k<=k2+1; k++) 
	{
        w[i][j1-1][k] = w[i][j2][k];
        w[i][j2+1][k] = w[i][j1][k];
    }

    /*** P, theta:  Z (0-gradient) Boundaries ***/
#pragma omp parallel for shared(p,theta_d) private(i,j)
    for (i=i1; i<=i2; i++)
    for (j=j1; j<=j2; j++)
	{
        p[i][j][k1-1] = p[i][j][k1];
        p[i][j][k2+1] = p[i][j][k2];
		theta_d[i][j][k1-1] = theta_d[i][j][k1];
        theta_d[i][j][k2+1] = theta_d[i][j][k2];
		theta_d[i][j][k1-2] = theta_d[i][j][k1];
        theta_d[i][j][k2+2] = theta_d[i][j][k2];
    }
    /*** P, theta: X (Symmetry) Boundaries ***/
#pragma omp parallel for shared(p,theta_d) private(j,k)
    for (j=j1; j<=j2; j++) 
    for (k=k1; k<=k2; k++) 
	{
        p[i1-1][j][k] = p[i1+1][j][k];
        p[i2+1][j][k] = p[i2-1][j][k];
		theta_d[i1-1][j][k] = theta_d[i1+1][j][k];
        theta_d[i2+1][j][k] = theta_d[i2-1][j][k];
		theta_d[i1-2][j][k] = theta_d[i1+2][j][k];
        theta_d[i2+2][j][k] = theta_d[i2-2][j][k];
    }
    /*** P, theta: Y (Periodic) Boundaries ***/
#pragma omp parallel for shared(p,theta_d) private(i,k)
    for (i=i1; i<=i2; i++)
    for (k=k1; k<=k2; k++) 
	{
        p[i][j1-1][k] = p[i][j2][k];
        p[i][j2+1][k] = p[i][j1][k];
        theta_d[i][j1-1][k] = theta_d[i][j2][k];
        theta_d[i][j2+1][k] = theta_d[i][j1][k];
        theta_d[i][j1-2][k] = theta_d[i][j2-1][k];
        theta_d[i][j2+2][k] = theta_d[i][j1+1][k];
    }

    /*** V: Z (0-Gradient) Boundaries ***/
#pragma omp parallel for shared(v) private(i,j)
    for (i=i1; i<=i2; i++)
    for (j=j1; j<=j2+1; j++) 
	{
        v[i][j][k1-1] = v[i][j][k1];
        v[i][j][k2+1] = v[i][j][k2];
    }
    /*** V: X (Symmetry) Boundaries ***/
#pragma omp parallel for shared(v) private(j,k)
    for (j=j1; j<=j2+1; j++)
    for (k=k1; k<=k2; k++) 
	{
        v[i1-1][j][k] = v[i1+1][j][k];
        v[i2+1][j][k] = v[i2-1][j][k];
    }
    /*** V: Y (Periodic) Boundaries ***/
#pragma omp parallel for shared(v) private(i,k)
    for (i=i1; i<=i2; i++)
    for (k=k1; k<=k2; k++)
	{
        v[i][j1-1][k] = v[i][j2][k];
        v[i][j2+1][k] = v[i][j1][k];
    }

	return;
}

