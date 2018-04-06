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

void stats(q2,i1,i2,j1,j2,q_av,qmax,qmin,n,nx,ny)
int i1,i2,j1,j2,n,nx,ny;
float q2[][ny+2],qmax[],qmin[],q_av[];
{
	int i,j;
	float q_max,q_min,q_average;
    q_max=q2[i1][j1];
    q_min=q2[i1][j1];
	q_average=0;
    for (i=i1 ;i<=i2 ;i++ )
	for (j=j1 ;j<=j2 ;j++ )
    {
	if (q2[i][j]>q_max) q_max=q2[i][j];
	if (q2[i][j]<q_min) q_min=q2[i][j];
	q_average=q_average+q2[i][j];
    }
    q_average=q_average/(i2-i1+1)/(j2-j1+1);

	qmax[n] = q_max;
    qmin[n] = q_min;
	q_av[n] = q_average;
	printf("Step %3d, Max = %8.5f\n",n,qmax[n]);
	printf("Step %3d, Min = %8.5f\n",n,qmin[n]);
	printf("Step %3d, Average = %8.5f\n",n,q_av[n]);   
	return;
}

