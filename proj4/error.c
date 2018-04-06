/*
 * ======================= advection ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
 *	c	real		true speed of wave
 *	dx	real		grid spacing
 *	dt	real		time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */
#include <stdio.h>
#include <stdlib.h>

void error(q1,qtrue,i1,i2,j1,j2,nx,ny,n,q_av,qtrue_av,cor_max,err_t,err_dissip,err_disper)
int i1,i2,j1,j2,nx,ny,n;
float q1[][ny],qtrue[][ny],q_av,qtrue_av,*cor_max,*err_t,*err_dissip,*err_disper;
{
	int i,j;
	float cor,err_sip,err_per;
	float upterm=0.0;
	float downterm_d=0.0;
	float downterm_T=0.0;
	float sigma_d=0.0;
	float sigma_T=0.0;
	float err=0.0;
	for (i=i1;i<=i2;i++)
	for (j=j1;j<=j2;j++)
	{
		upterm=upterm+(q1[i][j]-q_av)*(qtrue[i][j]-qtrue_av);
		downterm_d=downterm_d+(q1[i][j]-q_av)*(q1[i][j]-q_av);
		downterm_T=downterm_T+(qtrue[i][j]-qtrue_av)*(qtrue[i][j]-qtrue_av);
		err=err+(q1[i][j]-qtrue[i][j])*(q1[i][j]-qtrue[i][j]);
	}
	sigma_d=sqrt(downterm_d/(i2-i1+1)/(j2-j1+1));
	sigma_T=sqrt(downterm_T/(i2-i1+1)/(j2-j1+1));

	cor=upterm/(sqrt(downterm_d*downterm_T));
	err=err/(i2-i1+1)/(j2-j1+1);
	err_sip=(sigma_d-sigma_T)*(sigma_d-sigma_T)+(q_av-qtrue_av)*(q_av-qtrue_av);
	err_per=2*(1-cor)*sigma_d*sigma_T;

	*cor_max=upterm/(sqrt(downterm_d*downterm_T));
	*err_t=err;
	*err_dissip=err_sip;
	*err_disper=err_per;

	return;
}

