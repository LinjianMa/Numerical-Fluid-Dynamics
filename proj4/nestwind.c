/*
 *
 * ============================= nestwind ===========================
 *
 * ATMS 502/CSE 566   Brian Jewett, Univ. IL   Spring 2016
 * Set the nested grid wind fields
 *
 * Arguments:
 *	unest	real array	U field on nest, dimension (nx+1,ny)
 *	vnest	real array	V field on nest, dimension (nx,ny+1)
 *	nx,ny	integers	Used for u,v dimensions
 *	nestX1,	integers	Bottom left corner of nest,
 *	 nestY1			 in coarse grid coordinates
 *	dx	real		COARSE grid spacing
 *	ratio	integer		Nest refinement factor
 *	icase	integer		Flow: 1=rotation, 2=deformation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void nestwind(unest,vnest,nx,ny,nestX1,nestY1,dx,ratio,icase)
int nx,ny,nestX1,nestY1,ratio,icase;
float unest[nx+1][ny],vnest[nx][ny+1],dx;
{
	int i,j;
	float x0,y0,dxnest,x,y,k,pi;
/*
 * ... coordinates of bottom left corner of nest
 */
	dxnest = dx/(float)ratio;
	x0 = -0.5 + dx*(float)(nestX1-1);
	y0 = -0.5 + dx*(float)(nestY1-1);
/*
 * ... case 1: rotational flow
 */
	if (icase == 1) {

	  for (j=0; j<ny; j++) {
	    y = y0 + dxnest*(float)(j);
	    for (i=0; i<=nx; i++) {
	      x = x0 + dxnest*(float)(i) - 0.5*dxnest;
	      unest[i][j] = -2.*y;
	    }
	  }

	  for (j=0; j<=ny; j++) {
	    y = y0 + dxnest*(float)(j) - 0.5*dxnest;
	    for (i=0; i<nx; i++) {
	      x = x0 + dxnest*(float)(i);
	      vnest[i][j] = 2.*x;
	    }
	  }
/*
 * ... case 2: deformational flow
 */
	} else if (icase == 2) {

	  pi = 4.*atan(1.0);
	  k  = 6.0*pi;

	  for (j=0; j<ny; j++) {
	    y = y0 + dxnest*(float)(j);
	    for (i=0; i<=nx; i++) {
	      x = x0 + dxnest*(float)(i) - 0.5*dxnest;
	      unest[i][j] = sin(k*(x+.5)) * sin(k*(y+.5));
	    }
	  }

	  for (j=0; j<=ny; j++) {
	    y = y0 + dxnest*(float)(j) - 0.5*dxnest;
	    for (i=0; i<nx; i++) {
	      x = x0 + dxnest*(float)(i);
	      vnest[i][j] = cos(k*(x+.5)) * cos(k*(y+.5));
	    }
	  }
/*
 * ... check for invalid icase value
 */
	} else {
	  printf("nestwind: invalid icase=%d\n",icase);
	  exit(1);
	}

	return;
}

/*
 * Spring 2016
 */
