/*
 * ======================= Grid interpolation routine ======================
 *
 *  Routine dointerp - interpolation between coarse, nested grids
 *  ATMS 502/CSE 566   Brian Jewett, Univ. IL   Spring 2016
 *
 *  Arguments:
 * 	coarse	  real		in/out	coarse data array; NOTE: is (nxdim,nydim)
 * 	nest	  real		in/out	nest data array    NOTE: is (nxdim,nydim)
 * 	nx,ny     integer	in	dimensions of arrays without ghost points
 *	nxdim,	 integers	in	dimensions of coarse, nest
 *	 nydim
 *	i1,i2,	 integers	in	array bounds covering physical domain
 *	 j1,j2
 * 	nestX1,	  integers	in	nest start, end indices 
 * 	 nestX2,			  ... for I
 * 	 nestY1,
 * 	 nestY2				  ... for J
 * 	nstep	 integer	in	time step counter (for print messages)
 * 	ratio	 integer	in	nest refinement (time,space) ratio
 * 	flag	 integer	in	interpolation flag:
 * 					 -1=interp: coarse => nest, using
 * 					    old grid info where possible
 * 					    This is for a nest move.
 * 					  1=interp: coarse => nest
 * 					    This is for nest initialization.
 *					 10=interp: coarse => nest BCs
 *					    This sets nest BCs (boundary conditions)
 *					    only; nest interior is unchanged.
 * 					  2=interp: nest => coarse
 * 					    This is for feedback.
 * 	nestX1old, integers	in	If flag < 0, these are the OLD
 * 	 nestX2old,			  nested grid indices (coordinates
 * 	 nestY1old,			  on coarse grid) for interpolation.
 * 	 nestY1/2old  integers	in	  These are *ignored* if flag > 0.
 * 
 *  Internal:
 * 	VER	definition	local	if set = 1 below, verbose
 * 					(extra print statements done)
 */ 

#include <stdio.h>
#include <stdlib.h>

void dointerp(coarse,nest,nx,ny,nxdim,nydim,i1,i2,j1,j2,
          nestX1,nestX2,nestY1,nestY2,nstep,ratio,flag,
          nestX1old,nestX2old,nestY1old,nestY2old)
int nx,ny,nxdim,nydim,i1,i2,j1,j2,nestX1,nestX2,nestY1,nestY2,
    nstep,ratio,flag,nestX1old,nestX2old,nestY1old,nestY2old;
float coarse[][nydim],nest[][nydim];
{

#define VER 1

	int icoarse,jcoarse,bcwidth,ijnest,iNew1,iNew2,iOld1,iOld2,
      	  jNew1,jNew2,jOld1,jOld2,iold,jold,inew,jnew,inest,jnest,
          ncopy,istep;
	float *olddata;
	float ifraction,jfraction,offset,qBL,qBR,qTL,qTR,tmp1,tmp2;

/*
 *      If verbose option set, write out info now.
 */

	if (VER) {
	  printf("> Dointerp: flag %2d, nx/ny %d %d, nxdim/nydim %d %d, ij12 %d %d %d %d\n",
	    flag,nx,ny,nxdim,nydim,i1,i2,j1,j2);
	  if (flag == 1) {
	    printf("%s; step %3d, ratio %d; nest X %d:%d, Y %d:%d\n",
	      "  coarse to nest",nstep,ratio,nestX1,nestX2,nestY1,nestY2);
	  } else if (flag < 0) {
	    printf("%s; step %3d, ratio %d; from X,Y %d:%d %d:%d =to= %d:%d %d:%d\n",
	      "  coarse to nest w/nest copy",nstep,ratio,
	      nestX1old,nestX2old,nestY1old,nestY2old,nestX1,nestX2,nestY1,nestY2);
	  } else if (flag == 10) {
	    printf("%s; step %3d, ratio %d; nest X %d:%d, Y %d:%d\n",
	      "  setting nest BCs",nstep,ratio,nestX1,nestX2,nestY1,nestY2);
	  } else {
	    printf("%s; step %3d, ratio %d; nestX %d:%d, Y %d:%d\n",
	      "  nest to coarse",nstep,ratio,nestX1,nestX2,nestY1,nestY2);
	  }
	}
/*
 * 	Check for invalid arguments
 */
	if (nx < 10 || nx > 1000 || ny < 10 || ny > 1000) {
	  printf("ERROR: dointerp: bad nx or ny arguments: %d %d; Stop.\n",nx,ny);
	  exit(1);
	}
	if (nestX1 < i1 || nestX1 > i2  || 
      	    nestX2 < i1 || nestX2 > i2  || 
      	    nestY1 < j1 || nestY1 > j2  || 
      	    nestY2 < j1 || nestY2 > j2) {
	  printf("ERROR: dointerp: one or more invalid nest arguments.\n");
	  printf("i1,i2,j1,j2 = %d %d %d %d; nestX1,X2,Y1,Y2 = %d %d %d %d; Stop.\n",
	    i1,i2,j1,j2,nestX1,nestX2,nestY1,nestY2);
	  exit(1);
	}
	if (flag < 0) {
	  if (nestX1old < i1 || nestX1old > i2  || 
      	      nestX2old < i1 || nestX2old > i2  || 
      	      nestY1old < j1 || nestY1old > j2  || 
      	      nestY2old < j1 || nestY2old > j2) {
	  printf("ERROR: dointerp: one or more invalid nest arguments\n");
	  printf("i1,i2,j1,j2 = %d %d %d %d; nestX1old,X2old,Y1old,Y2old = %d %d %d %d; Stop.\n",
	    i1,i2,j1,j2,nestX1old,nestX2old,nestY1old,nestY2old);
	  exit(1);
	  }
	}
	if (ratio < 2 || ratio > 50) {
	  printf("ERROR: dointerp: bad nest refinement ratio = %d; Stop.\n",ratio);
	  exit(1);
	}
	if ((nx-1)%ratio != 0 ||
      	    (ny-1)%ratio != 0) {
	  printf("WARNING: dointerp: (nx-1) or (ny-1) are not even multiples of the grid nesting ratio: %d %d %d\n",
           nx,ny,ratio);
	}
	bcwidth = i1;
	if ( (i1 != j1) 			||
	     (nxdim != (nx+2*bcwidth)) 		||
	     (nydim != (ny+2*bcwidth)) ) {
	  printf("ERROR: dointerp finding inconsistent i1,j1,nx,ny, and/or nxdim,nydim.\n");
	  printf("i1,j1 = %d %d; nx,ny = %d %d; nxdim,nydim = %d %d\n",
	    i1,j1,nx,ny,nxdim,nydim);
	  exit(1);
	}

/*
 * ... FLAG < 0 or FLAG = 1 or 10:
 * ... Bilinear interpolation from coarse grid to nested grid
 * ... If flag<0, this is a grid move:  use old nested grid data
 * ... If flag=10, we set Only the nest boundary (border, here) points.
 */
	if (flag < 0 || flag == 1 || flag == 10) {
/*
 * . . . . If have old nest data, copy it to temporary space
 */
	  if (flag < 0) {
            olddata = (float *) malloc(nxdim*nydim*sizeof(float));
	    for (jnest=0; jnest<nydim; jnest++) {
	      for (inest=0; inest<nxdim; inest++) {
		olddata[inest + jnest*nydim] = nest[inest][jnest];
	      }
	    }
          }

/*
 * ....... Interpolate coarse grid to nested grid.
 * ....... nest() array is overwritten.
 */

/* . . . . Interior . . . . */

	  for (jnest=j1; jnest<j2; jnest++) {
	    jcoarse   = (jnest-j1)/ratio + nestY1;
	    jfraction = (float)((jnest-j1)%ratio) / (float)ratio;
	    if (flag == 10) {
	      if (jnest == j1) {
		istep = 1;
	      } else {
		istep = 9999;
	      }
	    } else {
	      istep = 1;
	    }
	    for (inest=i1; inest<i2; inest+=istep) {
	      icoarse   = (inest-i1)/ratio + nestX1;
	      ifraction = (float)((inest-i1)%ratio) / (float)ratio;
	      qBL  = coarse[icoarse  ][jcoarse  ];
	      qBR  = coarse[icoarse+1][jcoarse  ];
	      qTL  = coarse[icoarse  ][jcoarse+1];
	      qTR  = coarse[icoarse+1][jcoarse+1];
	      tmp1 = qBL + (qBR-qBL)*ifraction;
	      tmp2 = qTL + (qTR-qTL)*ifraction;
	      nest[inest][jnest] = tmp1 + (tmp2-tmp1)*jfraction;
	    }
	  }

/* . . . . Top edge . . . . */

	  for (inest=i1; inest<i2; inest++) {
	    icoarse   = (inest-i1)/ratio + nestX1;
	    ifraction = (float)((inest-i1)%ratio) / (float)ratio;
	    nest[inest][j2] = 
      	      coarse[icoarse][nestY2] + ifraction*
       	      (coarse[icoarse+1][nestY2]-coarse[icoarse][nestY2]);
	  }

/* . . . . Right edge . . . . */

	  for (jnest=j1; jnest<j2; jnest++) {
	    jcoarse   = (jnest-j1)/ratio + nestY1;
	    jfraction = ((jnest-j1)%ratio) / (float)ratio;
	    nest[i2][jnest] = 
      	      coarse[nestX2][jcoarse] + jfraction*
       	      (coarse[nestX2][jcoarse+1]-coarse[nestX2][jcoarse]);
	  }

/* . . . . Top corner . . . . */

	  nest[i2][j2] = coarse[nestX2][nestY2];
/*
 * . . . . If using old nest, copy overlapping data.
 */
	  if (flag < 0) {

/* . . . . . No overlap case . . . . */

	    if ( nestX2old < nestX1 || nestX1old > nestX2  || 
      	         nestY2old < nestY1 || nestY1old > nestY2) {
	      printf("Dointerp WARNING: *no overlap* between old, new nests.\n");
	      printf("Old nest location (%3d:%3d,%3d:%3d)\n",
		nestX1old,nestX2old,nestY1old,nestY2old);
	      printf("New nest location (%3d:%3d,%3d:%3d)\n",
		nestX1   ,nestX2   ,nestY1   ,nestY2);

	    } else {

/* . . . . . . X bounds for old, new nests . . . . */

	      if (nestX1 > nestX1old) {
		iOld1 = (nestX1-nestX1old)*ratio + i1;
		iOld2 = i2;
		iNew1 = i1;
		iNew2 = iNew1 + (iOld2-iOld1);
	      } else {
		iNew1 = (nestX1old-nestX1)*ratio + i1;
		iNew2 = i2;
		iOld1 = i1;
		iOld2 = iOld1 + (iNew2-iNew1);
	      }

/* . . . . . . Y bounds for old, new nests . . . . */

	      if (nestY1 > nestY1old) {
		jOld1 = (nestY1-nestY1old)*ratio + j1;
		jOld2 = j2;
		jNew1 = j1;
		jNew2 = jNew1 + (jOld2-jOld1);
	      } else {
		jNew1 = (nestY1old-nestY1)*ratio + j1;
		jNew2 = j2;
		jOld1 = j1;
		jOld2 = jOld1 + (jNew2-jNew1);
	      }

/* . . . . . . Summarize the planned nest data copy . . . . */

	      if (VER) {
	        ncopy = (iOld2-iOld1+1)*(jOld2-jOld1+1);
		printf("Dointerp: copying from old nest (%d:%d,%d:%d) to new (%d:%d,%d:%d); %5.1f %% overlap\n",
		  iOld1,iOld2,jOld1,jOld2,iNew1,iNew2,jNew1,jNew2,(100.*(float)ncopy/(float)(nx*ny)));
	      }

/*
 * . . . . . . Copy overlap data
 */

	      jnew = jNew1 - 1;
	      for (jold=jOld1; jold<=jOld2; jold++) {
		jnew = jnew + 1;
		inew = iNew1 - 1;
	        for (iold=iOld1; iold<=iOld2; iold++) {
		  inew = inew + 1;
		  nest[inew][jnew] = olddata[iold + jold*nydim];
		}
	      }

	    } /* end of non-zero area of nest to move */
	  } /* end of flag<0 case */
/*
 * 
 * ... FLAG = 2:
 * ... Interpolate from nested grid to coarse grid
 * ... Collocated point replacement here, rather
 * ... than averaging of neighboring nest values
 */

	} else if (flag == 2) {

	  for (jcoarse = (nestY1+1); jcoarse <= (nestY2-1); jcoarse++) {
	    jnest = (jcoarse-nestY1)*ratio + j1; /* WAS 1 */
	    for (icoarse = (nestX1+1); icoarse <= (nestX2-1); icoarse++) {
	      inest = (icoarse-nestX1)*ratio + i1; /* WAS 1 */
	      coarse[icoarse][jcoarse] = nest[inest][jnest];
	    }
	  }

/*
 * ... Otherwise: Error, unknown flag value.
 */

	} else {

	  printf("Interpolation flag=%d not recognized.\n",flag);
	  exit(1);

	}

/*
 * ... If we used temporary space, deallocate it now.
 */

	if (flag < 0) free(olddata);

/*
 * ... Done.
 */

	return;
}

/*
 * Spring 2016
 */
