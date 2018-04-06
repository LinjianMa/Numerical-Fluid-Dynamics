
/* -------------------------------------------------------
 * PUTFIELD - write a data field to the run history file.
 * -------------------------------------------------------
 * Arguments:
 *
 *   name	input   char *	  name of the field (only 4 chars used)
 *   datatime   input   float	  model integration time
 *   field	input	float	  array of data
 *   nx,ny,nz	input	int	  dimensions of field()
 *
 * If name = "*", the output file is closed (nothing else done)
 *
 * ==NOTE== make sure you pass a character string, e.g. "NAME"
 *	    rather than just a character, even if that string
 *	    is 1-char long.  A string is assumed as input.
 *
 * A header integer written to the file is set to 0 if the
 *   file is created from Fortran, and 1 from C.
 *
 * ATMS 502 / CSE 566, Spring 2016
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
void putfield(name,datatime,field,nx,ny,nz)
char *name;
float datatime,*field;
int nx,ny,nz;
{
	int i,sourcetype=1;
	int creat();
	static int fd;
	static int count=0;
	char tname[5];

/*
 * ... If name = '*' close file and return
 */
	if (name[0] == '*') {
	  printf("Output history file closed.\n");
	  close(fd);
	  count = -1;
	  return;
	}
	if (count < 0) {
	  printf("putfield: error: file has already been closed.\n");
	  exit(1);
	}

/*
 * ... First call: open file, write file data type as C.
 */
	if (count == 0) {
	  printf("Writing unformatted data file RunHistory.dat\n");
	  fd=creat("RunHistory.dat",0644);
	  if (fd<0) {
	    printf("putfield - error creating file. Stop.\n");
	    exit(1);
	  }
	  write(fd,&sourcetype,sizeof(int));	/* type=1 is flag for C data */
	}
	count++;

	strncpy(tname,name,4); tname[4]='\0';
	printf("Writing field %4d: %s(%3d,%3d,%3d) for T=%6.1f\n",
	  count,tname,nx,ny,nz,datatime);

/* ... Write header  */
	write(fd,&count,sizeof(int));
	write(fd,tname,4);
	write(fd,&datatime,sizeof(float));
	write(fd,&nx,sizeof(int));
	write(fd,&ny,sizeof(int));
	write(fd,&nz,sizeof(int));

/* ... Write array */
	write(fd,field,(nx*ny*nz)*sizeof(float));

	return;
}


