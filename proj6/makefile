# Makefile for parallel tests       ATMS 502/CSE 566, Spring 2016
#
# Lines starting with "#" are comments.
# Type "man make" for more information, or try
#
#          http://www.math.tau.ac.il/~danha/courses/software1/make-intro.html
#
# To compile this parallel test, just type "make".
#

CC	= icc
OPTIONS	= -O

default:
	@echo '> To build the p6 code, type  :  make p6'
	@echo '> To build the p6 parallel code, type  :  make parallel'
	@echo '> To remove the compiled programs from here:  make clean'

p6:
	@echo 'Compiling ...'
	$(CC) $(OPTIONS) -w p6.c advection.c bc.c ic.c stats.c advect1d.c PGF.c error.c diffusion.c putfield.c -o p6.1CPU
	@echo ' '
	@echo 'Ready to try.  From your login session, run it by typing:  p6.1CPU'

parallel:
	@echo 'Compiling parallel version ... '
	$(CC) $(OPTIONS) -openmp p6.c advection.c bc.c ic.c stats.c advect1d.c PGF.c error.c diffusion.c putfield.c -o p6.OMP
	@echo ' '
	@echo 'Ready to try.  From your login session, run it by typing:  p6.OMP'

clean:
	rm -f p6
