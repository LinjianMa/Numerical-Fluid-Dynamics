# Makefile for program #1             ATMS 502/CSE 566, Spring 2016
#
# Lines starting with "#" are comments.
# First line is default, so typing "make" makes executable named p2.
# Typing "make test_interp" compiles and creates executable test_interp.
#
# The executable has its name, then all dependencies (object files).
#   Beneath that is the statement to link them and create the program.
#
# The last statements say how to turn .f90 or .c files into .o files: compiling.
#
# Makefiles can get confusing, but they are quite useful in compiling programs.
# Type "man make" for more information, or try one of these help pages:
#
#          http://kiwi.atmos.colostate.edu/fortran/docs/fortran90-mar30.pdf
#	   http://genius2k.is-programmer.com/posts/40301.html
#	   https://www.sharcnet.ca/help/index.php/Make_utility_and_makefiles
#	   http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
#

CC	= ncargcc
OPTIONS	= -O
# use following instead, for debugging
#OPTIONS = -g -debug extended -traceback
# use following instead, for extensive debugging
#OPTIONS = -g -debug extended -traceback -Wuninitialized -Wcheck -check=uninit
OBJECTS = p2.o bc.o ic.o advection.o sfc.o stats.o contr.o advect1d.o error.o
PROGRAM	= p2
ARCHIVE = pgm2.tar

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make listing .... to make a printable listing ... you may want to adjust the page length
	@echo make archive .... to make an archive tar file you can transfer or submit

$(PROGRAM):	$(OBJECTS)
	$(CC) $(OPTIONS) -o $(PROGRAM) $(OBJECTS)

%.o:	%.c
	$(CC) -c $<

clean:	
	rm $(OBJECTS) $(PROGRAM) gmeta gmeta.zip *.gif

listing:
	@echo Creating code listing named listing.txt ...
	pr -F --length=58 --page-width=80 -n    \
           p2.c advection.c bc.c ic.c stats.c sfc.c contr.c advect1d.c error.c Makefile > listing.txt
	@echo listing.txt is `cat listing.txt | wc -l` lines long.
	@echo Your listing is ready.

archive:
	@echo Creating code archive tar-file $(ARCHIVE) ...
	tar cf $(ARCHIVE) p2.c advection.c bc.c ic.c stats.c sfc.c contr.c advect1d.c error.c Makefile
	@ls -l $(ARCHIVE)
	@echo Your archive is ready.
