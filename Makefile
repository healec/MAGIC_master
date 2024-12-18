#
#               Makefile for the mpiclaw code:
#                       (with MAGIC)
#
#       To make an MPI executable, type:          make xclawmpi
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawmpihdf
#       (from the application sub-directory)
#
#       To combine ASCII output files (one from
#       each processor at each time level) into
#       fort.qXXXX files for use with MATLAB, type:  make catfiles
#       (in directory with fort.qXXXX.YY files)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
FFLAGS = -O3
F90FLAGS = $(FFLAGS)
LFLAGS = $(FFLAGS)
F77    = ftn -c
F90    = ftn -c
LINK   = ftn
HDFLIBS = 

%.o: %.f
	$(F77) $(FFLAGS) $*.f -o $*.o

%.o: %.f90
	$(F90) $(F90FLAGS) $*.f90 -o $*.o

OBJECTS = \
  b4step3.o \
  bc3_mpi.o \
  qinit.o \
  rpn3euX.o \
  rpt3euX.o \
  rptt3euX.o \
  fcns.o \
  setaux.o \
  setprob.o \
  diffuse.o \
  conduct.o \
  ohsteady.o \
  src3.o

LIBOBJECTS = \
  mpilib/claw3ez_mpi_driver.o \
  mpilib/claw3ez_mpi.o \
  mpilib/claw3_mpi.o \
  lib/chkmth.o \
  lib/step3.o \
  lib/step3ds.o \
  lib/dimsp3.o \
  lib/flux3fw.o \
  lib/copyq3.o

OUTOBJECTS = mpilib/out3_mpi.o \
             mpilib/restart3_mpi.o

HDFOBJECTS = mpilib/out3_mpi_h5sep.o \
             mpilib/restart3_mpi_hdfsep.o

xclawmpi: $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclawmpi 

xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi \
	$(HDFLIBS)

# Type "make catfiles" to combine fort.qXXXX.YY output files from out3_mpi
# into fort.qXXXX files (for use with CLAWPACK's MATLAB graphics routines).
# This step is not necessary with the HDF output routines.
catfiles : fort.q0000

fort.q0000 : fort.q0000.00
	mpilib/catfiles

fort.q0000.00 :
	mpilib/catfiles

### DO NOT remove this line - make depends on it ###
