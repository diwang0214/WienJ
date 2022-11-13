OBJ =  module.o test1.o \
caljqso.o  matrixprod.o  readalm.o    readradwfd.o  ylm.o \
caljq.o       gpoint.o   readklist.o  readstruct.o \
jxjf.o     readalmd.o    readr2v.o    latgen.o

# compiler
F90  = mpiifort -fpp -DMPI

INCLUDE = -I${MKLROOT}/include
WFLAG = -nogen-interface
OFLAG = -O3 -g -traceback -static-intel
FFLAG = $(OFLAG) $(WFLAG)
LFLAG = $(OFLAG)

# blas and lapack libraries
# static linking
# LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
        ${MKLROOT}/lib/intel64/libmkl_sequential.a \
        ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

# dynamic linking
# LIBS = -L/${MKLROOT}/lib/intel64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lpthread
 
main : $(OBJ)
	$(F90) $(LFLAG) $(OBJ) -o wienj # $(LIBS) 
	cp wienj ../BIN   

.SUFFIXES: .o .f90

.f90.o :
	$(F90) $(FFLAG) $(INCLUDE) -c $*.f90

clean :
	rm -f *.o *.mod *~ wienj

