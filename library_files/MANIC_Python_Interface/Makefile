## Makefile for constructing the manic python modules.
## Options available are:


# choose compiler to use
#COMPILER = NAGF95
COMPILER = GNUFOR

## Compiler command & options
FC_NAGF95   = nagfor
FOPT_NAGF95 = -O0 -w=all -f90 -gline -C=all
FEXT_NAGF95 = -fpp -DNAGF95

FC_GNUFOR   = gfortran-mp-4.9
FOPT_NAGF95 = -O3
FOPT_DEBUG_NAGF95 = -O0 -w=all -f90 -gline -C=all
FEXT_GNUFOR = 

# select commands & options based on compiler
FC   = $(FC_$(COMPILER))
FOPT = $(FOPT_$(COMPILER))
FEXT = $(FEXT_$(COMPILER))


# define which libraries to use
LIBS =
LIBPATH = ../MANIC_SourceCode/
LOCLIBS = -L$(LIBPATH) -lmapmanic -I$(LIBPATH)modules/

# commands for f2py
FP = f2py-2.7
FPOPT = --f90exec=$(FC) -c -m
FPLIB = $(LIBPATH)libmapmanic.a -I$(LIBPATH)modules/



# aims

## all : compile all python modules
.PHONY : all
all: manic monitor parameters

## debug : compile all python modules with the fortran debugging options
.PHONY : debug
debug: 	FOPT = $(FOPT_DEBUG_$(COMPILER))
debug:  all	


## manic : compile the MANIC python module
.PHONY : manic
manic:  MANIC_Python.f90 MAP_Python_Initialize.o 
	$(FP) $(FPOPT) MANIC_Python MANIC_Python.f90 MAP_Python_Initialize.o $(FPLIB)

## monitor : compile the MAP_Monitor python module
.PHONY : monitor
monitor: $(LIBPATH)MAP_Monitor.f90
	$(FP) $(FPOPT) MAP_Monitor $(LIBPATH)MAP_Monitor.f90

## parameters : compile the MAP_Parameters python module
.PHONY : parameters
parameters: $(LIBPATH)MAP_Parameters.f90
	$(FP) $(FPOPT) MAP_Parameters $(LIBPATH)MAP_Parameters.f90 $(FPLIB)

## clean : remove all compiled python modules
.PHONY : clean
clean:
	rm -f *.o *.mod *.so
	rm -rf *.so.dSYM

## help : this help message
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<



# executables

MAP_Python_Initialize.o: MAP_Python_Initialize.f90
	$(FC) $(FOPT) $(LOCLIBS) -c $<


