#COMPILER (PARALLEL IF FPP_INEQ or FPP_ED ARE SET)
FC=mpif90
PLAT=intel
#PRECOMPILATION FLAG (leave blank for serial code)
#this is the flag for the inequivalent parallelism [MPI_INEQ]
FPP_INEQ=
#this is the flag for the ed parallelism [MPI]
FPP_ED=

.SUFFIXES: .f90

ifeq ($(PLAT),intel)
FFLAG +=-fpp -D_$(FPP_INEQ) -D_$(FPP_ED)
endif
ifeq ($(PLAT),gnu)
INCARGS=-I/opt/scifor/gnu/include -I/opt/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none -cpp -D_$(FPP_INEQ) -D_$(FPP_ED) $(INCARGS)
endif

#CHOOSE LINKING OPTIONS:
#if you intend to use mkl:
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) -lminpack -larpack -lparpack 
#ELSE:
ARGS=  -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack -larpack -lparpack

#REVISION SOFTWARE GIT:
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif

BRANCH=