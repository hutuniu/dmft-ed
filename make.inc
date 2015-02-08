#COMPILER (PARALLEL)
FC=gfortran
PLAT=gnu
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=



.SUFFIXES: .f90

ifeq ($(PLAT),intel)
FFLAG +=-fpp -D_$(FPP)
endif
ifeq ($(PLAT),gnu)
INCARGS=-I/opt/scifor/gnu/include -I/opt/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none -cpp -D_$(FPP) $(INCARGS)
endif

#CHOOSE LINKING OPTIONS:
#if you intend to use mkl:
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) -lminpack -larpack -lparpack 
#ELSE:
ARGS=  -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack -larpack -lparpack
