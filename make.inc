##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
##$ use mpif90 for parallel compiler
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu
##$ if using a parallel compiler (mpif90) check
##$ platform first using the cmd: $mpif90 -show
PLAT=gnu


##$ PRECOMPILATION FLAGS
##$ Leave blank for SERIAL CODE
##$
##$ parallelism on the inequivalent parallelism flag: MPI_INEQ
FPP_INEQ=MPI_INEQ

##$ parallelism on the ed code flag: MPI
FPP_ED=

##$ LOCATION OF THE scifor and dmft_tools DIRECTORIES
##$ is placed in different locations edit INCARGS here below
LIBDIR=$(HOME)/opt



##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
##$ if needed add yours to the list:
##$ --> HUBBARD MODELS:
#EXE=ed_hm_bethe
#EXE=ed_ahm_bethe
#EXE=ed_ahm_square
#EXE=ed_hm_2dsquare
#EXE=ed_hm_2b_cubic
#EXE=ed_hm_bethe_afm
#EXE=ed_hm_2bands_hyb_fcc3d
#EXE=ed_hm_mpitest

##$ --> PERIODIC ANDERSON & P-D MODELS
#EXE=ed_pam_1b
#EXE=ed_pam_2b
#EXE=ed_lda1b
#EXE=ed_lda
#EXE=ed_tddpam_lattice
#EXE=ed_tddpam_bethe

##$ --> B-H-Z MODELS
#EXE=ed_2x2bhz
#EXE=ed_bhz
#EXE=ed_bhz_afm
#EXE=ed_bhz_edge

##$ --> INHOMOGENEOUS (RDMFT)
#EXE=ed_ahm_disorder
#EXE=ed_hm_slab_hyb
#EXE=ed_ahm_stripe
#EXE=ed_ahm_finite_stripe
EXE=ed_nano
#EXE=ed_cdwhm_bethe


##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =drivers

##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin


ifeq ($(PLAT),intel)
FFLAG +=-fpp -D_$(FPP_INEQ) -D_$(FPP_ED) 
endif

ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBDIR)/scifor/gnu/include -I$(LIBDIR)/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none -cpp -D_$(FPP_INEQ) -D_$(FPP_ED) $(INCARGS)
endif




##$ CHOOSE LINKING OPTIONS:
##$ 
##$ If you intend to use mkl:
##$ 
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) -larpack -lparpack 
##$ 
##$ ELSE:
##$ 
ARGS=  -ldmftt -lscifor -lfftpack -lminpack -llapack -lblas  -larpack -lparpack




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
