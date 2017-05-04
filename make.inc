##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
##$ use mpif90 for parallel compiler
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu
##$ if using a parallel compiler (mpif90) check
##$ platform first using the cmd: $mpif90 -show
PLAT=gnu


##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
##$ if needed add yours to the list:
##$ --> HUBBARD MODELS:
#EXE=ed_hm_bethe
#EXE=ed_ahm_bethe
EXE=ed_ahm_square
#EXE=ed_hm_2dsquare
#EXE=ed_hm_2b_cubic
#EXE=ed_hm_bethe_afm
#EXE=ed_hm_2bands_hyb_fcc3d
#EXE=ed_hm_mpitest
#EXE=ed_hm_2bands_dos

##$ --> PERIODIC ANDERSON & P-D MODELS
#EXE=ed_pam_1b
#EXE=ed_pam_2b
#EXE=ed_lda1b
#EXE=ed_lda
#EXE=ed_tddpam_lattice
#EXE=ed_tddpam_bethe

##$ --> B-H-Z MODELS
#EXE=ed_bhz_2d
#EXE=ed_bhz_3d
#EXE=ed_bhz_2d_edge
#EXE=ed_bhz_afm_2d
#EXE=ed_bhz_afm_3d

#EXE=ed_bhz_2d_edge_afm
#EXE=ed_bhz_2d_afm4
#EXE=ed_bhz_2d_afm2

#EXE=ed_graphene
#EXE=ed_haldane
#EXE=ed_haldane_sc

##$ --> INHOMOGENEOUS (RDMFT)
#EXE=ed_ahm_disorder
#EXE=ed_hm_slab_hyb
#EXE=ed_ahm_stripe
#EXE=ed_ahm_finite_stripe
#EXE=ed_nano
#EXE=ed_nano_sc
#EXE=hartree_nano
#EXE=ed_cdwhm_bethe

##$ STO (nonSU2)
#EXE=ed_STO
#EXE=ed_SOC

EXE=ed_ahm_square_2bands

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =./drivers

##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin

INCARGS=$(shell pkg-config --cflags --libs-only-L dmft_tools scifor)

ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz $(INCARGS) 
OFLAG=-O3 -ftz -prof-gen -prof-use $(INCARGS) 
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created $(INCARGS) 
FPPFLAG =-fpp
endif

ifeq ($(PLAT),gnu)
# INCARGS=-I$(LIBDIR)/scifor/gnu/include -L$(LIBDIR)/scifor/gnu/lib
# INCARGS+=-I$(LIBDIR)/dmft_tools/gnu/include -L$(LIBDIR)/dmft_tools/gnu/lib
FFLAG = -O2 -ffree-line-length-none $(INCARGS) 
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none  $(INCARGS)
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none $(INCARGS)
FPPFLAG =-cpp
endif


##$ LINKING OPTIONS:
ARGS=  -ldmft_tools -lscifor


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
