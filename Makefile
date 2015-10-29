#include make.inc

#COMPILER (PARALLEL IF FPP_INEQ or FPP_ED ARE SET)
FC=mpif90
PLAT=gnu


#PRECOMPILATION FLAG (leave blank for serial code)
#this is the flag for the inequivalent parallelism [MPI_INEQ]
FPP_INEQ=
#this is the flag for the ed parallelism [MPI]
FPP_ED=MPI


.SUFFIXES: .f90

ifeq ($(PLAT),intel)
FFLAG +=-fpp -D_$(FPP_INEQ) -D_$(FPP_ED) 
endif

ifeq ($(PLAT),gnu)
#questi mi includono le libreirie customizzate (-ldmftt -lscifor), le altre (-lminpack -llapack -lblas  -larpack -lparpack) sono in /usr/lib
INCARGS=-I/home/xps13/Documenti/PhD/programs/libraries/scifor/gnu/include #-I/home/xps13/Documenti/PhD/thesis/adriano/dmft_tools/gnu/include
#FFLAG +=-O0 -p -g -Wall -fbacktrace -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -pedantic  -fbacktrace -fcheck=bounds
FFLAG +=-ffree-line-length-none -cpp -D_$(FPP_INEQ) -D_$(FPP_ED) $(INCARGS)
endif

#CHOOSE LINKING OPTIONS:
#if you intend to use mkl:
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) -larpack -lparpack 
#ELSE:
ARGS=  -ldmftt -lscifor -lminpack -llapack -lblas  -larpack -lparpack -lfftpack


#REVISION SOFTWARE GIT:
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif

#--> HUBBARD MODELS:
#EXE=ed_hm_bethe
#EXE=ed_ahm_bethe
#EXE=ed_ahm_square
#EXE=ed_hm_2dsquare
#EXE=ed_hm_2b_cubic
#EXE=ed_hm_bethe_afm
#EXE=ed_hm_2bands_hyb_fcc3d
#EXE=ed_hm_mpitest

#--> PERIODIC ANDERSON & P-D MODELS
#EXE=ed_pam_1b
#EXE=ed_pam_2b
#EXE=ed_lda1b
#EXE=ed_lda
#EXE=ed_tddpam_lattice
#EXE=ed_tddpam_bethe

#--> B-H-Z MODELS
#EXE=ed_2x2bhz
#EXE=ed_bhz
EXE=ed_STO
#EXE=ed_bhz_afm
#EXE=ed_bhz_edge

#INHOMO
#EXE=ed_ahm_disorder
#EXE=ed_hm_slab_hyb
#EXE=ed_ahm_stripe
#EXE=ed_ahm_finite_stripe
#EXE=ed_nano


#EXE=ed_cdwhm_bethe


DIR   =drivers
DIREXE=/home/xps13/Documenti/PhD/thesis/adriano/dmft_ed/TEST_merge_2


OBJS= MATRIX_SPARSE.o ED_BATH_TYPE.o ED_VARS_GLOBAL.o ED_INPUT_VARS.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_SETUP.o  ED_EIGENSPACE.o ED_BATH_DMFT.o ED_BATH_USER.o ED_BATH_FUNCTIONS.o ED_MATVEC.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o  ED_GLOC.o ED_WEISS.o ED_ENERGY.o ED_CHI2FIT.o ED_DIAG.o ED_MAIN.o  DMFT_ED.o

all: version compile completion
debug: debug version compile completion

debug: FFLAG=$(DFLAG)

compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FFLAG) -c $< 

completion:
	scifor_completion.sh $(DIR)/$(EXE).f90

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

