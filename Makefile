#=========================================================================
include sfmake.inc
#=========================================================================
#COMPILER (PARALLEL)
FC=$(MPIDIR)/bin/mpif90
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=MPI

#--> HUBBARD MODELS:
#EXE=ed_hm_bethe
#EXE=ed_ahm_bethe
#EXE=ed_hm_2dsquare
#EXE=ed_hm_bethe_afm

#--> PERIODIC ANDERSON & P-D MODELS
#EXE=ed_pam_1b
#EXE=ed_pam_2b
#EXE=ed_lda1b
#EXE=ed_lda
#EXE=ed_tddpam_lattice
#EXE=ed_tddpam_bethe

#--> B-H-Z MODELS
#EXE=ed_2x2bhz
EXE=ed_bhz
#EXE=ed_bhz_afm

DIR =drivers
DIREXE=$(HOME)/.bin

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

OBJS= MATRIX_SPARSE.o ED_BATH_TYPE.o ED_VARS_GLOBAL.o ED_INPUT_VARS.o  ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_EIGENSPACE.o ED_BATH.o ED_MATVEC.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD) -fpp -D_$(FPP)
all: INCLUDE=$(SFINCLUDE)
all: ARGS=$(SFARGS) $(MKLARGS) -larpack -lparpack
all:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB) -fpp -D_$(FPP)
debug: INCLUDE=$(SFINCLUDE_DEB)
debug: ARGS=$(SFARGS_DEB) $(MKLARGS) -larpack_debug -lparpack_debug
debug:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(INCLUDE) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

.f90.o:	
	$(FC) $(FLAG) -c $< $(INCLUDE) 


completion:
	$(SFCOMPLETION) $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

