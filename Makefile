#=========================================================================
include sfmake.inc
FC=ifort
#=========================================================================

#--> HUBBARD MODELS:
EXE=ed_hm_bethe
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
#EXE=ed_bhz
#EXE=ed_bhz_afm

DIR =drivers
DIREXE=$(HOME)/.bin

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#COMPILATION:
OBJS= MATRIX_SPARSE.o ED_EIGENSPACE.o ED_BATH_TYPE.o ED_INPUT_VARS.o ED_VARS_GLOBAL.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_BATH.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD) 
all: ARGS=$(SFLIBS)
all:compile 

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT)
opt: ARGS=$(SFLIBS)
opt:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS=$(SFLIBS_DEB)
debug:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 


completion:
	scifor_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc


version:
	@echo $(VER)

