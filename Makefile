include make.inc

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
#EXE=ed_bhz_afm
#EXE=ed_bhz_edge

#INHOMO
EXE=ed_ahm_disorder
#EXE=ed_hm_slab_hyb
#EXE=ed_nano

DIR =drivers
DIREXE=$(HOME)/.bin


OBJS= MATRIX_SPARSE.o ED_BATH_TYPE.o ED_VARS_GLOBAL.o ED_INPUT_VARS.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_EIGENSPACE.o ED_BATH.o ED_MATVEC.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o  ED_GLOC.o ED_WEISS.o ED_ENERGY.o ED_CHI2FIT.o ED_DIAG.o ED_MAIN.o ED_WRAP_AUX_FUNX.o ED_WRAP_MAIN.o ED_WRAP_GLOC.o ED_WRAP_WEISS.o ED_WRAP_ENERGY.o ED_WRAP_CHI2FIT.o DMFT_ED.o


all: version compile completion

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

