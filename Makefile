#=========================================================================
include sfmake.inc
#=========================================================================
#EXE=fulled_lda1b
#EXE=lanced_lda1b
#EXE =fulled_pam_2dsquare
EXE=fulled_pam_2imp
#EXE=fulled_hm_bethe
#EXE=lanced_hm_bethe
#EXE=testBATH
DIR =drivers
DIREXE=$(HOME)/.bin

#FC=$(SFMPI)/mpif90

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#COMPILATION:
OBJS= MATRIX_SPARSE.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_EIGENSPACE.o ED_VARS_GLOBAL.o ED_AUX_FUNX.o ED_BATH.o ED_GETH.o ED_GETGF.o ED_GETOBS.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

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

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)
