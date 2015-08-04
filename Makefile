include make.inc

OBJS= MATRIX_SPARSE.o ED_BATH_TYPE.o ED_VARS_GLOBAL.o ED_INPUT_VARS.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_EIGENSPACE.o ED_BATH.o ED_MATVEC.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o  ED_GLOC.o ED_WEISS.o ED_ENERGY.o ED_CHI2FIT.o ED_DIAG.o ED_MAIN.o  DMFT_ED.o

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

