include make.inc

OBJS= ED_SPARSE_MATRIX.o ED_INPUT_VARS.o ED_VARS_GLOBAL.o  ED_AUX_FUNX.o ED_IO.o ED_SETUP.o ED_EIGENSPACE.o ED_BATH.o ED_BATH_FUNCTIONS.o ED_MATVEC.o ED_HAMILTONIAN.o ED_DIAG.o ED_OBSERVABLES.o ED_GREENS_FUNCTIONS.o  ED_CHI2FIT.o  ED_MAIN.o DMFT_ED.o


#TRY TO USE THE DMFT_TOOLS ONES
#
#ED_ENERGY.o ==> DMFT_EKIN
#ED_GLOC.o ==> DMFT_GLOC
#ED_WEISS.o ==> DMFT_WEISS_FIELD

#ED_MPI2B.o

all: all version compile completion
mpi: mpi version compile completion
debug: debug version compile completion
debug_mpi:debug_mpi version compile completion

all: FPPFLAG+=-D_

mpi: FPPFLAG+=-D_MPI

debug: FFLAG=$(DFLAG)
debug: FPPFLAG+=-D_

debug_mpi: FFLAG=${DFLAG} 
debug_mpi: FPPFLAG+=-D_MPI


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FPPFLAG) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FPPFLAG) $(FFLAG) -c $< 

completion:
	scifor_completion.sh $(DIR)/$(EXE).f90

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

