include make.inc
OBJS= MATRIX_SPARSE.o ED_VARS_GLOBAL.o ED_INPUT_VARS.o ED_BATH_DMFT.o ED_BATH_USER.o ED_BATH_FUNCTIONS.o ED_AUX_FUNX.o ED_IO.o ED_SETUP.o ED_EIGENSPACE.o  ED_MATVEC.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o  ED_GLOC.o ED_WEISS.o ED_ENERGY.o ED_CHI2FIT.o ED_DIAG.o ED_MAIN.o DMFT_ED.o

all: serial version compile completion
mpi: mpi version compile completion
debug: debug version compile completion
debug_mpi:debug_mpi version compile completion


serial:FPPFLAG+=-D_

mpi:FPPFLAG+=-D_MPI

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

