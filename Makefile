include make.inc

#MPI FREE:
#+ED_SPARSE_MATRIX
#+ED_VARS_GLOBAL
#+ED_AUX_FUNX
#+ED_IO
#+ED_EIGENSPACE
#+ED_BATH
#+ED_BATH_FUNCTIONS

#+ED_SETUP: LOCAL MPI in *init_ed_structure*:


#MODULE DEFINED MPI:
#+ED_MATVEC: requires to call *ed_matvec_set/del_MPI* to set MPI_COMM
#+ED_HAMILTONIAN: requires to call *ed_hamiltonian_set/del_MPI* to set local MPI_COMM
## USED IN
#++> ED_DIAG: diagonalize_impurity (top level) has two versions: serial/parallel. It is in charge of setting and deleting MPI in related modules (MATVEC, HAMILTONIAN).
#             the rest of the module uses MPI implicitly as usual. Where requires pre-processing is used to separete serial and parallel channel
#
#++> ED_OBSERVABLES: almost MPI free. observables_impurity (top level) has optional MPI communicator. If MPI is defined communicator is used to set MPI_MASTER.
#
#

OBJS= ED_SPARSE_MATRIX.o ED_INPUT_VARS.o ED_VARS_GLOBAL.o  ED_AUX_FUNX.o ED_IO.o ED_SETUP.o ED_EIGENSPACE.o ED_BATH.o ED_BATH_FUNCTIONS.o ED_MATVEC.o ED_HAMILTONIAN.o ED_DIAG.o ED_OBSERVABLES.o ED_GREENS_FUNCTIONS.o  ED_CHI2FIT.o ED_MPI2B.o ED_GLOC.o ED_WEISS.o ED_ENERGY.o ED_MAIN.o DMFT_ED.o


#TRY TO USE THE DMFT_TOOLS ONES
#

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

