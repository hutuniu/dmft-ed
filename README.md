####TODOs LIST:
1. MPI: the whole code needs to be moved to a full and single MPI support. A set of global MPI variables are included in the `ED_VARS_GLOBAL.f90` module. The user pass the communicator to the few code interfaces (e.g. `ed_solve`) which will set the internal communicator, ids and size. The inequivalent sites case splits the user communicator into a set of new communicator, using the code input variable `mpi_colors`. Each color/group perform an ED calculation with of `MPI_SIZE/mpi_colors` size.

2. Nup/Ndw basis: reduce the calculation size moving the core part of the code to the Nup/Ndw basis.

3. direct HxV: restore the direct Matrix-Vector calculation to solve larger systems. Study strategies to speed-up the code such as:
	*  vector fragments compatible with non-zero columns in the Hamiltonian
	*  on-the-fly lookup to determine the out state after application of an operator (we use bisection now on the direct map)
 
 
##### Actual situation:
I removed `ARPACK_LANCZOS` and `PLAIN_LANCZOS`, now in `SF_SP_LINALG` (SciFor).

I changed all the interfaces to arpack and lanczos tridiagonalization. 
The code compiles after having changed the routines in `ED_MATVEC` to fulfill the new routines in SciFor.

todo: finish the MPI upgrade by hiding the internal MPI variables and force user to pass the MPI Communicator as input.