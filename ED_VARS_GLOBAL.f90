MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
#ifdef _MPI
  USE MPI
#endif
  implicit none


  !-------------------- ED  VARIABLES ----------------------!

  !SIZE OF THE PROBLEM
  !Ns       =              Number of levels per spin
  !Nlevels  = 2*Ns       = Total Number  of levels
  !Nhilbert = 2**Nlevels = Max size of the Hilbert space
  !Nsectors =              Number of sectors
  !=========================================================
  integer                                     :: Ns
  integer                                     :: Nlevels
  integer                                     :: Nhilbert
  integer                                     :: Nsectors
  integer                                     :: Nhel

  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable   :: impHloc           !local hamiltonian

  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:,:)          :: getsector
  integer,allocatable,dimension(:,:)          :: getCsector
  integer,allocatable,dimension(:,:)          :: getCDGsector
  integer,allocatable,dimension(:,:)          :: getBathStride
  integer,allocatable,dimension(:,:)          :: impIndex
  integer,allocatable,dimension(:)            :: getdim
  integer,allocatable,dimension(:)            :: getNup,getNdw
  integer,allocatable,dimension(:)            :: getSz
  integer,allocatable,dimension(:)            :: getN
  logical,allocatable,dimension(:)            :: twin_mask

  !Effective Bath used in the ED code (this is opaque to user)
  !PRIVATE
  !=========================================================
  type(effective_bath)                        :: dmft_bath

  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix)                         :: spH0
  integer,allocatable,dimension(:)            :: neigen_sector
  logical                                     :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                     :: zeta_function

  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSmats,impSAmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSreal,impSAreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGmats,impFmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGreal,impFreal

  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable            ::  ed_dens
  real(8),dimension(:),allocatable            ::  ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable            ::  ed_docc
  real(8),dimension(:),allocatable            ::  ed_phisc

  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  real(8)                                     :: ed_Ekin
  real(8)                                     :: ed_Epot
  real(8)                                     :: ed_Eint
  real(8)                                     :: ed_Ehartree
  real(8)                                     :: ed_Eknot
  real(8)                                     :: ed_Dust,ed_Dund,ed_Dse,ed_Dph

  !Impurity dennsity matrix
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:)  :: imp_density_matrix

  !MPI Parallel environment variables
  !PRIVATE
  !=========================================================
  integer                                     :: ED_GLOBAL_COMM



  !--------------- LATTICE WRAP VARIABLES -----------------!
  !Lattice size:
  !PUBLIC USE: (should be set by routine)
  !>DEBUG
  !THIS VARIABLE SHOULD BE REMOVED: THE CODE  WORKS WITH THE
  !# OF INEQ SITES SPECIFIED BY THE LENGHT OF THE BATH
  !<DEBUG
  !=========================================================
  integer                                     :: Nlat

  !Symmetry operations
  !=========================================================
  integer,allocatable,dimension(:)            :: indep_list
  integer,dimension(:),allocatable   :: map_lat2ind
  integer,dimension(:,:),allocatable :: map_ind2lat



  !OBSOLETE (to be removed associated to build_tight_binding_2dsquare)
  !Large matrices for Lattice Hamiltonian/GF
  !PUBLIC
  !=========================================================
  integer,dimension(:),allocatable            :: icol,irow
  integer,dimension(:,:),allocatable          :: ij2site

#ifdef _MPI

contains

  subroutine ED_MPI_Init_Print(comm)
    integer :: comm
    integer :: mpi_rank
    integer :: mpi_size
    integer :: mpi_ierr
    call MPI_Comm_rank(comm,mpi_rank,mpi_ierr)
    call MPI_Comm_size(comm,mpi_size,mpi_ierr)
    print*,"I am rank",mpi_rank,"of ",mpi_size
    call MPI_Barrier(comm,mpi_ierr)
    if(mpi_rank==0)print*,"----------------------------------"
    if(mpi_rank==0)print*,""
  end subroutine ED_MPI_Init_Print
  !
  function ED_MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function ED_MPI_Get_size
  !
  function ED_MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function ED_MPI_Get_rank
  !
  function ED_MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function ED_MPI_Get_master
#endif

END MODULE ED_VARS_GLOBAL
