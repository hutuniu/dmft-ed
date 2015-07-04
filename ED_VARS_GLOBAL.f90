MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
#ifdef _MPI_INEQ
  USE MPI
#endif
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

  !MPI Parallel environment variables
  !PUBLIC
  !=========================================================
  integer                                     :: ED_MPI_ID=0
  integer                                     :: ED_MPI_SIZE=1
  integer                                     :: ED_MPI_ERR




  !--------------- LATTICE WRAP VARIABLES -----------------!
  !Lattice size:
  !PUBLIC USE: (should be set by routine)
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
  real(8),dimension(:,:),allocatable          :: H0


END MODULE ED_VARS_GLOBAL
