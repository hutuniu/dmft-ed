MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
  !#ifdef _MPI
  USE MPI
  !#endif
  implicit none

  !SIZE OF THE PROBLEM
  !Ns   = # of levels per spin
  !Ntot = 2*Ns = total #  of levels
  !NN   = 2**Ntot = 2**(2Ns) max size of the Hilbert space
  !Nbo  =# number of bath sites (all sites - impurity sites)
  !Nsect=# of sectors
  !=========================================================
  integer                                     :: Ns,Ntot,NN
  integer                                     :: Nsect
  integer                                     :: Nbo


  !Lattice size:
  !=========================================================
  integer                                     :: Nlat

  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable   :: impHloc           !local hamiltonian


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)          :: getsector
  integer,allocatable,dimension(:,:)          :: getCsector
  integer,allocatable,dimension(:,:)          :: getCDGsector
  integer,allocatable,dimension(:,:)          :: getBathStride
  integer,allocatable,dimension(:,:)          :: impIndex
  integer,allocatable,dimension(:)            :: getdim,getnup,getndw,getsz
  logical,allocatable,dimension(:)            :: twin_mask


  !Effective Bath used in the ED code (this is opaque to user)
  !=========================================================
  type(effective_bath)                        :: dmft_bath

  !Variables for DIAGONALIZATION
  !=========================================================  
  type(sparse_matrix)                         :: spH0
  integer,allocatable,dimension(:)            :: neigen_sector
  logical                                     :: trim_state_list=.false.

  !Partition function
  !=========================================================
  real(8)                                     :: zeta_function


  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSmats,impSAmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSreal,impSAreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGmats,impFmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGreal,impFreal

  !Density and double occupancy
  !=========================================================
  real(8),dimension(:),allocatable            ::  ed_dens
  real(8),dimension(:),allocatable            ::  ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable            ::  ed_docc
  real(8),dimension(:),allocatable            ::  ed_phisc


  !Local energies and generalized double occupancies
  !=========================================================
  real(8)                                     :: ed_Ekin
  real(8)                                     :: ed_Epot
  real(8)                                     :: ed_Eint
  real(8)                                     :: ed_Ehartree
  real(8)                                     :: ed_Eknot
  real(8)                                     :: ed_Dust,ed_Dund,ed_Dse,ed_Dph


  !#ifdef _MPI
  !MPI Parallel environment variables
  !=========================================================
  integer                                     :: ED_MPI_ID=0
  integer                                     :: ED_MPI_SIZE=1
  integer                                     :: ED_MPI_ERR
  !#endif



  !Symmetry operations
  !=========================================================
  integer                                 :: Nsymm,Nindep
  integer,allocatable,dimension(:)        :: indep_list
  integer,dimension(:),allocatable        :: map_lat2ind
  integer,dimension(:,:),allocatable      :: map_ind2lat


  !OBSOLETE (to be removed associated to build_tight_binding_2dsquare)
  !Large matrices for Lattice Hamiltonian/GF
  !=========================================================
  integer,dimension(:),allocatable        :: icol,irow
  integer,dimension(:,:),allocatable      :: ij2site
  real(8),dimension(:,:),allocatable      :: H0
  complex(8),dimension(:,:,:),allocatable :: Hk,Hk_


END MODULE ED_VARS_GLOBAL
