MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
  implicit none

  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath
     real(8),dimension(:,:,:),allocatable        :: e     !local energies [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable        :: d     !SC amplitues   [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable        :: v     !spin-keep hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable        :: u     !spin-flip hyb. [Nspin][Norb][Nbath]
     complex(8),dimension(:),allocatable         :: vr    !diagonal hyb.  [Nbath]
     complex(8),dimension(:,:,:,:,:),allocatable :: h     !Replica hamilt [Nspin][Nspin][Norb][Norb][Nbath]
     logical(8),dimension(:,:,:,:,:),allocatable :: mask  !impHloc mask   [Nspin][Nspin][Norb][Norb][Re,Im]
     complex(8),dimension(:,:,:,:,:),allocatable :: LS    !Replica hamilt [Nspin][Nspin][Norb][Norb][LS,LSrot]
     logical                                     :: status=.false.
  end type effective_bath


  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable :: map
  end type sector_map

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate


  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  abstract interface

     !HAMILTONIAN CONSTUCTORS
     !DBLE
     subroutine d_build_hamiltonian(isector,Hmat)
       integer                         :: isector
       real(8),dimension(:,:),optional :: Hmat
     end subroutine d_build_hamiltonian
     !CMPLX
     subroutine c_build_hamiltonian(isector,Hmat)
       integer                            :: isector
       complex(8),dimension(:,:),optional :: Hmat
     end subroutine c_build_hamiltonian


     !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
     !dbleMat*dbleVec
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
     !dbleMat*cmplxVec
     subroutine dc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine dc_sparse_HxV
     !cmplxMat*cmplxVec
     subroutine cc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine cc_sparse_HxV

  end interface




  !-------------------------- ED  VARIABLES --------------------------!
  !SIZE OF THE PROBLEM
  !Ns       =              Number of levels per spin
  !Nlevels  = 2*Ns       = Total Number  of levels
  !Nsectors =              Number of sectors
  !=========================================================
  integer                                            :: Ns
  integer                                            :: Nlevels
  integer                                            :: Nsectors
  integer                                            :: Nhel

  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable          :: impHloc           !local hamiltonian

  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:,:)                 :: getsector
  integer,allocatable,dimension(:,:)                 :: getCsector
  integer,allocatable,dimension(:,:)                 :: getCDGsector
  integer,allocatable,dimension(:,:)                 :: getBathStride
  integer,allocatable,dimension(:,:)                 :: impIndex
  integer,allocatable,dimension(:)                   :: getDim,getDimUp,getDimDw
  integer,allocatable,dimension(:)                   :: getNup,getNdw
  integer,allocatable,dimension(:)                   :: getSz
  integer,allocatable,dimension(:)                   :: getN
  logical,allocatable,dimension(:)                   :: twin_mask

  !Effective Bath used in the ED code (this is opaque to user)
  !PRIVATE
  !=========================================================
  type(effective_bath)                               :: dmft_bath

  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix)                                :: spH0
  type(sparse_matrix)                                :: spH0up,spH0dw
  procedure(d_build_hamiltonian),pointer             :: ed_buildH_d=>null()
  procedure(c_build_hamiltonian),pointer             :: ed_buildH_c=>null()
  !
  procedure(dd_sparse_HxV),pointer                   :: spHtimesV_dd=>null()
  procedure(dc_sparse_HxV),pointer                   :: spHtimesV_dc=>null()
  procedure(cc_sparse_HxV),pointer                   :: spHtimesV_cc=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                   :: neigen_sector
  !--------------- LATTICE WRAP VARIABLES -----------------!  
  integer,allocatable,dimension(:,:)                 :: neigen_sectorii
  integer,allocatable,dimension(:)                   :: neigen_totalii
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function

  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSmats,impSAmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSreal,impSAreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGmats,impFmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGreal,impFreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0mats,impF0mats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0real,impF0real

  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: SAmatsii,SArealii        ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Fmatsii,Frealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]

  ! !Poles & Weights 
  ! !=========================================================
  ! real(8),allocatable,dimension(:,:,:,:,:,:)         :: GFpoles,GFweights


  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)                 :: spinChi_tau
  complex(8),allocatable,dimension(:,:)              :: spinChi_w
  complex(8),allocatable,dimension(:,:)              :: spinChi_iv


  !Diagonal/Off-diagonal charge-charge Susceptibilities
  !=========================================================  
  real(8),allocatable,dimension(:,:,:)               :: densChi_tau
  complex(8),allocatable,dimension(:,:,:)            :: densChi_w
  complex(8),allocatable,dimension(:,:,:)            :: densChi_iv

  !Mixed inter-orbital charge-charge Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)               :: densChi_mix_tau
  complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_w
  complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_iv

  !Total (orbital-sum) Density-density Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:)                   :: densChi_tot_tau
  complex(8),allocatable,dimension(:)                :: densChi_tot_w
  complex(8),allocatable,dimension(:)                :: densChi_tot_iv

  !Pair-Pair Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)                 :: pairChi_tau
  complex(8),allocatable,dimension(:,:)              :: pairChi_w
  complex(8),allocatable,dimension(:,:)              :: pairChi_iv



  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                   ::  ed_dens
  real(8),dimension(:),allocatable                   ::  ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable                   ::  ed_docc
  real(8),dimension(:),allocatable                   ::  ed_phisc

  !--------------- LATTICE WRAP VARIABLES -----------------!
  real(8),dimension(:,:),allocatable,save            :: nii,dii,mii,pii


  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  real(8)                                            :: ed_Ekin
  real(8)                                            :: ed_Epot
  real(8)                                            :: ed_Eint
  real(8)                                            :: ed_Ehartree
  real(8)                                            :: ed_Eknot
  real(8)                                            :: ed_Dust,ed_Dund,ed_Dse,ed_Dph
  !--------------- LATTICE WRAP VARIABLES -----------------!
  real(8),dimension(:,:),allocatable,save            :: ddii,eii


  !Impurity operators
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:)          :: imp_density_matrix
  complex(8),allocatable,dimension(:,:,:)            :: impStot
  complex(8),allocatable,dimension(:,:,:)            :: impLtot
  complex(8),allocatable,dimension(:)                :: impj_aplha
  complex(8),allocatable,dimension(:)                :: impj_aplha_sq
  complex(8)                                         :: impLdotS



contains

  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer :: N
    allocate(H%map(N))
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector


END MODULE ED_VARS_GLOBAL
