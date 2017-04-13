!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: free_unit,reg,free_units,txtfy,splot
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,inv_sym,inv_her
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO
  USE ED_EIGENSPACE
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_MATVEC
  USE ED_AUX_FUNX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  !
  implicit none
  private 



  public :: buildGf_impurity
  
  public :: buildChi_impurity

  public :: ed_greens_functions_set_MPI

  public :: ed_greens_functions_del_MPI



  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm

  !Auxiliary functions GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impDeltamats,impDeltareal
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpG0mats,invimpG0real
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpGmats,invimpGreal

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:)       :: auxGmats,auxGreal
  complex(8),allocatable,dimension(:,:,:)     :: auxGpoles,auxGweights


#ifdef _MPI
  integer                                     :: MpiComm=MPI_UNDEFINED
#endif
  logical                                     :: MpiStatus=.false.  
  integer                                     :: MPI_RANK=0
  integer                                     :: MPI_SIZE=1
  logical                                     :: MPI_MASTER=.true.
  integer                                     :: mpi_ierr



contains


  subroutine ed_greens_functions_set_MPI(comm)
#ifdef _MPI
    integer :: comm
    MpiComm  = comm
    MpiStatus = .true.
    MPI_RANK = get_Rank_MPI(MpiComm)
    MPI_SIZE = get_Size_MPI(MpiComm)
    MPI_MASTER=get_Master_MPI(MpiComm)
#else
    integer,optional :: comm
#endif
  end subroutine ed_greens_functions_set_MPI


  subroutine ed_greens_functions_del_MPI()
#ifdef _MPI
    MpiComm  = MPI_UNDEFINED
    MpiStatus = .false.
#endif
  end subroutine ed_greens_functions_del_MPI






  subroutine buildgf_impurity()
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    if(.not.allocated(impGmats))stop "buildgf_impurity: impGmats not allocated"
    if(.not.allocated(impGreal))stop "buildgf_impurity: impGreal not allocated"
    if(.not.allocated(impFmats))stop "buildgf_impurity: impFmats not allocated"
    if(.not.allocated(impFreal))stop "buildgf_impurity: impFreal not allocated"
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    if(.not.allocated(impSmats)) stop "buildgf_impurity: impSmats not allocated"
    if(.not.allocated(impSreal)) stop "buildgf_impurity: impSreal not allocated"
    if(.not.allocated(impSAmats))stop "buildgf_impurity: impSAmats not allocated"
    if(.not.allocated(impSAreal))stop "buildgf_impurity: impSAreal not allocated"    
    impSmats = zero
    impSreal = zero
    impSAmats = zero
    impSAreal = zero
    !
    if(.not.allocated(impG0mats)) stop "buildgf_impurity: impG0mats not allocated"
    if(.not.allocated(impG0real)) stop "buildgf_impurity: impG0real not allocated"
    if(.not.allocated(impF0mats)) stop "buildgf_impurity: impF0mats not allocated"
    if(.not.allocated(impF0real)) stop "buildgf_impurity: impF0real not allocated"
    impG0mats=zero
    impG0real=zero
    impF0mats=zero
    impF0real=zero
    !
    if(.not.allocated(GFpoles))   stop "buildgf_impurity: GFpoles not allocated"
    if(.not.allocated(GFweights)) stop "buildgf_impurity: GFweights not allocated"
    GFpoles=zero
    GFweights=zero
    !
    if(.not.allocated(impDeltamats)) allocate(impDeltamats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpG0mats)) allocate(invimpG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpGmats))  allocate( invimpGmats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impDeltareal)) allocate(impDeltareal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpG0real)) allocate(invimpG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpGreal))  allocate( invimpGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impDeltamats=zero
    invimpGmats=zero
    invimpG0mats=zero
    impDeltareal=zero
    invimpGreal=zero
    invimpG0real=zero
    !
    if(MPI_MASTER) write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default
       call build_gf_normal()
       call get_sigma_normal()
    case ("superc")
       call build_gf_superc()
       call get_sigma_superc()
    case ("nonsu2")
       call build_gf_nonsu2()
       call get_sigma_nonsu2()
    end select
    !
    if(mpi_master)then
       select case(ed_mode)
       case default
          call print_poles_weights_normal()
          call print_impSigma_normal()
          call print_impG_normal()
          call print_impG0_normal()
       case ("superc")
          call print_poles_weights_superc()
          call print_impSigma_superc()
          call print_impG_superc()
          call print_impG0_superc()
       case ("nonsu2")
          call print_poles_weights_nonsu2()
          call print_impSigma_nonsu2()
          call print_impG_nonsu2()
          call print_impG0_nonsu2()
       end select
    endif
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(invimpG0mats))deallocate(invimpG0mats)
    if(allocated(invimpGmats))deallocate(invimpGmats)
    if(allocated(impDeltamats))deallocate(impDeltamats)
    if(allocated(invimpG0real))deallocate(invimpG0real)
    if(allocated(invimpGreal))deallocate(invimpGreal)
    if(allocated(impDeltareal))deallocate(impDeltareal)
  end subroutine buildgf_impurity

  !+------------------------------------------------------------------+
  !                    GREEN'S FUNCTIONS 
  !+------------------------------------------------------------------+
  include 'ED_GREENS_FUNCTIONS/build_gf_normal.f90'
  include 'ED_GREENS_FUNCTIONS/build_gf_superc.f90'
  include 'ED_GREENS_FUNCTIONS/build_gf_nonsu2.f90'





  !+------------------------------------------------------------------+
  !                    SELF-ENERGY FUNCTIONS 
  !+------------------------------------------------------------------+
  include "ED_GREENS_FUNCTIONS/get_sigma_normal.f90"
  include "ED_GREENS_FUNCTIONS/get_sigma_superc.f90"
  include "ED_GREENS_FUNCTIONS/get_sigma_nonsu2.f90"


























  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    integer :: i
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    if(.not.allocated(spinChi_tau)) stop "buildChi_impurity: spinChi_tau not allocated"
    if(.not.allocated(spinChi_w))  stop "buildChi_impurity: spinChi_w not allocated"
    if(.not.allocated(spinChi_iv)) stop "buildChi_impurity: spinChi_iv not allocated"
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    call build_chi_spin()
    if(MPI_MASTER)call print_chi_spin()
    !
    !BUILD CHARGE SUSCEPTIBILITY
    if(.not.allocated(densChi_tau)) stop "buildChi_impurity: densChi_tau not allocated"
    if(.not.allocated(densChi_w))  stop "buildChi_impurity: densChi_w not allocated"
    if(.not.allocated(densChi_iv)) stop "buildChi_impurity: densChi_iv not allocated"
    if(.not.allocated(densChi_mix_tau))stop "buildChi_impurity: densChi_mix_tau not allocated"
    if(.not.allocated(densChi_mix_w))  stop "buildChi_impurity: densChi_mix_w not allocated"
    if(.not.allocated(densChi_mix_iv)) stop "buildChi_impurity: densChi_mix_iv not allocated"
    if(.not.allocated(densChi_tot_tau))stop "buildChi_impurity: densChi_tot_tau not allocated"
    if(.not.allocated(densChi_tot_w))  stop "buildChi_impurity: densChi_tot_w not allocated"
    if(.not.allocated(densChi_tot_iv)) stop "buildChi_impurity: densChi_tot_iv not allocated"
    densChi_tau=zero
    densChi_w=zero
    densChi_iv=zero
    densChi_mix_tau=zero
    densChi_mix_w=zero
    densChi_mix_iv=zero
    densChi_tot_tau=zero
    densChi_tot_w=zero
    densChi_tot_iv=zero
    call build_chi_dens()
    if(MPI_MASTER)call print_chi_dens()
    if(MPI_MASTER)call print_chi_dens_mix()
    if(MPI_MASTER)call print_chi_dens_tot()
    !
    !BUILD PAIR SUSCEPTIBILITY
    if(.not.allocated(pairChi_tau))stop "buildChi_impurity: pairChi_tau not allocated"
    if(.not.allocated(pairChi_w))  stop "buildChi_impurity: pairChi_w not allocated"
    if(.not.allocated(pairChi_iv)) stop "buildChi_impurity: pairChi_iv not allocated"
    pairChi_tau=zero
    pairChi_w=zero
    pairChi_iv=zero
    call build_chi_pair()
    if(MPI_MASTER)call print_chi_pair()
    !
    !
    call deallocate_grids
  end subroutine buildChi_impurity

  !+------------------------------------------------------------------+
  !                    SUSCEPTIBILITIES (SPIN, CHARGE, PAIR)
  !+------------------------------------------------------------------+
  include 'ED_GREENS_FUNCTIONS/build_chi_spin.f90'
  include 'ED_GREENS_FUNCTIONS/build_chi_dens.f90'
  include 'ED_GREENS_FUNCTIONS/build_chi_pair.f90'







  






  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids




end MODULE ED_GREENS_FUNCTIONS
