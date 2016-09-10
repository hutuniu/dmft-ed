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
  USE PLAIN_LANCZOS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_MATVEC
  USE ED_AUX_FUNX
  !
  implicit none
  private 

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer               :: state_vec
  complex(8),dimension(:),pointer            :: state_cvec
  real(8)                                    :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable           :: wm,tau,wr,vm

  !Non-interacting GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:):: impG0mats,impG0real
  complex(8),allocatable,dimension(:,:,:,:,:):: impF0mats,impF0real

  !Auxiliary functions GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:):: impDeltamats,impDeltareal
  complex(8),allocatable,dimension(:,:,:,:,:):: invimpG0mats,invimpG0real
  complex(8),allocatable,dimension(:,:,:,:,:):: invimpGmats,invimpGreal

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:)      :: auxGmats,auxGreal
  complex(8),allocatable,dimension(:,:,:)    :: auxGpoles,auxGweights

  !Poles & Weights 
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:,:,:) :: GFpoles,GFweights


  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)      :: spinChi_tau
  complex(8),allocatable,dimension(:,:)   :: spinChi_w
  complex(8),allocatable,dimension(:,:)   :: spinChi_iv


  !Diagonal/Off-diagonal charge-charge Susceptibilities
  !=========================================================  
  real(8),allocatable,dimension(:,:,:)    :: densChi_tau
  complex(8),allocatable,dimension(:,:,:) :: densChi_w
  complex(8),allocatable,dimension(:,:,:) :: densChi_iv

  !Mixed inter-orbital charge-charge Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)    :: densChi_mix_tau
  complex(8),allocatable,dimension(:,:,:) :: densChi_mix_w
  complex(8),allocatable,dimension(:,:,:) :: densChi_mix_iv

  !Total (orbital-sum) Density-density Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:)        :: densChi_tot_tau
  complex(8),allocatable,dimension(:)     :: densChi_tot_w
  complex(8),allocatable,dimension(:)     :: densChi_tot_iv

  !Pair-Pair Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)      :: pairChi_tau
  complex(8),allocatable,dimension(:,:)   :: pairChi_w
  complex(8),allocatable,dimension(:,:)   :: pairChi_iv


  ! <<<<<<< HEAD
  public :: buildgf_impurity
  public :: rebuildgf_impurity
  public :: buildchi_impurity

  public :: print_poles_weights_normal
  public :: print_sigma_normal
  public :: print_impG_normal
  public :: print_impG0_normal
  public :: print_poles_weights_superc
  public :: print_sigma_superc
  public :: print_impG_superc
  public :: print_impG0_superc
  public :: print_poles_weights_nonsu2
  public :: print_sigma_nonsu2
  public :: print_impG_nonsu2
  public :: print_impG0_nonsu2



  
contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Green's function calculation
  !+------------------------------------------------------------------+
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
    if(.not.allocated(impG0mats))    allocate(   impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impDeltamats)) allocate(impDeltamats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpG0mats)) allocate(invimpG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpGmats))  allocate( invimpGmats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impF0mats))    allocate(   impF0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impG0real))    allocate(   impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(impDeltareal)) allocate(impDeltareal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpG0real)) allocate(invimpG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpGreal))  allocate( invimpGreal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(impF0real))    allocate(   impF0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    !
    impDeltamats=zero
    invimpGmats=zero
    invimpG0mats=zero
    !
    impF0mats=zero
    impG0real=zero
    !
    impDeltareal=zero
    invimpGreal=zero
    invimpG0real=zero
    !
    impF0real=zero
    !
    if(.not.allocated(GFpoles))   allocate(GFpoles(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    if(.not.allocated(GFweights)) allocate(GFweights(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    GFpoles=zero
    GFweights=zero
    !
    if(ED_MPI_ID==0) write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default
       call build_gf_normal()
       call get_sigma_normal()
       call print_gf_normal()
    case ("superc")
       call build_gf_superc()
       call get_sigma_superc()
       call print_gf_superc()
    case ("nonsu2")
       call build_gf_nonsu2()
       call get_sigma_nonsu2()
       call print_gf_nonsu2()
    end select
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(GFpoles))deallocate(GFpoles)
    if(allocated(GFweights))deallocate(GFweights)
    if(allocated(impG0mats))deallocate(impG0mats)
    if(allocated(invimpG0mats))deallocate(invimpG0mats)
    if(allocated(invimpGmats))deallocate(invimpGmats)
    if(allocated(impDeltamats))deallocate(impDeltamats)
    if(allocated(impF0mats))deallocate(impF0mats)
    if(allocated(impG0real))deallocate(impG0real)
    if(allocated(invimpG0real))deallocate(invimpG0real)
    if(allocated(invimpGreal))deallocate(invimpGreal)
    if(allocated(impDeltareal))deallocate(impDeltareal)
    if(allocated(impF0real))deallocate(impF0real)
  end subroutine buildgf_impurity


  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Green's function calculation
  !+------------------------------------------------------------------+
  subroutine rebuildgf_impurity()
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
    if(.not.allocated(impG0mats)) allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impF0mats)) allocate(impF0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impG0real)) allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(impF0real)) allocate(impF0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impF0mats=zero
    impG0real=zero
    impF0real=zero
    !
    if(.not.allocated(GFpoles))   allocate(GFpoles(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    if(.not.allocated(GFweights)) allocate(GFweights(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    GFpoles=zero
    GFweights=zero
    !
    if(ED_MPI_ID==0) write(LOGfile,"(A)")"Rebuild impurity Greens functions:"
    select case(ed_mode)
    case default
       call rebuild_gf_normal()
       call get_sigma_normal()
    case ("superc")
       stop "ERROR: rebuild works for ed_mode=normal only"
    case ("nonsu2")
       stop "ERROR: rebuild works for ed_mode=normal only"
    end select
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(GFpoles))deallocate(GFpoles)
    if(allocated(GFweights))deallocate(GFweights)
    if(allocated(impG0mats))deallocate(impG0mats)
    if(allocated(impF0mats))deallocate(impF0mats)
    if(allocated(impG0real))deallocate(impG0real)
    if(allocated(impF0real))deallocate(impF0real)
  end subroutine rebuildgf_impurity





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
  !PURPOSE  : Print impurity Functions case:
  ! - Poles&Weights
  ! - impSigma
  ! - impG
  ! - impG0
  ! NORMAL case
  !+------------------------------------------------------------------+
  subroutine print_gf_normal
    call print_poles_weights_normal
    call print_sigma_normal
    call print_impG_normal
    call print_impG0_normal
  end subroutine print_gf_normal
  include "ED_GREENS_FUNCTIONS/print_gf_normal.f90"

  !SUPERConducting case
  !+------------------------------------------------------------------+
  subroutine print_gf_superc
    call print_poles_weights_superc
    call print_sigma_superc
    call print_impG_superc
    call print_impG0_superc
  end subroutine print_gf_superc
  include "ED_GREENS_FUNCTIONS/print_gf_superc.f90"

  !nonSU2 case
  !+------------------------------------------------------------------+
  subroutine print_gf_nonsu2
    call print_poles_weights_nonsu2
    call print_sigma_nonsu2
    call print_impG_nonsu2
    call print_impG0_nonsu2
  end subroutine print_gf_nonsu2
  include "ED_GREENS_FUNCTIONS/print_gf_nonsu2.f90"



  !+------------------------------------------------------------------+
  !PURPOSE  : Rebuild the impurity Green's functions
  ! - normal
  ! - superc: TODO
  ! - nonsu2: TODO
  !+------------------------------------------------------------------+
  subroutine rebuild_gf_normal
    integer                          :: i,ispin,isign,unit(1),iorb,jorb
    character(len=20)                :: suffix
    integer,dimension(:),allocatable :: getIorb,getJorb
    integer                          :: totNorb,l,j
    real(8)                          :: de,peso
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ("hybrid")             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select

    !Read the Poles&Weights => then it reconstructs the Gimp
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
       call open_units(reg(suffix))
       do isign=1,2
          do i=1,lanc_nGFiter
             read(unit(1),*)(GFpoles(ispin,ispin,iorb,jorb,isign,i),GFweights(ispin,ispin,iorb,jorb,isign,i),ispin=1,Nspin)
          enddo
       enddo
       call close_units
    enddo
    !
    impGmats=zero
    impGreal=zero
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          do isign=1,2
             do j=1,lanc_nGFiter
                de    = GFpoles(ispin,ispin,iorb,jorb,isign,j)
                peso  = GFweights(ispin,ispin,iorb,jorb,isign,j)
                do i=1,Lmats
                   impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(xi*wm(i)-de)
                enddo
                do i=1,Lreal
                   impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(dcmplx(wr(i),eps)-de)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  contains
    subroutine open_units(string)
      character(len=*) :: string
      unit=free_units(1)
      open(unit(1),file="Gpoles_weights"//string//reg(ed_file_suffix)//".ed")
    end subroutine open_units
    subroutine close_units()
      close(unit(1))
    end subroutine close_units
  end subroutine rebuild_gf_normal

















  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Susceptibility calculation
  !+------------------------------------------------------------------+
  subroutine buildchi_impurity()
    integer :: i
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    if(.not.allocated(spinChi_tau))allocate(spinChi_tau(Norb+1,0:Ltau))
    if(.not.allocated(spinChi_w))  allocate(spinChi_w(Norb+1,Lreal))
    if(.not.allocated(spinChi_iv)) allocate(spinChi_iv(Norb+1,0:Lmats))
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    call build_chi_spin()
    call print_chi_spin()
    if(allocated(spinChi_tau))deallocate(spinChi_tau)
    if(allocated(spinChi_w))deallocate(spinChi_w)
    if(allocated(spinChi_iv))deallocate(spinChi_iv)
    !
    !BUILD CHARGE SUSCEPTIBILITY
    if(.not.allocated(densChi_tau))allocate(densChi_tau(Norb,Norb,0:Ltau))
    if(.not.allocated(densChi_w))  allocate(densChi_w(Norb,Norb,Lreal))
    if(.not.allocated(densChi_iv)) allocate(densChi_iv(Norb,Norb,0:Lmats))
    if(.not.allocated(densChi_mix_tau))allocate(densChi_mix_tau(Norb,Norb,0:Ltau))
    if(.not.allocated(densChi_mix_w))  allocate(densChi_mix_w(Norb,Norb,Lreal))
    if(.not.allocated(densChi_mix_iv)) allocate(densChi_mix_iv(Norb,Norb,0:Lmats))
    if(.not.allocated(densChi_tot_tau))allocate(densChi_tot_tau(0:Ltau))
    if(.not.allocated(densChi_tot_w))  allocate(densChi_tot_w(Lreal))
    if(.not.allocated(densChi_tot_iv)) allocate(densChi_tot_iv(0:Lmats))
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
    call print_chi_dens()
    call print_chi_dens_mix()
    call print_chi_dens_tot()
    if(allocated(densChi_tau))deallocate(densChi_tau)
    if(allocated(densChi_w))deallocate(densChi_w)
    if(allocated(densChi_iv))deallocate(densChi_iv)
    if(allocated(densChi_mix_tau))deallocate(densChi_mix_tau)
    if(allocated(densChi_mix_w))deallocate(densChi_mix_w)
    if(allocated(densChi_mix_iv))deallocate(densChi_mix_iv)
    if(allocated(densChi_tot_tau))deallocate(densChi_tot_tau)
    if(allocated(densChi_tot_w))deallocate(densChi_tot_w)
    if(allocated(densChi_tot_iv))deallocate(densChi_tot_iv)
    !BUILD PAIR SUSCEPTIBILITY
    if(.not.allocated(pairChi_tau))allocate(pairChi_tau(Norb,0:Ltau))
    if(.not.allocated(pairChi_w))  allocate(pairChi_w(Norb,Lreal))
    if(.not.allocated(pairChi_iv)) allocate(pairChi_iv(Norb,0:Lmats))
    pairChi_tau=zero
    pairChi_w=zero
    pairChi_iv=zero
    call build_chi_pair()
    call print_chi_pair()
    if(allocated(pairChi_tau))deallocate(pairChi_tau)
    if(allocated(pairChi_w))deallocate(pairChi_w)
    if(allocated(pairChi_iv))deallocate(pairChi_iv)
    !
    call deallocate_grids
  end subroutine buildchi_impurity


  !+------------------------------------------------------------------+
  !                    SUSCEPTIBILITIES (SPIN, CHARGE, PAIR)
  !+------------------------------------------------------------------+
  include 'ED_GREENS_FUNCTIONS/build_chi_spin.f90'
  include 'ED_GREENS_FUNCTIONS/build_chi_dens.f90'
  include 'ED_GREENS_FUNCTIONS/build_chi_pair.f90'


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_chi_spin
    integer                               :: i,iorb
    integer                               :: unit(3)
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          unit=free_units(3)
          open(unit(1),file="spinChi_orb"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
          open(unit(2),file="spinChi_orb"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
          open(unit(3),file="spinChi_orb"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
          do i=0,Ltau
             write(unit(1),*)tau(i),spinChi_tau(iorb,i)
          enddo
          do i=1,Lreal
             if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(spinChi_w(iorb,i)),dreal(spinChi_w(iorb,i))
          enddo
          do i=0,Lmats
             write(unit(3),*)vm(i),dimag(spinChi_iv(iorb,i)),dreal(spinChi_iv(iorb,i))
          enddo
          do i=1,3
             close(unit(i))
          enddo
       enddo
       if(Norb>1)then
          iorb=Norb+1
          unit=free_units(3)
          open(unit(1),file="spinChi_tot"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
          open(unit(2),file="spinChi_tot"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
          open(unit(3),file="spinChi_tot"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
          do i=0,Ltau
             write(unit(1),*)tau(i),spinChi_tau(iorb,i)
          enddo
          do i=1,Lreal
             if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(spinChi_w(iorb,i)),dreal(spinChi_w(iorb,i))
          enddo
          do i=0,Lmats
             write(unit(3),*)vm(i),dimag(spinChi_iv(iorb,i)),dreal(spinChi_iv(iorb,i))
          enddo
          do i=1,3
             close(unit(i))
          enddo
       endif
    endif
  end subroutine print_chi_spin

  subroutine print_chi_dens
    integer                               :: i,j,iorb,jorb
    integer                               :: unit(3),unit_mix
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          do jorb=iorb,Norb
             unit=free_units(3)
             open(unit(1),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_tau"//reg(ed_file_suffix)//".ed")
             open(unit(2),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw"//reg(ed_file_suffix)//".ed")
             open(unit(3),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw"//reg(ed_file_suffix)//".ed")
             do i=0,Ltau
                write(unit(1),*)tau(i),densChi_tau(iorb,jorb,i)
             enddo
             do i=1,Lreal
                if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_w(iorb,jorb,i)),dreal(densChi_w(iorb,jorb,i))
             enddo
             do i=0,Lmats
                write(unit(3),*)vm(i),dimag(densChi_iv(iorb,jorb,i)),dreal(densChi_iv(iorb,jorb,i))
             enddo
             do i=1,3
                close(unit(i))
             enddo
          enddo
       enddo
    endif
  end subroutine print_chi_dens

  subroutine print_chi_dens_mix
    integer                               :: i,j,iorb,jorb
    integer                               :: unit(3),unit_mix
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          do jorb=1,Norb
             unit=free_units(3)
             open(unit(1),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_tau"//reg(ed_file_suffix)//".ed")
             open(unit(2),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw"//reg(ed_file_suffix)//".ed")
             open(unit(3),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw"//reg(ed_file_suffix)//".ed")
             do i=0,Ltau
                write(unit(1),*)tau(i),densChi_mix_tau(iorb,jorb,i)
             enddo
             do i=1,Lreal
                if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_mix_w(iorb,jorb,i)),dreal(densChi_mix_w(iorb,jorb,i))
             enddo
             do i=0,Lmats
                write(unit(3),*)vm(i),dimag(densChi_mix_iv(iorb,jorb,i)),dreal(densChi_mix_iv(iorb,jorb,i))
             enddo
             do i=1,3
                close(unit(i))
             enddo
          enddo
       enddo
    endif
  end subroutine print_chi_dens_mix

  subroutine print_chi_dens_tot
    integer                               :: i,j,iorb,jorb
    integer                               :: unit(3),unit_mix
    if(ED_MPI_ID==0)then
       unit=free_units(3)
       open(unit(1),file="densChi_tot_tau"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="densChi_tot_realw"//reg(ed_file_suffix)//".ed")
       open(unit(3),file="densChi_tot_iw"//reg(ed_file_suffix)//".ed")
       do i=0,Ltau
          write(unit(1),*)tau(i),densChi_tot_tau(i)
       enddo
       do i=1,Lreal
          if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_tot_w(i)),dreal(densChi_tot_w(i))
       enddo
       do i=0,Lmats
          write(unit(3),*)vm(i),dimag(densChi_tot_iv(i)),dreal(densChi_tot_iv(i))
       enddo
       do i=1,3
          close(unit(i))
       enddo
    endif
  end subroutine print_chi_dens_tot

  subroutine print_chi_pair
    integer                               :: i,iorb
    integer                               :: unit(3)
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          unit=free_units(3)
          open(unit(1),file="pairChi_orb"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
          open(unit(2),file="pairChi_orb"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
          open(unit(3),file="pairChi_orb"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
          do i=0,Ltau
             write(unit(1),*)tau(i),pairChi_tau(iorb,i)
          enddo
          do i=1,Lreal
             if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(pairChi_w(iorb,i)),dreal(pairChi_w(iorb,i))
          enddo
          do i=0,Lmats
             write(unit(3),*)vm(i),dimag(pairChi_iv(iorb,i)),dreal(pairChi_iv(iorb,i))
          enddo
          do i=1,3
             close(unit(i))
          enddo
       enddo
    endif
  end subroutine print_chi_pair










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
