!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
! chi2fit_interface
!
! fit_irred
!  + chi2_anderson/functions
!  + chi2_weiss/functions
!
! fit_irred_SC
!  + chi2_anderson_sc/functions
!  + chi2_weiss_sc/functinos
!
! fit_hybrd
!  + chi2_anderson/functions
!  + chi2_weiss/functions
!########################################################################
MODULE ED_CHI2FIT
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:matrix_inverse
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH

  implicit none
  private

  interface ed_chi2_fitgf
     module procedure chi2_fitgf_generic_normal
     module procedure chi2_fitgf_generic_normal_NOSPIN
     module procedure chi2_fitgf_generic_superc
     module procedure chi2_fitgf_generic_superc_NOSPIN
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf

  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb
  integer,dimension(:),allocatable      :: getIorb,getJorb
  integer                               :: Orb_indx,Spin_indx
  type(effective_bath)                  :: chi2_bath

contains


  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  !
  ! - CHI2_FITGF_GENERIC_NORMAL interface for the normal case 
  !
  !   * CHI2_FITGF_GENERIC_NORMAL_NOSPIN interface to fixed spin input
  !
  ! - CHI2_FITGF_GENERIC_SUPERC interface for the superconducting case 
  !
  !   * CHI2_FITGF_GENERIC_SUPERC_NOSPIN interface to fixed spin input
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_generic_normal_NOSPIN(fg,bath,ispin)
    complex(8),dimension(:,:,:)                      :: fg ![Norb][Norb][Niw]
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lfit) :: fg_
    real(8),dimension(:),intent(inout)               :: bath
    integer,optional                                 :: ispin
    integer                                          :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    if(size(fg,3)<Lfit)stop "chi2_fitgf_generic_normal_NOSPIN error: size[fg,3] < Lfit"
    fg_=zero
    fg_(ispin_,ispin_,:,:,1:Lfit) = fg(:,:,1:Lfit)
    call chi2_fitgf_generic_normal(fg_,bath,ispin_)
  end subroutine chi2_fitgf_generic_normal_NOSPIN

  subroutine chi2_fitgf_generic_normal(fg,bath,ispin)
    complex(8),dimension(:,:,:,:,:) :: fg ![Nspin][Nspin][Norb][Norb][Niw]
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin
    integer                         :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ED_MPI_ID==0)then
       if(size(fg,1)/=Nspin)stop"chi2_fitgf_generic_normal error: size[fg,1] != Nspin"
       if(size(fg,2)/=Nspin)stop"chi2_fitgf_generic_normal error: size[fg,2] != Nspin"
       if(size(fg,3)/=Norb)stop"chi2_fitgf_generic_normal error: size[fg,3] != Norb"
       if(size(fg,4)/=Norb)stop"chi2_fitgf_generic_normal error: size[fg,4] != Norb"
       select case(cg_method)
       case (0)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-nr"
       case (1)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-minimize"
       case(2)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-plus"
       case default
          stop "chi2_fitgf_generic_normal error: cg_method > 2"
       end select
    endif
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("superc")
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=superc but only NORMAL component provided"
          call sleep(1)
          call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("nonsu2")
          if(present(ispin))then
             write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but only ONE spin orientation required. disregarded"
             call sleep(1)
          endif
          call chi2_fitgf_normal_nonsu2(fg(:,:,:,:,:),bath)
          !
       end select
       !
    case ("hybrid")
       select case(ed_mode)
       case default
          call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("superc")
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=superc but only NORMAL component provided"
          call sleep(1)
          call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("nonsu2")
          if(present(ispin))then
             write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but only ONE spin orientation required. disregarded"
             call sleep(1)
          endif
          call chi2_fitgf_hybrid_nonsu2(fg(:,:,:,:,:),bath)
          !
       end select
    end select
#ifdef _MPI
    call MPI_BCAST(bath,size(bath),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_normal

  subroutine chi2_fitgf_generic_superc_NOSPIN(fg,bath,ispin)
    complex(8),dimension(:,:,:,:)                      :: fg ![2][Norb][Norb][Niw]
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lfit) :: fg_
    real(8),dimension(:),intent(inout)                 :: bath
    integer,optional                                   :: ispin
    integer                                            :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    if(size(fg,4)<Lfit)stop "chi2_fitgf_generic_superc_NOSPIN error: size[fg,4] < Lfit"
    fg_=zero
    fg_(:,ispin_,ispin_,:,:,1:Lfit) = fg(:,:,:,1:Lfit)
  end subroutine chi2_fitgf_generic_superc_NOSPIN

  subroutine chi2_fitgf_generic_superc(fg,bath,ispin)
    complex(8),dimension(:,:,:,:,:,:) :: fg ![2][Nspin][Nspin][Norb][Norb][Niw]
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin
    integer                         :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ED_MPI_ID==0)then
       if(size(fg,1)/=2)stop"chi2_fitgf_generic_normal error: size[fg,1] != 2"
       if(size(fg,2)/=Nspin)stop"chi2_fitgf_generic_normal error: size[fg,2] != Nspin"
       if(size(fg,3)/=Nspin)stop"chi2_fitgf_generic_normal error: size[fg,3] != Nspin"
       if(size(fg,4)/=Norb)stop"chi2_fitgf_generic_normal error: size[fg,4] != Norb"
       if(size(fg,5)/=Norb)stop"chi2_fitgf_generic_normal error: size[fg,5] != Norb"
       select case(cg_method)
       case (0)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-nr"
       case (1)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-minimize"
       case(2)
          if(ed_verbose<3)write(LOGfile,"(A)")"\Chi2 fit with CG-plus"
       case default
          stop "chi2_fitgf_generic_superc error: cg_method > 2"
       end select
    endif
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=normal but NORMAL & ANOMAL components provided. Fitting only the NORMAL."
          call sleep(1)
          call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("superc")
          call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("nonsu2")
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but NORMAL & ANOMAL components provided. Fitting only the NORMAL."
          call sleep(1)
          if(present(ispin))then
             write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but only ONE spin orientation required. Disregarded"
             call sleep(1)
          endif
          call chi2_fitgf_normal_nonsu2(fg(1,:,:,:,:,:),bath)
          !
       end select
       !
    case ("hybrid")
       select case(ed_mode)
       case default
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=normal but NORMAL & ANOMAL components provided. Fitting only NORMAL."
          call sleep(1)
          call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("superc")
          stop "chi2_fitgf ERROR: ed_mode=superc and bath_type=hybrid is not yet implemented."
          stop 'Error: Hybrid bath + SC is not implemented yet: ask the developer...'
          !
       case ("nonsu2")
          write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but NORMAL & ANOMAL components provided. Fitting only NORMAL."
          call sleep(1)
          if(present(ispin))then
             write(LOGfile,"(A)")"chi2_fitgf WARNING: ed_mode=nonsu2 but only ONE spin orientation required. Disregarded"
             call sleep(1)
          endif
          call chi2_fitgf_hybrid_nonsu2(fg(1,:,:,:,:,:),bath)
          !
       end select
    end select
#ifdef _MPI
    call MPI_BCAST(bath,size(bath),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_superc




  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !normal ED_bath
  include "ed_chi2_fitgf_normal_normal.f90"
  include "ed_chi2_fitgf_normal_superc.f90"
  include "ed_chi2_fitgf_normal_nonsu2.f90"

  !hybrid ED_bath
  include "ed_chi2_fitgf_hybrid_normal.f90"  
  !include "ed_chi2_fitgf_hybrid_superc.f90"  
  include "ed_chi2_fitgf_hybrid_nonsu2.f90"  
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

end MODULE ED_CHI2FIT
