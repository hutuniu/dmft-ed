module ED_WRAP_CHI2FIT
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_CHI2FIT
  USE SF_TIMER
  USE SF_LINALG
  implicit none
  private



contains


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: given a number of independent baths, evaluate N independent Delta/G0
  ! functions and fit them to update the effective baths for ED.
  ! As before we treat differently (for back-compatibility reasons) the two cases
  ! Hloc as input
  !+-----------------------------------------------------------------------------+!
  subroutine ed_fit_bath_sites_normal(bath,Delta,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary vars
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0d0
    do ilat=1+mpiID,Nsites,mpiSIZE
       bath_tmp(ilat,:)=bath(ilat,:)
       call set_Hloc(Hloc(ilat,:,:,:,:))
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(ilat,ispin_,ispin_,:,:,:),bath_tmp(ilat,:),ispin=ispin_)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(Delta(ilat,ispin_,ispin_,:,:,:),bath_tmp(ilat,:),ispin=ispin_)
          enddo
       end if
    end do
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    bath = bath_tmp
#endif
  end subroutine ed_fit_bath_sites_normal

  subroutine ed_fit_bath_sites_normal_1b(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    complex(8)               :: Delta_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    if(Norb>1)stop "ed_fit_bath_sites_hloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_fit_bath_sites_hloc_1b error: Nspin > 1 in 1-band routine" 
    Delta_(:,1,1,1,1,:) = Delta
    if(present(spin))then
       call ed_fit_bath_sites_normal(bath,Delta_,Hloc,spin)
    else
       call ed_fit_bath_sites_normal(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sites_normal_1b

  subroutine ed_fit_bath_sites_normal_mb(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    complex(8)               :: Delta_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    if(Nspin>1)stop "ed_fit_bath_sites_hloc_mb error: Nspin > 1 in M-band routine" 
    Delta_(:,1,1,:,:,:) = Delta
    if(present(spin))then
       call ed_fit_bath_sites_normal(bath,Delta_,Hloc,spin)
    else
       call ed_fit_bath_sites_normal(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sites_normal_mb













  !-------------------------------- SUPERCONDUCTING CASE -------------------------------------
  subroutine ed_fit_bath_sites_superc(bath,Delta,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    Nsites=size(bath,1)
    !
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0.d0
    !
    do ilat=1+mpiID,Nsites,mpiSIZE
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       !
       call set_superc(Hloc(ilat,:,:,:,:))
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(:,ilat,ispin_,ispin_,:,:,:),bath_tmp(ilat,:),ispin_)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(Delta(:,ilat,ispin_,ispin_,:,:,:),bath_tmp(ilat,:),ispin_)
          enddo
       endif
    end do
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    bath = bath_tmp
#endif
  end subroutine ed_fit_bath_sites_superc

  subroutine ed_fit_bath_sites_superc_1b(bath,Delta,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    complex(8)               :: Delta_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    if(Norb>1)stop "ed_fit_bath_sites_superc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_fit_bath_sites_superc_1b error: Nspin > 1 in 1-band routine" 
    Delta_(:,:,1,1,1,1,:) = Delta
    if(present(ispin))then
       call ed_fit_bath_sites_superc(bath,Delta_,Hloc,ispin)
    else
       call ed_fit_bath_sites_superc(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sites_superc_1b

  subroutine ed_fit_bath_sites_superc_mb(bath,Delta,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    complex(8)               :: Delta_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    if(Nspin>1)stop "ed_fit_bath_sites_superc_mb error: Nspin > 1 in M-band routine" 
    Delta_(:,:,1,1,:,:,:) = Delta
    if(present(ispin))then
       call ed_fit_bath_sites_superc(bath,Delta_,Hloc,ispin)
    else
       call ed_fit_bath_sites_superc(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sites_superc_mb








end module ED_WRAP_CHI2FIT





