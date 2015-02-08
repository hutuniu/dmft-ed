!###############################################################
! PROGRAM  : RDMFT_WRAP_ED
! PURPOSE  : Contains the main function performin RDMFT calculation
! using the ED solver
! AUTHORS  : Adriano Amaricci, Giacomo Mazza
!###############################################################
module ED_WRAP_CHI2FIT
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_CHI2FIT
  USE SF_TIMER
  USE SF_LINALG
  implicit none
  private


  interface ed_fit_bath
     module procedure &
          ed_fit_bath_sites_eloc,&
          ed_fit_bath_sites_hloc,&
          ed_fit_bath_sc_sites_eloc,&
          ed_fit_bath_sc_sites_hloc,&
          ed_fit_bath_sites_eloc_,&
          ed_fit_bath_sites_hloc_,&
          ed_fit_bath_sc_sites_eloc_,&
          ed_fit_bath_sc_sites_hloc_
  end interface ed_fit_bath


  public :: ed_fit_bath


contains


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: given a number of independent baths, evaluate N independent Delta/G0
  ! functions and fit them to update the effective baths for ED.
  ! As before we treat differently (for back-compatibility reasons) the two cases
  ! Eloc and Hloc as input
  !+-----------------------------------------------------------------------------+!
  subroutine ed_fit_bath_sites_eloc(bath,Delta,Eloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)
    integer,optional         :: spin
    !MPI auxiliary 
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2),size(bath,3))
    integer                  :: i,ilat,ispin,iorb,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    Nsites=size(bath,1)
    bath_tmp=0d0
    do ilat=1+mpiID,Nsites,mpiSIZE
       bath_tmp(ilat,:,:)=bath(ilat,:,:)
       !Set the local part of the Hamiltonian.
       Hloc=zero
       if(present(Eloc))then
          do ispin=1,Nspin
             do iorb=1,Norb
                i=iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                Hloc(ispin,ispin,iorb,iorb) = Eloc(i)
             enddo
          enddo
       endif
       call set_Hloc(Hloc)
       write(tmp_suffix,'(I4.4)')ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(spin))then
          ispin=spin
          if(ispin>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call chi2_fitgf(Delta(ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
       else
          do ispin=1,Nspin
             call chi2_fitgf(Delta(ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
          enddo
       end if
    end do
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_fit_bath_sites_eloc

  subroutine ed_fit_bath_sites_hloc(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    !MPI auxiliary vars
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2),size(bath,3))
    integer                  :: ilat,i,iorb,ispin
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    ! Check dimensions !
    Nsites=size(bath,1)
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    bath_tmp=0d0
    do ilat=1+mpiID,Nsites,mpiSIZE
       bath_tmp(ilat,:,:)=bath(ilat,:,:)
       call set_Hloc(Hloc(ilat,:,:,:,:))
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(spin))then
          ispin=spin
          if(ispin>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call chi2_fitgf(Delta(ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin=ispin)
       else
          do ispin=1,Nspin
             call chi2_fitgf(Delta(ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin=ispin)
          enddo
       end if
    end do
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_fit_bath_sites_hloc

  !-------------------------------- SUPERCONDUCTING CASE -------------------------------------
  subroutine ed_fit_bath_sc_sites_eloc(bath,Delta,Eloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)
    integer,optional         :: spin
    !MPI auxiliary
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2),size(bath,3))
    integer                  :: ilat,i,iorb,ispin,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    Nsites=size(bath,1)
    bath_tmp=0.d0
    !+- GET INDEPENDENT SITES HYBRIDIZATION FUNCTION AND FIT THE BATHS -+!
    do ilat=1+mpiID,Nsites,mpiSIZE
       bath_tmp(ilat,:,:) = bath(ilat,:,:)
       !Set the local part of the Hamiltonian.
       Hloc=zero
       if(present(Eloc))then
          do ispin=1,Nspin
             do iorb=1,Norb
                i=iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                Hloc(ispin,ispin,iorb,iorb) = Eloc(i)
             enddo
          enddo
       endif
       call set_Hloc(Hloc)
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(spin))then
          ispin=spin
          if(ispin>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call chi2_fitgf(Delta(:,ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
       else
          do ispin=1,Nspin
             call chi2_fitgf(Delta(:,ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
          enddo
       endif
    end do
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_fit_bath_sc_sites_eloc

  subroutine ed_fit_bath_sc_sites_hloc(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    !MPI auxiliary
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2),size(bath,3))
    integer                  :: ilat,i,iorb,ispin,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    Nsites=size(bath,1)
    bath_tmp=0.d0
    !+- GET INDEPENDENT SITES HYBRIDIZATION FUNCTION AND FIT THE BATHS -+!
    do ilat=1+mpiID,Nsites,mpiSIZE
       bath_tmp(ilat,:,:) = bath(ilat,:,:)
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(spin))then
          ispin=spin
          if(ispin>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call chi2_fitgf(Delta(:,ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
       else
          do ispin=1,Nspin
             call chi2_fitgf(Delta(:,ilat,ispin,ispin,:,:,:),bath_tmp(ilat,:,:),ispin)
          enddo
       endif
    end do
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_fit_bath_sc_sites_hloc




  !-------------------------------- 1-spin, 1-orb CASE -------------------------------------
  subroutine ed_fit_bath_sites_eloc_(bath,Delta,Eloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Lmats)
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)
    integer,optional         :: spin
    complex(8)               :: Delta_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    Delta_(:,1,1,1,1,:) = Delta
    if(present(Eloc).AND.present(spin))then
       call ed_fit_bath_sites_eloc(bath,Delta_,Eloc,spin)
    elseif(present(Eloc).AND..not.present(spin))then
       call ed_fit_bath_sites_eloc(bath,Delta_,Eloc)
    elseif(.not.present(Eloc).AND.present(spin))then
       call ed_fit_bath_sites_eloc(bath,Delta_,spin=spin)
    else
       call ed_fit_bath_sites_eloc(bath,Delta_)
    endif
  end subroutine ed_fit_bath_sites_eloc_

  subroutine ed_fit_bath_sites_hloc_(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    complex(8)               :: Delta_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    Delta_(:,1,1,1,1,:) = Delta
    if(present(spin))then
       call ed_fit_bath_sites_hloc(bath,Delta_,Hloc,spin)
    else
       call ed_fit_bath_sites_hloc(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sites_hloc_

  !-------------------------------- 1-spin, 1-orb CASE -------------------------------------
  !-------------------------------- SUPERCONDUCTING CASE -------------------------------------
  subroutine ed_fit_bath_sc_sites_eloc_(bath,Delta,Eloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Lmats)
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)
    integer,optional         :: spin
    complex(8)               :: Delta_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    Delta_(:,:,1,1,1,1,:) = Delta
    if(present(Eloc).AND.present(spin))then
       call ed_fit_bath_sc_sites_eloc(bath,Delta_,Eloc,spin)
    elseif(present(Eloc).AND..not.present(spin))then
       call ed_fit_bath_sc_sites_eloc(bath,Delta_,Eloc)
    elseif(.not.present(Eloc).AND.present(spin))then
       call ed_fit_bath_sc_sites_eloc(bath,Delta_,spin=spin)
    else
       call ed_fit_bath_sc_sites_eloc(bath,Delta_)
    endif
  end subroutine ed_fit_bath_sc_sites_eloc_

  subroutine ed_fit_bath_sc_sites_hloc_(bath,Delta,Hloc,spin)
    real(8),intent(inout)    :: bath(:,:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: spin
    complex(8)               :: Delta_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    Delta_(:,:,1,1,1,1,:) = Delta
    if(present(spin))then
       call ed_fit_bath_sc_sites_hloc(bath,Delta_,Hloc,spin)
    else
       call ed_fit_bath_sc_sites_hloc(bath,Delta_,Hloc)
    endif
  end subroutine ed_fit_bath_sc_sites_hloc_




end module ED_WRAP_CHI2FIT





