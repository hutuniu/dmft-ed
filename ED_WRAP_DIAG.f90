!###############################################################
! PROGRAM  : RDMFT_WRAP_ED
! PURPOSE  : Contains the main function performin RDMFT calculation
! using the ED solver
! AUTHORS  : Adriano Amaricci, Giacomo Mazza
!###############################################################
module ED_WRAP_DIAG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_DIAG
  USE ED_BATH
  USE SF_TIMER
  USE SF_LINALG
  implicit none
  private

  interface ed_solve_impurity
     module procedure &
          ed_solve_impurity_sites_eloc,&
          ed_solve_impurity_sites_hloc,&
          ed_solve_sc_impurity_sites_eloc,&
          ed_solve_sc_impurity_sites_hloc,&
          ed_solve_impurity_sites_eloc_,&
          ed_solve_impurity_sites_hloc_,&
          ed_solve_sc_impurity_sites_eloc_,&
          ed_solve_sc_impurity_sites_hloc_
  end interface ed_solve_impurity


  public :: ed_solve_impurity
  public :: ed_get_dens
  public :: ed_get_docc
  public :: ed_get_phisc
  public :: ed_get_epot
  !
  public :: init_lattice_baths

  real(8),dimension(:,:),allocatable,save :: nii,dii,pii
  real(8),dimension(:),allocatable,save   :: eii

contains






  !-------------------------------------------------------------------------------------------
  !PURPOSE: solve the impurity problems for each independent
  ! lattice site using ED. 
  !-------------------------------------------------------------------------------------------
  subroutine ed_solve_impurity_sites_eloc(bath,Smats,Sreal,Eloc,Uloc_ii,Ust_ii,Jh_ii)
    real(8)                  :: bath(:,:,:)                                     ![Nlat][Nb(1)][Nb(2)==Nspin]
    complex(8),intent(inout) :: Smats(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats) !
    complex(8),intent(inout) :: Sreal(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal) !
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)                   !A compact version of Hloc_ilat
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI 
    real(8)                  :: bath_tmp(size(bath,2),size(bath,3))
    complex(8)               :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1))
    !
    integer                  :: ilat,i,Nsites,iorb,ispin
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    Nsites=size(bath,1)
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(eii))deallocate(eii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(eii(Nsites))
    !Check the dimensions of the bath are ok:
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    nii   = 0d0 ; nii_tmp   = 0d0
    dii   = 0d0 ; dii_tmp   = 0d0
    eii   = 0d0 ; eii_tmp   = 0d0
    !
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
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
       !
       !Solve the impurity problem for the ilat-th site
       call ed_solver(bath(ilat,:,:))
       !<DEBUG this must be modified to use dedidacated routines rather than direct call
       Smats_tmp(ilat,:,:,:,:,:) = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:) = impSreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
       eii_tmp(ilat)             = ed_EPot+ed_Eknot
       !>DEBUG
       !
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_impurity_sites_eloc

  subroutine ed_solve_impurity_sites_hloc(bath,Smats,Sreal,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    !inputs
    real(8)                  :: bath(:,:,:) ![Nlat][Nb(1)][Nb(2) ==Nspin]
    complex(8),intent(inout) :: Smats(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout) :: Sreal(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI  auxiliary vars
    real(8)                  :: bath_tmp(size(bath,2),size(bath,3))
    complex(8)               :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1))
    ! 
    integer                  :: ilat,i
    integer                  :: Nsites
    integer                  :: Nb1,Nb2
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    ! Check dimensions !
    Nsites=size(bath,1)
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(eii))deallocate(eii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(eii(Nsites))
    !Check the dimensions of the bath are ok:
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    nii   = 0d0 ; nii_tmp   = 0d0
    dii   = 0d0 ; dii_tmp   = 0d0
    eii   = 0d0 ; eii_tmp   = 0d0
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       ! 
       !Solve the impurity problem for the ilat-th site
       call ed_solver(bath(ilat,:,:))
       !<DEBUG this must be modified to use dedidacated routines rather than direct call
       Smats_tmp(ilat,:,:,:,:,:) = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:) = impSreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
       eii_tmp(ilat)             = ed_EPot+ed_Eknot
       !>DEBUG
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_impurity_sites_hloc

  !-------------------------------- SUPERCONDUCTING CASE -------------------------------------
  subroutine ed_solve_sc_impurity_sites_eloc(bath,Smats,Sreal,Eloc,Uloc_ii,Ust_ii,Jh_ii)
    real(8)                  :: bath(:,:,:)
    complex(8),intent(inout) :: Smats(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats) !
    complex(8),intent(inout) :: Sreal(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal) !
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)                   !A compact version of Hloc_ilat
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI
    real(8)                  :: bath_tmp(size(bath,2),size(bath,3))
    complex(8)               :: Smats_tmp(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: pii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1))
    !
    integer                  :: ilat,i,iorb,ispin,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    Nsites=size(bath,1)
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites))
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    nii   = 0d0  ; nii_tmp   = 0d0
    dii   = 0d0  ; dii_tmp   = 0d0
    eii   = 0d0  ; eii_tmp   = 0d0
    pii   = 0d0  ; pii_tmp   = 0d0
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
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
       call ed_solver(bath(ilat,:,:))
       !
       !<DEBUG this must be modified to use dedidacated routines rather than direct call
       Smats_tmp(1,ilat,:,:,:,:,:) = impSmats(:,:,:,:,:)
       Smats_tmp(2,ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       Sreal_tmp(1,ilat,:,:,:,:,:) = impSreal(:,:,:,:,:)
       Sreal_tmp(2,ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
       pii_tmp(ilat,1:Norb)      = ed_phisc(1:Norb)
       eii_tmp(ilat)             = ed_EPot+ed_Eknot
       !>DEBUG
       !
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,2*Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,2*Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_sc_impurity_sites_eloc

  subroutine ed_solve_sc_impurity_sites_hloc(bath,Smats,Sreal,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    real(8)                  :: bath(:,:,:)
    complex(8),intent(inout) :: Smats(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats) !
    complex(8),intent(inout) :: Sreal(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal) !
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI
    real(8)                  :: bath_tmp(size(bath,2),size(bath,3))
    complex(8)               :: Smats_tmp(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: pii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1))
    !
    integer                  :: ilat,i,iorb,ispin,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    Nsites=size(bath,1)
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites))
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    nii   = 0d0  ; nii_tmp   = 0d0
    dii   = 0d0  ; dii_tmp   = 0d0
    eii   = 0d0  ; eii_tmp   = 0d0
    pii   = 0d0  ; pii_tmp   = 0d0
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       call ed_solver(bath(ilat,:,:))
       !
       !<DEBUG this must be modified to use dedidacated routines rather than direct call
       Smats_tmp(1,ilat,:,:,:,:,:) = impSmats(:,:,:,:,:)
       Smats_tmp(2,ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       Sreal_tmp(1,ilat,:,:,:,:,:) = impSreal(:,:,:,:,:)
       Sreal_tmp(2,ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
       pii_tmp(ilat,1:Norb)      = ed_phisc(1:Norb)
       eii_tmp(ilat)             = ed_EPot+ed_Eknot
       !>DEBUG
       !
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,2*Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,2*Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_sc_impurity_sites_hloc

  !-------------------------------- 1-spin, 1-orb CASE -------------------------------------
  subroutine ed_solve_impurity_sites_eloc_(bath,Smats,Sreal,Eloc,Uloc_ii)
    real(8)                  :: bath(:,:,:)                                     ![Nlat][Nb(1)][Nb(2)==Nspin]
    complex(8),intent(inout) :: Smats(size(bath,1),Lmats) !
    complex(8),intent(inout) :: Sreal(size(bath,1),Lreal) !
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)                   !A compact version of Hloc_ilat
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    complex(8)               :: Smats_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    Smats_(:,1,1,1,1,:) = Smats
    Sreal_(:,1,1,1,1,:) = Sreal
    if(present(Eloc).AND.present(Uloc_ii))then
       call ed_solve_impurity_sites_eloc(bath,Smats_,Sreal_,Eloc,Uloc_ii)
    elseif(present(Eloc).AND..not.present(Uloc_ii))then
       call ed_solve_impurity_sites_eloc(bath,Smats_,Sreal_,Eloc)
    elseif(.not.present(Eloc).AND.present(Uloc_ii))then
       call ed_solve_impurity_sites_eloc(bath,Smats_,Sreal_,Uloc_ii)
    else
       call ed_solve_impurity_sites_eloc(bath,Smats_,Sreal_)
    endif
    Smats = Smats_(:,1,1,1,1,:)
    Sreal = Sreal_(:,1,1,1,1,:)
  end subroutine ed_solve_impurity_sites_eloc_

  subroutine ed_solve_impurity_sites_hloc_(bath,Smats,Sreal,Hloc,Uloc_ii)
    real(8)                  :: bath(:,:,:)                                     ![Nlat][Nb(1)][Nb(2)==Nspin]
    complex(8),intent(inout) :: Smats(size(bath,1),Lmats) !
    complex(8),intent(inout) :: Sreal(size(bath,1),Lreal) !
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    complex(8)               :: Smats_(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    Smats_(:,1,1,1,1,:) = Smats
    Sreal_(:,1,1,1,1,:) = Sreal
    if(present(Uloc_ii))then
       call ed_solve_impurity_sites_hloc(bath,Smats_,Sreal_,Hloc,Uloc_ii)
    else
       call ed_solve_impurity_sites_hloc(bath,Smats_,Sreal_,Hloc)
    endif
    Smats = Smats_(:,1,1,1,1,:)
    Sreal = Sreal_(:,1,1,1,1,:)
  end subroutine ed_solve_impurity_sites_hloc_

  !-------------------------------- 1-spin, 1-orb CASE -------------------------------------
  !-------------------------------- SUPERCONDUCTING CASE -------------------------------------
  subroutine ed_solve_sc_impurity_sites_eloc_(bath,Smats,Sreal,Eloc,Uloc_ii)
    real(8)                  :: bath(:,:,:)
    complex(8),intent(inout) :: Smats(2,size(bath,1),Lmats)
    complex(8),intent(inout) :: Sreal(2,size(bath,1),Lreal)
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    complex(8)               :: Smats_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    Smats_(:,:,1,1,1,1,:) = Smats
    Sreal_(:,:,1,1,1,1,:) = Sreal
    if(present(Eloc).AND.present(Uloc_ii))then
       call ed_solve_sc_impurity_sites_eloc(bath,Smats_,Sreal_,Eloc,Uloc_ii)
    elseif(present(Eloc).AND..not.present(Uloc_ii))then
       call ed_solve_sc_impurity_sites_eloc(bath,Smats_,Sreal_,Eloc)
    elseif(.not.present(Eloc).AND.present(Uloc_ii))then
       call ed_solve_sc_impurity_sites_eloc(bath,Smats_,Sreal_,Uloc_ii)
    else
       call ed_solve_sc_impurity_sites_eloc(bath,Smats_,Sreal_)
    endif
    Smats = Smats_(:,:,1,1,1,1,:)
    Sreal = Sreal_(:,:,1,1,1,1,:)
  end subroutine ed_solve_sc_impurity_sites_eloc_

  subroutine ed_solve_sc_impurity_sites_hloc_(bath,Smats,Sreal,Hloc,Uloc_ii)
    real(8)                  :: bath(:,:,:)
    complex(8),intent(inout) :: Smats(2,size(bath,1),Lmats)
    complex(8),intent(inout) :: Sreal(2,size(bath,1),Lreal)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    complex(8)               :: Smats_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    Smats_(:,:,1,1,1,1,:) = Smats
    Sreal_(:,:,1,1,1,1,:) = Sreal
    if(present(Uloc_ii))then
       call ed_solve_sc_impurity_sites_hloc(bath,Smats_,Sreal_,Hloc,Uloc_ii)
    else
       call ed_solve_sc_impurity_sites_hloc(bath,Smats_,Sreal_,Hloc)
    endif
    Smats = Smats_(:,:,1,1,1,1,:)
    Sreal = Sreal_(:,:,1,1,1,1,:)
  end subroutine ed_solve_sc_impurity_sites_hloc_







  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of some local observables
  ! todo: find a more clever way to pass all the variables to the user other than
  ! writing and reading from file (or implement a routine to retrieve from file
  ! reading)
  !+-----------------------------------------------------------------------------+!
  function ed_get_dens(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    if(.not.allocated(nii))then
       yii=0d0
    else
       yii=nii
    endif
  end function ed_get_dens

  function ed_get_docc(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    if(.not.allocated(dii))then
       yii=0d0
    else
       yii=dii
    endif
  end function ed_get_docc

  function ed_get_phisc(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    if(.not.allocated(pii))then
       yii=0d0
    else
       yii=pii
    endif
  end function ed_get_phisc

  function ed_get_epot(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat)      :: yii
    if(.not.allocated(eii))then
       yii=0d0
    else
       yii=eii
    endif
  end function ed_get_epot





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize lattice baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine init_lattice_baths(bath)
    real(8),dimension(:,:,:) :: bath
    integer :: ilat
    logical :: check_dim
    character(len=5) :: tmp_suffix
    if(size(bath,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    do ilat=1,Nlat
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       call init_ed_solver(bath(ilat,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,mpiERR)
  end subroutine init_lattice_baths

end module ED_WRAP_DIAG





