!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN
  USE SF_CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !Get sparse sector Hamiltonian
  public :: build_H_normal_c
  public :: build_H_superc_c
  public :: build_H_nonsu2_c

  public :: ed_hamiltonian_set_MPI
  public :: ed_hamiltonian_del_MPI

#ifdef _MPI
  integer                      :: MpiComm=MPI_UNDEFINED
#else
  integer                      :: MpiComm=0
#endif
  logical                      :: MpiStatus=.false.
  integer                      :: MpiIerr
  integer                      :: MpiRank=0
  integer                      :: MpiSize=1
  integer                      :: mpiQ=1
  integer                      :: mpiR=0


contains


  subroutine ed_hamiltonian_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
    MpiStatus=.true.
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
#else
    integer,optional :: comm_
#endif
  end subroutine ed_hamiltonian_set_MPI


  subroutine ed_hamiltonian_del_MPI()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
    MpiStatus=.false.
#endif
  end subroutine ed_hamiltonian_del_MPI






  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix NORMAL CASE
  !+------------------------------------------------------------------+
  !DOUBLE COMPLEX
  subroutine build_H_normal_c(isector,Hmat)
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    integer                                :: isector
    type(sector_map)                       :: H,Hup,Hdw
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    !
    call build_sector(isector,H)
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !
    dim=getdim(isector)
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    call sp_init_matrix(spH0,mpiQ + mpiR)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_normal.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_c ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=zero
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_normal_c






  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix SUPERC CASE
  !+------------------------------------------------------------------+
  !DOUBLE COMPLEX
  subroutine build_H_superc_c(isector,Hmat)
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    integer                                :: isector
    type(sector_map)                       :: H,Hup,Hdw
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    !
    call build_sector(isector,H)
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !
    dim=getdim(isector)
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    call sp_init_matrix(spH0,mpiQ + mpiR)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_superc.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_superc_c ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=zero
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_superc_c



















  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix NONSU2 CASE
  !+------------------------------------------------------------------+
  !DOUBLE COMPLEX
  subroutine build_H_nonsu2_c(isector,Hmat)
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    integer                                :: isector
    type(sector_map)                       :: H,Hup,Hdw
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw    
    !
    !
    call build_sector(isector,H)
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !
    dim=getdim(isector)
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    call sp_init_matrix(spH0,mpiQ + mpiR)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_nonsu2.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_nonsu2_c ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=zero
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_nonsu2_c



end MODULE ED_HAMILTONIAN








! !DOUBLE PRECISION
! subroutine build_H_normal_d(isector,Hmat)
!   real(8),dimension(:,:),optional     :: Hmat
!   real(8),dimension(:,:),allocatable  :: Hredux
!   integer                             :: isector
!   type(sector_map)                    :: H,Hup,Hdw
!   integer,dimension(Nlevels)          :: ib
!   integer,dimension(Ns)               :: ibup,ibdw
!   integer                             :: dim,dimUp,dimDw
!   integer                             :: i,iup,idw
!   integer                             :: m,mup,mdw
!   integer                             :: ishift,ishift_up,ishift_dw
!   integer                             :: j,ms,impi
!   integer                             :: iorb,jorb,ispin,jspin,ibath
!   integer                             :: kp,k1,k2,k3,k4
!   integer                             :: alfa,beta
!   real(8)                             :: sg1,sg2,sg3,sg4
!   real(8),dimension(Norb)             :: nup,ndw
!   real(8)                             :: htmp,htmpup,htmpdw
!   real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!   logical                             :: Jcondition
!   integer                             :: first_state,last_state
!   integer                             :: first_state_up,last_state_up
!   integer                             :: first_state_dw,last_state_dw
!   !
!   !
!   call build_sector(isector,H)
!   !
!   if(spH0%status)call sp_delete_matrix(spH0) 
!   !
!   dim=getdim(isector)
!   mpiQ = dim/MpiSize
!   mpiR = 0
!   if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!   call sp_init_matrix(spH0,mpiQ + mpiR)
!   ishift      = MpiRank*mpiQ
!   first_state = MpiRank*mpiQ + 1
!   last_state  = (MpiRank+1)*mpiQ + mpiR
!   !
!   !
!   !Get diagonal hybridization
!   diag_hybr=0.0d0
!   if(bath_type/="replica")then
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
!            enddo
!         enddo
!      enddo
!   else
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
!            enddo
!         enddo
!      enddo
!   endif
!   !
!   !-----------------------------------------------!
!   include "ED_HAMILTONIAN/build_h_normal.f90"
!   !-----------------------------------------------!
!   !
!   deallocate(H%map)
!   !
!   if(present(Hmat))then
!      if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_d ERROR: size(Hmat) != dim**2"
!      if(MpiStatus)then
!         allocate(Hredux(dim,dim));Hredux=0d0
!         call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!         call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
! #endif
!      else
!         call sp_dump_matrix(spH0,Hmat)
!      endif
!   endif
!   !
! end subroutine build_H_normal_d





! !DOUBLE PRECISION
! subroutine build_H_superc_d(isector,Hmat)
!   real(8),dimension(:,:),optional     :: Hmat
!   real(8),dimension(:,:),allocatable  :: Hredux
!   integer                             :: isector
!   type(sector_map)                    :: H,Hup,Hdw
!   integer,dimension(Nlevels)          :: ib
!   integer,dimension(Ns)               :: ibup,ibdw
!   integer                             :: dim,dimUp,dimDw
!   integer                             :: i,iup,idw
!   integer                             :: m,mup,mdw
!   integer                             :: ishift,ishift_up,ishift_dw
!   integer                             :: j,ms,impi
!   integer                             :: iorb,jorb,ispin,jspin,ibath
!   integer                             :: kp,k1,k2,k3,k4
!   integer                             :: alfa,beta
!   real(8)                             :: sg1,sg2,sg3,sg4
!   real(8),dimension(Norb)             :: nup,ndw
!   real(8)                             :: htmp,htmpup,htmpdw
!   real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!   logical                             :: Jcondition
!   integer                             :: first_state,last_state
!   integer                             :: first_state_up,last_state_up
!   integer                             :: first_state_dw,last_state_dw
!   !
!   !
!   call build_sector(isector,H)
!   !
!   if(spH0%status)call sp_delete_matrix(spH0) 
!   !
!   dim=getdim(isector)
!   mpiQ = dim/MpiSize
!   mpiR = 0
!   if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!   call sp_init_matrix(spH0,mpiQ + mpiR)
!   ishift      = MpiRank*mpiQ
!   first_state = MpiRank*mpiQ + 1
!   last_state  = (MpiRank+1)*mpiQ + mpiR
!   !
!   !
!   !Get diagonal hybridization
!   diag_hybr=0.0d0
!   if(bath_type/="replica")then
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
!            enddo
!         enddo
!      enddo
!   else
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
!            enddo
!         enddo
!      enddo
!   endif
!   !
!   !-----------------------------------------------!
!   include "ED_HAMILTONIAN/build_h_superc.f90"
!   !-----------------------------------------------!
!   !
!   deallocate(H%map)
!   !
!   if(present(Hmat))then
!      if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_superc_d ERROR: size(Hmat) != dim**2"
!      if(MpiStatus)then
!         allocate(Hredux(dim,dim));Hredux=0d0
!         call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!         call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
! #endif
!      else
!         call sp_dump_matrix(spH0,Hmat)
!      endif
!   endif
!   !
! end subroutine build_H_superc_d





!   !DOUBLE PRECISION
! subroutine build_H_nonsu2_d(isector,Hmat)
!   real(8),dimension(:,:),optional     :: Hmat
!   real(8),dimension(:,:),allocatable  :: Hredux
!   integer                             :: isector
!   type(sector_map)                    :: H,Hup,Hdw
!   integer,dimension(Nlevels)          :: ib
!   integer,dimension(Ns)               :: ibup,ibdw
!   integer                             :: dim,dimUp,dimDw
!   integer                             :: i,iup,idw
!   integer                             :: m,mup,mdw
!   integer                             :: ishift,ishift_up,ishift_dw
!   integer                             :: j,ms,impi
!   integer                             :: iorb,jorb,ispin,jspin,ibath
!   integer                             :: kp,k1,k2,k3,k4
!   integer                             :: alfa,beta
!   real(8)                             :: sg1,sg2,sg3,sg4
!   real(8),dimension(Norb)             :: nup,ndw
!   real(8)                             :: htmp,htmpup,htmpdw
!   real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!   logical                             :: Jcondition
!   integer                             :: first_state,last_state
!   integer                             :: first_state_up,last_state_up
!   integer                             :: first_state_dw,last_state_dw    
!   !
!   !
!   call build_sector(isector,H)
!   !
!   if(spH0%status)call sp_delete_matrix(spH0) 
!   !
!   dim=getdim(isector)
!   mpiQ = dim/MpiSize
!   mpiR = 0
!   if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!   call sp_init_matrix(spH0,mpiQ + mpiR)
!   ishift      = MpiRank*mpiQ
!   first_state = MpiRank*mpiQ + 1
!   last_state  = (MpiRank+1)*mpiQ + mpiR
!   !
!   !
!   !Get diagonal hybridization
!   diag_hybr=0.0d0
!   if(bath_type/="replica")then
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
!            enddo
!         enddo
!      enddo
!   else
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do iorb=1,Norb
!               diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
!            enddo
!         enddo
!      enddo
!   endif
!   !
!   !-----------------------------------------------!
!   include "ED_HAMILTONIAN/build_h_nonsu2.f90"
!   !-----------------------------------------------!
!   !
!   deallocate(H%map)
!   !
!   if(present(Hmat))then
!      if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_nonsu2_d ERROR: size(Hmat) != dim**2"
!      if(MpiStatus)then
!         allocate(Hredux(dim,dim));Hredux=0d0
!         call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!         call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
! #endif
!      else
!         call sp_dump_matrix(spH0,Hmat)
!      endif
!   endif
!   !
! end subroutine build_H_nonsu2_d






!   !DOUBLE COMPLEX
!   subroutine build_H_normal_c(isector,Hmat)
!     complex(8),dimension(:,:),optional    :: Hmat
!     complex(8),dimension(:,:),allocatable :: Hredux
!     integer                               :: isector
!     type(sector_map)                      :: H,Hup,Hdw
!     integer,dimension(Nlevels)            :: ib
!     integer                               :: dim
!     integer                               :: i,j,m,ms,impi,ishift
!     integer                               :: iorb,jorb,ispin,jspin,ibath
!     integer                               :: kp,k1,k2,k3,k4
!     integer                               :: alfa,beta
!     real(8)                               :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)               :: nup,ndw
!     complex(8)                            :: htmp
!     complex(8),dimension(Nspin,Norb)      :: eloc
!     complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
!     logical                               :: Jcondition
!     integer                               :: first_state,last_state
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0)
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !Get diagonal part of Hloc
!     do ispin=1,Nspin
!        do iorb=1,Norb
!           eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
!        enddo
!     enddo
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_normal.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_normal_c






!   !+------------------------------------------------------------------+
!   !PURPOSE  : Build Hamiltonian sparse matrix SUPERC CASE
!   !+------------------------------------------------------------------+
!   !DOUBLE COMPLEX
!   subroutine build_H_superc_c(isector,Hmat)
!     complex(8),dimension(:,:),optional    :: Hmat
!     complex(8),dimension(:,:),allocatable :: Hredux
!     integer                               :: isector
!     type(sector_map)                      :: H,Hup,Hdw
!     integer,dimension(Nlevels)            :: ib
!     integer                               :: dim
!     integer                               :: i,j,m,ms,impi,ishift
!     integer                               :: iorb,jorb,ispin,jspin,ibath
!     integer                               :: kp,k1,k2,k3,k4
!     integer                               :: alfa,beta
!     real(8)                               :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)               :: nup,ndw
!     complex(8)                            :: htmp
!     complex(8),dimension(Nspin,Norb)      :: eloc
!     complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
!     logical                               :: Jcondition
!     integer                               :: first_state,last_state
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0)
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !Get diagonal part of Hloc
!     do ispin=1,Nspin
!        do iorb=1,Norb
!           eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
!        enddo
!     enddo
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_superc.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_superc_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_superc_c



















!   !+------------------------------------------------------------------+
!   !PURPOSE  : Build Hamiltonian sparse matrix NONSU2 CASE
!   !+------------------------------------------------------------------+
!   !DOUBLE COMPLEX
!   subroutine build_H_nonsu2_c(isector,Hmat)
!     complex(8),dimension(:,:),optional    :: Hmat
!     complex(8),dimension(:,:),allocatable :: Hredux
!     integer                               :: isector
!     type(sector_map)                      :: H,Hup,Hdw
!     integer,dimension(Nlevels)            :: ib
!     integer                               :: dim
!     integer                               :: i,j,m,ms,impi,ishift
!     integer                               :: iorb,jorb,ispin,jspin,ibath
!     integer                               :: kp,k1,k2,k3,k4
!     integer                               :: alfa,beta
!     real(8)                               :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)               :: nup,ndw
!     complex(8)                            :: htmp
!     complex(8),dimension(Nspin,Norb)      :: eloc
!     complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
!     logical                               :: Jcondition
!     integer                               :: first_state,last_state
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0)
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !Get diagonal part of Hloc
!     do ispin=1,Nspin
!        do iorb=1,Norb
!           eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
!        enddo
!     enddo
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_nonsu2.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_nonsu2_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_nonsu2_c
