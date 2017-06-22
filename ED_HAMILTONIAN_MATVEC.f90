!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN_MATVEC
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

  !>build sparse hamiltonian of the sector
  public  :: ed_buildH_c
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif
  !
  !
  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_cc
#endif
  !
  !
  !> Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector


  !> MPI local variables (shared)
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
  !
  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hup,Hdw






contains


  !####################################################################
  !                        AUXILIARY ROUTINES
  !####################################################################
  subroutine ed_hamiltonian_matvec_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
    MpiStatus=.true.
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
#else
    integer,optional :: comm_
#endif
  end subroutine ed_hamiltonian_matvec_set_MPI


  subroutine ed_hamiltonian_matvec_del_MPI()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
#else
    MpiComm = 0
#endif
    MpiStatus=.false.
    MpiRank=0
    MpiSize=1
    MpiQ=1
    MpiR=0
  end subroutine ed_hamiltonian_matvec_del_MPI


  subroutine setup_Hv_sector(isector)
    integer                   :: isector
    Hsector=isector
    Hstatus=.true.
    call build_sector(isector,H)
  end subroutine setup_Hv_sector


  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    Hsector=0
    Hstatus=.false.
  end subroutine delete_Hv_sector








  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(Hmat)
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    integer                                :: isector
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
    if(.not.Hstatus)stop "ed_buildH_c ERROR: Hsector NOT set"
    isector=Hsector
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
    include "ED_HAMILTONIAN_MATVEC/build_h.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "ed_buildH_c ERROR: size(Hmat) != dim**2"
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
  end subroutine ed_buildH_c














  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_cc(Nloc,v,Hv)
    integer                      :: Nloc
    complex(8),dimension(Nloc)   :: v
    complex(8),dimension(Nloc)   :: Hv
    integer                      :: i
    type(sparse_element),pointer :: c
    Hv=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_cc


#ifdef _MPI
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    type(sparse_element),pointer        :: c
    N=0
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    MpiSize = get_Size_MPI(MpiComm)
    mpiQ = get_Q_MPI(MpiComm,N)
    mpiR = get_R_MPI(MpiComm,N)
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
    vin                   = zero
    SendCounts(0:)        = mpiQ
    SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,MpiIerr)
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_mpi_cc
#endif












  !####################################################################
  !            SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
  !####################################################################
  subroutine directMatVec_cc(Nloc,vin,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: vin
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: isector
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
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    dim=getdim(isector)
    if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
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
    Hv=zero
    !-----------------------------------------------!
    include "ED_HAMILTONIAN_MATVEC/build_hxv.f90"
    !-----------------------------------------------!
    !
  end subroutine directMatVec_cc



#ifdef _MPI
  subroutine directMatVec_MPI_cc(Nloc,v,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: v
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: N
    complex(8),dimension(:),allocatable    :: vin
    integer,allocatable,dimension(:)       :: SendCounts,Displs
    integer                                :: isector
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
    if(.not.Hstatus)stop "directMatVec_MPI_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    dim=getdim(isector)
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
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    N=0
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
    if(N/=dim)stop "directMatVec_MPI_cc ERROR: N != dim(isector)"
    !
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
    vin                   = zero
    SendCounts(0:)        = mpiQ
    SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,MpiIerr)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN_MATVEC/build_hxv.f90"
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_cc
#endif








end MODULE ED_HAMILTONIAN_MATVEC







!   !+------------------------------------------------------------------+
!   !>NORMAL CASE
!   !+------------------------------------------------------------------+
!   subroutine build_H_normal_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
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
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
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
!   !>SUPERC CASE
!   !+------------------------------------------------------------------+
!   !DOUBLE COMPLEX
!   subroutine build_H_superc_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
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
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
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
!   !>NONSU2 CASE
!   !+------------------------------------------------------------------+
!   subroutine build_H_nonsu2_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw    
!     !
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
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
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
!   !+------------------------------------------------------------------+
!   !+------------------------------------------------------------------+
!   !+------------------------------------------------------------------+
