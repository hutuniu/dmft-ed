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

  !>build sparse hamiltonian of the sector
  public  :: build_H_normal_d
  public  :: build_H_superc_d
  public  :: build_H_nonsu2_d
  !
  public  :: build_H_normal_c
  public  :: build_H_superc_c
  public  :: build_H_nonsu2_c
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_dd
  public  :: spMatVec_dc
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_dd
  public  :: spMatVec_MPI_dc
  public  :: spMatVec_MPI_cc
#endif
  !
  !
  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_dd
  public  :: directMatVec_dc
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_dd
  public  :: directMatVec_MPI_dc
  public  :: directMatVec_MPI_cc
#endif
  !
  !
  !> Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI


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
  type(sector_map)             :: H,Hup,Hdw
  type(sparse_element),pointer :: c






contains











  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  !DOUBLE PRECISION:
  !
  !+------------------------------------------------------------------+
  !>NORMAL CASE
  !+------------------------------------------------------------------+
  subroutine build_H_normal_d(isector,Hmat)
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    integer                             :: isector
    integer,dimension(Nlevels)          :: ib
    integer,dimension(Ns)               :: ibup,ibdw
    integer                             :: dim,dimUp,dimDw
    integer                             :: i,iup,idw
    integer                             :: m,mup,mdw
    integer                             :: ishift,ishift_up,ishift_dw
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp,htmpup,htmpdw
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    integer                             :: first_state_up,last_state_up
    integer                             :: first_state_dw,last_state_dw
    !
    !
    call build_sector(isector,H)
    dim = getdim(isector)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !

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
    diag_hybr=0.0d0
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_normal.f90"
    !-----------------------------------------------!
    !
    call delete_sector(isector,H)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_d ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=0d0
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_normal_d


  !+------------------------------------------------------------------+
  !>SUPERC CASE
  !+------------------------------------------------------------------+
  subroutine build_H_superc_d(isector,Hmat)
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    integer                             :: isector
    integer,dimension(Nlevels)          :: ib
    integer,dimension(Ns)               :: ibup,ibdw
    integer                             :: dim,dimUp,dimDw
    integer                             :: i,iup,idw
    integer                             :: m,mup,mdw
    integer                             :: ishift,ishift_up,ishift_dw
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp,htmpup,htmpdw
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    integer                             :: first_state_up,last_state_up
    integer                             :: first_state_dw,last_state_dw
    !
    !
    call build_sector(isector,H)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
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
    diag_hybr=0.0d0
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_superc.f90"
    !-----------------------------------------------!
    !
    call delete_sector(isector,H)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_superc_d ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=0d0
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_superc_d


  !+------------------------------------------------------------------+
  !>NONSU2 CASE
  !+------------------------------------------------------------------+
  subroutine build_H_nonsu2_d(isector,Hmat)
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    integer                             :: isector
    integer,dimension(Nlevels)          :: ib
    integer,dimension(Ns)               :: ibup,ibdw
    integer                             :: dim,dimUp,dimDw
    integer                             :: i,iup,idw
    integer                             :: m,mup,mdw
    integer                             :: ishift,ishift_up,ishift_dw
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp,htmpup,htmpdw
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    integer                             :: first_state_up,last_state_up
    integer                             :: first_state_dw,last_state_dw    
    !
    !
    call build_sector(isector,H)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
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
    diag_hybr=0.0d0
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath))) 
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN/build_h_nonsu2.f90"
    !-----------------------------------------------!
    !
    call delete_sector(isector,H)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_nonsu2_d ERROR: size(Hmat) != dim**2"
       if(MpiStatus)then
          allocate(Hredux(dim,dim));Hredux=0d0
          call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
#ifdef _MPI
          call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MpiComm,MpiIerr)
#endif
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
    endif
    !
  end subroutine build_H_nonsu2_d
  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+



  !
  !>DOUBLE COMPLEX:
  !
  !+------------------------------------------------------------------+
  !>NORMAL CASE
  !+------------------------------------------------------------------+
  subroutine build_H_normal_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(Nlevels)            :: ib
    integer                               :: dim
    integer                               :: i,j,m,ms,impi,ishift
    integer                               :: iorb,jorb,ispin,jspin,ibath
    integer                               :: kp,k1,k2,k3,k4
    integer                               :: alfa,beta
    real(8)                               :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)               :: nup,ndw
    complex(8)                            :: htmp
    complex(8),dimension(Nspin,Norb)      :: eloc
    complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
    logical                               :: Jcondition
    integer                               :: first_state,last_state
    !
    call build_sector(isector,H)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
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
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
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
    call delete_sector(isector,H)
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
  !>SUPERC CASE
  !+------------------------------------------------------------------+
  subroutine build_H_superc_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(Nlevels)            :: ib
    integer                               :: dim
    integer                               :: i,j,m,ms,impi,ishift
    integer                               :: iorb,jorb,ispin,jspin,ibath
    integer                               :: kp,k1,k2,k3,k4
    integer                               :: alfa,beta
    real(8)                               :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)               :: nup,ndw
    complex(8)                            :: htmp
    complex(8),dimension(Nspin,Norb)      :: eloc
    complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
    logical                               :: Jcondition
    integer                               :: first_state,last_state
    !
    call build_sector(isector,H)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
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
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
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
    call delete_sector(isector,H)
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
  !>NONSU2 CASE
  !+------------------------------------------------------------------+
  subroutine build_H_nonsu2_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(Nlevels)            :: ib
    integer                               :: dim
    integer                               :: i,j,m,ms,impi,ishift
    integer                               :: iorb,jorb,ispin,jspin,ibath
    integer                               :: kp,k1,k2,k3,k4
    integer                               :: alfa,beta
    real(8)                               :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)               :: nup,ndw
    complex(8)                            :: htmp
    complex(8),dimension(Nspin,Norb)      :: eloc
    complex(8),dimension(Nspin,Norb,Nbath):: diag_hybr
    logical                               :: Jcondition
    integer                               :: first_state,last_state
    !
    call build_sector(isector,H)
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
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
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0.0d0)
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
    call delete_sector(isector,H)
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
  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+












  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! Lanczos algorithm using:
  ! - serial dble(H)*dble(V)
  ! - serial dble(H)*cmplx(V)
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI dble(H)*dble(V)
  ! - MPI dble(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_dd(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    integer                          :: i
    Hv=0d0
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_dd

  subroutine spMatVec_dc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    Hv=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_dc

  subroutine spMatVec_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
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
  subroutine spMatVec_mpi_dd(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    integer                          :: N
    real(8),dimension(:),allocatable :: vin
    integer                          :: i
    integer,allocatable,dimension(:) :: SendCounts,Displs
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dd ERRROR: MpiComm = MPI_UNDEFINED"
    MpiRank = get_rank_MPI(MpiComm)
    MpiSize = get_size_MPI(MpiComm)
    !
    !Get N_total by summing Nloc over all procs 
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)!call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,ierr)
    !
    !Get the chunks again (must equate Nloc for each proc)
    MpiQ  = get_Q_MPI(MpiComm,N)
    MpiR  = get_R_MPI(MpiComm,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1))
    allocate(displs(0:MpiSize-1))
    vin                   = 0d0
    SendCounts(0:)        = MpiQ
    SendCounts(MpiSize-1) = MpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*MpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Precision,vin,SendCounts,Displs,MPI_Double_Precision,MpiComm,MpiIerr)
    call Bcast_MPI(MpiComm,vin) !call MPI_Bcast(vin,N,MPI_Double_Precision,0,MpiComm,ierr)
    !
    !Perform the Mat*vec product. Each proc accumulates on the small vector Hv
    Hv=0d0
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_mpi_dd


  subroutine spMatVec_mpi_dc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer                             :: i
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dc ERRROR: MpiComm = MPI_UNDEFINED"
    MpiRank = get_rank_MPI(MpiComm)
    MpiSize = get_size_MPI(MpiComm)
    !
    !Get N_total by summing Nloc over all procs 
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Get the chunks again (must equate Nloc for each proc)
    MpiQ  = get_Q_MPI(MpiComm,N)
    MpiR  = get_R_MPI(MpiComm,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1))
    allocate(displs(0:MpiSize-1))
    vin                   = zero
    SendCounts(0:)        = MpiQ
    SendCounts(MpiSize-1) = MpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*MpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    call Bcast_MPI(MpiComm,vin)
    !
    !Perform the Mat*vec product. Each proc accumulates on the small vector Hv
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_mpi_dc

  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer                             :: i
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
    MpiRank = get_rank_MPI(MpiComm)
    MpiSize = get_size_MPI(MpiComm)
    !
    !Get N_total by summing Nloc over all procs 
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Get the chunks again (must equate Nloc for each proc)
    MpiQ  = get_Q_MPI(MpiComm,N)
    MpiR  = get_R_MPI(MpiComm,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(vin(N))
    allocate(SendCounts(0:MpiSize-1))
    allocate(displs(0:MpiSize-1))
    vin                   = zero
    SendCounts(0:)        = MpiQ
    SendCounts(MpiSize-1) = MpiQ+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*MpiQ
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    call Bcast_MPI(MpiComm,vin)
    !
    !Perform the Mat*vec product. Each proc accumulates on the small vector Hv
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








  !####################################################################
  !               RELATED AUXILIARY ROUTINES
  !####################################################################
  subroutine ed_hamiltonian_matvec_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
    MpiStatus=.true.
    ! MpiRank = get_Rank_MPI(MpiComm)
    ! MpiSize = get_Size_MPI(MpiComm)
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
    integer                              :: isector
    Hsector=isector
    call build_sector(isector,H)
  end subroutine setup_Hv_sector


  subroutine delete_Hv_sector(isector)
    Hsector=0
    call delete_sector(isector,H)
  end subroutine delete_Hv_sector

end MODULE ED_HAMILTONIAN
