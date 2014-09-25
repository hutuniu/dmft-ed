!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN
  USE CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private

  !Get sparse sector Hamiltonian
  public                       :: ed_buildH_d,ed_buildH_c


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_buildH_d(isector,Hmat)
    real(8),dimension(:,:),optional    :: Hmat
    real(8),dimension(:,:),allocatable :: Hredux
    integer                            :: isector
    integer,dimension(:),allocatable   :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Ntot)            :: ib
    integer                            :: mpiQ,mpiR                
    integer                            :: dim,iup,idw
    integer                            :: i,j,m,ms,iorb,jorb,ispin,impi,ishift
    integer                            :: kp,k1,k2,k3,k4
    real(8)                            :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)            :: nup,ndw
    real(8)                            :: htmp
    real(8),dimension(Nspin,Norb)      :: eloc
    logical                            :: Jcondition
    integer                            :: first_state,last_state
    !
    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    ishift     = ED_MPI_ID*mpiQ
    first_state= ED_MPI_ID*mpiQ+1
    last_state = (ED_MPI_ID+1)*mpiQ+mpiR
#else
    mpiQ=0
    mpiR=0
    call sp_init_matrix(spH0,dim)
    ishift     = 0
    first_state= 1
    last_state = dim
#endif
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=dreal(impHloc(ispin,ispin,iorb,iorb))
       enddo
    enddo
    !
    !-----------------------------------------------!
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ed_build_hamiltonian.f90"
    !-----------------------------------------------!
    !
    deallocate(Hmap)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "ED_HAMILTONIAN/ed_buildH_d: wrong dimensions in Hmat"
#ifdef _MPI
       allocate(Hredux(dim,dim))
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine ed_buildH_d



  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine ed_buildH_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(:),allocatable      :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Ntot)               :: ib
    integer                               :: mpiQ,mpiR                
    integer                               :: dim,iup,idw
    integer                               :: i,j,m,ms,iorb,jorb,ispin,impi,ishift
    integer                               :: kp,k1,k2,k3,k4
    real(8)                               :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)               :: nup,ndw
    complex(8)                            :: htmp
    complex(8),dimension(Nspin,Norb)      :: eloc
    logical                               :: Jcondition
    integer                               :: first_state,last_state
    !
    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
    !
    first_state= 1
    last_state = dim
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    ishift     = ED_MPI_ID*mpiQ
    first_state= ED_MPI_ID*mpiQ+1
    last_state = (ED_MPI_ID+1)*mpiQ+mpiR
#else
    mpiQ=0
    mpiR=0
    call sp_init_matrix(spH0,dim)
    ishift     = 0
    first_state= 1
    last_state = dim
#endif
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    !-----------------------------------------------!
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ed_build_hamiltonian.f90"
    !-----------------------------------------------!
    !
    deallocate(Hmap)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "ED_HAMILTONIAN/ed_buildH_d: wrong dimensions in Hmat"
#ifdef _MPI
       allocate(Hredux(dim,dim))
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine ed_buildH_c


end MODULE ED_HAMILTONIAN
