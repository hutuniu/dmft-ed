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
  implicit none
  private

  !Get sparse sector Hamiltonian
  public :: build_H_normal_d
  public :: build_H_normal_c
  public :: build_H_all_d
  public :: build_H_all_c


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine build_H_normal_d(isector,Hmat)
    real(8),dimension(:,:),optional     :: Hmat
#ifdef _MPI
    real(8),dimension(:,:),allocatable  :: Hredux
#endif
    integer                             :: isector
    type(sector_map)                    :: H,Hup,Hdw
    integer,dimension(Nlevels)          :: ib
    integer,dimension(Ns)               :: ibup,ibdw
    integer                             :: mpiQ,mpiR
    integer                             :: mpiQup,mpiRup
    integer                             :: mpiQdw,mpiRdw
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
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !
    dim=getdim(isector)
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    ishift      = ED_MPI_ID*mpiQ
    first_state = ED_MPI_ID*mpiQ+1
    last_state  = (ED_MPI_ID+1)*mpiQ+mpiR
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
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ED_HAMILTONIAN/build_h_normal.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    !<DEBUG
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_d ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_normal_d


  subroutine build_H_all_d(isector,Hmat)
    real(8),dimension(:,:),optional     :: Hmat
#ifdef _MPI
    real(8),dimension(:,:),allocatable  :: Hredux
#endif
    integer                             :: isector
    type(sector_map)                    :: H,Hup,Hdw
    integer,dimension(Nlevels)          :: ib
    integer,dimension(Ns)               :: ibup,ibdw
    integer                             :: mpiQ,mpiR
    integer                             :: mpiQup,mpiRup
    integer                             :: mpiQdw,mpiRdw
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
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
    !
    dim=getdim(isector)
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    ishift      = ED_MPI_ID*mpiQ
    first_state = ED_MPI_ID*mpiQ+1
    last_state  = (ED_MPI_ID+1)*mpiQ+mpiR
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
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ED_HAMILTONIAN/build_h_all.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    !<DEBUG
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_all_d ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_all_d










  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine build_H_normal_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
#ifdef _MPI
    complex(8),dimension(:,:),allocatable :: Hredux
#endif
    integer                               :: isector
    type(sector_map)                      :: H,Hup,Hdw
    integer,dimension(Nlevels)            :: ib
    integer                               :: mpiQ,mpiR                
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
    !
    dim=getdim(isector)
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
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ED_HAMILTONIAN/build_h_normal.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_c ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_normal_c

  subroutine build_H_all_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
#ifdef _MPI
    complex(8),dimension(:,:),allocatable :: Hredux
#endif
    integer                               :: isector
    type(sector_map)                      :: H,Hup,Hdw
    integer,dimension(Nlevels)            :: ib
    integer                               :: mpiQ,mpiR                
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
    !
    dim=getdim(isector)
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
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "ED_HAMILTONIAN/build_h_all.f90"
    !-----------------------------------------------!
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_all_c ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_all_c

end MODULE ED_HAMILTONIAN
