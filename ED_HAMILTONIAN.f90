MODULE ED_HAMILTONIAN
  USE SF_CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
  implicit none
  private

  !Get sparse sector Hamiltonian
  public :: build_H_normal_irrhyb_d
  public :: build_H_superc_irrhyb_d
  public :: build_H_nonsu2_irrhyb_d
  !
  public :: build_H_normal_replica_d
  public :: build_H_superc_replica_d
  public :: build_H_nonsu2_replica_d
  !
  public :: build_H_normal_irrhyb_c
  public :: build_H_superc_irrhyb_c
  public :: build_H_nonsu2_irrhyb_c
  !
  public :: build_H_normal_replica_c
  public :: build_H_superc_replica_c
  public :: build_H_nonsu2_replica_c

contains



  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                  IRREDUCIBLE & HYBRIDIZED BATH                        ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             NORMAL                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_normal_irrhyb_d '
  subroutine build_H_normal_irrhyb_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: dim
    integer                             :: i
    integer                             :: j
    integer                             :: m
    integer                             :: mpiQ,mpiR
    integer                             :: ishift
    integer                             :: first_state,last_state
    integer                             :: impi
    real(8)                             :: htmp        
    !
    type(sector_map)                    :: Hup,Hdw
    integer,dimension(Ns)               :: nup,ndw
    integer                             :: dimUp,dimDw
    integer                             :: iup,idw
    integer                             :: jup,jdw
    integer                             :: mup,mdw
    integer                             :: mpiQ_up,mpiR_up
    integer                             :: mpiQ_dw,mpiR_dw
    integer                             :: ishift_up,ishift_dw
    integer                             :: first_state_up,last_state_up
    integer                             :: first_state_dw,last_state_dw
    integer                             :: impi_up,impi_dw
    real(8)                             :: htmp_up,htmp_dw    
    !
    integer                             :: ms
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    !
    logical                             :: Jcondition
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_normal_irrhyb.f90"
    !
    deallocate(Hup%map,Hdw%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_normal_irrhyb_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_normal_irrhyb_c '
  subroutine build_H_normal_irrhyb_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: dim
    integer                                :: i
    integer                                :: j
    integer                                :: m
    integer                                :: mpiQ,mpiR
    integer                                :: ishift
    integer                                :: first_state,last_state
    integer                                :: impi
    real(8)                                :: htmp        
    !
    type(sector_map)                       :: Hup,Hdw
    integer,dimension(Ns)                  :: nup,ndw
    integer                                :: dimUp,dimDw
    integer                                :: iup,idw
    integer                                :: jup,jdw
    integer                                :: mup,mdw
    integer                                :: mpiQ_up,mpiR_up
    integer                                :: mpiQ_dw,mpiR_dw
    integer                                :: ishift_up,ishift_dw
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    integer                                :: impi_up,impi_dw
    real(8)                                :: htmp_up,htmp_dw    
    !
    integer                                :: ms
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    !
    logical                                :: Jcondition
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_normal_irrhyb.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_normal_irrhyb_c
#undef FNAME














  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                  IRREDUCIBLE & HYBRIDIZED BATH                        ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             SUPERC                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_superc_irrhyb_d '
  subroutine build_H_superc_irrhyb_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: mpiQ,mpiR
    integer                             :: dim
    integer                             :: i
    integer                             :: m
    integer                             :: ishift
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_superc_irrhyb.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_superc_irrhyb_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_superc_irrhyb_c '
  subroutine build_H_superc_irrhyb_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: mpiQ,mpiR
    integer                                :: dim
    integer                                :: i
    integer                                :: m
    integer                                :: ishift
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    real(8)                                :: htmp
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_superc_irrhyb.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_superc_irrhyb_c
#undef FNAME











  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                  IRREDUCIBLE & HYBRIDIZED BATH                        ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             NONSU2                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_nonsu2_irrhyb_d '
  subroutine build_H_nonsu2_irrhyb_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: mpiQ,mpiR
    integer                             :: dim
    integer                             :: i
    integer                             :: m
    integer                             :: ishift
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_nonsu2_irrhyb.f90"
    !
    deallocate(H%map)
    !
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_nonsu2_irrhyb_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_nonsu2_irrhyb_c '
  subroutine build_H_nonsu2_irrhyb_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: mpiQ,mpiR
    integer                                :: dim
    integer                                :: i
    integer                                :: m
    integer                                :: ishift
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    real(8)                                :: htmp
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=cmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/irrhyb/build_h_nonsu2_irrhyb.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_nonsu2_irrhyb_c
#undef FNAME




















  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                          REPLICA  BATH                                ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             NORMAL                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_normal_replica_d '
  subroutine build_H_normal_replica_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: mpiQ,mpiR
    integer                             :: dim
    integer                             :: i
    integer                             :: m
    integer                             :: ishift
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath)))!?? why ABS??
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_normal_replica.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_normal_replica_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_normal_replica_c '
  subroutine build_H_normal_replica_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: mpiQ,mpiR
    integer                                :: dim
    integer                                :: i
    integer                                :: m
    integer                                :: ishift
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    real(8)                                :: htmp
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_normal_replica.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_normal_replica_c
#undef FNAME













  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                          REPLICA  BATH                                ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             SUPERC                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_superc_replica_d '
  subroutine build_H_superc_replica_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: mpiQ,mpiR
    integer                             :: dim
    integer                             :: i
    integer                             :: m
    integer                             :: ishift
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath)))!?? why ABS??
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_superc_replica.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_superc_replica_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_superc_replica_c '
  subroutine build_H_superc_replica_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: mpiQ,mpiR
    integer                                :: dim
    integer                                :: i
    integer                                :: m
    integer                                :: ishift
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    real(8)                                :: htmp
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_superc_replica.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_superc_replica_c
#undef FNAME




















  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                          REPLICA  BATH                                ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !                             NONSU2                                    ! 
  !                                                                       !
  !***********************************************************************!
  !***********************************************************************!
  !                                                                       !
  !               DOUBLE PRECISION (SYMMETRIC) H                          ! 
  !                                                                       !
  !***********************************************************************!
#define FNAME 'build_H_nonsu2_replica_d '
  subroutine build_H_nonsu2_replica_d(isector,Hmat)
    integer                             :: isector
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Hredux
    real(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                    :: H
    integer,dimension(Nlevels)          :: ib
    integer                             :: mpiQ,mpiR
    integer                             :: dim
    integer                             :: i
    integer                             :: m
    integer                             :: ishift
    integer                             :: j,ms,impi
    integer                             :: iorb,jorb,ispin,jspin,ibath
    integer                             :: kp,k1,k2,k3,k4
    integer                             :: alfa,beta
    real(8)                             :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)             :: nup,ndw
    real(8)                             :: htmp
    logical                             :: Jcondition
    integer                             :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=0.0d0
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=abs(dreal(dmft_bath%vr(ibath)))!?? why ABS??
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_nonsu2_replica.f90"
    !
    deallocate(H%map)
    !
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=0d0
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine build_H_nonsu2_replica_d
#undef FNAME
  !
  !***********************************************************************!
  !                                                                       !
  !                 DOUBLE COMPLEX (HERMITIAN) H                          ! 
  !                                                                       !
  !***********************************************************************!
  !
#define FNAME 'build_H_nonsu2_replica_c '
  subroutine build_H_nonsu2_replica_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    complex(8),dimension(:,:),allocatable  :: Hredux
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    !
    type(sector_map)                       :: H
    integer,dimension(Nlevels)             :: ib
    integer                                :: mpiQ,mpiR
    integer                                :: dim
    integer                                :: i
    integer                                :: m
    integer                                :: ishift
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    real(8)                                :: htmp
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
          enddo
       enddo
    enddo
    !
    include "ED_HAMILTONIAN/replica/build_h_nonsu2_replica.f90"
    !
    deallocate(H%map)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop FNAME//" ERROR: size(Hmat) != dim**2"
#ifdef _MPI
       allocate(Hredux(dim,dim));Hredux=zero
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
  end subroutine build_H_nonsu2_replica_c
#undef FNAME


























































end MODULE ED_HAMILTONIAN
