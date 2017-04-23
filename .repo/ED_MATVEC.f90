! NOTE about the difference between *spHtimes_xx* and *lanc_spHtimes_yy* 
! routines: although similar the lanc_* routines perform a slightly different 
! procedure for the matrix-vector product.  
! In particular the in-out vectors here have full size (N), they are decomposed 
! internally and mat-vec product is performed on the reduced basis, then full 
! vectors are recollected.  in the spHtimes* routines, used in the P/ARPACK 
! routines vector are passed in reduced size (Nloc), so that the mat-vec procedure 
! is different.
MODULE ED_MATVEC
  USE SF_CONSTANTS, only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !Sparse Matrix-vector product using stored sparse matrix 
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif
  !
  public  :: ed_matvec_set_MPI
  public  :: ed_matvec_del_MPI


#ifdef _MPI
  integer                      :: MpiComm=MPI_UNDEFINED  
  integer                      :: ierr
  integer                      :: MpiSize
  integer                      :: Q
  integer                      :: R
#endif
  type(sparse_element),pointer :: c



contains



  subroutine ed_matvec_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
#else
    integer,optional :: comm_
#endif
  end subroutine ed_matvec_set_MPI


  subroutine ed_matvec_del_MPI()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
#endif
  end subroutine ed_matvec_del_MPI



  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
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
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    N=0
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,ierr)
    MpiSize = get_Size_MPI(MpiComm)
    Q = get_Q_MPI(MpiComm,N)
    R = get_R_MPI(MpiComm,N)
    allocate(vin(N))
    vin=dcmplx(0d0,0d0)
    allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
    SendCounts(0:)     = Q
    SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
    forall(i=0:MpiSize-1)Displs(i)=i*Q
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,ierr)
    call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,ierr)
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

END MODULE ED_MATVEC




! subroutine spMatVec_dd(Nloc,v,Hv)
!   integer                          :: Nloc
!   real(8),dimension(Nloc)          :: v
!   real(8),dimension(Nloc)          :: Hv
!   integer                          :: i
!   Hv=0d0
!   do i=1,Nloc
!      c => spH0%row(i)%root%next       
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%val*v(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
! end subroutine spMatVec_dd

! subroutine spMatVec_dc(Nloc,v,Hv)
!   integer                             :: Nloc
!   complex(8),dimension(Nloc)          :: v
!   complex(8),dimension(Nloc)          :: Hv
!   integer                             :: i
!   Hv=zero
!   do i=1,Nloc
!      c => spH0%row(i)%root%next       
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%val*v(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
! end subroutine spMatVec_dc

! #ifdef _MPI
! subroutine spMatVec_mpi_dd(Nloc,v,Hv)
!   integer                          :: Nloc
!   real(8),dimension(Nloc)          :: v
!   real(8),dimension(Nloc)          :: Hv
!   integer                          :: N
!   real(8),dimension(:),allocatable :: vin
!   integer                          :: i
!   integer,allocatable,dimension(:) :: SendCounts,Displs
!   N=0
!   if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dd ERRROR: MpiComm = MPI_UNDEFINED"
!   call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,ierr)
!   MpiSize = get_size_MPI(MpiComm)
!   Q = get_Q_MPI(MpiComm,N)
!   R = get_R_MPI(MpiComm,N)
!   allocate(vin(N))
!   vin=0d0
!   allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!   SendCounts(0:)     = Q
!   SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
!   forall(i=0:MpiSize-1)Displs(i)=i*Q
!   call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Precision,vin,SendCounts,Displs,MPI_Double_Precision,MpiComm,ierr)
!   call MPI_Bcast(vin,N,MPI_Double_Precision,0,MpiComm,ierr)
!   Hv=0d0
!   do i=1,Nloc
!      c => spH0%row(i)%root%next       
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%val*vin(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
! end subroutine spMatVec_mpi_dd

! subroutine spMatVec_mpi_dc(Nloc,v,Hv)
!   integer                             :: Nloc
!   complex(8),dimension(Nloc)          :: v
!   complex(8),dimension(Nloc)          :: Hv
!   integer                             :: i
!   integer                             :: N
!   complex(8),dimension(:),allocatable :: vin
!   integer,allocatable,dimension(:)    :: SendCounts,Displs
!   N=0
!   if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
!   call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,ierr)
!   MpiSize = get_Size_MPI(MpiComm)
!   Q = get_Q_MPI(MpiComm,N)
!   R = get_R_MPI(MpiComm,N)
!   allocate(vin(N))
!   vin=dcmplx(0d0,0d0)
!   allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!   SendCounts(0:)     = Q
!   SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
!   forall(i=0:MpiSize-1)Displs(i)=i*Q
!   call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,ierr)
!   call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,ierr)
!   Hv=zero
!   do i=1,Nloc                 !==spH0%Nrow
!      c => spH0%row(i)%root%next       
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%val*vin(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
! end subroutine spMatVec_mpi_dc
! #endif



!   !**********************************************************************
!   !**********************************************************************
!   !*                                                                   *!
!   !*                   PLAIN LANCZOS PRODUCTS                          *!
!   !*                                                                   *!
!   !**********************************************************************
!   !**********************************************************************
!   !+------------------------------------------------------------------+
!   !PURPOSE  : Perform the matrix-vector product H*v used in the
!   ! Plain Lanczos algorithm for GF using MPI
!   !+------------------------------------------------------------------+
!   subroutine lanc_spMatVec_dd(N,v,Hv)
!     integer                          :: N
!     real(8),dimension(N)             :: v
!     real(8),dimension(N)             :: Hv
!     integer                          :: i,j
!     Hv=0.d0
!     do i=1,N
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           Hv(i) = Hv(i) + c%val*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!   end subroutine lanc_spMatVec_dd

!   subroutine lanc_spMatVec_dc(N,v,Hv)
!     integer                             :: N
!     complex(8),dimension(N)             :: v
!     complex(8),dimension(N)             :: Hv
!     integer                             :: i,j
!     Hv=zero
!     do i=1,N
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           Hv(i) = Hv(i) + c%val*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!   end subroutine lanc_spMatVec_dc

!   subroutine lanc_spMatVec_cc(N,v,Hv)
!     integer                             :: N
!     complex(8),dimension(N)             :: v
!     complex(8),dimension(N)             :: Hv
!     integer                             :: i,j
!     Hv=zero
!     do i=1,N
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           Hv(i) = Hv(i) + c%cval*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!   end subroutine lanc_spMatVec_cc


! #ifdef _MPI
!   subroutine lanc_spMatVec_MPI_dd(N,v,Hv)
!     integer                          :: N
!     real(8),dimension(N)             :: v
!     real(8),dimension(N)             :: Hv
!     integer                          :: i,j
!     integer                          :: Nloc
!     real(8),dimension(:),allocatable :: vout
!     integer,allocatable,dimension(:) :: SendCounts,Displs
!     if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dd ERRROR: MpiComm = MPI_UNDEFINED"
!     MpiSize = get_Size_MPI(MpiComm)
!     Q = get_Q_MPI(MpiComm,N)
!     R = get_R_MPI(MpiComm,N)
!     Nloc = Q+R
!     allocate(vout(Nloc))
!     vout=0d0
!     do i=1,Nloc
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           vout(i) = vout(i) + c%val*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!     allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!     SendCounts(0:)     = Q
!     SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
!     forall(i=0:MpiSize-1)Displs(i)=i*Q
!     Hv=0d0
!     call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Precision,Hv,SendCounts,Displs,MPI_Double_Precision,MpiComm,ierr)
!     call MPI_Bcast(Hv,N,MPI_Double_Precision,0,MpiComm,ierr)
!   end subroutine lanc_spMatVec_MPI_dd


!   subroutine lanc_spMatVec_MPI_dc(N,v,Hv)
!     integer                             :: N
!     complex(8),dimension(N)             :: v
!     complex(8),dimension(N)             :: Hv
!     integer                             :: i,j
!     integer                             :: Nloc
!     complex(8),dimension(:),allocatable :: vout
!     integer,allocatable,dimension(:)    :: SendCounts,Displs
!     if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dd ERRROR: MpiComm = MPI_UNDEFINED"
!     MpiSize = get_Size_MPI(MpiComm)
!     Q = get_Q_MPI(MpiComm,N)
!     R = get_R_MPI(MpiComm,N)
!     Nloc = Q+R
!     allocate(vout(Nloc))
!     vout=zero
!     do i=1,Nloc
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           vout(i) = vout(i) + c%val*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!     allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!     SendCounts(0:)     = Q
!     SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
!     forall(i=0:MpiSize-1)Displs(i)=i*Q
!     Hv=zero
!     call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,MpiComm,ierr)
!     call MPI_Bcast(Hv,N,MPI_Double_Complex,0,MpiComm,ierr)
!   end subroutine lanc_spMatVec_MPI_dc


!   subroutine lanc_spMatVec_MPI_cc(N,v,Hv)
!     integer                             :: N
!     complex(8),dimension(N)             :: v
!     complex(8),dimension(N)             :: Hv
!     integer                             :: i,j
!     integer                             :: Nloc
!     complex(8),dimension(:),allocatable :: vout
!     integer,allocatable,dimension(:)    :: SendCounts,Displs
!     if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_dd ERRROR: MpiComm = MPI_UNDEFINED"
!     MpiSize = get_Size_MPI(MpiComm)
!     Q = get_Q_MPI(MpiComm,N)
!     R = get_R_MPI(MpiComm,N)
!     Nloc = Q+R
!     allocate(vout(Nloc))
!     vout=zero
!     do i=1,Nloc
!        c => spH0%row(i)%root%next       
!        matmul: do while(associated(c))
!           vout(i) = vout(i) + c%cval*v(c%col)
!           c => c%next
!        end do matmul
!     end do
!     nullify(c)
!     allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!     SendCounts(0:)     = Q
!     SendCounts(MpiSize-1) = Q+mod(N,MpiSize)
!     forall(i=0:MpiSize-1)Displs(i)=i*Q
!     Hv=zero
!     call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,MpiComm,ierr)
!     call MPI_Bcast(Hv,N,MPI_Double_Complex,0,MpiComm,ierr)
!   end subroutine lanc_spMatVec_MPI_cc
! #endif



