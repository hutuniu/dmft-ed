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
  USE SF_MPI
#endif
  implicit none
  private

  !Sparse Matrix-vector product using stored sparse matrix 
  public  :: spHtimesV_dd
  public  :: spHtimesV_cc
  !
  public  :: lanc_spHtimesV_dd
  public  :: lanc_spHtimesV_dc
  public  :: lanc_spHtimesV_cc

  integer                      :: ierr
  integer                      :: size
  logical                      :: master
  integer                      :: Q
  integer                      :: R

  type(sparse_element),pointer :: c



contains



  !**********************************************************************
  !**********************************************************************
  !*                                                                   *!
  !*                   ARPACK LANCZOS PRODUCTS                         *!
  !*                                                                   *!
  !**********************************************************************
  !**********************************************************************
  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! P/ARPACK Lanczos algorithm using serial double real, complex (MPI)
  !+------------------------------------------------------------------+
#ifdef _MPI
  subroutine spHtimesV_dd(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    integer                          :: N
    real(8),dimension(:),allocatable :: vin
    integer                          :: i
    integer,allocatable,dimension(:) :: SendCounts,Displs
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,ED_MPI_COMM,ierr)
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    allocate(vin(N))
    vin=0d0
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Precision,vin,SendCounts,Displs,MPI_Double_Precision,ED_MPI_COMM,ierr)
    call MPI_Bcast(vin,N,MPI_Double_Precision,0,ED_MPI_COMM,ierr)
    Hv=0d0
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spHtimesV_dd
  !
  !
  subroutine spHtimesV_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,ED_MPI_COMM,ierr)
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    allocate(vin(N))
    vin=dcmplx(0d0,0d0)
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(vin,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spHtimesV_cc
#else
  subroutine spHtimesV_dd(Nloc,v,Hv)
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
  end subroutine spHtimesV_dd
  !
  !
  subroutine spHtimesV_cc(Nloc,v,Hv)
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
  end subroutine spHtimesV_cc
#endif










  !**********************************************************************
  !**********************************************************************
  !*                                                                   *!
  !*                   PLAIN LANCZOS PRODUCTS                          *!
  !*                                                                   *!
  !**********************************************************************
  !**********************************************************************
  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! Plain Lanczos algorithm for GF using MPI
  !+------------------------------------------------------------------+
#ifdef _MPI
  subroutine lanc_spHtimesV_dd(N,v,Hv)
    integer                          :: N
    real(8),dimension(N)             :: v
    real(8),dimension(N)             :: Hv
    integer                          :: i,j
    integer                          :: Nloc
    real(8),dimension(:),allocatable :: vout
    integer,allocatable,dimension(:) :: SendCounts,Displs
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    allocate(vout(Nloc))
    vout=0d0
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          vout(i) = vout(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=0d0
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
  end subroutine lanc_spHtimesV_dd
  !
  subroutine lanc_spHtimesV_dc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
    integer                             :: i,j
    integer                             :: Nloc
    complex(8),dimension(:),allocatable :: vout
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    allocate(vout(Nloc))
    vout=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          vout(i) = vout(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=zero
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
  end subroutine lanc_spHtimesV_dc
  !
  subroutine lanc_spHtimesV_cc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
    integer                             :: i,j
    integer                             :: Nloc
    complex(8),dimension(:),allocatable :: vout
    integer                             :: i,j
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    allocate(vout(Nloc))
    vout=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          vout(i) = vout(i) + c%cval*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=zero
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
  end subroutine lanc_spHtimesV_cc
  !
  !
#else
  !
  !
  subroutine lanc_spHtimesV_dd(N,v,Hv)
    integer                          :: N
    real(8),dimension(N)             :: v
    real(8),dimension(N)             :: Hv
    integer                          :: i,j
    Hv=0.d0
    do i=1,N
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine lanc_spHtimesV_dd
  !
  subroutine lanc_spHtimesV_dc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
    integer                             :: i,j
    Hv=zero
    do i=1,N
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%val*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine lanc_spHtimesV_dc
  !
  subroutine lanc_spHtimesV_cc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
    integer                             :: i,j
    Hv=zero
    do i=1,N
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine lanc_spHtimesV_cc
  !
#endif









END MODULE ED_MATVEC
