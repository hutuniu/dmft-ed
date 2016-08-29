MODULE ED_MATVEC
  USE SF_CONSTANTS, only:zero
  USE MATRIX_SPARSE
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
#ifdef _MPI
  USE SF_MPI
#endif
  implicit none
  private


  ! IMPORTANT NOTE about the difference between *spHtimes_xx* and 
  ! *lanc_spHtimes_yy* routines:
  ! although similar the lanc_* routines perform a slightly 
  ! different procedure for the matrix-vector product. 
  ! In particular the in-out vectors here have full size (N), 
  ! they are decomposed internally and mat-vec product is 
  ! performed on the reduced basis, then full vectors are 
  ! recollected.
  ! in the spHtimes* routines, used in the P/ARPACK routines
  ! vector are passed in reduced size (Nloc), so that the 
  ! mat-vec procedure is different.

  !Sparse Matrix-vector product using stored sparse matrix 
  public                       :: spHtimesV_dd
  public                       :: spHtimesV_cc
  !
  public                       :: lanc_spHtimesV_dd
  public                       :: lanc_spHtimesV_dc
  public                       :: lanc_spHtimesV_cc


  integer                          :: ierr
  integer                          :: rank
  integer                          :: size
  logical                          :: master
  integer                          :: Q
  integer                          :: R


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! P/ARPACK Lanczos algorithm using serial double real, complex (MPI)
  !+------------------------------------------------------------------+
  subroutine spHtimesV_dd(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    integer                          :: N
    real(8),dimension(:),allocatable :: vin
    integer                          :: i
    integer,allocatable,dimension(:) :: SendCounts,Displs
#ifdef _MPI
    rank = MPI_Get_rank(ED_MPI_COMM)
    size = MPI_Get_size(ED_MPI_COMM)
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,ED_MPI_COMM,ierr)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    !
    allocate(vin(N));vin=0d0
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Precision,vin,SendCounts,Displs,MPI_Double_Precision,ED_MPI_COMM,ierr)
    call MPI_Bcast(vin,N,MPI_Double_Precision,0,ED_MPI_COMM,ierr)
    Hv=0d0
    call sp_matrix_vector_product_mpi_dd(spH0,N,vin,Nloc,Hv)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_dd


  subroutine spHtimesV_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer                             :: i
    integer,allocatable,dimension(:)    :: SendCounts,Displs
#ifdef _MPI
    rank = MPI_Get_rank(ED_MPI_COMM)
    size = MPI_Get_size(ED_MPI_COMM)
    N=0
    call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,ED_MPI_COMM,ierr)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    !
    allocate(vin(N))
    vin=dcmplx(0d0,0d0)
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(vin,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
    Hv=zero
    call sp_matrix_vector_product_mpi_cc(spH0,N,vin,Nloc,Hv)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_cc







  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! Plain Lanczos algorithm for GF using MPI
  !+------------------------------------------------------------------+
  subroutine lanc_spHtimesV_dd(N,v,Hv)
    integer                          :: N
    real(8),dimension(N)             :: v
    real(8),dimension(N)             :: Hv
#ifdef _MPI
    integer                          :: Nloc
    real(8),dimension(:),allocatable :: vout
    integer                          :: i,j
    integer,allocatable,dimension(:) :: SendCounts,Displs
    !
    rank = MPI_Get_rank(ED_MPI_COMM)
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    !
    allocate(vout(Nloc))
    vout=0d0
    call sp_matrix_vector_product_mpi_dd(spH0,N,v,Nloc,vout)
    !
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=0d0
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dd



  subroutine lanc_spHtimesV_dc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
#ifdef _MPI
    integer                             :: Nloc
    complex(8),dimension(:),allocatable :: vout
    integer                             :: i,j
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    !
    rank = MPI_Get_rank(ED_MPI_COMM)
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    !
    allocate(vout(Nloc))
    vout=zero
    call sp_matrix_vector_product_mpi_dc(spH0,N,v,Nloc,vout)
    !
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=zero
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
#else
    Hv=zero
    call sp_matrix_vector_product_dc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dc



  subroutine lanc_spHtimesV_cc(N,v,Hv)
    integer                             :: N
    complex(8),dimension(N)             :: v
    complex(8),dimension(N)             :: Hv
#ifdef _MPI
    integer                             :: Nloc
    complex(8),dimension(:),allocatable :: vout
    integer                             :: i,j
    integer,allocatable,dimension(:)    :: SendCounts,Displs
    !
    rank = MPI_Get_rank(ED_MPI_COMM)
    size = MPI_Get_size(ED_MPI_COMM)
    Q = MPI_Get_Q(ED_MPI_COMM,N)
    R = MPI_Get_R(ED_MPI_COMM,N)
    Nloc = Q+R
    !
    allocate(vout(Nloc))
    vout=zero
    call sp_matrix_vector_product_mpi_cc(spH0,N,v,Nloc,vout)
    !
    allocate(SendCounts(0:size-1),displs(0:size-1))
    SendCounts(0:)     = Q
    SendCounts(size-1) = Q+mod(N,size)
    forall(i=0:size-1)Displs(i)=i*Q
    Hv=zero
    call MPI_Allgatherv(vout(1:Nloc),Nloc,MPI_Double_Complex,Hv,SendCounts,Displs,MPI_Double_Complex,ED_MPI_COMM,ierr)
    call MPI_Bcast(Hv,N,MPI_Double_Complex,0,ED_MPI_COMM,ierr)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_cc









  !+------------------------------------------------------------------+
  !PURPOSE: given a vector vin, perform the matrix-vector multiplication
  ! H_sparse * vin and put the result in vout.
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dd(sparse,Ndim,vin,vout)
    integer                               :: Ndim
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(Ndim),intent(in)    :: vin
    real(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer          :: c
    integer                               :: i
    vout=0.d0
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_dd
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dc(sparse,Ndim,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0.d0,0.d0,8)
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_dc
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_csr(Nrow,sparse,vin,vout)
    integer                               :: Nrow
    type(sparse_matrix_csr),intent(in)    :: sparse
    real(8),dimension(Nrow),intent(in)    :: vin
    real(8),dimension(Nrow),intent(inout) :: vout
    integer                               :: i,pos
    vout=0.d0
    do i=1,Nrow
       do pos=sparse%rowIndex(i),sparse%rowIndex(i+1)-1
          vout(i) = vout(i) + sparse%values(pos)*vin(sparse%columns(pos))
       end do
    end do
  end subroutine sp_matrix_vector_product_csr

  subroutine sp_matrix_vector_product_cc(sparse,Ndim,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0.d0,0.d0,8)
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_cc



#ifdef _MPI
  subroutine sp_matrix_vector_product_mpi_dd(sparse,Ndim,vin,Nloc,vout)
    integer                               :: Ndim,Nloc
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(Ndim),intent(in)    :: vin
    real(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer          :: c
    integer                               :: i
    integer                               :: Nini,Nfin
    vout=0.d0
    do i=1,Nloc
       c => sparse%row(i)%root%next
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_dd
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_mpi_dc(sparse,Ndim,vin,Nloc,vout)
    integer                                  :: Ndim,Nloc
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    integer                                  :: Nini,Nfin
    vout=zero
    do i=1,Nloc
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_dc
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_mpi_cc(sparse,Ndim,vin,Nloc,vout)
    integer                                  :: Ndim,Nloc
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    integer                                  :: Nini,Nfin
    vout=zero
    do i=1,Nloc
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_cc
#endif

END MODULE ED_MATVEC
