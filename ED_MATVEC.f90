MODULE ED_MATVEC
  USE SF_CONSTANTS, only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
#ifdef _MPI
  USE SF_MPI
#endif
  implicit none
  private

  ! NOTE: FORMER CODE (NOT FULLFILLING SCIFOR RULES @ THE END OF THE MODULE)
  !
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
  public  :: spHtimesV_dd
  public  :: spHtimesV_cc
  !
  public  :: lanc_spHtimesV_dd
  public  :: lanc_spHtimesV_dc
  public  :: lanc_spHtimesV_cc

  public  :: sp_MatVec_Prod_dd !checked
  public  :: sp_MatVec_Prod_dc !checked
  public  :: sp_MatVec_Prod_cc !checked

  integer :: ierr
  integer :: rank
  integer :: size
  logical :: master
  integer :: Q
  integer :: R


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
  subroutine spHtimesV_dd(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
#ifdef _MPI
    integer                          :: N
    real(8),dimension(:),allocatable :: vin
    integer                          :: i
    integer,allocatable,dimension(:) :: SendCounts,Displs
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
    call sp_matrix_vector_product_dd(spH0,N,vin,Nloc,Hv)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,Nloc,v,Nloc,Hv)
#endif
  end subroutine spHtimesV_dd
  !
  !
  subroutine spHtimesV_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
#ifdef _MPI
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer                             :: i
    integer,allocatable,dimension(:)    :: SendCounts,Displs
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
    call sp_matrix_vector_product_cc(spH0,N,v,Nloc,Hv)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,Nloc,v,Nloc,Hv)
#endif
  end subroutine spHtimesV_cc










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
    call sp_matrix_vector_product_dd(spH0,N,v,Nloc,vout)
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
    call sp_matrix_vector_product_dd(spH0,N,v,N,Hv)
#endif
  end subroutine lanc_spHtimesV_dd
  !
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
    call sp_matrix_vector_product_dc(spH0,N,v,Nloc,vout)
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
    call sp_matrix_vector_product_dc(spH0,N,v,N,Hv)
#endif
  end subroutine lanc_spHtimesV_dc
  !
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
    call sp_matrix_vector_product_cc(spH0,N,v,Nloc,vout)
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
    call sp_matrix_vector_product_cc(spH0,N,v,N,Hv)
#endif
  end subroutine lanc_spHtimesV_cc








  !**********************************************************************
  !**********************************************************************
  !*                                                                   *!
  !*              SPARSE MATRIX - DENSE VECTOR PRODUCTS                *!
  !*                                                                   *!
  !**********************************************************************
  !**********************************************************************
  !+------------------------------------------------------------------+
  !PURPOSE: given a vector vin, perform the matrix-vector multiplication
  ! H_sparse * vin and put the result in vout.
  !+------------------------------------------------------------------+
  subroutine sp_MatVec_Prod_dd(sparse,Nin,vin,Nout,vout)
    type(sparse_matrix),intent(in)           :: sparse
    integer                                  :: Nin
    real(8),dimension(Nin),intent(in)        :: vin
    integer                                  :: Nout
    real(8),dimension(Nout),intent(inout)    :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=0.d0
    do i=1,Nout                 !Ndim
       c => sparse%row(i)%root%next       
       matmul: do               !can this be do while(associated(c))?
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_MatVec_Prod_dd
  !+------------------------------------------------------------------+
  subroutine sp_MatVec_Prod_dc(sparse,Nin,vin,Nout,vout)
    type(sparse_matrix),intent(in)           :: sparse
    integer                                  :: Nin
    complex(8),dimension(Nin),intent(in)     :: vin
    integer                                  :: Nout
    complex(8),dimension(Nout),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0d0,0d0,8)
    do i=1,Nout                 !Ndim
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_MatVec_Prod_dc
  !+------------------------------------------------------------------+
  subroutine sp_MatVec_Prod_cc(sparse,Nin,vin,Nout,vout)
    type(sparse_matrix),intent(in)           :: sparse
    integer                                  :: Nin
    complex(8),dimension(Nin),intent(in)     :: vin
    integer                                  :: Nout
    complex(8),dimension(Nout),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0d0,0d0,8)
    do i=1,Nout                 !Ndim
       c => sparse%row(i)%root%next
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_MatVec_Prod_cc



END MODULE ED_MATVEC




!+------------------------------------------------------------------+
!PURPOSE  : Perform the matrix-vector product H*v used in the
! P/ARPACK Lanczos algorithm using serial double real, complex (MPI)
!+------------------------------------------------------------------+
!   subroutine spHtimesV_dd(N,Nloc,v,Hv)
!     integer                    :: N,Nloc
!     real(8),dimension(Nloc)    :: v
!     real(8),dimension(Nloc)    :: Hv
! #ifdef _MPI
!     integer                    :: Q,R
!     real(8),dimension(N)       :: vin,vtmp
!     integer                    :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     vtmp=0.d0
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        vtmp(i)=v(i-ED_MPI_ID*Q)
!     enddo
!     call MPI_AllReduce(vtmp,vin,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
!     call sp_matrix_vector_product_dd(spH0,N,vin,Nloc,Hv)
! #else
!     call sp_matrix_vector_product_dd(spH0,Nloc,v,Nloc,Hv)
! #endif
!   end subroutine spHtimesV_dd
!
!   subroutine spHtimesV_dc(N,Nloc,v,Hv)
!     integer                    :: N,Nloc
!     complex(8),dimension(Nloc) :: v
!     complex(8),dimension(Nloc) :: Hv
! #ifdef _MPI
!     integer                    :: Q,R
!     complex(8),dimension(N)    :: vin,vtmp
!     integer                    :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     vtmp=zero
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        vtmp(i)=v(i-ED_MPI_ID*Q)
!     enddo
!     call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
!     call sp_matrix_vector_product_dc(spH0,N,vin,Nloc,Hv)
! #else
!     call sp_matrix_vector_product_dc(spH0,Nloc,v,Nloc,Hv)
! #endif
!   end subroutine spHtimesV_dc
!
!   subroutine spHtimesV_cc(N,Nloc,v,Hv)
!     integer                    :: N,Nloc
!     complex(8),dimension(Nloc) :: v
!     complex(8),dimension(Nloc) :: Hv
! #ifdef _MPI
!     integer                    :: Q,R
!     complex(8),dimension(N)    :: vin,vtmp
!     integer                    :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     vtmp=zero
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        vtmp(i)=v(i-ED_MPI_ID*Q)
!     enddo
!     call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
!     call sp_matrix_vector_product_cc(spH0,N,vin,Nloc,Hv)
! #else
!     call sp_matrix_vector_product_cc(spH0,Nloc,v,Nloc,Hv)
! #endif
!   end subroutine spHtimesV_cc





!+------------------------------------------------------------------+
!PURPOSE  : Perform the matrix-vector product H*v used in the
! Plain Lanczos algorithm for GF using MPI
!+------------------------------------------------------------------+
!   subroutine lanc_spHtimesV_dd(Nloc,N,v,Hv)
!     integer                 :: Nloc,N
!     real(8),dimension(N)    :: v
!     real(8),dimension(N)    :: Hv
! #ifdef _MPI
!     real(8),dimension(N)    :: Hvtmp
!     real(8),dimension(Nloc) :: vout
!     integer                 :: Q,R
!     integer                 :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     call sp_matrix_vector_product_dd(spH0,N,v,Nloc,vout)
!     Hvtmp=0.d0
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        Hvtmp(i)=vout(i-ED_MPI_ID*Q)
!     enddo
!     Hv=0d0
!     call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
! #else
!     call sp_matrix_vector_product_dd(spH0,N,v,N,Hv)
! #endif
!   end subroutine lanc_spHtimesV_dd
!
!   subroutine lanc_spHtimesV_dc(Nloc,N,v,Hv)
!     integer                    :: Nloc,N
!     complex(8),dimension(N)    :: v
!     complex(8),dimension(N)    :: Hv
! #ifdef _MPI
!     complex(8),dimension(N)    :: Hvtmp
!     complex(8),dimension(Nloc) :: vout
!     integer                    :: Q,R
!     integer                    :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     call sp_matrix_vector_product_dc(spH0,N,v,Nloc,vout)
!     Hvtmp=0.d0
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        Hvtmp(i)=vout(i-ED_MPI_ID*Q)
!     enddo
!     Hv=zero
!     call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
! #else
!     call sp_matrix_vector_product_dc(spH0,N,v,N,Hv)
! #endif
!   end subroutine lanc_spHtimesV_dc
!
!   subroutine lanc_spHtimesV_cc(Nloc,N,v,Hv)
!     integer                    :: Nloc,N
!     complex(8),dimension(N)    :: v
!     complex(8),dimension(N)    :: Hv
! #ifdef _MPI
!     complex(8),dimension(N)    :: Hvtmp
!     complex(8),dimension(Nloc) :: vout
!     integer                    :: Q,R
!     integer                    :: i
! #endif
!     !
! #ifdef _MPI
!     Q = N/ED_MPI_SIZE ; R = 0
!     if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
!     call sp_matrix_vector_product_cc(spH0,N,v,Nloc,vout)
!     Hvtmp=0.d0
!     do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
!        Hvtmp(i)=vout(i-ED_MPI_ID*Q)
!     enddo
!     Hv=zero
!     call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
! #else
!     call sp_matrix_vector_product_cc(spH0,N,v,N,Hv)
! #endif
!   end subroutine lanc_spHtimesV_cc
