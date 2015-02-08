MODULE ED_MATVEC
  USE SF_CONSTANTS, only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
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
  public                       :: spHtimesV_dc
  public                       :: spHtimesV_cc
  !
  public                       :: lanc_spHtimesV_dd
  public                       :: lanc_spHtimesV_dc
  public                       :: lanc_spHtimesV_cc



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! P/ARPACK Lanczos algorithm using serial double real, complex (MPI)
  !+------------------------------------------------------------------+
  subroutine spHtimesV_dd(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    real(8),dimension(Nloc)    :: v
    real(8),dimension(Nloc)    :: Hv
    integer                    :: Q,R
    real(8),dimension(N)       :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    vtmp=0.d0
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       vtmp(i)=v(i-ED_MPI_ID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
    call sp_matrix_vector_product_mpi_dd(spH0,N,vin,Nloc,Hv)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_dd

  subroutine spHtimesV_dc(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    integer                    :: Q,R
    complex(8),dimension(N)    :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    vtmp=zero
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       vtmp(i)=v(i-ED_MPI_ID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
    call sp_matrix_vector_product_mpi_dc(spH0,N,vin,Nloc,Hv)
#else
    Hv=zero
    call sp_matrix_vector_product_dc(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_dc

  subroutine spHtimesV_cc(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    integer                    :: Q,R
    complex(8),dimension(N)    :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    vtmp=zero
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       vtmp(i)=v(i-ED_MPI_ID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
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
  subroutine lanc_spHtimesV_dd(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    real(8),dimension(N)    :: v,Hv,Hvtmp
    real(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    call sp_matrix_vector_product_mpi_dd(spH0,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       Hvtmp(i)=vout(i-ED_MPI_ID*Q)
    enddo
    Hv=0.d0
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dd

  subroutine lanc_spHtimesV_dc(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    complex(8),dimension(N)    :: v,Hv,Hvtmp
    complex(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    call sp_matrix_vector_product_mpi_dc(spH0,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       Hvtmp(i)=vout(i-ED_MPI_ID*Q)
    enddo
    Hv=zero
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
    Hv=zero
    call sp_matrix_vector_product_dc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dc

  subroutine lanc_spHtimesV_cc(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    complex(8),dimension(N)    :: v,Hv,Hvtmp
    complex(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/ED_MPI_SIZE ; R = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))R=mod(N,ED_MPI_SIZE)
    call sp_matrix_vector_product_mpi_cc(spH0,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=ED_MPI_ID*Q+1,(ED_MPI_ID+1)*Q+R
       Hvtmp(i)=vout(i-ED_MPI_ID*Q)
    enddo
    Hv=zero
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_cc


END MODULE ED_MATVEC
