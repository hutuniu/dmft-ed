program testPP
  implicit none
  !USE -fpp -D_MPI (or -D_noD or nothing if not) to compile
  complex(8),dimension(2,2) :: Cmat


  Cmat = 0d0
  Cmat(1,1)=dcmplx(0.25d0,0d0)
  Cmat(2,2)=dcmplx(0d0,0.25d0)

#ifdef _MPI
#define MPI_IN MPI_COMM,Inumber,Darray,Carray
  print*,"with MPI"
  call foo(1,2,[1d0,0d0],Cmat)
#else
#define MPI_IN Inumber,Darray,Carray
  print*,"no MPI"
  call foo(2,[1d0,0d0],Cmat)
#endif

contains

  subroutine foo(MPI_IN)
#ifdef _MPI    
    integer :: MPI_COMM
#endif
    integer :: Inumber
    real(8) :: Darray(2)
    complex(8) :: Carray(2,2)

#ifdef _MPI
    print*,MPI_COMM
#endif
    print*,Inumber
    print*,Darray
    print*,Carray
  end subroutine foo

end program testPP
