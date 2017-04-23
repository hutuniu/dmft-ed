MODULE ED_MPI2B
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  !MPI Parallel environment variables
  !PUBLIC
  !=========================================================
  integer :: ED_MPI_COMM=MPI_UNDEFINED
  logical :: ED_MPI_STATUS=.false.
  integer :: ED_MPI_ID=0
  integer :: ED_MPI_SIZE=1
  integer :: ED_MPI_ERR
  logical :: ED_MPI_MASTER=.true.



contains

  subroutine ed_set_MPI(comm)
    integer :: comm
#ifdef _MPI
    ED_MPI_COMM = comm
    ED_MPI_STATUS = .true.
    ED_MPI_ID = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE=get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER=get_Master_MPI(ED_MPI_COMM)
#endif
  end subroutine ed_set_MPI


END MODULE ED_MPI2B
