MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable :: e ![Nspin][1/Norb][Nbath]
     real(8),dimension(:,:,:),allocatable :: v ![Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable :: w ![Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable :: d ![Nspin][Norb][Nbath]
     logical                              :: status=.false.
  end type effective_bath
END MODULE ED_BATH_TYPE
