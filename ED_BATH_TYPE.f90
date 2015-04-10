MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable   :: e  !local energies [Nspin][1/Norb][Nbath]
     real(8),dimension(:,:,:,:),allocatable :: v  !hybridizations [Nhel] [Nspin] [Norb] [Nbath]
     real(8),dimension(:,:,:),allocatable   :: d  !SC amplitues   [Nspin][Norb]  [Nbath]
     real(8),dimension(:,:,:,:),allocatable :: w  !nonSU2 hybrid. [Nhel] [Nhel]  [Norb] [Nbath]
     real(8),dimension(:,:,:),allocatable   :: h  !nonSU2 energies[Nhel] [Norb]  [Nbath]
     logical                                :: superc_status=.false.
     logical                                :: nonsu2_status=.false.
     logical                                :: status=.false.
  end type effective_bath

END MODULE ED_BATH_TYPE
