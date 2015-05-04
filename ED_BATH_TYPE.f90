MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable   :: e  !local energies [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable   :: d  !SC amplitues   [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable   :: v  !spin-keep hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable   :: vs !spin-flip hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:,:,:,:),allocatable :: w  !hybridizations [Nspin][Nspin][Norb][Nbath]
     logical                                :: status=.false.
  end type effective_bath

END MODULE ED_BATH_TYPE
