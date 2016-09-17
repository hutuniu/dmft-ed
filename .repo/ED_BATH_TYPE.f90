MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable          :: e     !local energies [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable          :: d     !SC amplitues   [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable          :: v     !spin-keep hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable          :: u     !spin-flip hyb. [Nspin][Norb][Nbath]
     complex(8),dimension(:),allocatable           :: vr    !diagonal hyb.  [Nbath]
     complex(8),dimension(:,:,:,:,:),allocatable   :: h     !Replica hamilt [Nspin][Nspin][Norb][Norb][Nbath]
     logical(8),dimension(:,:,:,:,:),allocatable   :: mask  !impHloc mask   [Nspin][Nspin][Norb][Norb][Re,Im]
     complex(8),dimension(:,:,:,:,:),allocatable   :: LS    !Replica hamilt [Nspin][Nspin][Norb][Norb][LS,LSrot]
     logical                                :: status=.false.
  end type effective_bath

END MODULE ED_BATH_TYPE
