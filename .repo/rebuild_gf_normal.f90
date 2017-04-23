!+------------------------------------------------------------------+
!PURPOSE  : Rebuild the impurity Green's functions
!+------------------------------------------------------------------+
subroutine rebuild_gf_normal
  integer                          :: i,ispin,isign,unit(1),iorb,jorb
  character(len=20)                :: suffix
  integer,dimension(:),allocatable :: getIorb,getJorb
  integer                          :: totNorb,l,j
  real(8)                          :: de,peso
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  select case(bath_type)
  case default                !Diagonal in both spin and orbital
     totNorb=Norb
     allocate(getIorb(totNorb),getJorb(totNorb))
     l=0
     do iorb=1,Norb
        L=l+1
        getIorb(l)=iorb
        getJorb(l)=iorb
     enddo
     totNorb=l
  case ("hybrid")             !Diagonal in spin only. Full Orbital structure
     totNorb=Norb*(Norb+1)/2
     allocate(getIorb(totNorb),getJorb(totNorb))
     l=0
     do iorb=1,Norb
        do jorb=iorb,Norb
           l=l+1
           getIorb(l)=iorb
           getJorb(l)=jorb
        enddo
     enddo
  end select
  !
  !Read the Poles&Weights => then it reconstructs the Gimp
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
     call open_units(reg(suffix))
     do isign=1,2
        do i=1,lanc_nGFiter
           read(unit(1),*)(GFpoles(ispin,ispin,iorb,jorb,isign,i),GFweights(ispin,ispin,iorb,jorb,isign,i),ispin=1,Nspin)
        enddo
     enddo
     call close_units
  enddo
  !
  impGmats=zero
  impGreal=zero
  do ispin=1,Nspin
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        do isign=1,2
           do j=1,lanc_nGFiter
              de    = GFpoles(ispin,ispin,iorb,jorb,isign,j)
              peso  = GFweights(ispin,ispin,iorb,jorb,isign,j)
              do i=1,Lmats
                 impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(xi*wm(i)-de)
              enddo
              do i=1,Lreal
                 impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(dcmplx(wr(i),eps)-de)
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
contains
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(1)
    open(unit(1),file="impGpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine rebuild_gf_normal
