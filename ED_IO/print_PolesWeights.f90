subroutine print_poles_weights_normal
  integer                                           :: i,ispin,isign,unit(1),iorb,jorb
  character(len=20)                                 :: suffix
  integer,dimension(:),allocatable                  :: getIorb,getJorb
  integer                                           :: totNorb,l
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
  case ('hybrid')             !Diagonal in spin only. Full Orbital structure
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
  if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
  !!
  !Print the impurity functions:
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
     call open_units(reg(suffix))
     do isign=1,2
        do i=1,lanc_nGFiter
           write(unit(1),"(6(F20.12,1x))")(GFpoles(ispin,ispin,iorb,jorb,isign,i),GFweights(ispin,ispin,iorb,jorb,isign,i),ispin=1,Nspin)
        enddo
     enddo
     call close_units()
  enddo
  !
contains
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(1)
    open(unit(1),file="impGpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  !
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine print_poles_weights_normal


subroutine print_poles_weights_superc
  integer                                               :: i,ispin,unit(1),iorb,jorb,isign
  character(len=20)                                     :: suffix
  integer,dimension(:),allocatable                      :: getIorb,getJorb
  integer                                               :: totNorb,l
  !
  select case(bath_type)
  case default
     totNorb=Norb
     allocate(getIorb(Norb),getJorb(Norb))
     l=0
     do iorb=1,Norb
        l=l+1
        getIorb(l)=iorb
        getJorb(l)=iorb
     enddo
  case ("hybrid")
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
  if(l/=totNorb)stop "print_gf_superc error counting the orbitals"
  !!
  !!PRINT OUT GF:
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
     call open_units(reg(suffix))
     do isign=1,2
        do i=1,lanc_nGFiter
           write(unit(1),"(6(F20.12,1x))")(GFpoles(ispin,ispin,iorb,jorb,isign,i),GFweights(ispin,ispin,iorb,jorb,isign,i),ispin=1,Nspin)
        enddo
     enddo
     call close_units
  enddo
  !
contains
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    open(unit(1),file="Gpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine print_poles_weights_superc




subroutine print_poles_weights_nonsu2
  integer                          :: i,isign,unit(1),iorb,jorb,ispin,jspin
  integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
  integer                          :: totNso,totNorb,totNspin,l,io,jo
  character(len=20)                :: suffix
  !
  select case(bath_type)
  case default
     totNorb =Norb
     totNspin=Nspin*(Nspin+1)/2
     totNso  =totNorb*totNspin
     allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
     l=0
     do iorb=1,Norb
        do ispin=1,Nspin
           do jspin=ispin,Nspin
              l=l+1
              getIorb(l)=iorb
              getIspin(l)=ispin
              getJorb(l)=iorb
              getJspin(l)=jspin
           enddo
        enddo
     enddo
  case ("hybrid")
     totNso  = (Norb*Nspin)**2
     allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
     l=0
     do iorb=1,Norb
        do jorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 l=l+1
                 getIorb(l)=iorb
                 getIspin(l)=ispin
                 getJorb(l)=jorb
                 getJspin(l)=jspin
              enddo
           enddo
        enddo
     enddo
  case ("replica")
     l=0
     do iorb=1,Norb
        do jorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 if(ed_verbose>=0.and.io<jo)cycle
                 if(dmft_bath%mask(ispin,jspin,iorb,jorb,1).or.dmft_bath%mask(ispin,jspin,iorb,jorb,2)) l=l+1
              enddo
           enddo
        enddo
     enddo
     totNso = l
     allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
     l=0
     do iorb=1,Norb
        do jorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 if(ed_verbose>=0.and.io<jo)cycle
                 if((.not.dmft_bath%mask(ispin,jspin,iorb,jorb,1)).and.(.not.dmft_bath%mask(ispin,jspin,iorb,jorb,2))) cycle
                 l=l+1
                 getIorb(l)=iorb
                 getIspin(l)=ispin
                 getJorb(l)=jorb
                 getJspin(l)=jspin
              enddo
           enddo
        enddo
     enddo
  end select
  if(l/=totNso)stop "print_gf_nonsu2 error counting the spin-orbitals"
  !!
  !!PRINT OUT GF:
  do l=1,totNso
     iorb=getIorb(l)
     jorb=getJorb(l)
     ispin=getIspin(l)
     jspin=getJspin(l)
     !
     suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
     call open_units(reg(suffix))
     !
     do isign=1,2
        do i=1,lanc_nGFiter
           write(unit(1),"(2(F20.12,1x))")GFpoles(ispin,jspin,iorb,iorb,isign,i),GFweights(ispin,jspin,iorb,iorb,isign,i)
        enddo
     enddo
     !
     call close_units()
  enddo
  !
contains
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    open(unit(1),file="Gpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine print_poles_weights_nonsu2
