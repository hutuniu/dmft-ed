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
  if(ED_MPI_ID==0)then
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
  endif
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






subroutine print_sigma_nonsu2
  integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
  integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
  integer                          :: totNso,totNorb,totNspin,l,io,jo
  character(len=20)                :: suffix
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  if(ED_MPI_ID==0)then
     do l=1,totNso
        iorb=getIorb(l)
        jorb=getJorb(l)
        ispin=getIspin(l)
        jspin=getJspin(l)
        !
        suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
        call open_units(reg(suffix))
        !
        if(ed_verbose<4)then
           do i=1,Lmats
              write(unit(1),"(F20.12,2(F20.12))")wm(i),dimag(impSmats(ispin,jspin,iorb,jorb,i)),dreal(impSmats(ispin,jspin,iorb,jorb,i))
           enddo
           do i=1,Lreal
              write(unit(1),"(F20.12,2(F20.12))")wr(i),dimag(impSreal(ispin,jspin,iorb,jorb,i)),dreal(impSreal(ispin,jspin,iorb,jorb,i))
           enddo
        endif
        !
        call close_units()
     enddo
  endif
  !
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
contains
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<4)then
       open(unit(1),file="impSigma"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impSigma"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  subroutine close_units()
    if(ed_verbose<4)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_sigma_nonsu2






subroutine print_impG_nonsu2
  integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
  integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
  integer                          :: totNso,totNorb,totNspin,l,io,jo
  character(len=20)                :: suffix
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  if(ED_MPI_ID==0)then
     do l=1,totNso
        iorb=getIorb(l)
        jorb=getJorb(l)
        ispin=getIspin(l)
        jspin=getJspin(l)
        !
        suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
        call open_units(reg(suffix))
        !
        if(ed_verbose<2)then
           do i=1,Lmats
              write(unit(1),"(F20.12,2(F20.12))")wm(i),dimag(impGmats(ispin,jspin,iorb,jorb,i)),dreal(impGmats(ispin,jspin,iorb,jorb,i))
           enddo
           do i=1,Lreal
              write(unit(2),"(F20.12,2(F20.12))")wr(i),dimag(impGreal(ispin,jspin,iorb,jorb,i)),dreal(impGreal(ispin,jspin,iorb,jorb,i))
           enddo
        endif
        !
        call close_units()
     enddo
  endif
  !
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
contains
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<2)then
       open(unit(1),file="impG"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impG"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  subroutine close_units()
    if(ed_verbose<2)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_impG_nonsu2







subroutine print_impG0_nonsu2
  integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
  integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
  integer                          :: totNso,totNorb,totNspin,l,io,jo
  character(len=20)                :: suffix
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  if(ED_MPI_ID==0)then
     do l=1,totNso
        iorb=getIorb(l)
        jorb=getJorb(l)
        ispin=getIspin(l)
        jspin=getJspin(l)
        !
        suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
        call open_units(reg(suffix))
        !
        if(ed_verbose<1)then
           do i=1,Lmats
              write(unit(1),"(F20.12,2(F20.12))")wm(i),dimag(impG0mats(ispin,jspin,iorb,jorb,i)),dreal(impG0mats(ispin,jspin,iorb,jorb,i))
           enddo
           do i=1,Lreal
              write(unit(2),"(F20.12,2(F20.12))")wr(i),dimag(impG0real(ispin,jspin,iorb,jorb,i)),dreal(impG0real(ispin,jspin,iorb,jorb,i))
           enddo
        endif
        call close_units()
     enddo
  endif
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
contains
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<1)then
       open(unit(1),file="impG0"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impG0"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  subroutine close_units()
    if(ed_verbose<1)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_impG0_nonsu2

