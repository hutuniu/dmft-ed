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
  if(ED_MPI_ID==0)then
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
end subroutine print_poles_weights_superc




subroutine print_sigma_superc
  integer                                               :: i,ispin,unit(4),iorb,jorb,isign
  character(len=20)                                     :: suffix
  integer,dimension(:),allocatable                      :: getIorb,getJorb
  integer                                               :: totNorb,l
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<4)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impSmats(ispin,ispin,iorb,jorb,i)),dreal(impSmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lmats
              write(unit(2),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impSAmats(ispin,ispin,iorb,jorb,i)),dreal(impSAmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(3),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impSreal(ispin,ispin,iorb,jorb,i)),dreal(impSreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(4),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impSAreal(ispin,ispin,iorb,jorb,i)),dreal(impSAreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
        endif
        call close_units
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
    if(ed_verbose<4)then
       open(unit(1),file="impSigma"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impSelf"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(3),file="impSigma"//string//"_realw"//reg(ed_file_suffix)//".ed")
       open(unit(4),file="impSelf"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  subroutine close_units()
    if(ed_verbose<4)then
       close(unit(1))
       close(unit(2))
       close(unit(3))
       close(unit(4))
    endif
  end subroutine close_units
end subroutine print_sigma_superc



subroutine print_impG_superc
  integer                                               :: i,ispin,unit(4),iorb,jorb,isign
  character(len=20)                                     :: suffix
  integer,dimension(:),allocatable                      :: getIorb,getJorb
  integer                                               :: totNorb,l
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<2)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impGmats(ispin,ispin,iorb,jorb,i)),dreal(impGmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lmats
              write(unit(2),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impFmats(ispin,ispin,iorb,jorb,i)),dreal(impFmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(3),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impGreal(ispin,ispin,iorb,jorb,i)),dreal(impGreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(4),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impFreal(ispin,ispin,iorb,jorb,i)),dreal(impFreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
        endif
        call close_units
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
       open(unit(2),file="impF"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(3),file="impG"//string//"_realw"//reg(ed_file_suffix)//".ed")
       open(unit(4),file="impF"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  subroutine close_units()
    if(ed_verbose<2)then
       close(unit(1))
       close(unit(2))
       close(unit(3))
       close(unit(4))
    endif
  end subroutine close_units
end subroutine print_impG_superc



subroutine print_impG0_superc
  integer                                               :: i,ispin,unit(4),iorb,jorb,isign
  character(len=20)                                     :: suffix
  integer,dimension(:),allocatable                      :: getIorb,getJorb
  integer                                               :: totNorb,l
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
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
  !!
  !!
  !!PRINT OUT GF:
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<1)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impG0mats(ispin,ispin,iorb,jorb,i)),dreal(impG0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lmats
              write(unit(2),"(F20.12,6(F20.12))")wm(i),&
                   (dimag(impF0mats(ispin,ispin,iorb,jorb,i)),dreal(impF0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(3),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impG0real(ispin,ispin,iorb,jorb,i)),dreal(impG0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(4),"(F20.12,6(F20.12))")wr(i),&
                   (dimag(impF0real(ispin,ispin,iorb,jorb,i)),dreal(impF0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
        endif
        call close_units
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
    if(ed_verbose<1)then
       open(unit(1),file="impG0"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impF0"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(3),file="impG0"//string//"_realw"//reg(ed_file_suffix)//".ed")
       open(unit(4),file="impF0"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  !
  subroutine close_units()
    if(ed_verbose<1)then
       close(unit(1))
       close(unit(2))
       close(unit(3))
       close(unit(4))
    endif
  end subroutine close_units
  !
end subroutine print_impG0_superc
