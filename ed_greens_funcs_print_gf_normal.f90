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
        call close_units()
     enddo
  endif
  !
contains
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_unit(1)
    open(unit(1),file="impGpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  !
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine print_poles_weights_normal




subroutine print_sigma_normal
  integer                                           :: i,ispin,isign,unit(2),iorb,jorb
  character(len=20)                                 :: suffix
  integer,dimension(:),allocatable                  :: getIorb,getJorb
  integer                                           :: totNorb,l
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
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<4)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),(dimag(impSmats(ispin,ispin,iorb,jorb,i)),dreal(impSmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(2),"(F20.12,6(F20.12))")wr(i),(dimag(impSreal(ispin,ispin,iorb,jorb,i)),dreal(impSreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
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
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<4)then
       open(unit(1),file="impSigma"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impSigma"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  !
  subroutine close_units()
    if(ed_verbose<4)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_sigma_normal





subroutine print_impg_normal
  integer                                           :: i,ispin,isign,unit(2),iorb,jorb
  character(len=20)                                 :: suffix
  integer,dimension(:),allocatable                  :: getIorb,getJorb
  integer                                           :: totNorb,l
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
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<2)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),(dimag(impGmats(ispin,ispin,iorb,jorb,i)),dreal(impGmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(2),"(F20.12,6(F20.12))")wr(i),(dimag(impGreal(ispin,ispin,iorb,jorb,i)),dreal(impGreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
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
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<2)then
       open(unit(1),file="impG"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impG"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  !
  subroutine close_units()
    if(ed_verbose<2)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_impg_normal





subroutine print_impg0_normal
  integer                                           :: i,ispin,isign,unit(2),iorb,jorb
  character(len=20)                                 :: suffix
  integer,dimension(:),allocatable                  :: getIorb,getJorb
  integer                                           :: totNorb,l
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
  if(ED_MPI_ID==0)then
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
        call open_units(reg(suffix))
        if(ed_verbose<1)then
           do i=1,Lmats
              write(unit(1),"(F20.12,6(F20.12))")wm(i),(dimag(impG0mats(ispin,ispin,iorb,jorb,i)),dreal(impG0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
           enddo
           do i=1,Lreal
              write(unit(2),"(F20.12,6(F20.12))")wr(i),(dimag(impG0real(ispin,ispin,iorb,jorb,i)),dreal(impG0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
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
  !
  subroutine open_units(string)
    character(len=*) :: string
    unit=free_units(size(unit))
    if(ed_verbose<1)then
       open(unit(1),file="impG0"//string//"_iw"//reg(ed_file_suffix)//".ed")
       open(unit(2),file="impG0"//string//"_realw"//reg(ed_file_suffix)//".ed")
    endif
  end subroutine open_units
  !
  subroutine close_units()
    if(ed_verbose<1)then
       close(unit(1))
       close(unit(2))
    endif
  end subroutine close_units
end subroutine print_impg0_normal
