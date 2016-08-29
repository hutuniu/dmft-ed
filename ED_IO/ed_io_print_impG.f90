subroutine print_impG_normal
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
end subroutine print_impG_normal









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










subroutine print_impG_nonsu2
  integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
  integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
  integer                          :: totNso,totNorb,totNspin,l
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
     totNorb =Norb*(Norb+1)/2
     totNspin=Nspin*(Nspin+1)/2
     totNso  =totNorb*totNspin
     allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
     l=0
     do iorb=1,Norb
        do jorb=iorb,Norb
           do ispin=1,Nspin
              do jspin=ispin,Nspin
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











subroutine ed_print_impG_normal_lattice(iprint)
  integer :: iprint
  integer :: ispin,jspin,iorb,jorb
  if(ED_MPI_MASTER)then
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     allocate(wm(Lmats))
     allocate(wr(Lreal))
     wm = pi/beta*(2*arange(1,Lmats)-1)
     wr = linspace(wini,wfin,Lreal)
     select case(iprint)
     case (0)
        write(LOGfile,*)"Gimp not written on file."
     case(1)                  !print only diagonal elements
        write(LOGfile,*)"write spin-orbital diagonal elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
              call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,iorb,:),wm)
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
              call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,iorb,:),wr)
           enddo
        enddo
     case(2)                  !print spin-diagonal, all orbitals 
        write(LOGfile,*)"write spin diagonal and all orbitals elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                 call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                 call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,jorb,:),wr)
              enddo
           enddo
        enddo
     case default                  !print all off-diagonals
        write(LOGfile,*)"write all elements:"
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                    call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,jspin,iorb,jorb,:),wm)
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                    call store_data("LGimp"//reg(suffix),Grealii(:,ispin,jspin,iorb,jorb,:),wr)
                 enddo
              enddo
           enddo
        enddo
     end select
  endif
end subroutine ed_print_impG_normal_lattice



subroutine ed_print_impG_superc_lattice(iprint)
  integer :: iprint
  integer :: ispin,jspin,iorb,jorb
  if(ED_MPI_MASTER)then
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     allocate(wm(Lmats))
     allocate(wr(Lreal))
     wm = pi/beta*(2*arange(1,Lmats)-1)
     wr = linspace(wini,wfin,Lreal)
     select case(iprint)
     case (0)
        write(LOGfile,*)"Gimp not written on file."
     case(1)                  !print only diagonal elements
        write(LOGfile,*)"write spin-orbital diagonal elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
              call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,iorb,:),wm)
              call store_data("LFimp"//reg(suffix),Fmatsii(:,ispin,ispin,iorb,iorb,:),wm)
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
              call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,iorb,:),wr)
              call store_data("LFimp"//reg(suffix),Frealii(:,ispin,ispin,iorb,iorb,:),wr)
           enddo
        enddo
     case(2)                  !print spin-diagonal, all orbitals 
        write(LOGfile,*)"write spin diagonal and all orbitals elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                 call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                 call store_data("LFimp"//reg(suffix),Fmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                 call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,jorb,:),wr)
                 call store_data("LFimp"//reg(suffix),Frealii(:,ispin,ispin,iorb,jorb,:),wr)
              enddo
           enddo
        enddo
     case default                  !print all off-diagonals
        write(LOGfile,*)"write all elements:"
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                    call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,jspin,iorb,jorb,:),wm)
                    call store_data("LFimp"//reg(suffix),Fmatsii(:,ispin,jspin,iorb,jorb,:),wm)
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                    call store_data("LGimp"//reg(suffix),Grealii(:,ispin,jspin,iorb,jorb,:),wr)
                    call store_data("LFimp"//reg(suffix),Frealii(:,ispin,jspin,iorb,jorb,:),wr)
                 enddo
              enddo
           enddo
        enddo
     end select
  endif
end subroutine ed_print_impG_superc_lattice



subroutine ed_print_impG_nonsu2_lattice(iprint)
  integer :: iprint
  integer :: ispin,jspin,iorb,jorb
  if(ED_MPI_MASTER)then
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     allocate(wm(Lmats))
     allocate(wr(Lreal))
     wm = pi/beta*(2*arange(1,Lmats)-1)
     wr = linspace(wini,wfin,Lreal)
     select case(iprint)
     case (0)
        write(LOGfile,*)"Gimp not written on file."
     case(1)                  !print only diagonal elements
        write(LOGfile,*)"write spin-orbital diagonal elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
              call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,iorb,:),wm)
              suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
              call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,iorb,:),wr)
           enddo
        enddo
     case(2)                  !print spin-diagonal, all orbitals 
        write(LOGfile,*)"write spin diagonal and all orbitals elements:"
        do ispin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                 call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                 suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                 call store_data("LGimp"//reg(suffix),Grealii(:,ispin,ispin,iorb,jorb,:),wr)
              enddo
           enddo
        enddo
     case default                  !print all off-diagonals
        write(LOGfile,*)"write all elements:"
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                    call store_data("LGimp"//reg(suffix),Gmatsii(:,ispin,jspin,iorb,jorb,:),wm)
                    suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                    call store_data("LGimp"//reg(suffix),Grealii(:,ispin,jspin,iorb,jorb,:),wr)
                 enddo
              enddo
           enddo
        enddo
     end select
  endif
end subroutine ed_print_impG_nonsu2_lattice
