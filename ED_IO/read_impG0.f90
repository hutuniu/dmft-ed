subroutine read_impG0_normal
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
  if(l/=totNorb)stop "read_gf_normal error counting the orbitals"
  !!
  !Print the impurity functions:
  do ispin=1,Nspin
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
        call sread("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impG0mats(ispin,ispin,iorb,jorb,:))
        call sread("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ispin,ispin,iorb,jorb,:))
     enddo
  enddo
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine read_impG0_normal






subroutine read_impG0_superc
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
  if(l/=totNorb)stop "read_gf_superc error counting the orbitals"
  !!
  !!
  !!
  !!PRINT OUT GF:
  do ispin=1,Nspin
     do l=1,totNorb
        iorb=getIorb(l)
        jorb=getJorb(l)
        suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
        call sread("impG0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"  ,wm,impG0mats(ispin,ispin,iorb,jorb,:))
        call sread("impF0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impF0mats(ispin,ispin,iorb,jorb,:))
        !
        call sread("impG0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ispin,ispin,iorb,jorb,:))
        call sread("impF0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed" ,wr,impF0real(ispin,ispin,iorb,jorb,:))
     enddo
  enddo
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine read_impG0_superc













subroutine read_impG0_nonsu2
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
                 if(io<jo)cycle
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
                 if(io<jo)cycle
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
  if(l/=totNso)stop "read_gf_nonsu2 error counting the spin-orbitals"
  !!
  !!PRINT OUT GF:
  do l=1,totNso
     iorb=getIorb(l)
     jorb=getJorb(l)
     ispin=getIspin(l)
     jspin=getJspin(l)
     !
     suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
     call sread("impG0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impG0mats(ispin,jspin,iorb,jorb,:))
     call sread("impG0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ispin,jspin,iorb,jorb,:))
  enddo
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine read_impG0_nonsu2



