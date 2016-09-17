!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_normal()
  integer :: iorb,jorb,ispin,i
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  !
  do ispin=1,Nspin
     do iorb=1,Norb
        if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))
        select case(ed_type)
        case default
           call lanc_build_gf_normal_d(iorb,ispin)
        case ('c')
           call lanc_build_gf_normal_c(iorb,ispin)
        end select
     enddo
  enddo
  !
  if(bath_type=="hybrid")then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")&
                   "Get G_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
              select case(ed_type)
              case default
                 call lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
              case ('c')
                 call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)                    
              end select
           enddo
        enddo
     enddo
     !Put here off-diagonal manipulation by symmetry:
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           impGmats(:,:,iorb,jorb,:) = 0.5d0*(impGmats(:,:,iorb,jorb,:) &
                - (one-xi)*impGmats(:,:,iorb,iorb,:) - (one-xi)*impGmats(:,:,jorb,jorb,:))
           impGreal(:,:,iorb,jorb,:) = 0.5d0*(impGreal(:,:,iorb,jorb,:) &
                - (one-xi)*impGreal(:,:,iorb,iorb,:) - (one-xi)*impGreal(:,:,jorb,jorb,:))
           impGmats(:,:,jorb,iorb,:) = impGmats(:,:,iorb,jorb,:)
           impGreal(:,:,jorb,iorb,:) = impGreal(:,:,iorb,jorb,:)
        enddo
     enddo
  endif
  !
  if(bath_type=="replica")then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if((dmft_bath%mask(ispin,ispin,iorb,jorb,1).eqv. .true.).or.(dmft_bath%mask(ispin,ispin,iorb,jorb,2).eqv. .true.))then
                 if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")&
                      "Get G_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
                 select case(ed_type)
                 case default
                    call lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
                 case ('c')
                    call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)                    
                 end select
              endif
           enddo
        enddo
     enddo
     !Put here off-diagonal manipulation by symmetry:
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if((dmft_bath%mask(ispin,ispin,iorb,jorb,1).eqv. .true.).or.(dmft_bath%mask(ispin,ispin,iorb,jorb,2).eqv. .true.))then
                 impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                      - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                 impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                      - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                 impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                 impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
              endif
           enddo
        enddo
     enddo
  endif
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
end subroutine build_gf_normal






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
    open(unit(1),file="Gpoles_weights"//string//reg(ed_file_suffix)//".ed")
  end subroutine open_units
  subroutine close_units()
    close(unit(1))
  end subroutine close_units
end subroutine rebuild_gf_normal





!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_d(iorb,ispin)
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,ispin,isite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map)                  :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)

     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_d


!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_c(iorb,ispin)
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: iorb,ispin,isite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim                     !loop over |gs> components m
           i=HI%map(m)                    !map m to Hilbert space state i
           ib = bdecomp(i,2*Ns)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)      !map r back to  jsector
              vvinit(j) = sgn*state_cvec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3,I15)")' del particle:',&
             getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_c






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_mix_d



!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_mix_c









!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_normal(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb,ispin)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb,jorb,ispin
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
  Egs = state_list%emin       !get the gs energy
  !
  if((finiteT).and.(beta*(Ei-Egs).lt.200))then
     pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
  elseif(.not.finiteT)then
     pesoBZ = vnorm2/zeta_function
  else
     pesoBZ=0.d0
  endif
  !
  !pesoBZ = vnorm2/zeta_function
  !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
  !
  itype=(3+isign)/2
  diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
  forall(i=1:Nlanc)Z(i,i)=1.d0
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  do j=1,nlanc
     de = diag(j)-Ei
     peso = pesoBZ*Z(1,j)*Z(1,j)
     do i=1,Lmats
        iw=xi*wm(i)
        impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     GFpoles(ispin,ispin,iorb,jorb,itype,j)   = isign*de
     GFweights(ispin,ispin,iorb,jorb,itype,j) = peso
  enddo


end subroutine add_to_lanczos_gf_normal
