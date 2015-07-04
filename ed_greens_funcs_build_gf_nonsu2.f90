!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_nonsu2()
  integer :: izero,iorb,jorb,ispin,jspin,i
  integer :: isect0,numstates
  real(8) :: norm0
  if(.not.allocated(impGmats))stop "build_gf_nonsu2: Gmats not allocated"
  if(.not.allocated(impGreal))stop "build_gf_nonsu2: Greal not allocated"
  impGmats=zero
  impGreal=zero
  write(LOGfile,"(A)")"Get impurity Greens functions:"

  !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
  do ispin=1,Nspin
     do iorb=1,Norb
        if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))
        select case(ed_type)
        case default
           call lanc_build_gf_nonsu2_diagOrb_diagSpin_d(iorb,ispin)
        case ('c')
           call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
        end select
     enddo
  enddo

  !Here we evaluate the same orbital, different spin GF: G_{aa}^{ss`}(z)
  do ispin=1,Nspin
     do jspin=ispin+1,Nspin
        do iorb=1,Norb
           if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
           select case(ed_type)
           case default
              call lanc_build_gf_nonsu2_diagOrb_mixSpin_d(iorb,ispin,jspin)
           case ('c')
              call lanc_build_gf_nonsu2_diagOrb_mixSpin_c(iorb,ispin,jspin)
           end select
        enddo
     enddo
  enddo
  !Put here off-diagonal manipulation by symmetry:
  do iorb=1,Norb
     do ispin=1,Nspin
        do jspin=ispin+1,Nspin
           !
           impGmats(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,iorb,:) & !this is Smix in our notation
                - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(jspin,jspin,iorb,iorb,:))
           !
           impGreal(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,iorb,:) & !this is Smix in our notation
                - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(jspin,jspin,iorb,iorb,:))
           !
           impGmats(jspin,ispin,iorb,iorb,:) = impGmats(ispin,jspin,iorb,iorb,:)
           impGreal(jspin,ispin,iorb,iorb,:) = impGreal(ispin,jspin,iorb,iorb,:)
        enddo
     enddo
  enddo



  if(bath_type=='hybrid')then
     !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))
              select case(ed_type)
              case default
                 call lanc_build_gf_nonsu2_mixOrb_diagSpin_d(iorb,jorb,ispin)
              case ('c')
                 call lanc_build_gf_nonsu2_mixOrb_diagSpin_c(iorb,jorb,ispin)                    
              end select
           enddo
        enddo
     enddo
     !Put here off-diagonal manipulation by symmetry:
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              !
              impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) & !this is Gmix in our notation
                   - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
              !
              impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                   - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
              !
              impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
              impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
           enddo
        enddo
     enddo


     !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss`}
     do ispin=1,Nspin
        do jspin=ispin+1,Nspin
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))
                 select case(ed_type)
                 case default
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_d(iorb,jorb,ispin,jspin)
                 case ('c')
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)                    
                 end select
              enddo
           enddo
        enddo
     enddo
     do ispin=1,Nspin
        do jspin=ispin+1,Nspin
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 !
                 impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) & !this is Fmix in our notation
                      - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(jspin,jspin,jorb,jorb,:))
                 !
                 impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                      - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(jspin,jspin,jorb,jorb,:))
                 !
                 impGmats(jspin,ispin,jorb,iorb,:) = impGmats(ispin,jspin,iorb,jorb,:)
                 impGreal(jspin,ispin,jorb,iorb,:) = impGreal(ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  endif
  !
end subroutine build_gf_nonsu2







!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
!PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_d(iorb,ispin)
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
  integer,allocatable,dimension(:) :: HImap,HJmap
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !APPLY c^+_{iorb,ispin}|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap) !note that here you are doing twice the map building...
        vvinit=0d0
        do m=1,idim                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsector
              vvinit(j) = sgn*state_vec(m)  !build the cdg_ispin|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY c_{iorb,ispin}|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_d


!PURPOSE: Evaluate the  different orbital IORB,JORB, same spin ISPIN impurity GF.
subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_d(iorb,jorb,ispin)
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
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,ispin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{jorb,ispin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,ispin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,ispin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_d


!PURPOSE: Evaluate the same same orbital IORB, different spin ISPIN,JSPIN impurity GF.
subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_d(iorb,ispin,jspin)
  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
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
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(iorb,jspin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector     =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{iorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{iorb,jspin})|gs>
     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{iorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,jspin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{iorb,jspin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_d


!PURPOSE: Evaluate the same different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_d(iorb,jorb,ispin,jspin)
  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
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
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,jspin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector     =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_d




















!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
!PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
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
  integer,allocatable,dimension(:) :: HImap,HJmap
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !ADD ONE PARTICLE with ISPIN:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap) !note that here you are doing twice the map building...
        vvinit=0d0
        do m=1,idim                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsector
              vvinit(j) = sgn*state_cvec(m)  !build the cdg_ispin|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !REMOVE ONE PARTICLE with ISPIN:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_c


!PURPOSE: Evaluate the same different orbital IORB,JORB, same spin ISPIN impurity GF.
subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_c(iorb,jorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector     =  es_return_sector(state_list,istate)
     state_e     =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,ispin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{jorb,ispin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,ispin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,ispin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_c


!PURPOSE: Evaluate the same same orbital IORB, different spin ISPIN,JSPIN impurity GF.
subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_c(iorb,ispin,jspin)
  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(iorb,jspin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector     =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{iorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{iorb,jspin})|gs>
     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{iorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{iorb,jspin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_c


!PURPOSE: Evaluate the same different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,jspin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do istate=1,numstates
     isector     =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' add particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,I3,I15)")' del particle:',getn(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=zero
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_cc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c













!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb,ispin,jspin)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb,jorb,ispin,jspin
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
  Egs = state_list%emin       !get the gs energy
  pesoBZ = vnorm2/zeta_function
  if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
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
        impGmats(ispin,jspin,iorb,jorb,i)=impGmats(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        impGreal(ispin,jspin,iorb,jorb,i)=impGreal(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     GFpoles(ispin,jspin,iorb,jorb,itype,j)   = isign*de
     GFweights(ispin,jspin,iorb,jorb,itype,j) = peso
  enddo


end subroutine add_to_lanczos_gf_nonsu2
