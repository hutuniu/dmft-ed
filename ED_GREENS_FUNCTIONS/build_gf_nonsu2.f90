!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_nonsu2()
  integer :: izero,iorb,jorb,ispin,jspin,i,io,jo,ifreq
  integer :: isect0,numstates
  real(8) :: norm0
  !
  if(.not.allocated(impGmats))stop "build_gf_nonsu2: Gmats not allocated"
  if(.not.allocated(impGreal))stop "build_gf_nonsu2: Greal not allocated"
  impGmats=zero
  impGreal=zero
  !
  select case(bath_type)
  case("normal")
     !
     !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
     do ispin=1,Nspin
        do iorb=1,Norb
           if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
           call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
        enddo
     enddo
     !same orbital, different spin GF: G_{aa}^{ss'}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
  case("hybrid")
     !
     !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
     do ispin=1,Nspin
        do iorb=1,Norb
           if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
           call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
        enddo
     enddo
     !
     !same orbital, different spin GF: G_{aa}^{ss'}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss'}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
  case("replica")
     !
     !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
     do ispin=1,Nspin
        do iorb=1,Norb
           if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
           call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
        enddo
     enddo
     !
     !same orbital, different spin GF: G_{aa}^{ss'}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss'}(z)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                    call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                 endif
              enddo
           enddo
        enddo
     enddo
     !Here we put the symmetry manipulation
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                    if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.).and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                    !
                    impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                    !
                    impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                         - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     if(ed_para)then
        call SOC_jz_symmetrize(impGmats)
        call SOC_jz_symmetrize(impGreal)
     endif
     !
  end select
end subroutine build_gf_nonsu2








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
  integer                          :: m,i,j,r
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map) :: HI,HJ
  !
  isite=impIndex(iorb,ispin)
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,state_list%size
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !ADD ONE PARTICLE with ISPIN:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ) !note that here you are doing twice the map building...
        vvinit=zero
        do m=1,idim                     !loop over |gs> components m
           i=HI%map(m)                    !map m to Hilbert space state i
           ib = bdecomp(i,2*Ns)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)      !map r back to  jsector
              vvinit(j) = sgn*state_cvec(m)  !build the cdg_ispin|gs> state
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,ispin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !REMOVE ONE PARTICLE with ISPIN:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
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
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_c





!PURPOSE: Evaluate the same different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
  !
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,jspin)  !orbital 2
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do istate=1,state_list%size
     isector     =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec  => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !
     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
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
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !APPLY (c_{iorb,ispin} + c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
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
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !
     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,jspin})|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
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
              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=1*xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,jspin})|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
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
              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        endif
        cnorm2=1*xi*norm2
        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
        deallocate(alfa_,beta_)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c













!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,jspin)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  integer                                    :: isign,iorb,jorb,ispin,jspin
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
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
  ! itype=(3+isign)/2
  diag             = 0.d0
  subdiag          = 0.d0
  Z                = eye(Nlanc)
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
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
  enddo
end subroutine add_to_lanczos_gf_nonsu2












! !+------------------------------------------------------------------+
! !PURPOSE  : DOUBLE PRECISION
! !+------------------------------------------------------------------+
! !PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
! subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_d(iorb,ispin)
!   real(8),allocatable              :: vvinit(:)
!   real(8),allocatable              :: alfa_(:),beta_(:)  
!   integer                          :: iorb,ispin,isite,isector,istate
!   integer                          :: idim,jsector
!   integer                          :: jdim
!   integer                          :: ib(Nlevels)
!   integer                          :: m,i,j,r,numstates
!   real(8)                          :: sgn,norm2,norm0
!   complex(8)                       :: cnorm2
!   integer                          :: Nitermax,Nlanc
!   type(sector_map) :: HI,HJ
!   !
!   Nitermax=lanc_nGFiter
!   allocate(alfa_(Nitermax),beta_(Nitermax))
!   !
!   isite=impIndex(iorb,ispin)
!   !
!   numstates=state_list%size
!   !   
!   if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!   do istate=1,numstates
!      isector    =  es_return_sector(state_list,istate)
!      state_e    =  es_return_energy(state_list,istate)
!      state_vec  => es_return_vector(state_list,istate)
!      norm0=sqrt(dot_product(state_vec,state_vec))
!      if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!      idim  = getdim(isector)
!      call build_sector(isector,HI)
!      !
!      !APPLY c^+_{iorb,ispin}|gs>
!      jsector = getCDGsector(ispin,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ) !note that here you are doing twice the map building...
!         vvinit=zero
!         do m=1,idim                     !loop over |gs> components m
!            i=HI%map(m)                    !map m to Hilbert space state i
!            ib = bdecomp(i,2*Ns)            !i into binary representation
!            if(ib(isite)==0)then          !if impurity is empty: proceed
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)      !map r back to  jsector
!               vvinit(j) = sgn*state_vec(m)  !build the cdg_ispin|gs> state
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,ispin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !APPLY c_{iorb,ispin}|gs>
!      jsector = getCsector(ispin,isector)
!      if(jsector/=0)then
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      nullify(state_vec)
!      deallocate(HI%map)
!      !
!   enddo
!   if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!   deallocate(alfa_,beta_)
! end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_d


! !PURPOSE: Evaluate the  different orbital IORB,JORB, same spin ISPIN impurity GF.
! subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_d(iorb,jorb,ispin)
!   integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
!   integer                          :: idim,jsector
!   integer                          :: jdim
!   integer                          :: ib(Nlevels)
!   integer                          :: m,i,j,r,numstates
!   real(8)                          :: sgn,norm2,norm0
!   complex(8)                       :: cnorm2
!   real(8),allocatable              :: vvinit(:)
!   complex(8),allocatable           :: cvinit(:)
!   real(8),allocatable              :: alfa_(:),beta_(:)
!   integer                          :: Nitermax,Nlanc
!   type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
!   !
!   Nitermax=lanc_nGFiter
!   allocate(alfa_(Nitermax),beta_(Nitermax))
!   isite=impIndex(iorb,ispin)  !orbital 1
!   jsite=impIndex(jorb,ispin)  !orbital 2
!   !
!   numstates=state_list%size
!   !   
!   if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!   do istate=1,numstates
!      isector    =  es_return_sector(state_list,istate)
!      state_e    =  es_return_energy(state_list,istate)
!      state_vec  => es_return_vector(state_list,istate)
!      norm0=sqrt(dot_product(state_vec,state_vec))
!      if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!      !
!      idim  = getdim(isector)
!      call build_sector(isector,HI)
!      !
!      !
!      !APPLY (c^+_{iorb,ispin} + c^+_{jorb,ispin})|gs>
!      jsector = getCDGsector(ispin,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,ispin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !APPLY (c_{iorb,ispin} + c_{jorb,ispin})|gs>
!      jsector = getCsector(ispin,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !
!      !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,ispin})|gs>
!      jsector = getCDGsector(ispin,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,ispin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,ispin})|gs>
!      jsector = getCsector(ispin,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      nullify(state_vec)
!      deallocate(HI%map)
!      !
!   enddo
!   if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!   deallocate(alfa_,beta_)
! end subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_d


! !PURPOSE: Evaluate the same same orbital IORB, different spin ISPIN,JSPIN impurity GF.
! subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_d(iorb,ispin,jspin)
!   integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
!   integer                          :: idim,jsector
!   integer                          :: jdim
!   integer                          :: ib(Nlevels)
!   integer                          :: m,i,j,r,numstates
!   real(8)                          :: sgn,norm2,norm0
!   complex(8)                       :: cnorm2
!   real(8),allocatable              :: vvinit(:)
!   complex(8),allocatable           :: cvinit(:)
!   real(8),allocatable              :: alfa_(:),beta_(:)
!   integer                          :: Nitermax,Nlanc
!   type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
!   !
!   Nitermax=lanc_nGFiter
!   allocate(alfa_(Nitermax),beta_(Nitermax))
!   isite=impIndex(iorb,ispin)  !orbital 1
!   jsite=impIndex(iorb,jspin)  !orbital 2
!   !
!   numstates=state_list%size
!   !   
!   if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!   do istate=1,numstates
!      isector     =  es_return_sector(state_list,istate)
!      state_e    =  es_return_energy(state_list,istate)
!      state_vec  => es_return_vector(state_list,istate)
!      norm0=sqrt(dot_product(state_vec,state_vec))
!      if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!      !
!      idim  = getdim(isector)
!      call build_sector(isector,HI)
!      !
!      !
!      !APPLY (c^+_{iorb,ispin} + c^+_{iorb,jspin})|gs>
!      jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,jspin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !APPLY (c_{iorb,ispin} + c_{iorb,jspin})|gs>
!      jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !
!      !EVALUATE (c^+_{iorb,ispin} + i*c^+_{iorb,jspin})|gs>
!      jsector = getCDGsector(ispin,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,jspin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !EVALUATE (c_{iorb,ispin} - xi*c_{iorb,jspin})|gs>
!      jsector = getCsector(ispin,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      nullify(state_vec)
!      deallocate(HI%map)
!      !
!   enddo
!   if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!   deallocate(alfa_,beta_)
! end subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_d


! !PURPOSE: Evaluate the same different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
! subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_d(iorb,jorb,ispin,jspin)
!   integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
!   integer                          :: idim,jsector
!   integer                          :: jdim
!   integer                          :: ib(Nlevels)
!   integer                          :: m,i,j,r,numstates
!   real(8)                          :: sgn,norm2,norm0
!   complex(8)                       :: cnorm2
!   real(8),allocatable              :: vvinit(:)
!   complex(8),allocatable           :: cvinit(:)
!   real(8),allocatable              :: alfa_(:),beta_(:)
!   integer                          :: Nitermax,Nlanc
!   type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
!   !
!   Nitermax=lanc_nGFiter
!   allocate(alfa_(Nitermax),beta_(Nitermax))
!   isite=impIndex(iorb,ispin)  !orbital 1
!   jsite=impIndex(jorb,jspin)  !orbital 2
!   !
!   numstates=state_list%size
!   !   
!   if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!   do istate=1,numstates
!      isector     =  es_return_sector(state_list,istate)
!      state_e    =  es_return_energy(state_list,istate)
!      state_vec  => es_return_vector(state_list,istate)
!      norm0=sqrt(dot_product(state_vec,state_vec))
!      if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!      !
!      idim  = getdim(isector)
!      call build_sector(isector,HI)
!      !
!      !
!      !APPLY (c^+_{iorb,ispin} + c^+_{jorb,jspin})|gs>
!      jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !APPLY (c_{iorb,ispin} + c_{jorb,jspin})|gs>
!      jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(vvinit(jdim))
!         call build_sector(jsector,HJ)
!         vvinit=0.d0
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               vvinit(j) = vvinit(j) + sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(vvinit,vvinit)
!         vvinit=vvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
!         endif
!         cnorm2=one*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
!         deallocate(vvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !
!      !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,jspin})|gs>
!      jsector = getCDGsector(ispin,isector)
!      if(jsector/=0)then 
!         jdim  = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==0)then
!               call cdg(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==0)then
!               call cdg(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,jspin})|gs>
!      jsector = getCsector(ispin,isector)
!      if(jsector/=0)then
!         jdim   = getdim(jsector)
!         if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!         allocate(cvinit(jdim))
!         call build_sector(jsector,HJ)
!         cvinit=zero
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(isite)==1)then
!               call c(isite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = sgn*state_vec(m)
!            endif
!         enddo
!         do m=1,idim
!            i=HI%map(m)
!            ib = bdecomp(i,2*Ns)
!            if(ib(jsite)==1)then
!               call c(jsite,i,r,sgn)
!               j=binary_search(HJ%map,r)
!               cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
!            endif
!         enddo
!         deallocate(HJ%map)
!         norm2=dot_product(cvinit,cvinit)
!         cvinit=cvinit/sqrt(norm2)
!         call ed_buildH_d(jsector)
!         nlanc=min(jdim,lanc_nGFiter)
!         allocate(alfa_(nlanc),beta_(nlanc))
!         if(MpiStatus)then
!            call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
!         else
!            call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
!         endif
!         cnorm2=+xi*norm2
!         call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
!         deallocate(cvinit)
!         if(spH0%status)call sp_delete_matrix(spH0)
!      endif
!      !
!      nullify(state_vec)
!      deallocate(HI%map)
!      !
!   enddo
!   if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!   deallocate(alfa_,beta_)
! end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_d



















!!PURPOSE: Evaluate the same different orbital IORB,JORB, same spin ISPIN impurity GF.
!subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_c(iorb,jorb,ispin)
!  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
!  integer                          :: idim,jsector
!  integer                          :: jdim
!  integer                          :: ib(Nlevels)
!  integer                          :: m,i,j,r
!  real(8)                          :: sgn,norm2,norm0
!  complex(8)                       :: cnorm2
!  complex(8),allocatable           :: vvinit(:)
!  real(8),allocatable              :: alfa_(:),beta_(:)
!  integer                          :: Nitermax,Nlanc
!  type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
!  !
!  isite=impIndex(iorb,ispin)  !orbital 1
!  jsite=impIndex(jorb,ispin)  !orbital 2
!  !
!  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!  !
!  do istate=1,state_list%size
!     isector     =  es_return_sector(state_list,istate)
!     state_e     =  es_return_energy(state_list,istate)
!     state_cvec  => es_return_cvector(state_list,istate)
!     norm0=sqrt(dot_product(state_cvec,state_cvec))
!     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!     !
!     idim  = getdim(isector)
!     call build_sector(isector,HI)
!     !
!     !
!     !APPLY (c^+_{iorb,ispin} + c^+_{jorb,ispin})|gs>
!     jsector = getCDGsector(ispin,isector)
!     if(jsector/=0)then 
!        jdim  = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==0)then
!              call cdg(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==0)then
!              call cdg(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=one*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,ispin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !APPLY (c_{iorb,ispin} + c_{jorb,ispin})|gs>
!     jsector = getCsector(ispin,isector)
!     if(jsector/=0)then
!        jdim   = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==1)then
!              call c(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==1)then
!              call c(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=one*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !
!     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{jorb,ispin})|gs>
!     jsector = getCDGsector(ispin,isector)
!     if(jsector/=0)then 
!        jdim  = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==0)then
!              call cdg(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==0)then
!              call cdg(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=+xi*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,ispin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !EVALUATE (c_{iorb,ispin} - xi*c_{jorb,ispin})|gs>
!     jsector = getCsector(ispin,isector)
!     if(jsector/=0)then
!        jdim   = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==1)then
!              call c(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==1)then
!              call c(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=+xi*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,ispin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     nullify(state_vec)
!     deallocate(HI%map)
!     !
!  enddo
!  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!  deallocate(alfa_,beta_)
!end subroutine lanc_build_gf_nonsu2_mixOrb_diagSpin_c



!!PURPOSE: Evaluate the same same orbital IORB, different spin ISPIN,JSPIN impurity GF.
!subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_c(iorb,ispin,jspin)
!  integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
!  integer                          :: idim,jsector
!  integer                          :: jdim
!  integer                          :: ib(Nlevels)
!  integer                          :: m,i,j,r
!  real(8)                          :: sgn,norm2,norm0
!  complex(8)                       :: cnorm2
!  complex(8),allocatable           :: vvinit(:)
!  real(8),allocatable              :: alfa_(:),beta_(:)
!  integer                          :: Nitermax,Nlanc
!  type(sector_map) :: HI,HJ    !map of the Sector S to Hilbert space H
!  !
!  isite=impIndex(iorb,ispin)  !orbital 1
!  jsite=impIndex(iorb,jspin)  !orbital 2
!  !
!  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
!  !
!  do istate=1,state_list%size
!     isector     =  es_return_sector(state_list,istate)
!     state_e    =  es_return_energy(state_list,istate)
!     state_cvec  => es_return_cvector(state_list,istate)
!     norm0=sqrt(dot_product(state_cvec,state_cvec))
!     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
!     !
!     idim  = getdim(isector)
!     call build_sector(isector,HI)
!     !
!     !
!     !APPLY (c^+_{iorb,ispin} + c^+_{iorb,jspin})|gs>
!     jsector = getCDGsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!     if(jsector/=0)then 
!        jdim  = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==0)then
!              call cdg(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==0)then
!              call cdg(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=one*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,jspin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !APPLY (c_{iorb,ispin} + c_{iorb,jspin})|gs>
!     jsector = getCsector(ispin,isector) !this is the same sector I'd get using getCDGsector(JSPIN,isector)
!     if(jsector/=0)then
!        jdim   = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==1)then
!              call c(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==1)then
!              call c(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=one*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !
!     !EVALUATE (c^+_{iorb,ispin} + i*c^+_{iorb,jspin})|gs>
!     jsector = getCDGsector(ispin,isector)
!     if(jsector/=0)then 
!        jdim  = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' add particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==0)then
!              call cdg(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==0)then
!              call cdg(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=+xi*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,jspin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     !EVALUATE (c_{iorb,ispin} - xi*c_{iorb,jspin})|gs>
!     jsector = getCsector(ispin,isector)
!     if(jsector/=0)then
!        jdim   = getdim(jsector)
!        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,I3)")' del particle:',getn(jsector)
!        allocate(vvinit(jdim))
!        call build_sector(jsector,HJ)
!        vvinit=zero
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(isite)==1)then
!              call c(isite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = sgn*state_cvec(m)
!           endif
!        enddo
!        do m=1,idim
!           i=HI%map(m)
!           ib = bdecomp(i,2*Ns)
!           if(ib(jsite)==1)then
!              call c(jsite,i,r,sgn)
!              j=binary_search(HJ%map,r)
!              vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
!           endif
!        enddo
!        deallocate(HJ%map)
!        norm2=dot_product(vvinit,vvinit)
!        vvinit=vvinit/sqrt(norm2)
!        !
!        call setup_Hv_sector(jsector)
!        if(ed_sparse_H)call ed_buildH_c()
!        !
!        nlanc=min(jdim,lanc_nGFiter)
!        allocate(alfa_(nlanc),beta_(nlanc))
!        if(MpiStatus)then
!           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
!        else
!           call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
!        endif
!        cnorm2=+xi*norm2
!        call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,jspin)
!        !
!        call delete_Hv_sector()
!        !
!        if(spH0%status)call sp_delete_matrix(spH0)
!     endif
!     !
!     nullify(state_cvec)
!     deallocate(HI%map)
!     !
!  enddo
!  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
!  deallocate(alfa_,beta_)
!end subroutine lanc_build_gf_nonsu2_diagOrb_mixSpin_c
