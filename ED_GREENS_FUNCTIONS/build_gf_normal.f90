!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_normal()
  integer :: iorb,jorb,ispin,i
  logical :: MaskBool
  !
  !
  !NORMAL: (default)
  do ispin=1,Nspin
     do iorb=1,Norb
        if(MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
        call lanc_build_gf_normal_c(iorb,ispin)
     enddo
  enddo
  !
  !
  !HYBRID:
  if(bath_type/="normal")then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if(MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
              !if(hybrid)always T; if(replica)T iff following condition is T
              MaskBool=.true.   
              if(bath_type=="replica")MaskBool=(dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
              if(.not.MaskBool)cycle
              call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
           enddo
        enddo
     enddo
     !
     !
     !Put here off-diagonal manipulation by symmetry:
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              !
              !if(hybrid)always T; if(replica)T iff following condition is T
              MaskBool=.true.   
              if(bath_type=="replica")MaskBool=(dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
              !
              impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                   - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
              impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                   - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
              !>>ACTHUNG: this relation might not be true, it depends on the value of the impHloc_ij
              ! if impHloc_ij is REAL then it is true. if CMPLX hermiticity must be ensured
              impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
              impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
           enddo
        enddo
     enddo
  endif
  !
end subroutine build_gf_normal











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
  isite=impIndex(iorb,ispin)
  !
  if(MPI_MASTER)call start_timer
  !
  do istate=1,state_list%size
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
        if(ed_verbose==3.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
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
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        if(ed_verbose==3.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)        
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
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(MPI_MASTER)call stop_timer
end subroutine lanc_build_gf_normal_c








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
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  if(MPI_MASTER)call start_timer
  !
  do istate=1,state_list%size
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
        if(ed_verbose==3.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
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
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_iorb + c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose==3.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose==3.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
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
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,cvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
        endif
        call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(cvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose==3.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
        !
        call setup_Hv_sector(jsector)
        if(ed_sparse_H)call ed_buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_cc,cvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
        endif
        call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
        !
        call delete_Hv_sector()
        !
        deallocate(cvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  !
  if(MPI_MASTER)call stop_timer
  !
end subroutine lanc_build_gf_normal_mix_c


















!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  integer                                    :: isign,iorb,jorb,ispin
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
        impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
  enddo
end subroutine add_to_lanczos_gf_normal








