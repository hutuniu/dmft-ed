!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_superc()
  integer :: izero,iorb,jorb,ispin,i
  integer :: isect0,numstates
  real(8) :: norm0
  logical :: verbose
  verbose=.false.;if(ed_verbose<1)verbose=.true.
  if(.not.allocated(Gaux_mats))allocate(Gaux_mats(3,Lmats))
  if(.not.allocated(Gaux_real))allocate(Gaux_real(3,Lreal))
  Gaux_mats=zero
  Gaux_real=zero
  !
  ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
  do iorb=1,Norb
     if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G&F_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))
     call lanc_build_gf_superc_d(iorb,verbose)
  enddo
  do iorb=1,Norb
     impGmats(ispin,ispin,iorb,iorb,:) = Gaux_mats(1,:) !this is the up,up function
     impGreal(ispin,ispin,iorb,iorb,:) = Gaux_real(1,:)
     impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(Gaux_mats(3,:)-Gaux_mats(2,:)-Gaux_mats(1,:))
     impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(Gaux_real(3,:)-Gaux_real(2,:)-Gaux_real(1,:))
  enddo
  !
  !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
  if(bath_type=='hybrid')then
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G_l"//&
                reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
           call lanc_build_gf_superc_mix_d(iorb,jorb,verbose)
        enddo
     enddo

  endif

  deallocate(Gaux_mats,Gaux_real)
end subroutine build_gf_superc







!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_superc_d(iorb,iverbose)
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0,isz0,jsz0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap,HJmap
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)

     !APPLY CDG_UP
     jsect0 = getCDGsector(1,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply cdg_up:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb)==0)then           !if impurity is empty: proceed
              call cdg(iorb,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ichan=1)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !APPLY C_DW
     jsect0 = getCsector(2,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply c_dw:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb+Ns)==1)then           !if impurity is empty: proceed
              call c(iorb+Ns,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ichan=2)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !APPLY CDG_UP + C_DW
     isz0 = getsz(isect0)
     if(isz0<Ns)then
        jsz0   = isz0+1
        jsect0 = getsector(jsz0,1)
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply c_dw:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb)==0)then           !if impurity is empty: proceed
              call cdg(iorb,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb+Ns)==1)then           !if impurity is empty: proceed
              call c(iorb+Ns,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = vvinit(j) + sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ichan=3)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif





     !APPLY C_UP
     jsect0 = getCsector(1,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply cdg_up:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb)==1)then           !if impurity is empty: proceed
              call c(iorb,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ichan=1)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !APPLY CDG_DW
     jsect0 = getCDGsector(2,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply c_dw:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb+Ns)==0)then           !if impurity is empty: proceed
              call cdg(iorb+Ns,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ichan=2)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !APPLY C_UP + CDG_DW
     isz0 = getsz(isect0)
     if(isz0>-Ns)then
        jsz0   = isz0-1
        jsect0 = getsector(jsz0,1)
        jdim0  = getdim(jsect0)
        if(iverbose_.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A,I3,I15)")'apply c_dw:',getsz(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb)==1)then           !if impurity is empty: proceed
              call c(iorb,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb+Ns)==0)then           !if impurity is empty: proceed
              call cdg(iorb+Ns,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = vvinit(j) + sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call ed_buildH_d(jsect0)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ichan=3)
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
end subroutine lanc_build_gf_superc_d






!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_superc(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb,ichan)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb,jorb,ichan
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
  Egs = state_list%emin       !get the gs energy
  pesoBZ = vnorm2/zeta_function
  if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
  !
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
        Gaux_mats(ichan,i)=Gaux_mats(ichan,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        Gaux_real(ichan,i)=Gaux_real(ichan,i) + peso/(iw-isign*de)
     enddo
  enddo
end subroutine add_to_lanczos_gf_superc
