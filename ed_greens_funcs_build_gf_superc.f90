!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_superc()
  integer    :: iorb,jorb,ispin,i,isign
  complex(8) :: barGmats(Norb,Lmats),barGreal(Norb,Lreal)
  !
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  if(.not.allocated(auxGmats))allocate(auxGmats(3,Lmats))
  if(.not.allocated(auxGreal))allocate(auxGreal(3,Lreal))
  if(.not.allocated(auxGpoles))allocate(auxGpoles(2,3,lanc_nGFiter))
  if(.not.allocated(auxGweights))allocate(auxGweights(2,3,lanc_nGFiter))
  auxGmats=zero
  auxGreal=zero
  auxGpoles=zero
  auxGweights=zero
  barGmats=zero
  barGreal=zero
  !
  ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
  do iorb=1,Norb
     auxGmats=zero
     auxGreal=zero
     if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get G&F_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))
     call lanc_build_gf_superc_d(iorb)
     !
     impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{iorb,iorb} = G_{up,up;iorb,iorb}
     barGmats(                 iorb,:) = auxGmats(2,:) !this is \bar{G}_{iorb,iorb} = \bar{G}_{dw,dw;iorb,iorb}
     impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)
     barGreal(                 iorb,:) = auxGreal(2,:)
     impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-auxGmats(1,:)-auxGmats(2,:))
     impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-auxGreal(1,:)-auxGreal(2,:))

     !
     GFpoles(ispin,ispin,iorb,iorb,:,:)   = auxGpoles(:,1,:)
     GFweights(ispin,ispin,iorb,iorb,:,:) = auxGweights(:,1,:)
     !

     ! UNCOMMENT THIS AND FOLLOWING LINES MARKED WITH >ANOMAL TO USE THE MORE GENERAL ALGORITHM
     ! FOR THE EVALUATION OF THE ANOMALOUS GF
     ! >ANOMAL
     ! impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
     ! impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
     ! <ANOMAL
  enddo
  !
  !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
  if(bath_type=='hybrid')then
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")&
                "Get G_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
           call lanc_build_gf_superc_mix_d(iorb,jorb)
           impGmats(ispin,ispin,iorb,jorb,:) = auxGmats(3,:)
           impGreal(ispin,ispin,iorb,jorb,:) = auxGreal(3,:)
        enddo
     enddo
     !
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
           impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
        enddo
     enddo
  endif
  deallocate(auxGmats,auxGreal)
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
end subroutine build_gf_superc







!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_superc_d(iorb)
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim,isz,jsz
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
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  !
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
     !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
     jsector = getCDGsector(1,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c^+_up:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb)==0)then
              call cdg(iorb,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=1)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
     jsector = getCsector(1,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c_up:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb)==1)then
              call c(iorb,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=1)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif




     !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
     jsector = getCsector(2,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A23,I3,I15)")'apply c_dw:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb+Ns)==1)then
              call c(iorb+Ns,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=2)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
     jsector = getCDGsector(2,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A23,I3,I15)")'apply c^+_dw:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(iorb+Ns)==0)then           !if impurity is empty: proceed
              call cdg(iorb+Ns,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsector
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
        cnorm2=one*norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=2)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif




     !EVALUATE [c^+_{up,iorb} + c_{dw,iorb}]|gs> --> A_{iorb,iorb}
     isz = getsz(isector)
     if(isz<Ns)then
        jsz   = isz+1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c^+_up + c_dw:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb)==0)then
              call cdg(iorb,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb+Ns)==1)then
              call c(iorb+Ns,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=3)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE [c_{up,iorb} + c^+_{dw,iorb}]|gs>  --> A_{iorb,iorb}
     isz = getsz(isector)
     if(isz>-Ns)then
        jsz   = isz-1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A23,I3,I15)")'apply c_up + c^+_dw:',getsz(jsector),jdim
        allocate(HJmap(jdim),vvinit(jdim))
        call build_sector(jsector,HJmap)
        vvinit=0.d0
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb)==1)then
              call c(iorb,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(iorb+Ns)==0)then
              call cdg(iorb+Ns,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=3)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif



     ! ! <ANOMAL
     ! !EVALUATE [c^+_{up,iorb} + xi*c_{dw,iorb}]|gs> --> -xi*B_{iorb,iorb}
     ! isz = getsz(isector)
     ! if(isz<Ns)then
     !    jsz   = isz+1
     !    jsector = getsector(jsz,1)
     !    jdim  = getdim(jsector)
     !    if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c^+_up + xi*c_dw:',getsz(jsector),jdim
     !    allocate(HJmap(jdim),cvinit(jdim))
     !    call build_sector(jsector,HJmap)
     !    cvinit=0.d0
     !    do m=1,idim
     !       i=HImap(m)
     !       call bdecomp(i,ib)
     !       if(ib(iorb)==0)then
     !          call cdg(iorb,i,r,sgn)
     !          j=binary_search(HJmap,r)
     !          cvinit(j) = sgn*state_vec(m)
     !       endif
     !    enddo
     !    do m=1,idim
     !       i=HImap(m)
     !       call bdecomp(i,ib)
     !       if(ib(iorb+Ns)==1)then
     !          call c(iorb+Ns,i,r,sgn)
     !          j=binary_search(HJmap,r)
     !          cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
     !       endif
     !    enddo
     !    deallocate(HJmap)
     !    norm2=dot_product(cvinit,cvinit)
     !    cvinit=cvinit/sqrt(norm2)
     !    alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
     !    call ed_buildH_d(jsector)
     !    call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
     !    cnorm2=-xi*norm2
     !    call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=3)
     !    deallocate(cvinit)
     !    if(spH0%status)call sp_delete_matrix(spH0)
     ! endif
     ! !
     ! !EVALUATE [c_{up,iorb} - xi*c^+_{dw,iorb}]|gs> --> -xi*B_{iorb,iorb}
     ! isz = getsz(isector)
     ! if(isz>-Ns)then
     !    jsz   = isz-1
     !    jsector = getsector(jsz,1)
     !    jdim  = getdim(jsector)
     !    if(ed_verbose<1.AND.ED_MPI_ID==0)&
     !         write(LOGfile,"(A23,I3,I15)")'apply c_up - xi*c^+_dw:',getsz(jsector),jdim
     !    allocate(HJmap(jdim),cvinit(jdim))
     !    call build_sector(jsector,HJmap)
     !    cvinit=0.d0
     !    do m=1,idim
     !       i=HImap(m)
     !       call bdecomp(i,ib)
     !       if(ib(iorb)==1)then
     !          call c(iorb,i,r,sgn)
     !          j=binary_search(HJmap,r)
     !          cvinit(j) = sgn*state_vec(m)
     !       endif
     !    enddo
     !    do m=1,idim
     !       i=HImap(m)
     !       call bdecomp(i,ib)
     !       if(ib(iorb+Ns)==0)then
     !          call cdg(iorb+Ns,i,r,sgn)
     !          j=binary_search(HJmap,r)
     !          cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
     !       endif
     !    enddo
     !    deallocate(HJmap)
     !    norm2=dot_product(cvinit,cvinit)
     !    cvinit=cvinit/sqrt(norm2)
     !    alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
     !    call ed_buildH_d(jsector)
     !    call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc,lanc_spHtimesV_dc)
     !    cnorm2=-xi**norm2
     !    call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=3)
     !    deallocate(cvinit)
     !    if(spH0%status)call sp_delete_matrix(spH0)
     ! endif
     ! ! <ANOMAL
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_superc_d





!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_superc_mix_d(iorb,jorb)
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,jorb,isector,istate
  integer                          :: idim,jsector,isite
  integer                          :: jdim,isz,jsz,jsite
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  integer,allocatable,dimension(:) :: HImap,HJmap
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,1)  !orbital alfa_up
  jsite=impIndex(jorb,2)  !orbital beta_dw
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  !
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)


     !EVALUATE [c^+_{up,iorb} + c_{dw,jorb}]|gs> --> A_{iorb,jorb}
     isz = getsz(isector)
     if(isz<Ns)then
        jsz   = isz+1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c^+_{up,iorb} + c_{dw,jorb}:',getsz(jsector),jdim
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=3)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE [c_{up,iorb} + c^+_{dw,jorb}]|gs>  --> A_{iorb,jorb}
     isz = getsz(isector)
     if(isz>-Ns)then
        jsz   = isz-1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A23,I3,I15)")'apply c_{up,iorb} + c^+_{dw,jorb}:',getsz(jsector),jdim
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=3)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif


     !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
     isz = getsz(isector)
     if(isz<Ns)then
        jsz   = isz+1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A23,I3,I15)")'apply c^+_{up,iorb} + xi*c_{dw,horb}:',getsz(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=0.d0
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
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
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
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,1,ichan=3)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
     isz = getsz(isector)
     if(isz>-Ns)then
        jsz   = isz-1
        jsector = getsector(jsz,1)
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.ED_MPI_ID==0)&
             write(LOGfile,"(A23,I3,I15)")'apply c_{up,iorb} - xi*c^+_{dw,jorb}:',getsz(jsector),jdim
        allocate(HJmap(jdim),cvinit(jdim))
        call build_sector(jsector,HJmap)
        cvinit=0.d0
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
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
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
        cnorm2=-xi**norm2
        call add_to_lanczos_gf_superc(cnorm2,state_e,nlanc,alfa_,beta_,-1,ichan=3)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_superc_mix_d




!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_superc(vnorm2,Ei,nlanc,alanc,blanc,isign,ichan)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,ichan
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
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
        auxGmats(ichan,i)=auxGmats(ichan,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        auxGreal(ichan,i)=auxGreal(ichan,i) + peso/(iw-isign*de)
     enddo
     ! GFpoles(ispin,ispin,iorb,jorb,itype,j)   = isign*de
     ! GFweights(ispin,ispin,iorb,jorb,itype,j) = peso
     auxGpoles(itype,ichan,j)   = isign*de
     auxGweights(itype,ichan,j) = peso
  enddo
end subroutine add_to_lanczos_gf_superc
