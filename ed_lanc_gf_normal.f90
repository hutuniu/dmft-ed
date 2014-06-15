!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine lanc_ed_getgf_normal()
  integer :: izero,iorb,jorb,ispin,i
  integer :: isect0,numstates
  real(8) :: norm0
  complex(8) :: Greal(Norb,Lreal)
  call allocate_grids

  if(.not.allocated(impGmats))allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
  if(.not.allocated(impGreal))allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
  impGmats=zero
  impGreal=zero

  do ispin=1,Nspin
     do iorb=1,Norb
        if(ed_verbose<3)write(LOGfile,"(A)")"Get G_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))
        select case(ed_type)
        case default
           call lanc_ed_buildgf_d(iorb,ispin,.false.)
        case ('c')
           call lanc_ed_buildgf_c(iorb,ispin,.false.)
        end select
     enddo
  enddo

  if(bath_type=='hybrid')then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if(ed_verbose<3)write(LOGfile,"(A)")"Get G_l"//&
                   reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
              select case(ed_type)
              case default
                 call lanc_ed_buildgf_mix_d(iorb,jorb,ispin,.false.)
              case ('c')
                 call lanc_ed_buildgf_mix_c(iorb,jorb,ispin,.false.)                    
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
  !Print impurity functions:
  call print_imp_gf

  !<DEBUG
  call cv_ed_buildgf_d(1,1)
  !>DEBUG

  deallocate(wm,wr,tau,vm)
  deallocate(impGmats,impGreal)
end subroutine lanc_ed_getgf_normal




!<DEBUG
!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine cv_ed_buildgf_d(iorb,ispin)
  complex(8)                       :: Greal(Lreal)
  complex(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: xvec(:)
  integer                          :: iorb,ispin,isite,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0,err,Egs
  complex(8)                       :: cnorm2,omega
  integer                          :: Nitermax,Nlanc,iter,iomega
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap,HJmap
  complex(8),dimension(:),allocatable :: Aw
  allocate(Aw(Lreal));Aw=zero
  !
  iverbose_=.false.
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !   
  if(ed_verbose<2)call start_progress
  do izero=1,numstates
     if(ed_verbose<1)call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)

     !ADD ONE PARTICLE:
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'add particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=zero
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        !norm2=dot_product(vvinit,vvinit)
        !vvinit=vvinit/sqrt(norm2)
        deallocate(HJmap)
        call ed_buildH_d(jsect0)
        Egs = state_list%emin
        allocate(xvec(jdim0))
        do iomega=1,Lreal
           omega=dcmplx(wm(iomega),eps)
           xvec=one
           call bcg_solve(atimes_,asolve_,vvinit,xvec,1,1.d-9,100,iter,err,.true.)
           Aw(iomega) = -dimag(dot_product(vvinit,xvec))/pi
        enddo
        deallocate(xvec)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'del particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=zero
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJmap)
        !norm2=dot_product(vvinit,vvinit)
        !vvinit=vvinit/sqrt(norm2)
        call ed_buildH_d(jsect0)
        Egs = state_list%emin        
        allocate(xvec(jdim0))
        do iomega=1,Lreal
           omega=dcmplx(wm(iomega),eps)
           xvec=one
           call bcg_solve(atimes_,asolve_,vvinit,xvec,1,1.d-9,100,iter,err,.false.)
           Aw(iomega) = Aw(iomega)-dimag(dot_product(vvinit,xvec))/pi
        enddo
        deallocate(xvec)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<2)call stop_progress


  do i=1,Lreal
     write(876,*)wr(i),dreal(Aw(i)),dimag(Aw(i))
  enddo


contains

  subroutine atimes_(vin,vout)
    complex(8),dimension(:)                  :: vin
    complex(8),dimension(size(vin))          :: vout
    if(size(vin)/=jdim0.OR.size(vin)/=jdim0)stop "atimes: error in dimensions"
    call sp_matrix_vector_product_dc(spH0,jdim0,vin,vout)
    vout = (omega+egs)*vin - vout
  end subroutine atimes_


  subroutine asolve_(vin,vout)
    complex(8),dimension(:)                  :: vin
    complex(8),dimension(size(vin))          :: vout
    real(8),dimension(size(vin))          :: diag
    if(size(vin)/=jdim0.OR.size(vin)/=jdim0)stop "asolve: error in dimensions"
    call sp_get_diagonal(spH0,diag)
    do i=1,size(vin)
       if((omega+egs-diag(i))/=0.d0)vout(i)=vin(i)/(omega+egs-diag(i))
    enddo
  end subroutine asolve_



end subroutine cv_ed_buildgf_d

!>DEBUG






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_ed_buildgf_d(iorb,ispin,iverbose)
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,ispin,isite,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0
  integer                          :: ib(Ntot)
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
  isite=impIndex(iorb,ispin)
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !   
  if(ed_verbose<2)call start_progress
  do izero=1,numstates
     if(ed_verbose<1)call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)

     !ADD ONE PARTICLE:
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'add particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        call lanczos_plain_set_htimesv_d(lanc_spHtimesV_dd)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'del particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_d(lanc_spHtimesV_dd)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<2)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildgf_d


!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_ed_buildgf_c(iorb,ispin,iverbose)
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: iorb,ispin,isite,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0
  integer                          :: ib(Ntot)
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
  isite=impIndex(iorb,ispin)
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !   
  if(ed_verbose<2)call start_progress
  do izero=1,numstates
     if(ed_verbose<1)call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)

     !ADD ONE PARTICLE:
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'add particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim0                     !loop over |gs> components m
           i=HImap(m)                    !map m to Hilbert space state i
           call bdecomp(i,ib)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)      !map r back to  jsect0
              vvinit(j) = sgn*state_cvec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0  = getdim(jsect0)
        if(iverbose_)write(LOGfile,"(A,2I3,I15)")'del particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<2)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildgf_c






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_ed_buildgf_mix_d(iorb,jorb,ispin,iverbose)
  integer                          :: iorb,jorb,ispin,isite,jsite,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !   
  if(ed_verbose<2)call start_progress
  do izero=1,numstates
     if(ed_verbose<1)call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)write(*,"(A,2I3,I15)")'add particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_d(lanc_spHtimesV_dd)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0   = getdim(jsect0)
        if(iverbose_)write(*,"(A,2I3,I15)")'del particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_d(lanc_spHtimesV_dd)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)write(*,"(A,2I3,I15)")'add particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),cvinit(jdim0))
        call build_sector(jsect0,HJmap)
        cvinit=zero
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_dc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=-xi*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0   = getdim(jsect0)
        if(iverbose_)write(*,"(A,2I3,I15)")'del particle:',&
             getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),cvinit(jdim0))
        call build_sector(jsect0,HJmap)
        cvinit=zero
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_dc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_d(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=-xi*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<2)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildgf_mix_d



!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_ed_buildgf_mix_c(iorb,jorb,ispin,iverbose)
  integer                          :: iorb,jorb,ispin,isite,jsite,isect0,izero
  integer                          :: idim0,jsect0
  integer                          :: jdim0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  !
  nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !   
  if(ed_verbose<2)call start_progress
  do izero=1,numstates
     if(ed_verbose<1)call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim0  = getdim(isect0)
     allocate(HImap(idim0))
     call build_sector(isect0,HImap)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)&
             write(*,"(A,2I3,I15)")'add particle:',getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0   = getdim(jsect0)
        if(iverbose_)&
             write(*,"(A,2I3,I15)")'del particle:',getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),vvinit(jdim0))
        call build_sector(jsect0,HJmap)
        vvinit=0.d0
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim0
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
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=one*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsect0 = getCDGsector(ispin,isect0)
     if(jsect0/=0)then 
        jdim0  = getdim(jsect0)
        if(iverbose_)&
             write(*,"(A,2I3,I15)")'add particle:',getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),cvinit(jdim0))
        call build_sector(jsect0,HJmap)
        cvinit=zero
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=-xi*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsect0 = getCsector(ispin,isect0)
     if(jsect0/=0)then
        jdim0   = getdim(jsect0)
        if(iverbose_)&
             write(*,"(A,2I3,I15)")'del particle:',getnup(jsect0),getndw(jsect0),jdim0
        allocate(HJmap(jdim0),cvinit(jdim0))
        call build_sector(jsect0,HJmap)
        cvinit=zero
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim0
           i=HImap(m)
           call bdecomp(i,ib)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJmap,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJmap)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
        !call setup_Hv_sector(jsect0)
        call ed_buildH_c(jsect0)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
        call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
        !call delete_Hv_sector()
        call lanczos_plain_delete_htimesv
        cnorm2=-xi*norm2
        call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<2)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildgf_mix_c









!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb,ispin)
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb,jorb,ispin
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
        impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
     enddo
  enddo
end subroutine add_to_lanczos_gf
