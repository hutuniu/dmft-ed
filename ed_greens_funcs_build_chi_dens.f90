!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Spin Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_chi_dens()
  integer :: iorb,jorb,ispin
  logical :: verbose
  verbose=.false.;if(ed_verbose<1)verbose=.true. 
  write(LOGfile,"(A)")"Get impurity dens Chi:"
  do iorb=1,Norb
     select case(ed_type)
     case default
        call lanc_ed_build_densChi_d(iorb,verbose)
     case ('c')
        call lanc_ed_build_densChi_c(iorb,verbose)
     end select
  enddo
  select case(ed_type)
  case default
     call lanc_ed_build_densChi_tot_d(verbose)
  case ('c')
     stop
     call lanc_ed_build_densChi_tot_c(verbose)
  end select
  denschi_tau = Denschi_tau/zeta_function
  denschi_w   = denschi_w/zeta_function
  denschi_iv  = denschi_iv/zeta_function
end subroutine build_chi_dens



!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Spin Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_chi_dens_mb()
  integer :: iorb,jorb,ispin
  logical :: verbose
  verbose=.false.;if(ed_verbose<1)verbose=.true. 
  write(LOGfile,"(A)")"Get impurity dens Chi:"
  do iorb=1,Norb
     select case(ed_type)
     case default
        call lanc_ed_build_densChi_d(iorb,verbose)
     case ('c')
        call lanc_ed_build_densChi_c(iorb,verbose)
     end select
  enddo

  do iorb=1,Norb
     do jorb=iorb+1,Norb
        call lanc_ed_build_densChi_mix_d(iorb,jorb,verbose)
     end do
  end do
  
  do iorb=1,Norb
     do jorb=iorb+1,Norb
        denschi_w(iorb,jorb,:) = 0.5d0*(denschi_w(iorb,jorb,:) &
             -(one+xi)*denschi_w(iorb,iorb,:) - (one+xi)*denschi_w(jorb,jorb,:))
     enddo
  enddo


  do iorb=1,Norb
     do jorb=1,Norb
!        do ispin=1,Nspin
        ispin=1
        call lanc_ed_build_Chi_mix_d(iorb,jorb,ispin)
!        end do
     end do
  end do


  select case(ed_type)
  case default
     call lanc_ed_build_densChi_tot_d(verbose)
  case ('c')
     call lanc_ed_build_densChi_tot_c(verbose)
  end select
  denschi_tau = Denschi_tau/zeta_function
  denschi_w   = denschi_w/zeta_function
  denschi_iv  = denschi_iv/zeta_function
end subroutine build_chi_dens_mb





!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
! single orbital: \chi = <S_a(\tau)S_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_d(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero,isign
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  complex(8)                       :: cnorm2
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Evaluating dens Chi_Orb"//reg(txtfy(iorb))//":"
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
     allocate(HImap(idim0),vvinit(idim0))
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")'Apply N:',getnup(isect0),getndw(isect0),idim0
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HImap)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)


     ! vvinit=0.d0
     ! do m=1,idim0                     !loop over |gs> components m
     !    i=HImap(m)
     !    call bdecomp(i,ib)
     !    sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
     !    vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     ! enddo
     ! deallocate(HImap)
     ! norm0=dot_product(vvinit,vvinit)
     ! vvinit=vvinit/sqrt(norm0)
     ! alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     ! call ed_buildH_d(isect0)
     ! call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     ! cnorm2=one*norm0
     ! isign=1
     ! call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)     
     ! isign=-1
     ! call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_d

subroutine lanc_ed_build_densChi_c(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero,isign
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  complex(8)                       :: cnorm2
  real(8),allocatable              :: alfa_(:),beta_(:)
  complex(8),allocatable           :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
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
     idim0      =  getdim(isect0)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0),vvinit(idim0))
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")'Apply N:',getnup(isect0),getndw(isect0),idim0
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_cvec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HImap)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isect0)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_c




subroutine lanc_ed_build_densChi_mix_d(iorb,jorb,iverbose)
  integer                          :: iorb,jorb,isite,isect0,izero,isign
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  complex(8)                       :: cnorm2
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable              :: cvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Evaluating dens Chi_Orb"//reg(txtfy(iorb))//":"
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     ! properties of the ground states
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     ! apply N_{iorb} + N_{jorb}



     allocate(HImap(idim0),vvinit(idim0),cvinit(idim0))
     call build_sector(isect0,HImap)

     !build the (N_iorb+N_jorb)|gs> state
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A)")'Apply N_{iorb} + N_{jorb}:'
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        vvinit(m) = vvinit(m) + sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,jorb)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,jorb)

     !build the (N_iorb-xi*N_jorb)|gs> state
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")'Apply N_iorb + xi*N_jorb:'
     cvinit=zero
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        cvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        cvinit(m) = cvinit(m) - xi*sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(cvinit,cvinit)
     cvinit=cvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dc)
     cnorm2=xi*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,jorb)

     !apply N_{iorb} + xi*N_{jorb}
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A)")'Apply N_iorb + xi*N_jorb:'
     cvinit=zero
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        cvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        cvinit(m) = cvinit(m) + xi*sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(cvinit,cvinit)
     cvinit=cvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dc)
     cnorm2=xi*norm0
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,iorb,jorb)


     deallocate(vvinit)
     deallocate(HImap)

     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_mix_d





!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the total Spin susceptibility \Chi_spin for a 
! single orbital: \chi = \sum_a <S_a(\tau)S_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_tot_d(iverbose)
  integer                          :: iorb,isite,isect0,izero,isign
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  complex(8)                       :: cnorm2
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
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
     allocate(HImap(idim0),vvinit(idim0))
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")'Apply N:',getnup(isect0),getndw(isect0),idim0
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0  
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = sum(dble(ib(1:Norb)))+sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = sgn*state_vec(m) 
     enddo
     deallocate(HImap)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,Norb+1,Norb+1)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,Norb+1,Norb+1)     
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_d

subroutine lanc_ed_build_densChi_tot_c(iverbose)
  integer                          :: iorb,isite,isect0,izero,isign
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  complex(8)                       :: cnorm2
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  complex(8),allocatable           :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
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
     idim0      =  getdim(isect0)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0),vvinit(idim0))
     if(iverbose_.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")'Apply N:',getnup(isect0),getndw(isect0),idim0
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = sum(dble(ib(1:Norb)))+sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = sgn*state_cvec(m) 
     enddo
     deallocate(HImap)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isect0)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,Norb+1,Norb+1)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,nitermax,alfa_,beta_,isign,Norb+1,Norb+1)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_c



subroutine lanc_ed_build_Chi_mix_d(iorb,jorb,ispin)
  real(8),allocatable              :: vvinit(:),vvinit_(:),tmp_vect
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,jorb,ispin,isite,isector,istate
  integer                          :: idim,jsector,jsector_
  integer                          :: jdim,jdim_
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  integer,allocatable,dimension(:) :: HImap,HJmap,HJmap_,HJmap_tmp

  write(*,*) 'LANC_CHI_MIX',iorb,jorb

  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
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

     !+- Apply c^dg_jorb c_iorb -+!
     ispin=1
     isite=impIndex(iorb,ispin)
     jsector = getCsector(ispin,isector)     
     if(jsector/=0)then
        jdim  = getdim(jsector)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
           end if
        enddo
     endif
     !     
     jsector_ = getCDGsector(ispin,jsector)
     if(jsector_/=0) then       
        jdim_  = getdim(jsector_)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(HJmap_(jdim_),vvinit_(jdim_))
        call build_sector(jsector_,HJmap_)
        vvinit_=0.d0
        do m=1,jdim
           i=HJmap(m)
           call bdecomp(i,ib)
           isite=impIndex(jorb,ispin)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap_,r)
              vvinit_(j) = sgn*vvinit(m)
           endif
        enddo
     end if
     
     !<DEBUG
     ! deallocate(vvinit); allocate(vvinit(idim))
     ! vvinit=0.d0
     ! do m=1,idim                     !loop over |gs> components m
     !    i=HImap(m)
     !    call bdecomp(i,ib)
     !    sgn = dble(ib(iorb))
     !    vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     ! enddo     
     ! do m=1,idim
     !    write(776,*) m,vvinit_(m),vvinit(m)
     ! end do
     !DEBUG>

     !
     deallocate(HJmap,vvinit,HJmap_)
     !
     ispin=2
     isite=impIndex(iorb,ispin)
     jsector = getCsector(ispin,isector)     
     if(jsector/=0)then
        jdim  = getdim(jsector)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
           end if
        enddo
     endif
     !     
     jsector_ = getCDGsector(ispin,jsector)
     if(jsector_/=0) then       
        jdim_  = getdim(jsector_)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(HJmap_(jdim_))
        call build_sector(jsector_,HJmap_)
        !vvinit_=0.d0
        do m=1,jdim
           i=HJmap(m)
           call bdecomp(i,ib)
           isite=impIndex(jorb,ispin)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap_,r)
              vvinit_(j) = vvinit_(j) + sgn*vvinit(m)
           endif
        enddo
     end if
     

     norm2=dot_product(vvinit_,vvinit_)
     vvinit_=vvinit_/sqrt(norm2)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
     call ed_buildH_d(jsector_)
     call lanczos_plain_tridiag_d(vvinit_,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
     cnorm2=one*norm2
     call add_to_lanczos_mixChi(cnorm2,state_e,nlanc,alfa_,beta_,+1,iorb,jorb)


     deallocate(HJmap,HJmap_)
     deallocate(vvinit,vvinit_)
        

     
     !+- Apply c^dg_iorb c_jorb -+!
     ispin=1
     isite=impIndex(jorb,ispin)
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
     endif
     !
     isite=impIndex(iorb,ispin)
     jsector_ = getCDGsector(ispin,jsector)
     if(jsector_/=0) then       
        jdim_  = getdim(jsector_)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(HJmap_(jdim_),vvinit_(jdim_))
        call build_sector(jsector_,HJmap_)
        vvinit_=0.d0
        do m=1,jdim
           i=HJmap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap_,r)
              vvinit_(j) = sgn*vvinit(m)
           endif
        enddo
     end if
     !
     deallocate(Hjmap,Hjmap_,vvinit)
     !
     ispin=2
     isite=impIndex(jorb,ispin)
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
     endif
     !
     isite=impIndex(iorb,ispin)
     jsector_ = getCDGsector(ispin,jsector)
     if(jsector_/=0) then       
        jdim_  = getdim(jsector_)
        !if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(HJmap_(jdim_))
        call build_sector(jsector_,HJmap_)
        !vvinit_=0.d0
        do m=1,jdim
           i=HJmap(m)
           call bdecomp(i,ib)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJmap_,r)
              vvinit_(j) = vvinit_(j) +  sgn*vvinit(m)
           endif
        enddo
     end if


     !<DEBUG
     ! deallocate(vvinit); allocate(vvinit(idim))
     ! vvinit=0.d0
     ! do m=1,idim                     !loop over |gs> components m
     !    i=HImap(m)
     !    call bdecomp(i,ib)
     !    sgn = dble(ib(iorb+Ns))+dble(ib(iorb))
     !    vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     ! enddo     
     ! do m=1,idim
     !    write(777,*) m,vvinit_(m),vvinit(m)
     ! end do
     ! stop
     !DEBUG>



     norm2=dot_product(vvinit_,vvinit_)
     vvinit_=vvinit_/sqrt(norm2)
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
     call ed_buildH_d(jsector_)
     call lanczos_plain_tridiag_d(vvinit_,alfa_,beta_,nlanc,lanc_spHtimesV_dd)
     cnorm2=one*norm2
     call add_to_lanczos_mixChi(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb)
     deallocate(vvinit,vvinit_)
     deallocate(HJmap,HJmap_)     
     !
     nullify(state_vec)
     deallocate(HImap)
     !
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_Chi_mix_d








!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_densChi(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb)
  complex(8) :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
  real(8)                                    :: Ei,Ej,Egs,de
  integer                                    :: nlanc,isign
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: iorb,jorb
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs    = state_list%emin
  pesoF  = vnorm2/zeta_function 
  pesoBZ = 1d0
  if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
  !
  diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
  forall(i=1:Nlanc)Z(i,i)=1.d0
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  do j=1,nlanc
     Ej     = diag(j)
     dE     = Ej-Ei
     pesoAB = Z(1,j)*Z(1,j)
     peso   = pesoF*pesoAB*pesoBZ

     !Only real freq calculations
     if(iorb.le.Norb.and.jorb.le.Norb) then
        do i=1,Lreal
           !denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
           denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + isign*peso/(dcmplx(wr(i),eps)-isign*dE)           
        enddo
     else
        do i=1,Lreal
           denschi_tot_w(i)=denschi_tot_w(i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
           !denschi_tot_w(i)=denschi_tot_w(i) + isign*peso/(dcmplx(wr(i),eps)-isign*dE)
        end do
     end if



     !Matsubara:
     !treat separately the first bosonic Matsubara freq.
     !    if(iorb.le.Norb.and.jorb.le.Norb) then
     !       if(beta*dE < 1)then
     !          denschi_iv(iorb,jorb,0)=denschi_iv(iorb,jorb,0) + peso*2*beta
     !       else
     !          denschi_iv(iorb,jorb,0)=denschi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE !there is a factor 2 we do not know
     !       endif
     !       do i=1,Lmats
     !          denschi_iv(iorb,jorb,i)=denschi_iv(iorb,jorb,i) + peso*2*dE/(vm(i)**2+dE**2)
     !       enddo
     !       !Imag. time:
     !       do i=0,Ltau
     !          denschi_tau(iorb,jorb,i)=denschi_tau(iorb,jorb,i) + peso*(exp(-tau(i)*de)+exp(-(beta-tau(i))*de))
     !       enddo
     !       !Real freq.: misses a factor 2
     !       do i=1,Lreal
     !          !denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
     !          !denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + peso*2.d0*dE/(dcmplx(wr(i),eps)**2.d0-dE**2.d0)           
     !          denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + isign*peso/(dcmplx(wr(i),eps)-isign*dE)           
     !       enddo
     !    else
     !       if(beta*dE < 1)then
     !          denschi_tot_iv(0)=denschi_tot_iv(0) + peso*2*beta
     !       else
     !          denschi_tot_iv(0)=denschi_tot_iv(0) + peso*2*(1d0-exp(-beta*dE))/dE !there is a factor 2 we do not know
     !       endif
     !       do i=1,Lmats
     !          denschi_tot_iv(i)=denschi_tot_iv(i) + peso*2*dE/(vm(i)**2+dE**2)
     !       enddo
     !       !Imag. time:
     !       do i=0,Ltau
     !          denschi_tot_tau(i)=denschi_tot_tau(i) + peso*(exp(-tau(i)*de)+exp(-(beta-tau(i))*de))
     !       enddo
     !       !Real freq.: misses a factor 2
     !       do i=1,Lreal
     !          denschi_tot_w(i)=denschi_tot_w(i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
     !       enddo
     !    end if
  enddo

end subroutine add_to_lanczos_densChi




subroutine add_to_lanczos_mixChi(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb)
  complex(8) :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
  real(8)                                    :: Ei,Ej,Egs,de
  integer                                    :: nlanc,isign
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: iorb,jorb
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs    = state_list%emin
  pesoF  = vnorm2/zeta_function 
  pesoBZ = 1d0
  if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
  !
  diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
  forall(i=1:Nlanc)Z(i,i)=1.d0
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  do j=1,nlanc
     Ej     = diag(j)
     dE     = Ej-Ei
     pesoAB = Z(1,j)*Z(1,j)
     peso   = pesoF*pesoAB*pesoBZ
     do i=1,Lreal
        !denschi_w(iorb,jorb,i)=denschi_w(iorb,jorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
        denschi_mix_w(iorb,jorb,i)=denschi_mix_w(iorb,jorb,i) + isign*peso/(dcmplx(wr(i),eps)-isign*dE)           
     enddo
     if(iorb.ne.jorb) then
        write(777,*) j,dE
     else
        write(776,*) j,dE
     end if
  enddo
  if(iorb.ne.jorb) then
     write(777,*) 
  else
     write(776,*) 
  end if


end subroutine add_to_lanczos_mixChi

