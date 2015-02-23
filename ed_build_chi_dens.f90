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
     call lanc_ed_build_densChi_tot_c(verbose)
  end select
  denschi_tau = Denschi_tau/zeta_function
  denschi_w   = denschi_w/zeta_function
  denschi_iv  = denschi_iv/zeta_function

end subroutine build_chi_dens




!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
! single orbital: \chi = <S_a(\tau)S_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_d(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
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
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     call add_to_lanczos_densChi(norm0,state_e,nitermax,alfa_,beta_,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_d

subroutine lanc_ed_build_densChi_c(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
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
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_cvec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HImap)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isect0)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     call add_to_lanczos_densChi(norm0,state_e,nitermax,alfa_,beta_,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_c






!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the total Spin susceptibility \Chi_spin for a 
! single orbital: \chi = \sum_a <S_a(\tau)S_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_tot_d(iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
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
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     call add_to_lanczos_densChi(norm0,state_e,nitermax,alfa_,beta_,Norb+1)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_d

subroutine lanc_ed_build_densChi_tot_c(iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
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
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isect0)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     call add_to_lanczos_densChi(norm0,state_e,nitermax,alfa_,beta_,Norb+1)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_c





!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_densChi(vnorm,Ei,nlanc,alanc,blanc,iorb)
  real(8)                                    :: vnorm,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
  integer                                    :: nlanc
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs    = state_list%emin
  pesoF  = vnorm**2/zeta_function 
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
     !Matsubara:
     !treat separately the first bosonic Matsubara freq.
     if(beta*dE < 1)then
        denschi_iv(iorb,0)=denschi_iv(iorb,0) + peso*2*beta
     else
        denschi_iv(iorb,0)=denschi_iv(iorb,0) + peso*2*(1d0-exp(-beta*dE))/dE !there is a factor 2 we do not know
     endif
     do i=1,Lmats
        denschi_iv(iorb,i)=denschi_iv(iorb,i) + peso*2*dE/(vm(i)**2+dE**2)
     enddo
     !Imag. time:
     do i=0,Ltau
        denschi_tau(iorb,i)=denschi_tau(iorb,i) + peso*(exp(-tau(i)*de)+exp(-(beta-tau(i))*de))
     enddo
     !Real freq.: misses a factor 2
     do i=1,Lreal
        denschi_w(iorb,i)=denschi_w(iorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
     enddo
  enddo
end subroutine add_to_lanczos_densChi
