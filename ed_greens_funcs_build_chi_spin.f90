!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Spin Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_chi_spin()
  integer :: iorb
  write(LOGfile,"(A)")"Get impurity spin Chi:"
  do iorb=1,Norb
     if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
     select case(ed_type)
     case default
        call lanc_ed_build_spinChi_d(iorb)
     case ('c')
        call lanc_ed_build_spinChi_c(iorb)
     end select
  enddo
  if(Norb>1)then
     if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Get Chi_spin_tot"
     select case(ed_type)
     case default
        call lanc_ed_build_spinChi_tot_d()
     case ('c')
        call lanc_ed_build_spinChi_tot_c()
     end select
  endif
  spinChi_tau = SpinChi_tau/zeta_function
  spinChi_w   = spinChi_w/zeta_function
  spinChi_iv  = spinChi_iv/zeta_function
end subroutine build_chi_spin




!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
! single orbital: \chi = <S_a(\tau)S_a(0)>
! note: as S_a is hermitian particle and holes contributions (isign=1,-1)
! are identical so work out only one lanczos tridiag. work out the 
! reduction for both values of isign in the same call.
!+------------------------------------------------------------------+
subroutine lanc_ed_build_spinChi_d(iorb)
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  integer                          :: Nitermax
  type(sector_map) :: HI    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     isector     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))-dble(ib(iorb+Ns))
        vvinit(m) = 0.5d0*sgn*state_vec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HI%map)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isector)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     !particles
     isign=1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     !holes
     isign=-1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_spinChi_d

subroutine lanc_ed_build_spinChi_c(iorb)
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  complex(8),allocatable           :: vvinit(:)
  integer                          :: Nitermax
  type(sector_map) :: HI    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     isector     =  es_return_sector(state_list,izero)
     idim      =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))-dble(ib(iorb+Ns))
        vvinit(m) = 0.5d0*sgn*state_cvec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HI%map)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isector)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     !particles
     isign=1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     !holes
     isign=-1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_spinChi_c






!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the total Spin susceptibility \Chi_spin for a 
! single orbital: \chi =  <[\sum_a S_a(\tau)][\sum_a S_a(0)]>
! note: as S_a is hermitian particle and holes contributions (isign=1,-1)
! are identical so work out only one lanczos tridiag. work out the 
! reduction for both values of isign in the same call.
!+------------------------------------------------------------------+
subroutine lanc_ed_build_spinChi_tot_d()
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  integer                          :: Nitermax
  type(sector_map) :: HI    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     isector     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim  
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = sum(dble(ib(1:Norb)))-sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = 0.5d0*sgn*state_vec(m) 
     enddo
     deallocate(HI%map)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isector)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_dd)
     !particles
     isign=1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,Norb+1)
     !holes
     isign=-1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,Norb+1)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_spinChi_tot_d

subroutine lanc_ed_build_spinChi_tot_c()
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  complex(8),allocatable           :: vvinit(:)
  integer                          :: Nitermax
  type(sector_map) :: HI    !map of the Sector S to Hilbert space H
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=state_list%size
  !
  if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer
  do izero=1,numstates
     isector     =  es_return_sector(state_list,izero)
     idim       =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.ED_MPI_ID==0)write(LOGfile,"(A,2I3)")'Apply Sz:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = sum(dble(ib(1:Norb)))-sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = 0.5d0*sgn*state_cvec(m) 
     enddo
     deallocate(HI%map)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isector)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax,lanc_spHtimesV_cc)
     !particles
     isign=1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,Norb+1)
     !holes
     isign=-1
     call add_to_lanczos_spinChi(norm0,state_e,nitermax,alfa_,beta_,isign,Norb+1)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_spinChi_tot_c






!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_spinChi(vnorm,Ei,nlanc,alanc,blanc,isign,iorb)
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
  select case(isign)
  case (1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*beta
        else
           spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           spinChi_iv(iorb,i)=spinChi_iv(iorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
        enddo
        do i=0,Ltau
           spinChi_tau(iorb,i)=spinChi_tau(iorb,i) + peso*exp(-tau(i)*de)
        enddo
        do i=1,Lreal
           spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
        enddo
     enddo
  case (-1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*beta
        else
           spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           spinChi_iv(iorb,i)=spinChi_iv(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
        enddo
        do i=0,Ltau
           spinChi_tau(iorb,i)=spinChi_tau(iorb,i) + peso*exp(-(beta-tau(i))*dE)
        enddo
        do i=1,Lreal
           spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
        enddo
     enddo
  case default
     stop "add_to_lanczos_spinChi: isign not in {-1,1}"
  end select


  ! do j=1,nlanc
  !    Ej     = diag(j)
  !    dE     = Ej-Ei
  !    pesoAB = Z(1,j)*Z(1,j)
  !    peso   = pesoF*pesoAB*pesoBZ
  ! !Matsubara: treat separately the first bosonic Matsubara freq.
  ! if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
  !    spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*2*beta
  ! else
  !    spinChi_iv(iorb,0)=spinChi_iv(iorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
  ! endif
  ! do i=1,Lmats
  !    spinChi_iv(iorb,i)=spinChi_iv(iorb,i) + peso*2*dE/(vm(i)**2+dE**2)
  ! enddo
  ! !Imag. time:
  ! do i=0,Ltau
  !    spinChi_tau(iorb,i)=spinChi_tau(iorb,i) + peso*(exp(-tau(i)*de)+exp(-(beta-tau(i))*de))
  ! enddo
  ! !Real freq.: misses a factor 2
  ! ![ (exp(-beta*DeltaE)-1)/(w+xi*eta - DeltaE) + (1-exp(-beta*DeltaE))/(w+xi*eta + DeltaE) ]
  ! ! The second term is missing: does it contribute only to w<0? I guess so...
  ! do i=1,Lreal
  !    spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
  !    !spinChi_w(iorb,i)=spinChi_w(iorb,i) + peso*(1.d0-exp(-beta*de))/(dcmplx(wr(i),eps)+de)
  ! enddo
  ! enddo
end subroutine add_to_lanczos_spinChi
