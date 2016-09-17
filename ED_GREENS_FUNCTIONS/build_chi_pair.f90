!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Pair Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_chi_pair()
  integer :: iorb
  write(LOGfile,"(A)")"Get impurity pair Chi:"
  do iorb=1,Norb
     if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
     select case(ed_type)
     case default
        call lanc_ed_build_pairChi_d(iorb)
     case ('c')
        call lanc_ed_build_pairChi_c(iorb)
     end select
  enddo
  pairChi_tau = PairChi_tau/zeta_function
  pairChi_w   = pairChi_w/zeta_function
  pairChi_iv  = pairChi_iv/zeta_function
end subroutine build_chi_pair




!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Pair susceptibility \Chi_pair for a 
! single orbital: \chi = <Phi_a(\tau)Phi_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_pairChi_d(iorb)
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,i1,i2,j,r
  real(8)                          :: norm0,sgn,sgn1,sgn2
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
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do izero=1,numstates
     isector    =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     
     call build_sector(isector,HI)
     !Build the C_{iorb,up}C_{iorb,dw}|eigvec> 
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply C_{iorb,up}C_{iorb,dw}:',getsz(isector)
     allocate(vvinit(idim))
     vvinit=0.d0
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        if(ib(iorb+Ns)==0.OR.ib(iorb)==0)cycle
        call c(iorb+Ns,i,i1,sgn1)
        call c(iorb,i1,i2,sgn2)
        j = binary_search(HI%map,i2)
        vvinit(j) = sgn1*sgn2*state_vec(m)
     enddo
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isector)
     call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nitermax)
     isign=-1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     deallocate(vvinit)
     !Build the CDG_{iorb,dw}CDG_{iorb,up}|eigvec> 
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply CDG_{iorb,dw}CDG_{iorb,up}:',getsz(isector)
     allocate(vvinit(idim))
     vvinit=0.d0
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        if(ib(iorb+Ns)==1.OR.ib(iorb)==1)cycle
        call cdg(iorb,i,i1,sgn1)
        call cdg(iorb+Ns,i1,i2,sgn2)
        j = binary_search(HI%map,i2)
        vvinit(j) = sgn1*sgn2*state_vec(m)
     enddo
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isector)
     call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nitermax)
     isign=1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     if(spH0%status)call sp_delete_matrix(spH0)
     deallocate(vvinit)
     deallocate(HI%map)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_pairChi_d

subroutine lanc_ed_build_pairChi_c(iorb)
  integer                          :: iorb,isite,isector,izero
  integer                          :: numstates
  integer                          :: nlanc,idim
  integer                          :: iup0,idw0,isign
  integer                          :: ib(Nlevels)
  integer                          :: m,i,i1,i2,j,r
  real(8)                          :: norm0,sgn,sgn1,sgn2
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
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do izero=1,numstates
     isector    =  es_return_sector(state_list,izero)
     idim       =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)

     call build_sector(isector,HI)
     !Build the C_{iorb,up}C_{iorb,dw}|eigvec> 
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply C_{iorb,up}C_{iorb,dw}:',getsz(isector)
     allocate(vvinit(idim))
     vvinit=0.d0
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        if(ib(iorb+Ns)==0.OR.ib(iorb)==0)cycle
        call c(iorb+Ns,i,i1,sgn1)
        call c(iorb,i1,i2,sgn2)
        j = binary_search(HI%map,i2)
        vvinit(j) = sgn1*sgn2*state_cvec(m)
     enddo
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isector)
     call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nitermax)
     isign=-1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     deallocate(vvinit)
     !Build the CDG_{iorb,dw}CDG_{iorb,up}|eigvec> 
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply CDG_{iorb,dw}CDG_{iorb,up}:',getsz(isector)
     allocate(vvinit(idim))
     vvinit=0.d0
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        if(ib(iorb+Ns)==1.OR.ib(iorb)==1)cycle
        call cdg(iorb,i,i1,sgn1)
        call cdg(iorb+Ns,i1,i2,sgn2)
        j = binary_search(HI%map,i2)
        vvinit(j) = sgn1*sgn2*state_cvec(m)
     enddo
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_c(isector)
     call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nitermax)
     isign=1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,nitermax,alfa_,beta_,isign,iorb)
     if(spH0%status)call sp_delete_matrix(spH0)
     deallocate(vvinit)
     deallocate(HI%map)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_pairChi_c









!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_pairChi(vnorm,Ei,nlanc,alanc,blanc,isign,iorb)
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
           pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*beta
        else
           pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           pairChi_iv(iorb,i)=pairChi_iv(iorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
        enddo
        do i=0,Ltau
           pairChi_tau(iorb,i)=pairChi_tau(iorb,i) + peso*exp(-tau(i)*de)
        enddo
        do i=1,Lreal
           pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
        enddo
     enddo
  case (-1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*beta
        else
           pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           pairChi_iv(iorb,i)=pairChi_iv(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
        enddo
        do i=0,Ltau
           pairChi_tau(iorb,i)=pairChi_tau(iorb,i) + peso*exp(-(beta-tau(i))*dE)
        enddo
        do i=1,Lreal
           pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
        enddo
     enddo
  case default
     stop "add_to_lanczos_pairChi: isign not in {-1,1}"
  end select


  ! do j=1,nlanc
  !    Ej     = diag(j)
  !    dE     = Ej-Ei
  !    pesoAB = Z(1,j)*Z(1,j)
  !    peso   = pesoF*pesoAB*pesoBZ
  !    !Matsubara:
  !    !treat separately the first bosonic Matsubara freq.
  !    if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
  !       pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*2*beta
  !    else
  !       pairChi_iv(iorb,0)=pairChi_iv(iorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
  !    endif
  !    do i=1,Lmats
  !       pairChi_iv(iorb,i)=pairChi_iv(iorb,i) + peso*2*dE/(vm(i)**2+dE**2)
  !    enddo
  !    !Imag. time:
  !    do i=0,Ltau
  !       pairChi_tau(iorb,i)=pairChi_tau(iorb,i) + peso*(exp(-tau(i)*de)+exp(-(beta-tau(i))*de))
  !    enddo
  !    !Real freq.: misses a factor 2
  !    ![ (exp(-beta*DeltaE)-1)/(w+xi*eta - DeltaE) + (1-exp(-beta*DeltaE))/(w+xi*eta + DeltaE) ]
  !    ! The second term is missing: does it contribute only to w<0? I guess so...
  !    do i=1,Lreal
  !       pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(exp(-beta*de)-1.d0)/(dcmplx(wr(i),eps)-de)
  !       !pairChi_w(iorb,i)=pairChi_w(iorb,i) + peso*(1.d0-exp(-beta*de))/(dcmplx(wr(i),eps)+de)
  !    enddo
  ! enddo



end subroutine add_to_lanczos_pairChi
