!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Pair Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_chi_pair()
  integer :: iorb
  write(LOGfile,"(A)")"Get impurity pair Chi:"
  do iorb=1,Norb
     write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
     call lanc_ed_build_pairChi_c(iorb)
  enddo
  pairChi_tau = PairChi_tau/zeta_function
  pairChi_w   = pairChi_w/zeta_function
  pairChi_iv  = pairChi_iv/zeta_function
end subroutine build_chi_pair




!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Pair susceptibility \Chi_pair for a 
! single orbital: \chi = <Phi_a(\tau)Phi_a(0)>
!+------------------------------------------------------------------+
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
  call start_timer
  !
  do izero=1,state_list%size
     isector    =  es_return_sector(state_list,izero)
     idim       =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     call build_sector(isector,HI)
     !
     !Build the C_{iorb,up}C_{iorb,dw}|eigvec> 
     if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply C_{iorb,up}C_{iorb,dw}:',getsz(isector)
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
     !
     call setup_Hv_sector(isector)
     if(ed_sparse_H)call ed_buildH_c()
     !
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
     endif
     isign=-1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,alfa_,beta_,isign,iorb)
     !
     call delete_Hv_sector()
     !
     if(spH0%status)call sp_delete_matrix(spH0)
     deallocate(vvinit,alfa_,beta_)
     !
     !Build the CDG_{iorb,dw}CDG_{iorb,up}|eigvec> 
     if(ed_verbose==3)write(LOGfile,"(A,2I3)")'Apply CDG_{iorb,dw}CDG_{iorb,up}:',getsz(isector)
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
     !
     call setup_Hv_sector(isector)
     if(ed_sparse_H)call ed_buildH_c()
     !
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
     endif
     isign=1 !<== ACTHUNG!!!! check this is the correct value of isign
     call add_to_lanczos_pairChi(norm0,state_e,alfa_,beta_,isign,iorb)
     !
     call delete_Hv_sector()
     !
     if(spH0%status)call sp_delete_matrix(spH0)
     deallocate(vvinit,alfa_,beta_)
     deallocate(HI%map)
     nullify(state_cvec)
  enddo
  call stop_timer(LOGfile)
end subroutine lanc_ed_build_pairChi_c










!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_pairChi(vnorm,Ei,alanc,blanc,isign,iorb)
  real(8)                                    :: vnorm,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
  integer                                    :: nlanc
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  integer                                    :: isign,iorb
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
  !
  pesoF  = vnorm**2/zeta_function 
  pesoBZ = 1d0
  if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
  !
  diag             = 0.d0
  subdiag          = 0.d0
  Z                = eye(Nlanc)
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
end subroutine add_to_lanczos_pairChi







