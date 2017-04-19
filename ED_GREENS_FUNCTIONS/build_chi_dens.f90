!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Charge-Charge Susceptibility <n_a(tau)n_b(0)>
!+------------------------------------------------------------------+
subroutine build_chi_dens()
  integer :: iorb,jorb
  write(LOGfile,"(A)")"Get impurity dens Chi:"
  do iorb=1,Norb
     if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get Chi_dens_diag_l"//reg(txtfy(iorb))
     select case(ed_type)
     case default
        call lanc_ed_build_densChi_diag_d(iorb)
     case ('c')
        call lanc_ed_build_densChi_diag_c(iorb)
     end select
  enddo
  !
  !
  if(Norb>1)then
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
           select case(ed_type)
           case default
              call lanc_ed_build_densChi_offdiag_d(iorb,jorb)
           case ('c')
              stop "ed_greens_funcs_build_chi_dens: lanc_ed_build_densChi_offdiac_C is not implemented. Sorry."
           end select
        end do
     end do
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           denschi_w(iorb,jorb,:) = 0.5d0*( denschi_w(iorb,jorb,:) - (one+xi)*denschi_w(iorb,iorb,:) - (one+xi)*denschi_w(jorb,jorb,:))
        enddo
     enddo
     !
     do iorb=1,Norb
        do jorb=1,Norb
           if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
           select case(ed_type)
           case default
              call lanc_ed_build_densChi_mix_d(iorb,jorb)
           case ('c')
              stop "ed_greens_funcs_build_chi_dens: lanc_ed_build_densChi_mix_C is not implemented. Sorry."
           end select
        end do
     end do
     !
     if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get Chi_dens_tot"
     select case(ed_type)
     case default
        call lanc_ed_build_densChi_tot_d()
     case ('c')
        call lanc_ed_build_densChi_tot_c()     
     end select
  endif
  !
  denschi_tau = Denschi_tau/zeta_function
  denschi_w   = denschi_w/zeta_function
  denschi_iv  = denschi_iv/zeta_function
  !
end subroutine build_chi_dens










!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Charge-Charge susceptibility \Chi_dens for  
! the orbital diagonal case: \chi_dens_aa = <N_a(\tau)N_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_diag_d(iorb)
  integer             :: iorb,isite,isector,izero
  integer             :: numstates
  integer             :: nlanc,idim
  integer             :: iup0,idw0,isign
  integer             :: ib(Nlevels)
  integer             :: m,i,j,r
  real(8)             :: norm0,sgn
  complex(8)          :: cnorm2
  real(8),allocatable :: alfa_(:),beta_(:)
  real(8),allocatable :: vvinit(:)
  integer             :: Nitermax
  type(sector_map)    :: HI    !map of the Sector S to Hilbert space H
  !
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do izero=1,state_list%size
     isector     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HI%map)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     call ed_buildH_d(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
     endif
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_diag_d

subroutine lanc_ed_build_densChi_diag_c(iorb)
  integer                :: iorb,isite,isector,izero
  integer                :: numstates
  integer                :: nlanc,idim
  integer                :: iup0,idw0,isign
  integer                :: ib(Nlevels)
  integer                :: m,i,j,r
  real(8)                :: norm0,sgn
  complex(8)             :: cnorm2
  real(8),allocatable    :: alfa_(:),beta_(:)
  complex(8),allocatable :: vvinit(:)
  integer                :: Nitermax
  type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
  !
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do izero=1,state_list%size
     isector     =  es_return_sector(state_list,izero)
     idim      =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_cvec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HI%map)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     call ed_buildH_c(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
     endif
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_diag_c



!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Charge-Charge susceptibility \Chi_dens for
! the orbital off-diagonal case: \chi_dens_ab = <N_a(\tau)N_b(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_offdiag_d(iorb,jorb)
  integer                :: iorb,jorb,isite,isector,izero,isign
  integer                :: numstates
  integer                :: nlanc,idim
  integer                :: iup0,idw0
  integer                :: ib(Nlevels)
  integer                :: m,i,j,r
  complex(8)             :: cnorm2
  real(8)                :: norm0,sgn
  real(8),allocatable    :: alfa_(:),beta_(:)
  real(8),allocatable    :: vvinit(:)
  complex(8),allocatable :: cvinit(:)
  integer                :: Nitermax
  type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do izero=1,state_list%size
     ! properties of the ground states
     isector     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim),cvinit(idim))
     call build_sector(isector,HI)
     !
     !build the (N_iorb+N_jorb)|gs> state
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A)")'Apply N_iorb + N_jorb:'
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        vvinit(m) = vvinit(m) + sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     call ed_buildH_d(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
     endif
     cnorm2=one*norm0
     !particle and holes excitations all at once
     isign=1                    !<---
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
     isign=-1                   !<---
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
     !
     !build the (N_iorb - xi*N_jorb)|gs> state
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A)")'Apply N_iorb + xi*N_jorb:'
     cvinit=zero
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        cvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        cvinit(m) = cvinit(m) - xi*sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(cvinit,cvinit)
     cvinit=cvinit/sqrt(norm0)
     call ed_buildH_d(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
     endif
     cnorm2=xi*norm0
     isign=1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
     !
     !build the (N_iorb + xi*N_jorb)|gs> state
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A)")'Apply N_iorb + xi*N_jorb:'
     cvinit=zero
     do m=1,idim
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = dble(ib(iorb))+dble(ib(iorb+Ns))
        cvinit(m) = sgn*state_vec(m)   
        !
        sgn = dble(ib(jorb))+dble(ib(jorb+Ns))
        cvinit(m) = cvinit(m) + xi*sgn*state_vec(m)   
        !
     enddo
     norm0=dot_product(cvinit,cvinit)
     cvinit=cvinit/sqrt(norm0)
     call ed_buildH_d(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_dc,cvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_dc,cvinit,alfa_,beta_)
     endif
     cnorm2=xi*norm0
     isign=-1
     call add_to_lanczos_densChi(cnorm2,state_e,alfa_,beta_,isign,iorb,jorb)
     deallocate(vvinit)
     deallocate(HI%map)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_offdiag_d






!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the TOTAL Charge-Charge susceptibility \Chi_dens  
! \chi_dens_tot = <N(\tau)N(0)>, N=sum_a N_a
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_tot_d()
  integer             :: iorb,isite,isector,izero
  integer             :: numstates
  integer             :: nlanc,idim
  integer             :: iup0,idw0,isign
  integer             :: ib(Nlevels)
  integer             :: m,i,j,r
  complex(8)          :: cnorm2
  real(8)             :: norm0,sgn
  real(8),allocatable :: alfa_(:),beta_(:)
  real(8),allocatable :: vvinit(:)
  integer             :: Nitermax
  type(sector_map)    :: HI    !map of the Sector S to Hilbert space H
  !
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do izero=1,state_list%size
     isector     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim  
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = sum(dble(ib(1:Norb)))+sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = sgn*state_vec(m) 
     enddo
     deallocate(HI%map)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     call ed_buildH_d(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
     endif
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)
     isign=-1
     call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)     
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_d

subroutine lanc_ed_build_densChi_tot_c()
  integer                :: iorb,isite,isector,izero
  integer                :: numstates
  integer                :: nlanc,idim
  integer                :: iup0,idw0,isign
  integer                :: ib(Nlevels)
  integer                :: m,i,j,r
  complex(8)             :: cnorm2
  real(8)                :: norm0,sgn
  real(8),allocatable    :: alfa_(:),beta_(:)
  complex(8),allocatable :: vvinit(:)
  integer                :: Nitermax
  type(sector_map)       :: HI    !map of the Sector S to Hilbert space H
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do izero=1,state_list%size
     isector     =  es_return_sector(state_list,izero)
     idim      =  getdim(isector)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     allocate(vvinit(idim))
     if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")'Apply N:',getnup(isector),getndw(isector)
     call build_sector(isector,HI)
     vvinit=0.d0
     do m=1,idim                     !loop over |gs> components m
        i=HI%map(m)
        ib = bdecomp(i,2*Ns)
        sgn = sum(dble(ib(1:Norb)))+sum(dble(ib(Ns+1:Ns+Norb)))
        vvinit(m) = sgn*state_cvec(m) 
     enddo
     deallocate(HI%map)
     norm0=dot_product(vvinit,vvinit)
     vvinit=vvinit/sqrt(norm0)
     call ed_buildH_c(isector)
     nlanc=min(idim,lanc_nGFiter)
     allocate(alfa_(nlanc),beta_(nlanc))
     if(MpiStatus)then
        call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvinit,alfa_,beta_)
     else
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
     endif
     cnorm2=one*norm0
     isign=1
     call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)
     isign=-1
     call add_to_lanczos_densChi_tot(cnorm2,state_e,alfa_,beta_,isign)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_tot_c















!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the inter-orbital charge susceptibility \Chi_mix 
! \chi_mix = <C^+_a(\tau)N_a(0)>
!+------------------------------------------------------------------+
subroutine lanc_ed_build_densChi_mix_d(iorb,jorb)
  integer             :: iorb,jorb,ispin
  real(8),allocatable :: vvinit(:),vvinit_tmp(:)
  real(8),allocatable :: alfa_(:),beta_(:)
  integer             :: isite,jsite,istate
  integer             :: isector,jsector,ksector
  integer             :: idim,jdim,kdim
  type(sector_map)    :: HI,HJ,HK
  integer             :: ib(Nlevels)
  integer             :: m,i,j,r,numstates
  real(8)             :: sgn,norm2,norm0
  complex(8)          :: cnorm2
  integer             :: Nitermax,Nlanc
  !
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do istate=1,state_list%size
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !+- Apply Sum_ispin c^dg_{jorb,ispin} c_{iorb,ispin} -+!
     do ispin=1,Nspin
        isite=impIndex(iorb,ispin)
        jsector = getCsector(ispin,isector)
        if(jsector/=0)then
           jdim  = getdim(jsector)
           allocate(vvinit_tmp(jdim))
           call build_sector(jsector,HJ)
           vvinit_tmp=0d0
           do m=1,idim
              i=HI%map(m)
              ib = bdecomp(i,2*Ns)
              if(ib(isite)==1)then
                 call c(isite,i,r,sgn)
                 j=binary_search(HJ%map,r)
                 vvinit_tmp(j) = sgn*state_vec(m)
              end if
           enddo
        endif
        jsite = impIndex(jorb,ispin)
        ksector = getCDGsector(ispin,jsector)
        if(ksector/=0) then       
           kdim  = getdim(ksector)
           allocate(vvinit(kdim)) !<==== ACTHUNG! 
           call build_sector(ksector,HK)
           vvinit=0d0              !<==== ACTHUNG! 
           do m=1,jdim
              i=HJ%map(m)
              ib = bdecomp(i,2*Ns)
              if(ib(jsite)==0)then
                 call cdg(jsite,i,r,sgn)
                 j=binary_search(HK%map,r)
                 vvinit(j) = sgn*vvinit_tmp(m)
              endif
           enddo
        end if
        deallocate(HJ%map,HK%map,vvinit_tmp)
        !
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        call ed_buildH_d(ksector)
        nlanc=min(kdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_densChi_mix(cnorm2,state_e,alfa_,beta_,1,iorb,jorb)
        deallocate(vvinit)
     enddo
     !
     !
     !+- Apply Sum_ispin c^dg_{iorb,ispin} c_{jorb,ispin} -+!
     do ispin=1,Nspin
        jsite=impIndex(jorb,ispin)
        jsector = getCsector(ispin,isector)
        if(jsector/=0)then
           jdim  = getdim(jsector)
           allocate(vvinit_tmp(jdim))
           call build_sector(jsector,HJ)
           vvinit_tmp=0d0
           do m=1,idim
              i=HI%map(m)
              ib = bdecomp(i,2*Ns)
              if(ib(jsite)==1)then
                 call c(jsite,i,r,sgn)
                 j=binary_search(HJ%map,r)
                 vvinit_tmp(j) = sgn*state_vec(m)
              endif
           enddo
        endif
        isite = impIndex(iorb,ispin)
        ksector = getCDGsector(ispin,jsector)
        if(ksector/=0) then       
           kdim  = getdim(ksector)
           allocate(vvinit(kdim))
           call build_sector(ksector,HK)
           vvinit=0d0
           do m=1,jdim
              i=HJ%map(m)
              ib = bdecomp(i,2*Ns)
              if(ib(isite)==0)then
                 call cdg(isite,i,r,sgn)
                 j=binary_search(HK%map,r)
                 vvinit(j) = sgn*vvinit_tmp(m)
              endif
           enddo
        end if
        deallocate(HJ%map,HK%map,vvinit_tmp)
        !
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        call ed_buildH_d(ksector)
        nlanc=min(kdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        if(MpiStatus)then
           call sp_lanc_tridiag(MpiComm,spHtimesV_dd,vvinit,alfa_,beta_)
        else
           call sp_lanc_tridiag(spHtimesV_dd,vvinit,alfa_,beta_)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_densChi_mix(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb)
        deallocate(vvinit)
     enddo
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_ed_build_densChi_mix_d











!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,isign,iorb,jorb)
  integer                                    :: iorb,jorb,isign
  complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
  real(8)                                    :: Ei,Ej,Egs,de
  integer                                    :: nlanc
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
  !
  pesoF  = vnorm2/zeta_function 
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
           densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) - peso*beta
        else
           densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*(exp(-beta*dE)-1d0)/dE 
        endif
        do i=1,Lmats
           densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
        enddo
        do i=0,Ltau
           densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + peso*exp(-tau(i)*de)
        enddo
        do i=1,Lreal
           densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
        enddo
     enddo
  case (-1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*beta
        else
           densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
        enddo
        do i=0,Ltau
           densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + peso*exp(-(beta-tau(i))*dE)
        enddo
        do i=1,Lreal
           densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
        enddo
     enddo
  case default
     stop "add_to_lanczos_densChi: isign not in {-1,1}"
  end select
end subroutine add_to_lanczos_densChi




subroutine add_to_lanczos_densChi_mix(vnorm2,Ei,alanc,blanc,isign,iorb,jorb)
  integer                                    :: iorb,jorb,isign
  complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
  real(8)                                    :: Ei,Ej,Egs,de
  integer                                    :: nlanc
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
  !
  pesoF  = vnorm2/zeta_function 
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
           densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) - peso*beta
        else
           densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*(exp(-beta*dE)-1d0)/dE 
        endif
        do i=1,Lmats
           densChi_mix_iv(iorb,jorb,i)=densChi_mix_iv(iorb,jorb,i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
        enddo
        do i=0,Ltau
           densChi_mix_tau(iorb,jorb,i)=densChi_mix_tau(iorb,jorb,i) + peso*exp(-tau(i)*de)
        enddo
        do i=1,Lreal
           densChi_mix_w(iorb,jorb,i)=densChi_mix_w(iorb,jorb,i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
        enddo
     enddo
  case (-1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*beta
        else
           densChi_mix_iv(iorb,jorb,0)=densChi_mix_iv(iorb,jorb,0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           densChi_mix_iv(iorb,jorb,i)=densChi_mix_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
        enddo
        do i=0,Ltau
           densChi_mix_tau(iorb,jorb,i)=densChi_mix_tau(iorb,jorb,i) + peso*exp(-(beta-tau(i))*dE)
        enddo
        do i=1,Lreal
           densChi_mix_w(iorb,jorb,i)=densChi_mix_w(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
        enddo
     enddo
  case default
     stop "add_to_lanczos_densChi_mix: isign not in {-1,1}"
  end select
end subroutine add_to_lanczos_densChi_mix



subroutine add_to_lanczos_densChi_tot(vnorm2,Ei,alanc,blanc,isign)
  complex(8)                                 :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
  real(8)                                    :: Ei,Ej,Egs,de
  integer                                    :: nlanc,isign
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
  !
  pesoF  = vnorm2/zeta_function 
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
           densChi_tot_iv(0)=densChi_tot_iv(0) - peso*beta
        else
           densChi_tot_iv(0)=densChi_tot_iv(0) + peso*(exp(-beta*dE)-1d0)/dE 
        endif
        do i=1,Lmats
           densChi_tot_iv(i)=densChi_tot_iv(i) + peso*(exp(-beta*dE)-1d0)/(dcmplx(0d0,vm(i)) - dE)
        enddo
        do i=0,Ltau
           densChi_tot_tau(i)=densChi_tot_tau(i) + peso*exp(-tau(i)*de)
        enddo
        do i=1,Lreal
           densChi_tot_w(i)=densChi_tot_w(i) + peso*(exp(-beta*dE)-1.d0)/(dcmplx(wr(i),eps) - dE)
        enddo
     enddo
  case (-1)
     do j=1,nlanc
        Ej     = diag(j)
        dE     = Ej-Ei
        pesoAB = Z(1,j)*Z(1,j)
        peso   = pesoF*pesoAB*pesoBZ
        if(beta*dE < 1d-1)then     !abs(X - (1-exp(-X)) is about 5*10^-3 for X<10^-1 this is a satisfactory bound
           densChi_tot_iv(0)=densChi_tot_iv(0) + peso*beta
        else
           densChi_tot_iv(0)=densChi_tot_iv(0) + peso*(1d0-exp(-beta*dE))/dE 
        endif
        do i=1,Lmats
           densChi_tot_iv(i)=densChi_tot_iv(i) + peso*(1d0-exp(-beta*dE))/(dcmplx(0d0,vm(i)) + dE)
        enddo
        do i=0,Ltau
           densChi_tot_tau(i)=densChi_tot_tau(i) + peso*exp(-(beta-tau(i))*dE)
        enddo
        do i=1,Lreal
           densChi_tot_w(i)=densChi_tot_w(i) + peso*(1d0-exp(-beta*dE))/(dcmplx(wr(i),eps) + dE)
        enddo
     enddo
  case default
     stop "add_to_lanczos_densChi_tot: isign not in {-1,1}"
  end select

end subroutine add_to_lanczos_densChi_tot

