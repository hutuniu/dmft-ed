!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_superc(fg,bath_,ispin)
  complex(8),dimension(:,:,:,:)        :: fg
  real(8),dimension(:,:),intent(inout) :: bath_
  integer                              :: ispin
  real(8),dimension(3*Nbath)           :: a
  integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
  real(8)                              :: chi
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  real(8)                              :: w
  character(len=20)                    :: suffix
  integer                              :: unit
  if(size(fg,1)/=2)stop"CHI2_FITGF: wrong dimension 1 in chi2_input"
  if(size(fg,2)/=Norb)stop"CHI2_FITGF: wrong dimension 2 in chi2_input"
  if(size(fg,3)/=Norb)stop"CHI2_FITGF: wrong dimension 3 in chi2_input"
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
  Ldelta = Lfit
  if(Ldelta>size(fg,4))Ldelta=size(fg,4)
  !
  allocate(Fdelta(2,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
  !
  select case(Cg_weight)
  case default
     Wdelta=dble(Ldelta)
  case(1)
     Wdelta=1.d0
  case(2)
     Wdelta=(/(real(i,8),i=1,Ldelta)/)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
     Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
     a(1:Nbath)         = dmft_bath%e(ispin,iorb,1:Nbath) 
     a(Nbath+1:2*Nbath) = dmft_bath%d(ispin,iorb,1:Nbath)
     a(2*Nbath+1:3*Nbath) = dmft_bath%v(ispin,iorb,1:Nbath)
     if(cg_method==0)then
        if(cg_scheme=='weiss')then
           call fmin_cg(a,chi2_weiss_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        else
           call fmin_cg(a,chi2_delta_irred_sc,dchi2_delta_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.,istop=cg_stop,eps=cg_eps)
        endif
     elseif(cg_method==1)then
        if(cg_scheme=='weiss')then
           call fmin_cgminimize(a,chi2_weiss_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        else
           call fmin_cgminimize(a,chi2_delta_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        endif
     elseif(cg_method==2)then
        if(cg_scheme=='weiss')then
           call fmin_cgplus(a,chi2_weiss_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        else
           call fmin_cgplus(a,chi2_delta_irred_sc,dchi2_delta_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        endif
     else
        stop "ED_CHI2FIT: error cg_method > 2"
     end if
     if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)") 'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,"  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
     if(ed_verbose<2)then
        suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
        unit=free_unit()
        open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
        write(unit,"(ES18.9,1x,I5)") chi,iter
        close(unit)
     endif
     dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
     dmft_bath%d(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
     dmft_bath%v(ispin,iorb,1:Nbath) = a(2*Nbath+1:3*Nbath)
  enddo
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  unit=free_unit()
  open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  call write_bath(dmft_bath,unit)
  close(unit)
  if(ed_verbose<3)call write_fit_result(ispin)
  call copy_bath(dmft_bath,bath_)
  call deallocate_bath(dmft_bath)
  deallocate(Fdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(2),det
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
       Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
       fgand=zero
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          w = Xdelta(i)
          if(cg_scheme=='weiss')then
             fgand(1) = xi*w + xmu - impHloc(ispin,ispin,iorb,iorb) - delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
             fgand(2) = -fdelta_bath_mats(ispin,iorb,xi*w,dmft_bath)
             det     =  abs(fgand(1))**2 + (fgand(2))**2
             fgand(1) = conjg(fgand(1))/det
             fgand(2) = fgand(2)/det
          else
             fgand(1) = delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
             fgand(2) = fdelta_bath_mats(ispin,iorb,xi*w,dmft_bath)
          endif
          write(unit,"(10F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand(1)),dreal(Fdelta(1,i)),dreal(fgand(1)), &
               dimag(Fdelta(2,i)),dimag(fgand(2)),dreal(Fdelta(2,i)),dreal(fgand(2))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_superc



!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
! The gradient is not evaluated.
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_delta_irred_sc(a) result(chi2_sc)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0
  real(8)                        ::  chi2_sc
  integer                        ::  i,iorb
  chi2_sc = 0.d0 
  iorb=Orb_indx
  do i=1,Ldelta
     g0(:,i)   = fg_delta_irred_sc(xdelta(i),iorb,a)
  enddo
  chi2_sc=sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:)) 
end function chi2_delta_irred_sc
! the analytic GRADIENT of \chi^2
!+-------------------------------------------------------------+
function dchi2_delta_irred_sc(a) result(dchi2_sc)
  real(8),dimension(:)                   ::  a
  complex(8),dimension(2,Ldelta)         ::  g0
  real(8),dimension(size(a))             ::  dchi2_sc,df
  complex(8),dimension(2,Ldelta,size(a)) ::  dg0
  integer                                ::  i,j,iorb
  dchi2_sc = 0.d0 
  df = 0.d0
  iorb=Orb_indx
  do i=1,Ldelta
     g0(:,i)   = fg_delta_irred_sc(xdelta(i),iorb,a)
     dg0(:,i,:)= grad_fg_delta_irred_sc(xdelta(i),iorb,a)
  enddo
  do j=1,size(a)
     df(j) = &
          sum( dreal(Fdelta(1,:)-g0(1,:))*dreal(dg0(1,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(1,:)-g0(1,:))*dimag(dg0(1,:,j))/Wdelta(:) ) +&
          sum( dreal(Fdelta(2,:)-g0(2,:))*dreal(dg0(2,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(2,:)-g0(2,:))*dimag(dg0(2,:,j))/Wdelta(:) )
  enddo
  dchi2_sc=-2.d0*df
end function dchi2_delta_irred_sc

! the \Delta_Anderson function used in \chi^2 and d\chi^2
! \Delta = \sum_l V_l^2/(iw-e_l)
!+-------------------------------------------------------------+
function fg_delta_irred_sc(w,iorb,a) result(gg)
  real(8)                      :: w
  real(8),dimension(:)         :: a
  real(8),dimension(size(a)/3) :: eps,vps,dps
  complex(8)                   :: gg(2),delta(2),x
  integer                      :: i,iorb,Nb
  Nb=size(a)/3
  eps=a(1:Nb)
  dps=a(Nb+1:2*Nb)
  vps=a(2*Nb+1:3*Nb)
  gg = zero
  x = xi*w
  gg(1) = -sum(vps(:)**2*(x+eps(:))/(dimag(x)**2+eps(:)**2+dps(:)**2))
  gg(2) =  sum(dps(:)*vps(:)**2/(dimag(x)**2+eps(:)**2+dps(:)**2))
end function fg_delta_irred_sc
! the gradient d\Delta_Anderson function used in d\chi^2
! d\Delta = \grad_{e_k,d_k,V_k}Delta(e_k,d_k,V_k)
!+-------------------------------------------------------------+
function grad_fg_delta_irred_sc(w,iorb,a) result(dgz)
  real(8)                         :: w
  real(8),dimension(:)            :: a
  real(8),dimension(size(a)/3)    :: eps,vps,dps
  complex(8),dimension(2,size(a)) :: dgz
  complex(8)                      :: x,den
  integer                         :: i,iorb,Nb
  Nb=size(a)/3
  eps=a(1:Nb)
  dps=a(Nb+1:2*Nb)
  vps=a(2*Nb+1:3*Nb)
  x = xi*w
  do i=1,Nb
     den = dimag(x)**2+eps(i)**2+dps(i)**2
     dgz(1,i) = -vps(i)*vps(i)*(1.d0/den - 2.d0*eps(i)*(x+eps(i))/den**2)
     dgz(1,i+Nb) = 2.d0*vps(i)*vps(i)*dps(i)*(x+eps(i))/den**2
     dgz(1,i+2*Nb) = -2.d0*vps(i)*(x+eps(i))/den
     !
     dgz(2,i) = -2.d0*vps(i)*vps(i)*dps(i)*eps(i)/den**2
     dgz(2,i+Nb) = vps(i)*vps(i)*(1.d0/den - 2.d0*dps(i)*dps(i)/den**2)
     dgz(2,i+2*Nb) = 2.d0*vps(i)*dps(i)/den
  enddo
end function grad_fg_delta_irred_sc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0 function 
!         in the SUPERCONDUCTING case.
! The gradient is not evaluated.
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_weiss_irred_sc(a) result(chi2_sc)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0
  real(8)                        ::  chi2_sc
  integer                        ::  i,iorb
  chi2_sc = 0.d0 
  iorb=Orb_indx
  do i=1,Ldelta
     g0(:,i)   = fg_weiss_irred_sc(xdelta(i),iorb,a)
  enddo
  chi2_sc = sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) 
  chi2_sc = chi2_sc + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:)) 
end function chi2_weiss_irred_sc
! the non interacting GF (~ inverse weiss field) 
!+-------------------------------------------------------------+
function fg_weiss_irred_sc(w,iorb,a) result(gg)
  real(8)                      :: w
  real(8),dimension(:)         :: a
  complex(8)                   :: gg(2),g0(2),delta(2),det
  integer                      :: i,iorb,ispin
  ispin = Spin_indx
  delta = fg_delta_irred_sc(w,iorb,a)
  g0(1) = xi*w + xmu - impHloc(ispin,ispin,iorb,iorb) - delta(1)
  g0(2) = -delta(2)
  det   = abs(g0(1))**2 + (g0(2))**2
  gg(1) = conjg(g0(1))/det
  gg(2) = g0(2)/det
end function fg_weiss_irred_sc
