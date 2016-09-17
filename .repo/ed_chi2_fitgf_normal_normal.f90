!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)        :: fg
  real(8),dimension(:),intent(inout) :: bath_
  integer                            :: ispin
  ! real(8),dimension(2*Nbath)         :: a
  real(8),dimension(:),allocatable   :: a
  integer                            :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  complex(8)                         :: fgand
  real(8)                            :: w
  character(len=20)                  :: suffix
  integer                            :: unit
  if(size(fg,1)/=Norb)stop"chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(Fdelta(1,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !forall(i=1:Ldelta)Xdelta(i)=pi/beta*dble(2*i-1)
  !
  select case(cg_weight)
  case default
     Wdelta=1d0*Ldelta
  case(1)
     Wdelta=1d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)!(/(dble(i),i=1,Ldelta)/)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_bath(dmft_bath)
  call allocate_bath(chi2_bath)
  call set_bath(bath_,dmft_bath)
  Asize = get_array_bath_dimension()
  allocate(a(Asize))
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     Fdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
     !
     call dmft_bath2array(dmft_bath,a,ispin,iorb)
     ! a(1:Nbath)         = dmft_bath%e(ispin,iorb,1:Nbath)
     ! a(Nbath+1:2*Nbath) = dmft_bath%v(ispin,iorb,1:Nbath)
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(a,chi2_weiss_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case ("delta")
           call fmin_cg(a,chi2_delta_irred,dchi2_delta_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(a,chi2_weiss_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(a,chi2_delta_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(a,chi2_weiss_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(a,chi2_delta_irred,dchi2_delta_irred,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
     end select
     if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)")&
          "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
     if(ed_verbose<2)then
        suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
        unit=free_unit()
        open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
        write(unit,"(ES18.9,1x,I5)") chi,iter
        close(unit)
     endif
     dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
     dmft_bath%v(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
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
    complex(8)        :: fgand
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Fdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
       fgand=zero
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          w = Xdelta(i)
          if(cg_scheme=='weiss')then
             fgand = xi*w + xmu - impHloc(ispin,ispin,iorb,iorb) - delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
             fgand = one/fgand
          else
             fgand = delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
          endif
          write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand),dreal(Fdelta(1,i)),dreal(fgand)
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_normal



!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
! The gradient is evaluated analytically d\chi^2.
! \Delta_Anderson  and its gradient d\Delta_Anderson are evaluated 
! as separate functions.  
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_delta_irred(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0
  real(8)                      ::  chi2
  integer                      ::  i,iorb,ispin
  type(effective_bath)         ::  dmft_bath
  chi2 = 0.d0 
  iorb=Orb_indx
  ispin=Spin_indx
  !<DEBUG   
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     g0(i)   = delta_bath_mats(ispin,ispin,iorb,iorb,xdelta(i),chi2_bath)
  enddo
  !>DEBUG
  ! do i=1,Ldelta
  !    g0(i)   = fg_delta_irred(xdelta(i),iorb,a)
  ! enddo
  chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
end function chi2_delta_irred
! the analytic GRADIENT of \chi^2
!+-------------------------------------------------------------+
function dchi2_delta_irred(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(size(a))           :: df
  complex(8),dimension(Ldelta)         :: g0
  complex(8),dimension(Ldelta,size(a)) :: dg0
  integer                              :: i,j,iorb
  df=0.d0
  iorb=Orb_indx
  do i=1,Ldelta
     g0(i)    = fg_delta_irred(xdelta(i),iorb,a)
     dg0(i,:) = grad_fg_delta_irred(xdelta(i),iorb,a)
  enddo
  do j=1,size(a)
     df(j)=sum( dreal(Fdelta(1,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
          sum(  dimag(Fdelta(1,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
  enddo
  dchi2 = -2.d0*df
end function dchi2_delta_irred
! ! the \Delta_Anderson function used in \chi^2 and d\chi^2
! ! \Delta = \sum_l V_l^2/(iw-e_l)
! !+-------------------------------------------------------------+
! function fg_delta_irred(w,iorb,a) result(gg)
!   real(8)                      :: w
!   real(8),dimension(:)         :: a
!   real(8),dimension(size(a)/2) :: eps,vps
!   complex(8)                   :: gg,delta
!   integer                      :: i,iorb,Nb
!   Nb=size(a)/2
!   eps=a(1:Nb)
!   vps=a(Nb+1:2*Nb)
!   delta=sum(vps(:)**2/(xi*w-eps(:)))
!   gg=delta
! end function fg_delta_irred
! the gradient d\Delta_Anderson function used in d\chi^2
! d\Delta = \grad_{V_k,e_k}\sum_l V_l^2/(iw-e_l)
!+-------------------------------------------------------------+
function grad_fg_delta_irred(w,iorb,a) result(dgz)
  real(8)                         :: w,sgn
  real(8),dimension(:)            :: a
  real(8),dimension(size(a)/2)    :: eps,vps
  complex(8),dimension(size(a))   :: dgz
  complex(8)                      :: gg,delta
  integer                         :: i,iorb,Nb
  dgz=zero
  Nb=size(a)/2
  eps=a(1:Nb)
  vps=a(Nb+1:2*Nb)
  do i=1,Nb
     dgz(i)    = vps(i)*vps(i)/(xi*w-eps(i))**2
     dgz(i+Nb) = 2.d0*vps(i)/(xi*w-eps(i))
  enddo
end function grad_fg_delta_irred


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distanec of G_0  function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
! G_0 is evaluated in a function below
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_weiss_irred(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0
  real(8)                      ::  chi2
  integer                      ::  i,iorb
  chi2 = 0.d0 
  iorb=Orb_indx
  do i=1,Ldelta   !Number of freq. in common to the module
     g0(i)   = fg_weiss_irred(xdelta(i),iorb,a)
  enddo
  chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
end function chi2_weiss_irred
! the inverse non-interacting GF (~inverse Weiss Field) 
! used in \chi^2(\caG_0 - G_0) 
! \G_0 = [iw_n + xmu - H_loc - \Delta]^-1
!+-------------------------------------------------------------+
function fg_weiss_irred(w,iorb,a) result(gg)
  real(8)                      :: w
  real(8),dimension(:)         :: a
  complex(8)                   :: gg,delta
  integer                      :: i,iorb,ispin
  ispin=Spin_indx
  delta=fg_delta_irred(w,iorb,a)
  gg = one/(xi*w+xmu-impHloc(ispin,ispin,iorb,iorb)-delta)
end function fg_weiss_irred
