!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Hybrid bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_hybrid_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)          :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:,:),intent(inout) :: bath_
  integer                              :: ispin
  real(8),dimension(:),allocatable     :: array_bath
  integer                              :: iter,stride,ifirst,ilast,i,j,corb,l,Asize
  integer                              :: iorb,jorb
  real(8)                              :: chi
  complex(8),dimension(size(fg,3))     :: g0
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  real(8)                              :: w
  character(len=20)                    :: suffix
  integer                              :: unit
  !
  if(size(fg,1)/=Norb)stop"chi2_fitgf_hybrid_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_hybrid_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_hybrid_normal error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
  corb=0
  do iorb=1,Norb
     do jorb=iorb,Norb
        corb=corb+1
        getIorb(corb)=iorb
        getJorb(corb)=jorb
     enddo
  enddo
  totNorb=corb
  if(totNorb/=(Norb*(Norb+1)/2))stop "chi2_fitgf_hybrid_normal error counting the orbitals"
  !
  allocate(Fdelta(totNorb,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  do i=1,totNorb
     Fdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  select case(cg_weight)
  case default
     Wdelta=1d0*Ldelta
  case(1)
     Wdelta=1d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
  !
  call allocate_bath(chi2_bath)
  call set_bath(bath_,chi2_bath)
  !
  Asize = get_array_bath_dimension()
  allocate(array_bath(Asize))
  !
  Spin_indx=ispin
  !
  call dmft_bath2array(dmft_bath,array_bath,ispin)
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,&
             chi2_weiss_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol  ,&
             istop=cg_stop ,&
             eps=cg_eps)
     case ("delta")
        call fmin_cg(array_bath,&
             chi2_delta_hybrid_normal,&
             grad_chi2_delta_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol  ,&
             istop=cg_stop ,&
             eps=cg_eps)
     case default
        stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,&
             chi2_weiss_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case ("delta")
        call fmin_cgminimize(array_bath,&
             chi2_delta_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  case (2)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgplus(array_bath,&
             chi2_weiss_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case ("delta")
        call fmin_cgplus(array_bath,&
             chi2_delta_hybrid_normal,&
             grad_chi2_delta_hybrid_normal,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5)")&
       'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
       "  <--  all Orbs, Spin"//reg(txtfy(ispin))
  !
  if(ed_verbose<2)then
     suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
  endif
  !
  call array2dmft_bath(array_bath,dmft_bath,ispin)
  !
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  unit=free_unit()
  open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  call write_bath(dmft_bath,unit)
  close(unit)
  if(ed_verbose<3)call write_fit_result(ispin)
  call deallocate_bath(dmft_bath)
  deallocate(Fdelta,Xdelta,Wdelta)
  deallocate(getIorb,getJorb)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer                              :: i,j,l,m,iorb,jorb,ispin,jspin
    real(8)                              :: w
    character(len=20)                    :: suffix
    integer                              :: unit
    complex(8),dimension(Norb,Norb)      :: gwf
    complex(8),dimension(totNorb,Ldelta) :: fgand
    do i=1,Ldelta
       w=Xdelta(i)
       if(cg_scheme=='weiss')then
          do l=1,Norb
             do m=1,Norb
                gwf(l,m) = weiss_bath_mats(ispin,ispin,l,m,xi*w,dmft_bath)
             enddo
          enddo
       else
          do l=1,Norb
             do m=1,Norb
                gwf(l,m) = delta_bath_mats(ispin,ispin,l,m,xi*w,dmft_bath)
             enddo
          enddo
       endif
       !
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          fgand(l,i)=gwf(iorb,jorb)
       enddo
       !
    enddo
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(l,i)),dimag(fgand(l,i)),&
               dreal(Fdelta(l,i)),dreal(fgand(l,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_hybrid_normal




!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
! The gradient is evaluated analytically d\chi^2.
! \Delta_Anderson  and its gradient d\Delta_Anderson are evaluated 
! as separate functions.  
! HYBRID bath. 
! SPIN  DIAGONAL ; ORBITAL NON-DIAGONAL
!+-------------------------------------------------------------+
function chi2_delta_hybrid_normal(a) result(chi2)
  real(8),dimension(:)         :: a
  real(8),dimension(totNorb)   :: chi_orb
  complex(8),dimension(Ldelta) :: Delta_orb
  real(8)                      :: chi2,w
  integer                      :: i,l,iorb,jorb,ispin
  chi_orb = 0.d0
  ispin=Spin_indx
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        w = xdelta(i)
        Delta_orb(i) = delta_bath_mats(ispin,ispin,iorb,jorb,xi*w,chi2_bath)
     enddo
     chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi_orb)
  !
end function chi2_delta_hybrid_normal
! the analytic GRADIENT of \chi^2
!+-------------------------------------------------------------+
function dchi2_delta_hybrid_normal(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(totNorb,size(a))   :: df
  complex(8),dimension(Ldelta)         :: g0
  complex(8),dimension(Ldelta,size(a)) :: dg0
  integer                              :: i,j,l,iorb,jorb,ispin
  real(8)                              :: w
  df=0.d0
  ispin=Spin_indx
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        w = xdelta(i)
        g0(i)    = delta_bath_mats(ispin,ispin,iorb,jorb,xi*w,chi2_bath)
        dg0(i,:) = grad_delta_bath_mats(ispin,ispin,iorb,jorb,xi*w,chi2_bath,size(a))
     enddo
     !
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Fdelta(l,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
             sum( dimag(Fdelta(l,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
     enddo
  enddo
  !
  dchi2 = -2.d0*sum(df,1)     !sum over all orbital indices
  !
end function dchi2_delta_hybrid_normal


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distanec of G_0  function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
! G_0 is evaluated in a function below
! HYBRID bath.
! SPIN DIAGONAL & ORBITAL OFF-DIAGONAL
!+-------------------------------------------------------------+
function chi2_weiss_hybrd(a) result(chi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(totNorb)           :: chi_orb
  complex(8),dimension(Ldelta)         :: Delta_orb
  real(8)                              :: chi2,w
  integer                              :: i,l,iorb,jorb,ispin
  chi_orb = 0.d0
  ispin=Spin_indx
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        w      = Xdelta(i)
        Delta_orb(i)= weiss_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     enddo
     chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi_orb)
  !
end function chi2_weiss_hybrd
!
