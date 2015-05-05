!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_hybrid_superc(fg,bath_,ispin)
  complex(8),dimension(:,:,:,:)        :: fg ![2][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)   :: bath_
  integer                              :: ispin
  real(8),dimension(:),allocatable     :: array_bath
  integer                              :: iter,stride,stride_spin,stride_orb,i,j,corb,l,Asize
  integer                              :: iorb,jorb
  real(8)                              :: chi
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  character(len=20)                    :: suffix
  integer                              :: unit
  !
  if(size(fg,1)/=2)stop"chi2_fitgf_normal_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
  if(size(fg,3)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_hybrid_superc: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
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
  if(totNorb/=(Norb*(Norb+1)/2))stop "chi2_fitgf_hybrid_superc error counting the orbitals"
  !
  allocate(Gdelta(totNorb,Ldelta))
  allocate(Fdelta(totNorb,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  do i=1,totNorb
     Gdelta(i,1:Ldelta) = fg(1,getIorb(i),getJorb(i),1:Ldelta)
     Fdelta(i,1:Ldelta) = fg(2,getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  select case(Cg_weight)
  case default
     Wdelta=dble(Ldelta)
  case(1)
     Wdelta=1.d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  !
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
  !
  call allocate_bath(chi2_bath)
  call set_bath(bath_,chi2_bath)
  !
  Asize = get_chi2_bath_size()
  allocate(array_bath(Asize))
  !
  Spin_indx=ispin
  !
  call dmft_bath2chi2_bath(dmft_bath,array_bath,ispin)
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,&
             chi2_weiss_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol  ,&
             istop=cg_stop ,&
             eps=cg_eps)
     case ("delta")
        call fmin_cg(array_bath,&
             chi2_delta_hybrid_superc,&
             grad_chi2_delta_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol  ,&
             istop=cg_stop ,&
             eps=cg_eps)
     case default
        stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,&
             chi2_weiss_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case ("delta")
        call fmin_cgminimize(array_bath,&
             chi2_delta_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
     end select
     !
  case (2)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgplus(array_bath,&
             chi2_weiss_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case ("delta")
        call fmin_cgplus(array_bath,&
             chi2_delta_hybrid_superc,&
             grad_chi2_delta_hybrid_superc,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
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
  call chi2_bath2dmft_bath(array_bath,dmft_bath,ispin)
  !
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  !
  ! unit=free_unit()
  ! open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  ! call write_bath(dmft_bath,unit)
  ! close(unit)
  call save_bath(dmft_bath)
  !
  if(ed_verbose<3)call write_fit_result(ispin)
  !
  call copy_bath(dmft_bath,bath_)
  call deallocate_bath(dmft_bath)
  deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
  deallocate(getIorb,getJorb)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer                              :: i,j,l,m,iorb,jorb,ispin,jspin
    real(8)                              :: w
    complex(8),dimension(2,Norb,Norb)    :: gwf
    complex(8),dimension(totNorb,Ldelta) :: fgand,ffand
    do i=1,Ldelta
       w=Xdelta(i)
       if(cg_scheme=='weiss')then
          gwf(1,:,:) = g0and_bath_mats(ispin,ispin,xi*w,dmft_bath)
          gwf(2,:,:) = f0and_bath_mats(ispin,ispin,xi*w,dmft_bath)
       else
          do l=1,Norb
             do m=1,Norb
                gwf(1,l,m) = delta_bath_mats(ispin,ispin,l,m,xi*w,dmft_bath)
                gwf(2,l,m) = fdelta_bath_mats(ispin,ispin,l,m,xi*w,dmft_bath)
             enddo
          enddo
       endif
       !
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          fgand(l,i)=gwf(1,iorb,jorb)
          ffand(l,i)=gwf(2,iorb,jorb)
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
          write(unit,"(10F24.15)")Xdelta(i),&
               dimag(Gdelta(l,i)),dimag(fgand(l,i)),&
               dreal(Gdelta(l,i)),dreal(fgand(l,i)),&
               dimag(Fdelta(l,i)),dimag(ffand(l,i)),&
               dreal(Fdelta(l,i)),dreal(ffand(l,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_hybrid_superc



!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_delta_hybrid_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  real(8),dimension(totNorb)     :: chi_orb
  complex(8),dimension(2,Ldelta) ::  g0
  real(8)                        ::  chi2,w
  integer                        ::  i,l,iorb,jorb,ispin
  ispin=Spin_indx
  call chi2_bath2dmft_bath(a,chi2_bath,ispin)
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        w = xdelta(i)
        g0(1,i) = delta_bath_mats(ispin,ispin,iorb,jorb,xi*w,chi2_bath)
        g0(2,i) = fdelta_bath_mats(ispin,ispin,iorb,jorb,xi*w,chi2_bath)
     enddo
     chi_orb(l) = &
          sum(abs(Gdelta(l,:)-g0(1,:))**2/Wdelta(:)) + &
          sum(abs(Fdelta(l,:)-g0(2,:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi_orb)
  !
end function chi2_delta_hybrid_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson 
! function in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function grad_chi2_delta_hybrid_superc(a) result(dchi2)
  real(8),dimension(:)                   ::  a
  real(8),dimension(size(a))             ::  dchi2
  real(8),dimension(totNorb,size(a))     ::  df
  complex(8),dimension(2,Ldelta)         ::  g0
  complex(8),dimension(2,Ldelta,size(a)) ::  dg0
  integer                                ::  i,j,l,iorb,jorb,ispin
  real(8)                                ::  w
  ispin=Spin_indx
  call chi2_bath2dmft_bath(a,chi2_bath,ispin)
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        w  = Xdelta(i)
        g0(1,i)  = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
        g0(2,i)  = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
        dg0(1,i,:) = grad_delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath,size(a))
        dg0(2,i,:) = grad_fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath,size(a))
     enddo
     !
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Gdelta(l,:)-g0(1,:))*dreal(dg0(1,:,j))/Wdelta(:) ) + &
             sum( dimag(Gdelta(l,:)-g0(1,:))*dimag(dg0(1,:,j))/Wdelta(:) ) + &
             sum( dreal(Fdelta(l,:)-g0(2,:))*dreal(dg0(2,:,j))/Wdelta(:) ) + &
             sum( dimag(Fdelta(l,:)-g0(2,:))*dimag(dg0(2,:,j))/Wdelta(:) )
     enddo
     !
  enddo
  !
  dchi2=-2.d0*sum(df,dim=1)
  !
end function grad_chi2_delta_hybrid_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance for G_0 function 
! in the SUPERCONDUCTING case.
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_hybrid_superc(a) result(chi2)
  real(8),dimension(:)                     ::  a
  real(8),dimension(totNorb)               ::  chi2_orb
  complex(8),dimension(2,Norb,Norb,Ldelta) ::  g0
  real(8)                                  ::  chi2,w
  integer                                  ::  i,l,iorb,jorb,ispin
  ispin=Spin_indx
  call chi2_bath2dmft_bath(a,chi2_bath,ispin)
  do i=1,Ldelta
     w   = Xdelta(i)
     g0(1,:,:,i)  = g0and_bath_mats(ispin,ispin,xi*w,chi2_bath)
     g0(2,:,:,i)  = f0and_bath_mats(ispin,ispin,xi*w,chi2_bath)
  enddo
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     chi2_orb(l) = sum(abs(Gdelta(l,:)-g0(1,iorb,jorb,:))**2/Wdelta(:)) &
          + sum(abs(Fdelta(l,:)-g0(2,iorb,jorb,:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi2_orb)
  !
end function chi2_weiss_hybrid_superc
