!##################################################################
! THE CALCULATION OF THE \chi^2 FUNCTIONS USE PROCEDURES FURTHER 
! BELOW TO EVALUATE INDEPENDENTLY THE ANDERSON MODEL:
!  - DELTA, 
!  -\GRAD DELTA
!  - G0
! THE LATTER ARE ADAPTED FROM THE PROCEDURES:
! DELTA_BATH_MATS
! GRAD_DELTA_BATH_MATS
! G0 BATH_MATS
! FOR, YOU NEED TO DECOMPOSE THE a INPUT ARRAY INTO ELEMENTS.
!##################################################################


!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)        :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  integer                            :: ispin
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: iter,stride,iorb,i,io,j,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=20)                  :: suffix
  integer                            :: unit
  !
  if(size(fg,1)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(Gdelta(1,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
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
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !Asize = get_chi2_bath_size()
  !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  Asize = Nbath + Nbath
  allocate(array_bath(Asize))
  !
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     !
     Gdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
     !
     !Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%e(ispin,iorb,i)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%v(ispin,iorb,i)
     enddo
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_normal_normal,grad_chi2_delta_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,chi2_weiss_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,chi2_delta_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(array_bath,chi2_weiss_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,chi2_delta_normal_normal,grad_chi2_delta_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     !
     write(LOGfile,"(A,ES18.9,A,I5,A)")&
          "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
     !
     suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
     !
     !Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        dmft_bath%e(ispin,iorb,i) = array_bath(io) 
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%v(ispin,iorb,i) = array_bath(io)
     enddo
     !
  enddo
  !
  call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  call write_fit_result(ispin)
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(Ldelta)
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Gdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
          fgand = g0and_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
          fgand = delta_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),dimag(Gdelta(1,i)),dimag(fgand(i)),dreal(Gdelta(1,i)),dreal(fgand(i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_normal








!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  real(8)                      ::  chi2
  complex(8),dimension(Ldelta) ::  Delta
  !
  Delta = delta_normal_normal(a)
  !
  chi2=sum(abs(Gdelta(1,:)-Delta(:))**2/Wdelta(:))
  !
end function chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_normal(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(size(a))           :: df
  complex(8),dimension(Ldelta)         :: Delta
  complex(8),dimension(Ldelta,size(a)) :: dDelta
  integer                              :: j
  !
  Delta   = delta_normal_normal(a)
  dDelta  = grad_delta_normal_normal(a)
  !
  do j=1,size(a)
     df(j)=sum( dreal(Gdelta(1,:)-Delta(:))*dreal(dDelta(:,j))/Wdelta(:) ) + &
          sum(  dimag(Gdelta(1,:)-Delta(:))*dimag(dDelta(:,j))/Wdelta(:) )
  enddo
  !
  dchi2 = -2.d0*df
  !
end function grad_chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0and
  real(8)                      ::  chi2,w
  !
  g0and  = g0and_normal_normal(a)
  !
  chi2=sum(abs(Gdelta(1,:)-g0and(:))**2/Wdelta(:))
  !
end function chi2_weiss_normal_normal





!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################

function delta_normal_normal(a) result(Delta)
  real(8),dimension(:)         :: a
  complex(8),dimension(Ldelta) :: Delta
  integer                      :: i,io,stride
  real(8),dimension(Nbath)     :: eps,vps
  !
  !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io) 
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  do i=1,Ldelta
     Delta(i) = sum( vps(:)*vps(:)/(xi*Xdelta(i) - eps(:)) )
  enddo
  !
end function delta_normal_normal


function grad_delta_normal_normal(a) result(dDelta)
  real(8),dimension(:)                 :: a
  complex(8),dimension(Ldelta,size(a)) :: dDelta
  integer                              :: i,k,ik,io,stride
  real(8),dimension(Nbath)             :: eps,vps
  complex(8)                           :: iw
  !
  !
  !\grad_{E_{a}(k)} \Delta_{bb}^{rr} = [ V_{a}(k)*V_{a}(k) / ( iw_n - E_{a}(k) )**2 ]
  !
  !\grad_{V_{a}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}(k) / ( iw_n - E_{a}(k) ) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io) 
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(:,ik) = vps(k)*vps(k)/(xi*Xdelta(:) - eps(k))**2
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(:,ik) = 2d0*vps(k)/(xi*Xdelta(:) - eps(k))
  enddo
  !
end function grad_delta_normal_normal

function g0and_normal_normal(a) result(G0and)
  real(8),dimension(:)         :: a
  complex(8),dimension(Ldelta) :: G0and,Delta
  integer                      :: i,io,iorb,ispin
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta(:) = delta_normal_normal(a)
  G0and(:) = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(:)
  G0and(:) = one/G0and(:)
  !
end function g0and_normal_normal























