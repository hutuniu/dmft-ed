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
!PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_superc(fg,bath_,ispin)
  complex(8),dimension(:,:,:,:)        :: fg ![2][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)   :: bath_
  integer                              :: ispin
  real(8),dimension(:),allocatable     :: array_bath
  integer                              :: iter,stride,i,io,j,iorb,Asize
  real(8)                              :: chi
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  character(len=20)                    :: suffix
  integer                              :: unit
  !
  if(size(fg,1)/=2)stop "chi2_fitgf_normal_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_superc: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
  !
  allocate(Gdelta(1,Ldelta))
  allocate(Fdelta(1,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*dble(2*arange(1,Ldelta)-1)
  !
  select case(Cg_weight)
  case default
     Wdelta=dble(Ldelta)
  case(1)
     Wdelta=1.d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     wdelta=Xdelta
  end select
  !
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !D_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  Asize = Nbath + Nbath + Nbath
  allocate(array_bath(Asize))
  !
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     !
     Gdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
     Fdelta(1,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
     !
     !3*Nbath == Nbath + Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%e(ispin,iorb,i)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%d(ispin,iorb,i)
     enddo
     stride = 2*Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%v(ispin,iorb,i)
     enddo
     !
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_normal_superc,grad_chi2_delta_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,chi2_weiss_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,chi2_delta_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(array_bath,chi2_weiss_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,chi2_delta_normal_superc,grad_chi2_delta_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     write(LOGfile,"(A,ES18.9,A,I5,A)")&
          'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
     !
     if(ed_verbose<2)then
        suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
        unit=free_unit()
        open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
        write(unit,"(ES18.9,1x,I5)") chi,iter
        close(unit)
     endif
     !
     stride = 0
     do i=1,Nbath
        io = stride + i
        dmft_bath%e(ispin,iorb,i) = array_bath(io)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%d(ispin,iorb,i) = array_bath(io) 
     enddo
     stride = 2*Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%v(ispin,iorb,i) = array_bath(io)
     enddo
     !
  enddo
  if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  if(ed_verbose<3)call write_fit_result(ispin)
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(2,Ldelta)
    integer           :: i,j,iorb,ispin,gunit,funit
    real(8)           :: w
    do iorb=1,Norb
       !
       Gdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
       Fdelta(1,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
       !
       if(cg_scheme=='weiss')then
          fgand(1,:) = g0and_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
          fgand(2,:) = f0and_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
       else
          fgand(1,:) = delta_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
          fgand(2,:) = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*Xdelta(:),dmft_bath)
       endif
       !
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       if(cg_scheme=='weiss')then
          gunit=free_unit()
          open(gunit,file="fit_weiss"//reg(suffix)//".ed")
          funit=free_unit()
          open(funit,file="fit_fweiss"//reg(suffix)//".ed")
       else
          gunit=free_unit()
          open(gunit,file="fit_delta"//reg(suffix)//".ed")
          funit=free_unit()
          open(funit,file="fit_fdelta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(gunit,"(10F24.15)")Xdelta(i),dimag(Gdelta(1,i)),dimag(fgand(1,i)),dreal(Gdelta(1,i)),dreal(fgand(1,i))
          write(funit,"(10F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand(2,i)),dreal(Fdelta(1,i)),dreal(fgand(2,i))
       enddo
       close(gunit)
       close(funit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_superc







! !##################################################################
! ! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
! !##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_delta_normal_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  real(8)                        ::  chi2
  complex(8),dimension(2,Ldelta) ::  Delta
  !
  Delta(:,:) = delta_normal_superc(a)
  !
  chi2 = sum( abs(Gdelta(1,:)-Delta(1,:))**2/Wdelta(:) ) + sum( abs(Fdelta(1,:)-Delta(2,:))**2/Wdelta(:) )
  !
end function chi2_delta_normal_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson 
! function in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_superc(a) result(dchi2)
  real(8),dimension(:)                   ::  a
  real(8),dimension(size(a))             ::  dchi2
  real(8),dimension(size(a))             ::  df
  complex(8),dimension(2,Ldelta)         ::  Delta
  complex(8),dimension(2,Ldelta,size(a)) ::  dDelta
  integer                                ::  j
  !
  Delta(:,:)    = delta_normal_superc(a)
  dDelta(:,:,:) = grad_delta_normal_superc(a)
  !
  do j=1,size(a)
     df(j) = &
          sum( dreal(Gdelta(1,:)-Delta(1,:))*dreal(dDelta(1,:,j))/Wdelta(:) ) + &
          sum( dimag(Gdelta(1,:)-Delta(1,:))*dimag(dDelta(1,:,j))/Wdelta(:) ) + &
          sum( dreal(Fdelta(1,:)-Delta(2,:))*dreal(dDelta(2,:,j))/Wdelta(:) ) + &
          sum( dimag(Fdelta(1,:)-Delta(2,:))*dimag(dDelta(2,:,j))/Wdelta(:) )
  enddo
  !
  dchi2=-2.d0*df
  !
end function grad_chi2_delta_normal_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance for G_0 function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_weiss_normal_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0and
  real(8)                        ::  chi2
  chi2=0d0
  !
  g0and(:,:)  = g0and_normal_superc(a)
  !
  chi2 =        sum(abs(Gdelta(1,:)-g0and(1,:))**2/Wdelta(:))
  chi2 = chi2 + sum(abs(Fdelta(1,:)-g0and(2,:))**2/Wdelta(:))
  !
end function chi2_weiss_normal_superc








!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_normal_superc(a) result(Delta)
  real(8),dimension(:)            :: a
  complex(8),dimension(2,Ldelta)  :: Delta
  integer                         :: i,k,io,stride
  real(8),dimension(Nbath)        :: eps,vps,dps
  real(8),dimension(Nbath)        :: Den
  !
  !\Delta_{aa} = - \sum_k [ V_{a}(k) * V_{a}(k) * (iw_n + E_{a}(k)) / Den(k) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io)
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     dps(i) = a(io) 
  enddo
  stride = 2*Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  do i=1,Ldelta
     Delta(1,i) = -sum( vps(:)*vps(:)*( xi*Xdelta(i) + eps(:) )/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
     Delta(2,i) =  sum( dps(:)*vps(:)*vps(:)/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
  enddo
  !
end function delta_normal_superc

function grad_delta_normal_superc(a) result(dDelta)
  real(8),dimension(:)                   :: a
  complex(8),dimension(2,Ldelta,size(a)) :: dDelta
  integer                                :: i,k,ik,io,stride
  real(8),dimension(Nbath)               :: eps,vps,dps
  real(8),dimension(Ldelta,Nbath)        :: Den
  !
  !\grad_{E_{a}(k)} \Delta_{bb} = -V_{a}(k)*V_{a}(k)*[ 1/den(k) - 2*E_{a}(k)*(iw_n + E_{a}(k))/den(k)**2 ]
  !
  !\grad_{\D_{a}(k)} \Delta_{bb} = V_{a}(k)*V_{a}(k)*\D_{a}(k)*(iw_n + E_{a}(k)) /den(k)**2
  !
  !\grad_{ V_{a}(k)} \Delta_{bb} =  2*V_{a}(k)*(iw_n + E_{a}(k))/den(k)
  !
  !
  !
  !\grad_{E_{a}(k)} \FDelta_{aa} = -2 * V_{a}(k) * V_{a}(k) * E_{a}(k) * \Delta_{a}(k) / Den**2
  !
  !\grad_{\Delta_{a}(k)} \FDelta_{aa} = V_{a}(k) * V_{a}(k) * [ 1/den - 2* \Delta_{a}(k)*\Delta_{a}(k)/den**2 ]
  !
  !\grad_{ V_{a}(k)} \FDelta_{aa} =  2 * V_{a}(k) * \Delta_{a}(k) / den
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io)
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     dps(i) = a(io) 
  enddo
  stride = 2*Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  !Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
  forall(i=1:Ldelta,k=1:Nbath)Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = -vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*eps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2)
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = 2d0*vps(k)*vps(k)*dps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2
  enddo
  stride = 2*Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = -2d0*vps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)
  enddo
  !
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = -2d0*vps(k)*vps(k)*eps(k)*dps(k)/Den(:,k)**2
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*dps(k)*dps(k)/Den(:,k)**2)
  enddo
  stride = 2*Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = 2d0*vps(k)*dps(k)/Den(:,k)
  enddo
  !
end function grad_delta_normal_superc

function g0and_normal_superc(a) result(G0and)
  real(8),dimension(:)            :: a
  complex(8),dimension(2,Ldelta)  :: G0and,Delta
  real(8),dimension(Ldelta)       :: det
  complex(8),dimension(Ldelta)    :: fg,ff
  integer                         :: iorb,ispin
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta    = delta_normal_superc(a)
  !
  fg(:)    = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(1,:)
  ff(:)    =                                                     -  Delta(2,:)
  det(:)   = abs(fg(:))**2 + ff(:)**2
  G0and(1,:) = conjg(fg(:))/det(:)
  G0and(2,:) = ff(:)/det(:)
  !
end function g0and_normal_superc


! !<DEBUG 
! !+-------------------------------------------------------------+
! !PURPOSE: Evaluate the \chi^2 distance of G_0 function 
! !         in the SUPERCONDUCTING case.
! ! The gradient is not evaluated.
! ! NORMAL bath.
! ! SPIN & ORBITAL DIAGONAL
! !+-------------------------------------------------------------+
! function chi2_weiss_normal_superc(a) result(chi2)
!   real(8),dimension(:)           ::  a
!   complex(8),dimension(2,Ldelta) ::  g0
!   real(8)                        ::  chi2
!   integer                        ::  i
!   do i=1,Ldelta
!      g0(:,i)   = fg_weiss_irred_sc(xdelta(i),a)
!   enddo
!   chi2 =        sum(abs(Gdelta(1,:)-g0(1,:))**2/Wdelta(:))
!   chi2 = chi2 + sum(abs(Fdelta(1,:)-g0(2,:))**2/Wdelta(:))
! end function chi2_weiss_normal_superc
! !
! ! the \Delta_Anderson function used in \chi^2 and d\chi^2
! ! \Delta = \sum_l V_l^2/(iw-e_l)
! !+-------------------------------------------------------------+
! function fg_delta_irred_sc(w,a) result(gg)
!   real(8)                      :: w
!   real(8),dimension(:)         :: a
!   real(8),dimension(size(a)/3) :: eps,vps,dps
!   complex(8)                   :: gg(2),delta(2),x
!   integer                      :: i,Nb
!   Nb=size(a)/3
!   eps=a(1:Nb)
!   dps=a(Nb+1:2*Nb)
!   vps=a(2*Nb+1:3*Nb)
!   gg = zero
!   x = xi*w
!   gg(1) = -sum(vps(:)**2*(x+eps(:))/(dimag(x)**2+eps(:)**2+dps(:)**2))
!   gg(2) =  sum(dps(:)*vps(:)**2/(dimag(x)**2+eps(:)**2+dps(:)**2))
! end function fg_delta_irred_sc
! !
! ! the non interacting GF (~ inverse weiss field) 
! !+-------------------------------------------------------------+
! function fg_weiss_irred_sc(w,a) result(gg)
!   real(8)                      :: w
!   real(8),dimension(:)         :: a
!   complex(8)                   :: gg(2),g0(2),delta(2),det
!   integer                      :: i,iorb,ispin
!   ispin = Spin_indx
!   iorb  = Orb_indx
!   delta = fg_delta_irred_sc(w,a)
!   g0(1) = xi*w + xmu - impHloc(ispin,ispin,iorb,iorb) - delta(1)
!   g0(2) = -delta(2)
!   det   = abs(g0(1))**2 + (g0(2))**2
!   gg(1) = conjg(g0(1))/det
!   gg(2) = g0(2)/det
! end function fg_weiss_irred_sc
! !>DEBUG



