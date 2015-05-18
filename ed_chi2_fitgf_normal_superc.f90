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
  if(size(fg,1)/=2)stop"chi2_fitgf_normal_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
  if(size(fg,3)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
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
     wdelta=Xdelta
  end select
  !
  !
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
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
           call fmin_cg(array_bath,&
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,&
                chi2_delta_normal_superc,&
                grad_chi2_delta_normal_superc,&
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
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,&
                chi2_delta_normal_superc,&
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
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,&
                chi2_delta_normal_superc,&
                grad_chi2_delta_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)")&
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
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  !
  call save_bath(dmft_bath)
  !
  if(ed_verbose<3)call write_fit_result(ispin)
  call copy_bath(dmft_bath,bath_)
  call deallocate_bath(dmft_bath)
  deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(2),det
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Gdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
       Fdelta(1,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
       fgand=zero
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          w = Xdelta(i)
          if(cg_scheme=='weiss')then
             fgand(1) = g0and_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
             fgand(2) = f0and_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          else
             fgand(1) = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
             fgand(2) = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          endif
          write(unit,"(10F24.15)")Xdelta(i),dimag(Gdelta(1,i)),dimag(fgand(1)),dreal(Gdelta(1,i)),dreal(fgand(1)), &
               dimag(Fdelta(1,i)),dimag(fgand(2)),dreal(Fdelta(1,i)),dreal(fgand(2))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_superc







!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_delta_normal_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  real(8)                        ::  chi2
  complex(8),dimension(2,Ldelta) ::  Delta
  integer                        ::  i
  !
  Delta(:,:) = delta_normal_superc(a)
  !
  chi2=sum( abs(Gdelta(1,:)-Delta(1,:))**2/Wdelta(:) ) + sum( abs(Fdelta(1,:)-Delta(2,:))**2/Wdelta(:) )
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
  integer                                ::  i,j
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
     Den(:)     =  Xdelta(i)**2 + eps(:)**2 + dps(:)**2
     Delta(1,i) = -sum( vps(:)*vps(:)*( xi*Xdelta(i) + eps(:) )/den(:) )
     Delta(2,i) =  sum( dps(:)*vps(:)*vps(:)/den(:) )
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
  complex(8),dimension(Ldelta)    :: fg,ff,det
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
