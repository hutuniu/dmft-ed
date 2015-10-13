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
subroutine chi2_fitgf_hybrid_superc(fg,bath_,ispin)
  complex(8),dimension(:,:,:,:,:,:)    :: fg ![2][Nspin][Nspin][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)   :: bath_
  integer                              :: ispin
  real(8),dimension(:),allocatable     :: array_bath
  integer                              :: iter,stride,i,j,io,corb,l,Asize
  integer                              :: iorb,jorb
  real(8)                              :: chi
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  character(len=20)                    :: suffix
  integer                              :: unit
  !
  if(size(fg,1)/=2)stop"chi2_fitgf_hybrid_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Nspin)stop"chi2_fitgf_hybrid_superc error: size[fg,2]!=Nspin"
  if(size(fg,3)/=Nspin)stop"chi2_fitgf_hybrid_superc error: size[fg,3]!=Nspin"
  if(size(fg,4)/=Norb)stop"chi2_fitgf_hybrid_superc error: size[fg,4]!=Norb"
  if(size(fg,5)/=Norb)stop"chi2_fitgf_hybrid_superc error: size[fg,5]!=Norb"
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
     Gdelta(i,1:Ldelta) = fg(1,ispin,ispin,getIorb(i),getJorb(i),1:Ldelta)
     Fdelta(i,1:Ldelta) = fg(2,ispin,ispin,getIorb(i),getJorb(i),1:Ldelta)
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
  Spin_indx=ispin 
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !D_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
  Asize = Nbath + Nbath + Norb*Nbath
  allocate(array_bath(Asize))
  !
  !Nbath + Nbath + Norb*Nbath
  stride = 0
  do i=1,Nbath 
     io = stride + i
     array_bath(io)   = dmft_bath%e(ispin,1,i)
  enddo
  stride = Nbath
  do i=1,Nbath 
     io = stride + i
     array_bath(io)    = dmft_bath%d(ispin,1,i)
  enddo
  stride = Nbath + Nbath
  do iorb=1,Norb
     do i=1,Nbath
        io = stride + i + (iorb-1)*Nbath
        array_bath(io) = dmft_bath%v(ispin,iorb,i)
     enddo
  enddo
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,chi2_weiss_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case ("delta")
        call fmin_cg(array_bath,chi2_delta_hybrid_superc,grad_chi2_delta_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case default
        stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
     end select
     !
  case (2)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgplus(array_bath,chi2_weiss_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgplus(array_bath,chi2_delta_hybrid_superc,grad_chi2_delta_hybrid_superc,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  write(LOGfile,"(A,ES18.9,A,I5)")&
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
  stride = 0
  do i=1,Nbath 
     io = stride + i
     dmft_bath%e(ispin,1,i) = array_bath(io)   
  enddo
  stride = Nbath
  do i=1,Nbath 
     io = stride + i
     dmft_bath%d(ispin,1,i)   = array_bath(io)
  enddo
  stride = 2*Nbath
  do iorb=1,Norb
     do i=1,Nbath
        io = stride + i + (iorb-1)*Nbath
        dmft_bath%v(ispin,iorb,i) = array_bath(io)
     enddo
  enddo
  !
  if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  if(ed_verbose<3)call write_fit_result(ispin)
  !
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
  deallocate(getIorb,getJorb)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer                                :: i,j,l,m,iorb,jorb,ispin,jspin
    integer                                :: gunit,funit
    real(8)                                :: w
    complex(8),dimension(Norb,Norb,Ldelta) :: fgand,ffand
    if(cg_scheme=='weiss')then
       fgand(:,:,:) = g0and_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
       ffand(:,:,:) = f0and_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
    else
       fgand(:,:,:) = delta_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
       ffand(:,:,:) = fdelta_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
    endif
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
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
          write(gunit,"(10F24.15)")Xdelta(i),&
               dimag(Gdelta(l,i)),dimag(fgand(iorb,jorb,i)),&
               dreal(Gdelta(l,i)),dreal(fgand(iorb,jorb,i))
          write(funit,"(10F24.15)")Xdelta(i),&
               dimag(Fdelta(l,i)),dimag(ffand(iorb,jorb,i)),&
               dreal(Fdelta(l,i)),dreal(ffand(iorb,jorb,i))
       enddo
       close(gunit)
       close(funit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_hybrid_superc





!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_delta_hybrid_superc(a) result(chi2)
  real(8),dimension(:)                     :: a
  real(8)                                  :: chi2
  real(8),dimension(totNorb)               :: chi_orb
  complex(8),dimension(2,Norb,Norb,Ldelta) :: Delta
  integer                                  ::  i,l,iorb,jorb
  !  
  Delta = delta_hybrid_superc(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     chi_orb(l) = &
          sum(abs(Gdelta(l,:)-Delta(1,iorb,jorb,:))**2/Wdelta(:)) + &
          sum(abs(Fdelta(l,:)-Delta(2,iorb,jorb,:))**2/Wdelta(:))
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
  real(8),dimension(:)                             ::  a
  real(8),dimension(size(a))                       ::  dchi2
  real(8),dimension(totNorb,size(a))               ::  df
  complex(8),dimension(2,Norb,Norb,Ldelta)         ::  Delta
  complex(8),dimension(2,Norb,Norb,Ldelta,size(a)) ::  dDelta
  integer                                          ::  i,j,l,iorb,jorb
  !
  Delta  = delta_hybrid_superc(a)
  dDelta = grad_delta_hybrid_superc(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     !
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Gdelta(l,:)-Delta(1,iorb,jorb,:))*dreal(dDelta(1,iorb,jorb,:,j))/Wdelta(:) ) + &
             sum( dimag(Gdelta(l,:)-Delta(1,iorb,jorb,:))*dimag(dDelta(1,iorb,jorb,:,j))/Wdelta(:) ) + &
             sum( dreal(Fdelta(l,:)-Delta(2,iorb,jorb,:))*dreal(dDelta(2,iorb,jorb,:,j))/Wdelta(:) ) + &
             sum( dimag(Fdelta(l,:)-Delta(2,iorb,jorb,:))*dimag(dDelta(2,iorb,jorb,:,j))/Wdelta(:) )
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
  complex(8),dimension(2,Norb,Norb,Ldelta) ::  g0and
  real(8)                                  ::  chi2
  integer                                  ::  i,l,iorb,jorb
  !
  g0and = g0and_hybrid_superc(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     chi2_orb(l) = &
          sum(abs(Gdelta(l,:)-g0and(1,iorb,jorb,:))**2/Wdelta(:))  + &
          sum(abs(Fdelta(l,:)-g0and(2,iorb,jorb,:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi2_orb)
  !
end function chi2_weiss_hybrid_superc







!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_hybrid_superc(a) result(Delta)
  real(8),dimension(:)                     :: a
  complex(8),dimension(2,Norb,Norb,Ldelta) :: Delta
  integer                                  :: iorb,jorb
  integer                                  :: i,io,k,l,stride
  real(8),dimension(Nbath)                 :: eps,dps
  real(8),dimension(Norb,Nbath)            :: vops
  real(8),dimension(Ldelta,Nbath)          :: Den
  !
  !\Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (iw_n + E(k)) / Den(k) ]
  !
  !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
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
  do l=1,Norb
     do i=1,Nbath
        io = stride + i + (l-1)*Nbath
        vops(l,i) = a(io)
     enddo
  enddo
  !
  forall(i=1:Ldelta,k=1:Nbath)& !den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
       Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
  !
  do i=1,Ldelta
     do iorb=1,Norb
        Delta(1,iorb,iorb,i) = -sum( vops(iorb,:)*vops(iorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
        !
        Delta(2,iorb,iorb,i) = -sum( dps(:)*vops(iorb,:)*vops(iorb,:)/Den(i,:) )
        do jorb=iorb+1,Norb
           Delta(1,iorb,jorb,i) = -sum( vops(iorb,:)*vops(jorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
           Delta(1,jorb,iorb,i) = -sum( vops(jorb,:)*vops(iorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
           !
           Delta(2,iorb,jorb,i) = -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/Den(i,:) )
           Delta(2,jorb,iorb,i) = -sum( dps(:)*vops(jorb,:)*vops(iorb,:)/Den(i,:) )
        enddo
     enddo
  enddo
  !
end function delta_hybrid_superc

function grad_delta_hybrid_superc(a) result(dDelta)
  real(8),dimension(:)                             :: a
  complex(8),dimension(2,Norb,Norb,Ldelta,size(a)) :: dDelta
  integer                                          :: iorb,jorb
  integer                                          :: i,k,io,ik,l,stride
  real(8),dimension(Nbath)                         :: eps,dps
  real(8),dimension(Norb,Nbath)                    :: vops
  real(8),dimension(Ldelta,Nbath)                  :: Den
  real(8),dimension(Norb,Norb)                     :: delta_orb
  !
  !\grad_{E_{1}(k)} \Delta_{ab} = -V_{a}(k)*V_{b}(k)*[ 1/den(k) - 2*E(k)*(iw_n + E(k))/den(k)**2 ]
  !
  !\grad_{\D_{1}(k)} \Delta_{ab} = V_{a}(k)*V_{b}(k)*\D(k)*(iw_n + E(k)) /den(k)**2
  !
  !\grad_{ V_{g}(k)} \Delta_{ab} =  [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ]*(iw_n + E_{a}(k))/den(k)
  !
  !
  !\grad_{E(k)} \FDelta_{ab} = -2 * V_{a}(k) * V_{b}(k) * E(k) * \D(k) / Den**2
  !
  !\grad_{\D(k)} \FDelta_{ab} = V_{a}(k) * V_{b}(k) * [ 1/den - 2* \D(k)*\D(k)/den**2 ]
  !
  !\grad_{ V_{g}(k)} \FDelta_{ab} =  [delta(g,a)*V_{b}(k) + delta(g,b)*V_{a}(k)] * \D(k) / den
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
  do l=1,Norb
     do i=1,Nbath
        io = stride + i + (l-1)*Nbath
        vops(l,i) = a(io)
     enddo
  enddo
  !
  forall(i=1:Ldelta,k=1:Nbath)& !den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
       Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
  !
  delta_orb = zeye(Norb)
  !
  do iorb=1,Norb
     do jorb=1,Norb
        !
        stride=0
        do k=1,Nbath
           ik = stride + k
           dDelta(1,iorb,jorb,:,ik) = -vops(iorb,k)*vops(jorb,k)*(1d0/Den(:,k) - 2d0*eps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2)
           dDelta(2,iorb,jorb,:,ik) = -2d0*vops(iorb,k)*vops(jorb,k)*eps(k)*dps(k)/Den(:,k)**2
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dDelta(1,iorb,jorb,:,ik) = 2d0*vops(iorb,k)*vops(jorb,k)*dps(k)*(xi*Xdelta(:) + eps(k))/Den(:,k)**2
           dDelta(2,iorb,jorb,:,ik) = vops(iorb,k)*vops(jorb,k)*( 1d0/Den(:,k) - 2d0*dps(k)*dps(k)/Den(:,k)**2 )
        enddo
        stride=2*Nbath
        do l=1,Norb
           do k=1,Nbath
              ik = stride + k + (l-1)*Nbath
              dDelta(1,iorb,jorb,:,ik) = (delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))*(xi*Xdelta(:) + eps(k))/Den(:,k)
              dDelta(2,iorb,jorb,:,ik) = (delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))*dps(k)/Den(:,k)
           enddo
        enddo
        !
     enddo
  enddo
  !
end function grad_delta_hybrid_superc

function g0and_hybrid_superc(a) result(G0and)
  real(8),dimension(:)                     :: a
  complex(8),dimension(2,Norb,Norb,Ldelta) :: G0and,Delta
  complex(8),dimension(2,Norb,Norb)        :: zeta
  complex(8),dimension(2*Norb,2*Norb)      :: fgorb
  integer                                  :: i,k,ispin
  !
  ispin  = Spin_indx
  !
  Delta = delta_hybrid_superc(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta(1,:,:) = (xi*Xdelta(i)+xmu)*zeye(Norb)
     zeta(2,:,:) = (xi*Xdelta(i)-xmu)*zeye(Norb)
     fgorb(1:Norb,1:Norb)                     = zeta(1,:,:) - impHloc(ispin,ispin,:,:)  - Delta(1,:,:,i)
     fgorb(1:Norb,Norb+1:Norb+Norb)           =                                         - Delta(2,:,:,i)
     fgorb(Norb+1:Norb+Norb,1:Norb)           =                                         - Delta(2,:,:,i)
     fgorb(Norb+1:Norb+Norb,Norb+1:Norb+Norb) = zeta(2,:,:) + impHloc(ispin,ispin,:,:)  + conjg( Delta(1,:,:,i) )
     call inv(fgorb)
     G0and(1,:,:,i) = fgorb(1:Norb,1:Norb)
     G0and(2,:,:,i) = fgorb(1:Norb,Norb+1:Norb+Norb)
  enddo
  !
end function g0and_hybrid_superc
