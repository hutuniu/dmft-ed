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
!PURPOSE  : Chi^2 interface for
!+-------------------------------------------------------------+
subroutine chi2_fitgf_replica(fg,bath_)
  complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: i,j,iorb,jorb,ispin,jspin,io,jo
  integer                            :: iter,stride,count,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=20)                  :: suffix
  integer                            :: unit
  !
  if(size(fg,1)/=Nspin)stop "chi2_fitgf_replica error: size[fg,1]!=Nspin"
  if(size(fg,2)/=Nspin)stop "chi2_fitgf_replica error: size[fg,2]!=Nspin"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_replica error: size[fg,3]!=Norb"
  if(size(fg,4)/=Norb)stop "chi2_fitgf_replica error: size[fg,4]!=Norb"
  !
  check= check_bath_dimension(bath_,impHloc)
  if(.not.check)stop "chi2_fitgf_replica error: wrong bath dimensions"
  !
  call allocate_dmft_bath(dmft_bath)
  call init_dmft_bath_mask(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  array_bath=bath_
  !
  count=0
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              if ((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv..true.).or.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv..true.))then
                 count=count+1
                 getIspin(count) = ispin
                 getIorb(count)  = iorb
                 getJspin(count) = jspin
                 getJorb(count)  = jorb
              endif
           enddo
        enddo
     enddo
  enddo
  totNso=count
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
  !
  allocate(Gdelta(totNso,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  allocate(array_bath(size(bath_)))
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
  do i=1,totNso
     Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,chi2_weiss_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case ("delta")
        call fmin_cg(array_bath,chi2_delta_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case default
        stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
     end select
     !
  case (2)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgplus(array_bath,chi2_weiss_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgplus(array_bath,chi2_delta_replica,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  write(LOGfile,"(A,ES18.9,A,I5,A)")&
       "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
       "  <--  All Orbs, All Spins"
  !
  if(ed_verbose<2)then
     suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
  endif
  !
  call set_dmft_bath(array_bath,dmft_bath)           ! *** array_bath --> dmft_bath ***    (per write fit result)
  !
  if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  !if(ed_verbose<3)call write_fit_result(ispin)
  !
  call get_dmft_bath(dmft_bath,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
  !bath_=array_bath                                  ! *** array_bath --> bath_ ***    (analogo della riga sopra)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(Nspin,Nspin,Norb,Norb,Ldelta)
    integer           :: i,j,s,l,iorb,jorb,ispin,jspin
    real(8)           :: w
    if(cg_scheme=='weiss')then
       fgand(:,:,:,:,:) = g0and_bath_mats(xi*Xdelta(:),dmft_bath)
    else
       fgand(:,:,:,:,:) = delta_bath_mats(xi*Xdelta(:),dmft_bath)
    endif
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       suffix="_l"//reg(txtfy(iorb))//&
            "_m"//reg(txtfy(jorb))//&
            "_s"//reg(txtfy(ispin))//&
            "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          w = Xdelta(i)
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
               dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result

end subroutine chi2_fitgf_replica


!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for
!+-------------------------------------------------------------+
subroutine chi2_fitgf_replica_normal(fg,bath_,ispin_)
  complex(8),dimension(:,:,:)        :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: i,j,iorb,jorb,ispin,jspin,io,jo
  integer                            :: iter,stride,count,Asize,ispin_
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=20)                  :: suffix
  integer                            :: unit
  !
  if(size(fg,1)/=Norb)stop "chi2_fitgf_replica_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_replica_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_,impHloc)
  if(.not.check)stop "chi2_fitgf_replica_normal error: wrong bath dimensions"
  !
  call allocate_dmft_bath(dmft_bath)
  call init_dmft_bath_mask(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  count=0
  Asize=0
  do iorb=1,Norb
     do jorb=1,Norb
        if ((dmft_bath%mask(ispin_,ispin_,iorb,jorb,1).eqv..true.).or.(dmft_bath%mask(ispin_,ispin_,iorb,jorb,2).eqv..true.))then
           count=count+1
           getIorb(count)  = iorb
           getJorb(count)  = jorb
           if (dmft_bath%mask(ispin_,ispin_,iorb,jorb,1).eqv..true.)Asize=Asize+1
           if (dmft_bath%mask(ispin_,ispin_,iorb,jorb,2).eqv..true.)Asize=Asize+1
        endif
     enddo
  enddo
  totNso=count
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(Gdelta(totNso,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  allocate(array_bath(Asize))
  !
  if(mod(size(bath_),2).ne.0)stop"something wrong bath size is not even"
  if(Asize.ne.int(size(bath_)/2))stop"something wrong Asize is not equal to size(bath_)/2"
  !
  if (ispin_==1) then
     array_bath=bath_(1:Asize)
     Spin_mask=1
  elseif (ispin_==2) then
     array_bath=bath_(1+Asize:2*Asize)
     Spin_mask=2
  endif
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
  do i=1,totNso
     Gdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,chi2_weiss_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case ("delta")
        call fmin_cg(array_bath,chi2_delta_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     case default
        stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  case (2)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgplus(array_bath,chi2_weiss_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case ("delta")
        call fmin_cgplus(array_bath,chi2_delta_replica_normal,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     case default
        stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  write(LOGfile,"(A,ES18.9,A,I5,A)")&
       "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
       "  <--  All Orbs, All Spins"
  !
  if(ed_verbose<2)then
     suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
  endif
  !
  bath_=0.0d0
  bath_(1:Asize)=array_bath
  bath_(1+Asize:2*Asize)=array_bath
  !
  call set_dmft_bath(bath_,dmft_bath)
  !
  if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  !if(ed_verbose<3)call write_fit_result(ispin)
  !
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer                                :: i,j,l,m,iorb,jorb,ispin,jspin
    real(8)                                :: w
    complex(8),dimension(Norb,Norb,Ldelta) :: fgand

    if(cg_scheme=='weiss')then
       fgand(:,:,:) = g0and_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
    else
       fgand(:,:,:) = delta_bath_mats(ispin,ispin,xi*Xdelta(:),dmft_bath)
    endif
    !

    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(Gdelta(l,i)),dimag(fgand(iorb,jorb,i)),&
               dreal(Gdelta(l,i)),dreal(fgand(iorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result

end subroutine chi2_fitgf_replica_normal



!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_replica(a) result(chi2)
  real(8),dimension(:)                               ::  a
  real(8)                                            ::  chi2
  real(8),dimension(totNso)                          ::  chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) ::  Delta
  integer                                            ::  i,l,iorb,jorb,ispin,jspin
  !
  Delta = delta_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     chi2_so(l) = sum( abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))**2/Wdelta(:) )
  enddo
  !
  chi2=sum(chi2_so)
  !
end function chi2_delta_replica

function chi2_delta_replica_normal(a) result(chi2)
  real(8),dimension(:)                               ::  a
  real(8)                                            ::  chi2
  real(8),dimension(totNso)                          ::  chi2_so
  complex(8),dimension(Norb,Norb,Ldelta)             ::  Delta
  integer                                            ::  i,l,iorb,jorb
  !
  Delta = delta_replica_normal(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     chi2_so(l) = sum( abs(Gdelta(l,:)-Delta(iorb,jorb,:))**2/Wdelta(:) )
  enddo
  !
  chi2=sum(chi2_so)
  !
end function chi2_delta_replica_normal


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_replica(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(totNso)                          :: chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
  integer                                            :: i,l,iorb,jorb,ispin,jspin
  !
  g0and = g0and_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     chi2_so(l) = sum( abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))**2/Wdelta(:) )
  enddo
  !
  chi2=sum(chi2_so)
  !
end function chi2_weiss_replica

function chi2_weiss_replica_normal(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(totNso)                          :: chi2_so
  complex(8),dimension(Norb,Norb,Ldelta)             :: g0and
  integer                                            :: i,l,iorb,jorb
  !
  g0and = g0and_replica_normal(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     chi2_so(l) = sum( abs(Gdelta(l,:)-g0and(iorb,jorb,:))**2/Wdelta(:) )
  enddo
  !
  chi2=sum(chi2_so)
  !
end function chi2_weiss_replica_normal






!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_replica(a) result(Delta)
  real(8),dimension(:)                                :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)  :: Delta
  integer                                             :: ispin,jspin,iorb,jorb,ibath
  integer                                             :: i,io,jo,ndx
  real(8),dimension(Nspin*Norb,Nspin*Norb)            :: V_k
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)  :: invH_k
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)  :: Delta_so
  type(effective_bath)                                :: dmft_bath_tmp

  call allocate_dmft_bath(dmft_bath_tmp)
  call set_dmft_bath(a,dmft_bath_tmp) 

  Delta=zero
  Delta_so=zero
  invH_k=zero
  do i=1,Ldelta
     !
     ndx=0
     do ibath=1,Nbath
        !
        V_k=0.0d0
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1) * Norb
                    jo = jorb + (jspin-1) * Norb
                    ndx=ndx+1
                    if ((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv..true.).or.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv..true.))then
                       invH_k(io,jo,i)=a(ndx)
                       V_k(io,io)=a(ndx)
                    endif
                 enddo
              enddo
           enddo
           !
           invH_k(:,:,i) = eye(Nspin*Norb) * xi * Xdelta(i) - invH_k(:,:,i)
           call inv(invH_k(:,:,i))
           !
           Delta_so(:,:,i)=Delta_so(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
           !
        enddo
        !
        Delta(:,:,:,:,i)=so2nn_reshape(Delta_so(:,:,i),Nspin,Norb)
        !
     enddo
  enddo
  call deallocate_dmft_bath(dmft_bath_tmp)
  !
end function delta_replica

function delta_replica_normal(a) result(Delta)
  real(8),dimension(:)                                :: a
  complex(8),dimension(Norb,Norb,Ldelta)              :: Delta
  integer                                             :: iorb,jorb,ibath
  integer                                             :: i,io,jo,ndx
  real(8),dimension(Norb,Norb)                        :: V_k
  complex(8),dimension(Norb,Norb,Ldelta)              :: invH_k
  type(effective_bath)                                :: dmft_bath_tmp

  call allocate_dmft_bath(dmft_bath_tmp)
  call set_dmft_bath(a,dmft_bath_tmp) 

  Delta=zero
  invH_k=zero
  do i=1,Ldelta
     !
     ndx=0
     do ibath=1,Nbath
        !
        V_k=0.0d0
        do iorb=1,Norb
           do jorb=1,Norb
              ndx=ndx+1
              if ((dmft_bath%mask(Spin_mask,Spin_mask,iorb,jorb,1).eqv..true.).or.(dmft_bath%mask(Spin_mask,Spin_mask,iorb,jorb,2).eqv..true.))then
                 V_k(iorb,iorb)=a(ndx)
                 invH_k(iorb,jorb,i)=a(ndx)
              endif
           enddo
        enddo
        !
        invH_k(:,:,i) = eye(Norb) * xi * Xdelta(i) - invH_k(:,:,i)
        call inv(invH_k(:,:,i))
        !
        Delta(:,:,i)=Delta(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
        !
     enddo
  enddo
  call deallocate_dmft_bath(dmft_bath_tmp)
  !
end function delta_replica_normal






function g0and_replica(a) result(G0and)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
  complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: zeta,fgorb
  integer                                            :: i,Nso
  integer                                            :: iorb,jorb,ispin,jspin,io,jo
  !
  Nso = Norb*Nspin
  !
  Delta = delta_replica(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 fgorb(io,jo)   = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
              enddo
           enddo
        enddo
     enddo
     call inv(fgorb)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
              enddo
           enddo
        enddo
     enddo
  enddo
  !
end function g0and_replica


function g0and_replica_normal(a) result(G0and)
  real(8),dimension(:)                        :: a
  complex(8),dimension(Norb,Norb,Ldelta)      :: G0and,Delta
  complex(8),dimension(Norb,Norb)             :: zeta,fgorb
  integer                                     :: i,Nso
  integer                                     :: iorb,jorb
  !
  Nso = Norb
  !
  Delta = delta_replica_normal(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
     do iorb=1,Norb
        do jorb=1,Norb
           fgorb(iorb,jorb)   = zeta(iorb,jorb) - impHloc(Spin_mask,Spin_mask,iorb,jorb) - Delta(iorb,jorb,i)
        enddo
     enddo
     call inv(fgorb)
     !
     G0and(:,:,i) = fgorb
     !
  enddo
  !
end function g0and_replica_normal

