!NORMAL:
function grad_delta_bath_mats_main(x,dmft_bath_,N) result(dDelta)
  complex(8),intent(in)                         :: x
  type(effective_bath)                          :: dmft_bath_
  integer,intent(in)                            :: N
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dDelta
  integer                                       :: iorb,jorb,ispin,jspin
  real(8),dimension(Nbath)                      :: eps,dps,vps,den
  real(8),dimension(Norb,Nbath)                 :: vops
  real(8),dimension(Norb,Norb)                  :: delta_orb
  integer                                       :: ik,k,stride,gorb
  !
  dDelta=zero
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              !
              !\grad_{E_{a}(k)} \Delta_{bb}^{rr} = [ V_{a}(k)*V_{a}(k) / ( iw_n - E_{a}(k) )**2 ]
              !
              !\grad_{V_{a}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}(k) / ( iw_n - E_{a}(k) ) ]
              !
              stride=0
              do k=1,Nbath
                 ik = stride + k
                 dDelta(ispin,ispin,iorb,iorb,ik) = vps(k)*vps(k)/(x - eps(k))**2
              enddo
              stride=Nbath
              do k=1,Nbath
                 ik = stride + k
                 dDelta(ispin,ispin,iorb,iorb,ik) = 2d0*vps(k)/(x - eps(k))
              enddo
           enddo
        enddo
        !
     case ("superc")
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              ! Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
              forall(k=1:Nbath)den(k) = dimag(x)**2 + eps(k)**2 + dps(k)**2
              !
              !\grad_{E_{a}(k)} \Delta_{bb} = -V_{a}(k)*V_{a}(k)*[ 1/den(k) - 2*E_{a}(k)*(iw_n + E_{a}(k))/den(k)**2 ]
              !
              !\grad_{\D_{a}(k)} \Delta_{bb} = V_{a}(k)*V_{a}(k)*\D_{a}(k)*(iw_n + E_{a}(k)) /den(k)**2
              !
              !\grad_{ V_{a}(k)} \Delta_{bb} =  2*V_{a}(k)*(iw_n + E_{a}(k))/den(k)
              !
              stride=0
              do k=1,Nbath
                 ik = stride + k
                 dDelta(ispin,ispin,iorb,iorb,ik) = -vps(k)*vps(k)*(1d0/den(k) - 2d0*eps(k)*( x + eps(k))/den(k)**2)
              enddo
              stride=Nbath
              do k=1,Nbath
                 ik = stride + k
                 dDelta(ispin,ispin,iorb,iorb,ik) = 2d0*vps(k)*vps(k)*dps(k)*(x + eps(k))/den(k)**2
              enddo
              stride=2*Nbath
              do k=1,Nbath
                 ik = stride + k
                 dDelta(ispin,ispin,iorb,iorb,ik) = 2d0*vps(k)*(x + eps(k))/den(k)
              enddo
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           delta_orb = eye(Norb)
           !
           !\grad_{E_{1}(k)} \Delta_{ab} = [ V_{a}(k)*V_{b}(k) / ( iw_n - E_{1}(k) )**2 ]
           !
           !\grad_{V_{g}(k)} \Delta_{ab} = [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ] / ( iw_n - E_{1}(k) )
           !
           do iorb=1,Norb
              do jorb=1,Norb
                 stride=0
                 do k=1,Nbath
                    ik = stride + k
                    dDelta(ispin,ispin,iorb,jorb,ik) = vops(iorb,k)*vops(jorb,k)/(x - eps(k))**2
                 enddo
                 stride=Nbath
                 do gorb=1,Norb
                    do k=1,Nbath
                       ik = stride + k + (gorb-1)*Nbath
                       dDelta(ispin,ispin,iorb,jorb,ik) = (delta_orb(gorb,iorb)*vops(jorb,k) + delta_orb(gorb,jorb)*vops(iorb,k))/(x - eps(k))
                    enddo
                 enddo
                 !
              enddo
           enddo
        enddo
        !
     case ("superc")
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           delta_orb = eye(Norb)
           ! Den(k) = (w_n**2 + E(k)**2 + \D(k)**2
           forall(k=1:Nbath)den(k) = dimag(x)**2 + eps(k)**2 + dps(k)**2
           !
           !\grad_{E_{1}(k)} \Delta_{ab} = -V_{a}(k)*V_{b}(k)*[ 1/den(k) - 2*E(k)*(iw_n + E(k))/den(k)**2 ]
           !
           !\grad_{\D_{1}(k)} \Delta_{ab} = V_{a}(k)*V_{b}(k)*\D(k)*(iw_n + E(k)) /den(k)**2
           !
           !\grad_{ V_{g}(k)} \Delta_{ab} =  [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ]*(iw_n + E_{a}(k))/den(k)
           !
           do iorb=1,Norb
              do jorb=1,Norb
                 stride=0
                 do k=1,Nbath
                    ik = stride + k
                    dDelta(ispin,ispin,iorb,jorb,ik) = -vops(iorb,k)*vops(jorb,k)*(1d0/den(k) - 2d0*eps(k)*( x + eps(k))/den(k)**2)
                 enddo
                 stride=Nbath
                 do k=1,Nbath
                    ik = stride + k
                    dDelta(ispin,ispin,iorb,jorb,ik) = 2d0*vops(iorb,k)*vops(jorb,k)*dps(k)*(x + eps(k))/den(k)**2
                 enddo
                 stride=2*Nbath
                 do gorb=1,Norb
                    do k=1,Nbath
                       ik = stride + k + (gorb-1)*Nbath
                       dDelta(ispin,ispin,iorb,jorb,ik) = (delta_orb(gorb,iorb)*vops(jorb,k) + delta_orb(gorb,jorb)*vops(iorb,k))*(x + eps(k))/den(k)
                    enddo
                 enddo
                 !
              enddo
           enddo
        enddo
        !
     end select
  end select
end function grad_delta_bath_mats_main

function grad_delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_,N) result(dG0out)
  integer,intent(in)                            :: ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  complex(8),intent(in)                         :: x
  complex(8),dimension(Norb,Norb,N)             :: dG0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dDelta
  dDelta = grad_delta_bath_mats_main(x,dmft_bath_,N)
  dG0out = dDelta(ispin,jspin,:,:,:)
end function grad_delta_bath_mats_ispin_jspin


function grad_delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(dG0out)
  integer,intent(in)                            :: iorb,jorb,ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  complex(8),intent(in)                         :: x
  complex(8),dimension(N)                       :: dG0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dDelta
  dDelta = grad_delta_bath_mats_main(x,dmft_bath_,N)
  dG0out = dDelta(ispin,jspin,iorb,jorb,:)
end function grad_delta_bath_mats_ispin_jspin_iorb_jorb

function grad_delta_bath_mats_main_(x,bath_,N) result(dDelta)
  complex(8),intent(in)                         :: x
  real(8),dimension(:)                          :: bath_
  integer,intent(in)                            :: N
  type(effective_bath)                          :: dmft_bath_
  logical                                       :: check
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dDelta
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dDelta = grad_delta_bath_mats_main(x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_delta_bath_mats_main_

function grad_delta_bath_mats_ispin_jspin_(ispin,jspin,x,bath_,N) result(dG0out)
  complex(8),intent(in)                         :: x
  real(8),dimension(:)                          :: bath_
  integer,intent(in)                            :: ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  logical                                       :: check
  complex(8),dimension(Norb,Norb,N)             :: dG0out
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dG0out = grad_delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_delta_bath_mats_ispin_jspin_

function grad_delta_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_,N) result(dG0out)
  complex(8),intent(in)   :: x
  real(8),dimension(:)    :: bath_
  integer,intent(in)      :: iorb,jorb,ispin,jspin,N
  type(effective_bath)    :: dmft_bath_
  logical                 :: check
  complex(8),dimension(N) :: dG0out
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dG0out = grad_delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_delta_bath_mats_ispin_jspin_iorb_jorb_









!ANOMALous:
function grad_fdelta_bath_mats_main(x,dmft_bath_,N) result(dFdelta)
  complex(8),intent(in)                         :: x
  type(effective_bath)                          :: dmft_bath_
  integer,intent(in)                            :: N
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dFdelta
  integer                                       :: iorb,jorb,ispin,jspin
  real(8),dimension(Nbath)                      :: eps,dps,vps,den
  real(8),dimension(Norb,Nbath)                 :: vops
  real(8),dimension(Norb,Norb)                  :: delta_orb
  integer                                       :: ik,k,stride,gorb
  !
  dFdelta=zero
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              ! Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
              forall(k=1:Nbath)den(k) = dimag(x)**2 + eps(k)**2 + dps(k)**2
              !
              !\grad_{E_{a}(k)} \FDelta_{aa} = -2 * V_{a}(k) * V_{a}(k) * E_{a}(k) * \Delta_{a}(k) / Den**2
              !
              !\grad_{\Delta_{a}(k)} \FDelta_{aa} = V_{a}(k) * V_{a}(k) * [ 1/den - 2* \Delta_{a}(k)*\Delta_{a}(k)/den**2 ]
              !
              !\grad_{ V_{a}(k)} \FDelta_{aa} =  2 * V_{a}(k) * \Delta_{a}(k) / den
              stride=0
              do k=1,Nbath
                 ik = stride + k
                 dFdelta(ispin,ispin,iorb,iorb,ik) = -2d0*vps(k)*vps(k)*eps(k)*dps(k)/den(k)**2
              enddo
              stride=Nbath
              do k=1,Nbath
                 ik = stride + k
                 dFdelta(ispin,ispin,iorb,iorb,ik) = vps(k)*vps(k)*(1d0/den(k) - 2d0*dps(k)*dps(k)/den(k)**2)
              enddo
              stride=2*Nbath
              do k=1,Nbath
                 ik = stride + k
                 dFdelta(ispin,ispin,iorb,iorb,ik) = 2d0*vps(k)*dps(k)/den(k)
              enddo
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
        !
     case ("superc")
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           delta_orb = eye(Norb)
           ! Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
           forall(k=1:Nbath)den(k) = dimag(x)**2 + eps(k)**2 + dps(k)**2
           !
           !\grad_{E(k)} \FDelta_{ab} = -2 * V_{a}(k) * V_{b}(k) * E(k) * \D(k) / Den**2
           !
           !\grad_{\D(k)} \FDelta_{ab} = V_{a}(k) * V_{b}(k) * [ 1/den - 2* \D(k)*\D(k)/den**2 ]
           !
           !\grad_{ V_{g}(k)} \FDelta_{ab} =  [delta(g,a)*V_{b}(k) + delta(g,b)*V_{a}(k)] * \D(k) / den
           !
           do iorb=1,Norb
              do jorb=1,Norb
                 stride=0
                 do k=1,Nbath
                    ik = stride + k
                    dFdelta(ispin,ispin,iorb,jorb,ik) = -2d0*vops(iorb,k)*vops(jorb,k)*eps(k)*dps(k)/den(k)**2
                 enddo
                 stride=Nbath
                 do k=1,Nbath
                    ik = stride + k
                    dFdelta(ispin,ispin,iorb,jorb,ik) = vops(iorb,k)*vops(jorb,k)*(1d0/den(k) - 2d0*dps(k)*dps(k)/den(k)**2)
                 enddo
                 stride=2*Nbath
                 do gorb=1,Norb
                    do k=1,Nbath
                       ik = stride + k + (gorb-1)*Nbath
                       dFdelta(ispin,ispin,iorb,jorb,ik) = (delta_orb(gorb,iorb)*vops(jorb,k) + delta_orb(gorb,jorb)*vops(iorb,k))*dps(k)/den(k)
                    enddo
                 enddo
                 !
              enddo
           enddo
        enddo
        !
     end select
     !
  end select
end function grad_fdelta_bath_mats_main

function grad_fdelta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_,N) result(dG0out)
  integer,intent(in)                            :: ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  complex(8),intent(in)                         :: x
  complex(8),dimension(Norb,Norb,N)             :: dG0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dFdelta
  dFdelta = grad_fdelta_bath_mats_main(x,dmft_bath_,N)
  dG0out = dFdelta(ispin,jspin,:,:,:)
end function grad_fdelta_bath_mats_ispin_jspin


function grad_fdelta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(dG0out)
  integer,intent(in)                            :: iorb,jorb,ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  complex(8),intent(in)                         :: x
  complex(8),dimension(N)                       :: dG0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dFdelta
  dFdelta = grad_fdelta_bath_mats_main(x,dmft_bath_,N)
  dG0out = dFdelta(ispin,jspin,iorb,jorb,:)
end function grad_fdelta_bath_mats_ispin_jspin_iorb_jorb

function grad_fdelta_bath_mats_main_(x,bath_,N) result(dFdelta)
  complex(8),intent(in)                         :: x
  real(8),dimension(:)                          :: bath_
  integer,intent(in)                            :: N
  type(effective_bath)                          :: dmft_bath_
  logical                                       :: check
  complex(8),dimension(Nspin,Nspin,Norb,Norb,N) :: dFdelta
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dFdelta = grad_fdelta_bath_mats_main(x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_fdelta_bath_mats_main_

function grad_fdelta_bath_mats_ispin_jspin_(ispin,jspin,x,bath_,N) result(dG0out)
  complex(8),intent(in)                         :: x
  real(8),dimension(:)                          :: bath_
  integer,intent(in)                            :: ispin,jspin,N
  type(effective_bath)                          :: dmft_bath_
  logical                                       :: check
  complex(8),dimension(Norb,Norb,N)             :: dG0out
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dG0out = grad_fdelta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_fdelta_bath_mats_ispin_jspin_

function grad_fdelta_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_,N) result(dG0out)
  complex(8),intent(in)   :: x
  real(8),dimension(:)    :: bath_
  integer,intent(in)      :: iorb,jorb,ispin,jspin,N
  type(effective_bath)    :: dmft_bath_
  logical                 :: check
  complex(8),dimension(N) :: dG0out
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dG0out = grad_fdelta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_fdelta_bath_mats_ispin_jspin_iorb_jorb_
