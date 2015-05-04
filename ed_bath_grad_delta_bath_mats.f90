!NORMAL:
function grad_delta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(dfg)
  integer,intent(in)                        :: iorb,jorb,ispin,jspin
  type(effective_bath)                      :: dmft_bath_
  complex(8),intent(in)                     :: x
  integer                                   :: ik,ih,k,N,stride,korb
  complex(8)                                :: dfg(N)
  real(8),dimension(Nbath)                  :: eps,dps,vps,den
  real(8),dimension(Norb,Nbath)             :: vops
  real(8),dimension(Nspin,Nbath)            :: hps
  real(8),dimension(Nspin,Nspin,Nbath)      :: wps
  real(8),dimension(Nspin,Nspin,Norb,Nbath) :: wops
  real(8),dimension(Norb,Norb)              :: delta_orb
  real(8),dimension(2,2)                    :: delta_spin=reshape([1d0,0d0,0d0,1d0],[2,2]) !column order
  !
  if(N/=get_array_bath_dimension())stop "grad_delta_bath_mats_1 error: N != get_array_bath_dimension() "
  dfg=zero
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default
        eps = dmft_bath_%e(ispin,iorb,1:Nbath)
        vps = dmft_bath_%v(ispin,iorb,1:Nbath)
        !\grad_{E_{a}(k)} \Delta_{bb}^{rr} = [ V_{a}(k)*V_{a}(k) / ( iw_n - E_{a}(k) )**2 ]
        !
        !\grad_{V_{a}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}(k) / ( iw_n - E_{a}(k) ) ]
        !
        stride=0
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = vps(k)*vps(k)/(x - eps(k))**2
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = 2d0*vps(k)/(x - eps(k))
        enddo
        !
     case ("superc")
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
           dfg(ik) = -vps(k)*vps(k)*(1d0/den(k) - 2d0*eps(k)*( x + eps(k))/den(k)**2)
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = 2d0*vps(k)*vps(k)*dps(k)*(x + eps(k))/den(k)**2
        enddo
        stride=2*Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = 2d0*vps(k)*(x + eps(k))/den(k)
        enddo
        !
     case ("nonsu2")
        !we assume that the %e and %w are correctly popolated in this channel
        hps = dmft_bath_%e(        1:Nspin,iorb,1:Nbath)
        wps = dmft_bath_%w(1:Nspin,1:Nspin,iorb,1:Nbath)
        !
        !\grad_H_{a}^{h}(k) \Delta_{bb}^{ss`} =  W_{a}^{sh}(k) * W_{a}^{hs`}(k) / (iw_n - H_{a}^{h}(k))**2
        !
        !\grad_W_{a}^{hh`}(k) \Delta_{bb}^{ss`} =  delta(h,s) W_{a}^{h`s`}(k)/(iw_n - H_{a}^{h`}(k)) + delta(h`,s`) W_{a}^{sh}(k)/(iw_n - H_{a}^{h}(k))
        !
        stride=0
        do ih=1,Nspin
           do k=1,Nbath
              ik = stride + k + (ih-1)*Nbath
              dfg(ik) = wps(ispin,ih,k)*wps(ih,jspin,k)/(x - hps(ih,k))**2
           enddo
        enddo
        stride=Nspin*Nbath
        do ih=1,Nspin
           do jh=1,Nspin
              do k=1,Nbath
                 ik = stride + k + (jh-1)*Nbath + (ih-1)*Nspin*Nbath
                 dfg(ik) = &
                      delta_spin(ih,ispin)*wps(jh,jspin,k)/(x - hps(jh,k)) + &
                      delta_spin(jh,jspin)*wps(ispin,ih,k)/(x - hps(ih,k))
              enddo
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")
     !
     select case(ed_mode)
     case default
        eps  = dmft_bath_%e(  ispin,1     ,1:Nbath)
        vops = dmft_bath_%v(1,ispin,1:Norb,1:Nbath)
        delta_orb = eye(Norb)
        !
        !\grad_{E_{1}(k)} \Delta_{ab} = [ V_{a}(k)*V_{b}(k) / ( iw_n - E_{1}(k) )**2 ]
        !
        !\grad_{V_{g}(k)} \Delta_{ab} = [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ] / ( iw_n - E_{1}(k) )
        !
        stride=0
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = vops(iorb,k)*vops(jorb,k)/(x - eps(k))**2
        enddo
        stride=Nbath
        do korb=1,Norb
           do k=1,Nbath
              ik = stride + k + (korb-1)*Nbath
              dfg(ik) = (delta_orb(korb,iorb)*vops(jorb,k) + delta_orb(korb,jorb)*vops(iorb,k))/(x - eps(k))
           enddo
        enddo
        !
     case ("superc")
        eps  = dmft_bath_%e(  ispin,1     ,1:Nbath)
        dps  = dmft_bath_%d(  ispin,1     ,1:Nbath)
        vops = dmft_bath_%v(1,ispin,1:Norb,1:Nbath)
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
        stride=0
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = -vops(iorb,k)*vops(jorb,k)*(1d0/den(k) - 2d0*eps(k)*( x + eps(k))/den(k)**2)
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = 2d0*vops(iorb,k)*vops(jorb,k)*dps(k)*(x + eps(k))/den(k)**2
        enddo
        stride=2*Nbath
        do korb=1,Norb
           do k=1,Nbath
              ik = stride + k + (korb-1)*Nbath
              dfg(ik) = (delta_orb(korb,iorb)*vops(jorb,k) + delta_orb(korb,jorb)*vops(iorb,k))*(x + eps(k))/den(k)
           enddo
        enddo
        !
     case ("nonsu2")
        !we assume that the %e and %w are correctly popolated in this channel
        hps  = dmft_bath_%e(       1:Nspin,1     ,1:Nbath)
        wops = dmft_bath_%w(1:Nspin,1:Nspin,1:Norb,1:Nbath)
        delta_orb = eye(Norb)
        !
        !\grad_H^{h}(k) \Delta_{ab}^{ss`} =  W_{a}^{sh}(k) * W_{b}^{hs`}(k) / (iw_n - H^{h}(k))**2
        !
        !\grad_W_{g}^{hh`}(k) \Delta_{bb}^{ss`} =  delta(h,s) W_{a}^{h`s`}(k)/(iw_n - H_{a}^{h`}(k)) + delta(h`,s`) W_{a}^{sh}(k)/(iw_n - H_{a}^{h}(k))
        !
        stride=0
        do ih=1,Nspin
           do k=1,Nbath
              ik = stride + k + (ih-1)*Nbath
              dfg(ik) = wops(ispin,ih,iorb,k)*wops(ih,jspin,jorb,k)/(x - hps(ih,k))**2
           enddo
        enddo
        stride=Nspin*Nbath
        do ih=1,Nspin
           do jh=1,Nspin
              do korb=1,Norb
                 do k=1,Nbath
                    ik = stride + k + (korb-1)*Nbath + (jh-1)*Norb*Nbath + (ih-1)*Nspin*Norb*Nbath
                    dfg(ik) = &
                         delta_orb(korb,iorb)*delta_spin(ih,ispin)*wops(jh,jspin,jorb,k)/(x - hps(jh,k)) + &
                         delta_orb(korb,jorb)*delta_spin(jh,jspin)*wops(ispin,ih,iorb,k)/(x - hps(ih,k))
                 enddo
              enddo
           enddo
        enddo
        !
     end select
  end select
end function grad_delta_bath_mats_1

function grad_delta_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_,N) result(dfg)
  integer,intent(in)                      :: iorb,jorb,ispin,jspin,N
  type(effective_bath)                    :: dmft_bath_
  real(8),dimension(:)                    :: bath_
  logical                                 :: check
  complex(8),intent(in)                   :: x
  complex(8)                              :: dfg(N)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dfg = grad_delta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_delta_bath_mats_2






!ANOMALous:
function grad_fdelta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(dfg)
  integer,intent(in)                      :: iorb,ispin,jorb,jspin
  type(effective_bath)                    :: dmft_bath_
  complex(8),intent(in)                   :: x
  integer                                 :: ik,k,N,stride,ih,korb
  complex(8)                              :: dfg(N)
  real(8),dimension(Norb,Norb)            :: delta_orb
  real(8),dimension(Nbath)                :: eps,dps,vps,den
  real(8),dimension(Norb,Nbath)           :: vops
  dfg=zero
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        eps = dmft_bath_%e(  ispin,iorb,1:Nbath)
        dps = dmft_bath_%d(  ispin,iorb,1:Nbath)
        vps = dmft_bath_%v(1,ispin,iorb,1:Nbath)
        ! Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
        forall(k=1:Nbath)den(k) = dimag(x)**2 + eps(k)**2 + dps(k)**2
        !
        !\grad_{E_{a}(k)} \FDelta_{aa} = -2 * V_{a}(k) * V_{a}(k) * E_{a}(k) * \Delta_{a}(k) / Den**2
        !
        !\grad_{\Delta_{a}(k)} \FDelta_{aa} = V_{a}(k) * V_{a}(k) * [ 1/den - 2* \Delta_{a}(k)*\Delta_{a}(k)/den**2 ]
        !
        !\grad_{ V_{a}(k)} \FDelta_{aa} =  2 * V_{a}(k) * \Delta_{a}(k) / den
        !
        stride=0
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = -2d0*vps(k)*vps(k)*eps(k)*dps(k)/den(k)**2
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = vps(k)*vps(k)*(1d0/den(k) - 2d0*dps(k)*dps(k)/den(k)**2)
        enddo
        stride=2*Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = 2d0*vps(k)*dps(k)/den(k)
        enddo
        !
     case ("nonsu2")
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=nonsu2, bath_type=normal"
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
        eps  = dmft_bath_%e(  ispin,1     ,1:Nbath)
        dps  = dmft_bath_%d(  ispin,1     ,1:Nbath)
        vops = dmft_bath_%v(1,ispin,1:Norb,1:Nbath)
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
        stride=0
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = -2d0*vops(iorb,k)*vops(jorb,k)*eps(k)*dps(k)/den(k)**2
        enddo
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           dfg(ik) = vops(iorb,k)*vops(jorb,k)*(1d0/den(k) - 2d0*dps(k)*dps(k)/den(k)**2)
        enddo
        stride=2*Nbath
        do korb=1,Norb
           do k=1,Nbath
              ik = stride + k + (korb-1)*Nbath
              dfg(ik) = (delta_orb(korb,iorb)*vops(jorb,k) + delta_orb(korb,jorb)*vops(iorb,k))*dps(k)/den(k)
           enddo
        enddo
        !
     case ("nonsu2")
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=nonsu2, bath_type=hybrid"
        !
     end select
     !
  end select
end function grad_fdelta_bath_mats_1

function grad_fdelta_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_,N) result(dfg)
  integer,intent(in)            :: iorb,ispin,jorb,jspin,N
  type(effective_bath)          :: dmft_bath_
  real(8),dimension(:)          :: bath_
  logical                       :: check
  complex(8),intent(in)         :: x
  complex(8)                    :: dfg(N)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  dfg =  grad_fdelta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N)
  call deallocate_bath(dmft_bath_)
end function grad_fdelta_bath_mats_2


