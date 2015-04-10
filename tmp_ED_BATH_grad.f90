!+-------------------------------------------------------------------+
!PURPOSE  : compute the hybridization function for a given spin and 
! orbital indices ispin and iorb at a point x from the input 
! effective_bath type dmft_bath
!+-------------------------------------------------------------------+
!NORMAL:
!Matsubara:
function grad_delta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(fg)
  integer,intent(in)                 :: iorb,jorb,ispin,jspin
  type(effective_bath)               :: dmft_bath_
  complex(8),intent(in)              :: x
  integer                            :: ik,k,N,stride
  complex(8)                         :: fg(N),Den
  real(8),dimension(2,2)             :: delta_hel=reshape([1d0,0d0,0d0,1d0],[2,2]) !column order
  real(8),dimension(Norb,Norb)       :: delta_orb
  real(8),dimension(Nbath)           :: eps,dps,vps
  real(8),dimension(Norb,Nbath)      :: vops
  real(8),dimension(Nhel,Nbath)      :: hps
  real(8),dimension(Nhel,Nspin,Nbath):: vhps
  real(8),dimension(Nhel,Nhel,Nbath) :: wps

  if(N/=get_array_bath_dimension())stop "grad_delta_bath_mats_1 error: N != get_array_bath_dimension() "
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default

        !\grad_{E_{a}^{s}(k)} \Delta_{bb}^{rr} = [ V_{a}^{s}(k)*V_{a}^{s}(k) / ( iw_n - E_{a}^{s}(k) )**2 ]
        stride=0
        do k=1,Nbath
           ik = stride + k
           fg(ik) = dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,iorb,k)/(x - dmft_bath_%e(ispin,iorb,k))**2
        enddo
        !\grad_{V_{a}^{s}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}^{s}(k) / ( iw_n - E_{a}^{s}(k) ) ]
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           fg(ik) = 2d0*dmft_bath_%v(1,ispin,iorb,k)/(x - dmft_bath_%e(ispin,iorb,k))
        enddo
        !
        !
     case ("superc")
        ! Den = (w_n**2 + E_{a}^{s}(k)**2 + \Delta_{a}^{s}(k)**2
        !\grad_{E_{a}^{s}(k)} \Delta_{bb}^{rr} = -V_{a}^{s}(k)*V_{a}^{s}(k)*[ 1/den - 2*E_{a}^{s}(k)*(iw_n + E_{a}^{s}(k))/den**2 ]
        stride=0
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = -dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,iorb,k)*(1d0/den - 2d0*dmft_bath_%e(ispin,iorb,k)*( x + dmft_bath_%e(ispin,iorb,k))/den**2)
        enddo
        !\grad_{\Delta_{a}^{s}(k)} \Delta_{bb}^{rr} = V_{a}^{s}(k)*V_{a}^{s}(k)*\Delta_{a}^{s}(k)*(iw_n + E_{a}^{s}(k)) /den**2
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = 2d0*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%d(ispin,iorb,k)*(x + dmft_bath_%e(ispin,iorb,k))/den**2
        enddo
        !\grad_{ V_{a}^{s}(k)} \Delta_{bb}^{rr} =  2*V_{a}^{s}(k)*(iw_n + E_{a}^{s}(k))/den
        stride=2*Nbath
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = 2d0*dmft_bath_%v(1,ispin,iorb,k)*(x + dmft_bath_%e(ispin,iorb,k))/den
        enddo
        !
        !
     case ("nonsu2")
        call pull_nonsu2_components(dmft_bath_)
        !\grad_H_{a}^{h}(k) \Delta_{bb}^{ss`} =  W_{a}^{sh}(k) * W_{a}^{hs`}(k) / (iw_n - H_{a}^{h}(k))**2
        stride=0
        do ihel=1,Nhel
           do k=1,Nbath
              ik = stride + k + (ihel-1)*Nbath
              fg(ik) = dmft_bath%w(ispin,ihel,iorb,k)*dmft_bath%w(ihel,jspin,iorb,k)/(x - dmft_bath_%h(ihel,iorb,k))**2
           enddo
        enddo
        !\grad_W_{a}^{hh`}(k) \Delta_{bb}^{ss`} =  delta(h,s) W_{a}^{h`s`}(k)/(iw_n - H_{a}^{h`}(k)) + delta(h`,s`) W_{a}^{sh}(k)/(iw_n - H_{a}^{h}(k))
        stride=Nhel*Nbath
        do ihel=1,Nhel
           do jhel=1,Nhel
              do k=1,Nbath
                 ik = stride + k + (jhel-1)*Nbath + (ihel-1)*Nhel*Nbath
                 fg(ik) = &
                      delta_hel(ihel,ispin)*[ dmft_bath%w(jhel,jspin,iorb,k)/(x - dmft_bath_%h(jhel,iorb,k)) ] + &
                      delta_hel(jhel,jspin)*[ dmft_bath%w(ispin,ihel,iorb,k)/(x - dmft_bath_%h(ihel,iorb,k)) ]
              enddo
           enddo
        enddo
        !
        !
     end select
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        !\grad_{E_{1}^{s}(k)} \Delta_{ab}^{rr} = [ V_{a}^{s}(k)*V_{b}^{s}(k) / ( iw_n - E_{1}^{s}(k) )**2 ]
        stride=0
        do k=1,Nbath
           ik = stride + k
           fg(ik) = dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,jorb,k)/(x - dmft_bath_%e(ispin,1,k))**2
        enddo
        !\grad_{V_{a}^{s}(k)} \Delta_{bc}^{rr} = [ 2*V_{a}^{s}(k) / ( iw_n - E_{a}^{s}(k) ) ]
        stride=Nbath
        do k=1,Nbath
           if( iorb == jorb )then
              fg(iorb*stride+k) = 2d0*dmft_bath_%v(1,ispin,iorb,k)/(x - dmft_bath_%e(ispin,1,k))
           else
              fg(iorb*stride+k) = dmft_bath_%v(1,ispin,jorb,k)/(x - dmft_bath_%e(ispin,1,k))
              fg(jorb*stride+k) = dmft_bath_%v(1,ispin,iorb,k)/(x - dmft_bath_%e(ispin,1,k))
           endif
        enddo
        !
        !
     case ("superc")
        ! should be:
        !\Delta_{ab}^{ss} = - \sum_k [ V_{a}^{s}(k) * V_{b}^{s}(k) * (iw_n + E^{s}(k)) /
        !                             (w_n**2 + E^{s}(k)**2 + \Delta^{s}(k)**2) ]
        stop "Delta_bath_mats error: called with ed_mode=superc, bath_type=hybrid. THIS IS NOT YET CHECKED"
        fg = -sum( dmft_bath_%v(1,ispin,iorb,:)*dmft_bath_%v(1,ispin,jorb,:)*(x + dmft_bath_%e(ispin,1,:))/&
             (dimag(x)**2 + dmft_bath_%e(ispin,1,:)**2 + dmft_bath_%d(ispin,1,:)**2))
        !
        !
     case ("nonsu2")
        !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{hs`}(k)/(iw_n - EH^{h}(k))]
        call set_nonsu2_components(dmft_bath_)
        fg = zero
        do ihel=1,Nhel
           fg = fg + sum( dmft_bath_%w(ispin,ihel,iorb,:)*dmft_bath_%w(ihel,jspin,jorb,:)/(x - dmft_bath_%h(ihel,1,:))
        enddo
        !
        !
     end select
  end select
end function grad_delta_bath_mats_1
!
!
!
!ANOMALous:
!Matsubara:
function grad_fdelta_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_,N) result(fg)
  integer,intent(in)    :: iorb,ispin,jorb,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  integer               :: ik,k,N,stride
  complex(8)            :: fg(N),Den
  if(N/=get_array_bath_dimension())stop "grad_fdelta_bath_mats_1 error: N != get_array_bath_dimension() "
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=normal, bath_type=normal"
     case ("superc")
        !Den = (w_n**2 + E_{a}^{s}(k)**2 + \Delta_{a}^{s}(k)**2
        !
        !\grad_{E_{a}^{s}(k)} \FDelta_{aa}^{ss} = -2 * V_{a}^{s}(k) * V_{a}^{s}(k) * E_{a}^{s}(k) * \Delta_{a}^{s}(k) / Den**2
        stride=0
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = -2d0*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%e(ispin,iorb,k)*dmft_bath_%d(ispin,iorb,k)/den**2
        enddo
        !\grad_{\Delta_{a}^{s}(k)} \FDelta_{aa}^{ss} = V_{a}^{s}(k) * V_{a}^{s}(k) * [ 1/den - 2* \Delta_{a}^{s}(k)*\Delta_{a}^{s}(k)/den**2 ]
        stride=Nbath
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = 2d0*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%d(ispin,iorb,k)*(x + dmft_bath_%e(ispin,iorb,k))/den**2
        enddo
        !\grad_{ V_{a}^{s}(k)} \FDelta_{aa}^{ss} =  2 * V_{a}^{s}(k) * \Delta_{a}^{s}(k) / den
        stride=2*Nbath
        do k=1,Nbath
           ik = stride + k
           Den    = dimag(x)**2 + dmft_bath_%e(ispin,iorb,k)**2 + dmft_bath_%d(ispin,iorb,k)**2
           fg(ik) = 2d0*dmft_bath_%v(1,ispin,iorb,k)*dmft_bath_%d(ispin,iorb,k)/den
        enddo
        !
     case ("nonsu2")
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=nonsu2, bath_type=normal"
     end select
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
     case ("superc")
        ! should be:
        !\FDelta_{ab}^{ss} = - \sum_k [ \Delta^{s}(k) * V_{a}^{s}(k) * V_{b}^{s}(k) /
        !                             (w_n**2 + E^{s}(k)**2 + \Delta^{s}(k)**2) ]
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=superc, bath_type=hybrid. THIS IS NOT YET CHECKED"
        fg = -sum( dmft_bath_d%(ispin, 1,:)*dmft_bath_%v(1,ispin,iorb,:)*dmft_bath_%v(1,ispin,jorb,:)/&
             ( dimag(x)**2 + dmft_bath_%e(ispin,   1,:)**2 + dmft_bath_%d(ispin,   1,:)**2 ) )
        !
        !
     case ("nonsu2")
        stop "Grad_Fdelta_bath_mats error: called with ed_mode=nonsu2, bath_type=hybrid"
     end select
  end select
end function grad_fdelta_bath_mats_1
!
