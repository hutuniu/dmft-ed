  Breplica: if(bath_type/="replica") then
     !
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     htmp=zero
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(alfa)        !UP
           htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(alfa+Ns) !DW
        enddo
     enddo
     !
     hv(impi) = hv(impi) + htmp*vin(i)
     !
  else
     !
     !
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     htmp=zero
     do kp=1,Nbath
        do iorb=1,Norb
           alfa = getBathStride(iorb,kp)
           htmp = htmp + dmft_bath%h(1,1,iorb,iorb,kp)*ib(alfa)              !UP
           htmp = htmp + dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ib(alfa+Ns)   !DW
        enddo
     enddo
     !
     hv(impi) = hv(impi) + htmp*vin(i)
     !
     !off-diagonal elements
     !
     !1. same spin:
     do kp=1,Nbath
        do iorb=1,Norb
           do jorb=1,Norb
              !this loop considers only the orbital off-diagonal terms
              !because iorb=jorb can not have simultaneously
              !occupation 0 and 1, as required by this if Jcondition:
              !UP
              alfa = getBathStride(iorb,kp)
              beta = getBathStride(jorb,kp)
              Jcondition = &
                   (dmft_bath%h(1,1,iorb,jorb,kp)/=zero)               .AND. &
                   (ib(beta)==1)                                       .AND. &
                   (ib(alfa)==0)
              if (Jcondition)then
                 call c(beta,m,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                 !
                 hv(impi) = hv(impi) + htmp*vin(j)
                 !
              endif
              !DW
              alfa = getBathStride(iorb,kp) + Ns
              beta = getBathStride(jorb,kp) + Ns
              Jcondition = &
                   (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero)       .AND. &
                   (ib(beta)==1)                                       .AND. &
                   (ib(alfa)==0)
              if (Jcondition)then
                 call c(beta,m,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                 !
                 hv(impi) = hv(impi) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     enddo
     !
     !2. spin-flip part (only for the nonSU2 channel!)
     if(ed_mode=="nonsu2")then
        do kp=1,Nbath
           do ispin=1,Nspin
              jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
              do iorb=1,Norb
                 do jorb=1,Norb
                    alfa = getBathStride(iorb,kp) + (ispin-1)*Ns
                    beta = getBathStride(jorb,kp) + (jspin-1)*Ns
                    Jcondition=&
                         (dmft_bath%h(ispin,jspin,iorb,jorb,kp)/=zero) .AND. &
                         (ib(beta)==1)                                 .AND. &
                         (ib(alfa)==0)
                    if(Jcondition)then
                       call c(beta,m,k1,sg1)
                       call cdg(alfa,k1,k2,sg2)
                       j = binary_search(H%map,k2)
                       htmp = dmft_bath%h(ispin,jspin,iorb,jorb,kp)*sg1*sg2
                       !
                       hv(impi) = hv(impi) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
  endif Breplica



  !anomalous pair-creation/destruction
  if(ed_mode=="superc")then
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           !\Delta_l c_{\up,ms} c_{\dw,ms}
           if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==1) .AND. (ib(ms+Ns)==1) )then
              call c(ms,m,k1,sg1)
              call c(ms+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2
              !
              hv(impi) = hv(impi) + htmp*vin(j)
              !
           endif
           !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
           if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==0) .AND. (ib(ms+Ns)==0) )then
              call cdg(ms+Ns,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2 !
              !
              hv(impi) = hv(impi) + htmp*vin(j)
              !
           endif
        enddo
     enddo
  endif



  ! ! CLUSTER REPLICA-HAMILTONIAN - no inter-cluster couplings
  ! do kp=1,Nbath
  !    do iorb=1,Norb
  !       do jorb=1,Norb
  !          do ispin=1,Nspin
  !             do jspin=1,Nspin
  !                !
  !                if(dmft_bath%h(ispin,jspin,iorb,jorb,kp)/=zero) then
  !                   !
  !                   alfa = iorb + kp*Norb + (ispin-1)*Ns
  !                   beta = jorb + kp*Norb + (jspin-1)*Ns
  !                   !
  !                   !diagonal elements
  !                   if ((ispin==jspin).and.(iorb==jorb)) then
  !                      htmp = dmft_bath%h(ispin,jspin,iorb,jorb,kp)*real(ib(alfa),8)
  !                      !
  !                      hv(impi) = hv(impi) + htmp*vin(i)
  !                      !
  !                   endif
  !                   !
  !                   !off-diagonal elements
  !                   if ((ib(beta)==1) .AND. (ib(alfa)==0)) then
  !                      call c(beta,m,k1,sg1)
  !                      call cdg(alfa,k1,k2,sg2)
  !                      j = binary_search(H%map,k2)
  !                      htmp = dmft_bath%h(ispin,jspin,iorb,jorb,kp)*sg1*sg2
  !                      !
  !                      hv(impi) = hv(impi) + htmp*vin(j)
  !                      !
  !                   endif
  !                   !
  !                endif
  !                !
  !             enddo
  !          enddo
  !       enddo
  !    enddo
  ! enddo
