  !diagonal, spin conserving:
  do iorb=1,Norb
     do kp=1,Nbath
        ms=getBathStride(iorb,kp)
        !
        ! IMP UP <--> BATH UP
        if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==1) .AND. (ib(ms)==0) )then
           call c(iorb,m,k1,sg1)
           call cdg(ms,k1,k2,sg2)
           j = binary_search(H%map,k2)
           htmp = diag_hybr(1,iorb,kp)*sg1*sg2
           hv(impi) = hv(impi) + htmp*vin(j)
        endif
        if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(ms)==1) )then
           call c(ms,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j=binary_search(H%map,k2)
           htmp = diag_hybr(1,iorb,kp)*sg1*sg2
           hv(impi) = hv(impi) + htmp*vin(j)
        endif
        !
        !IMP DW <--> BATH DW
        if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==1) .AND. (ib(ms+Ns)==0) )then
           call c(iorb+Ns,m,k1,sg1)
           call cdg(ms+Ns,k1,k2,sg2)
           j=binary_search(H%map,k2)
           htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
           hv(impi) = hv(impi) + htmp*vin(j)
        endif
        if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==0) .AND. (ib(ms+Ns)==1) )then
           call c(ms+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j=binary_search(H%map,k2)
           htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
           hv(impi) = hv(impi) + htmp*vin(j)
        endif
     enddo
  enddo


  !Off-diagonal, spin-flipping: (only nonsu2 & !replica bath)
  if((ed_mode=="nonsu2").and.(bath_type/="replica"))then
     do iorb=1,Norb
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           !
           ! IMP UP <--> BATH DW
           if( (ib(iorb)==1) .AND. (ib(ms+Ns)==0) )then
              call c(iorb,m,k1,sg1)
              call cdg(ms+Ns,k1,k2,sg2)
              j = binary_search(H%map,k2)
              htmp = dmft_bath%u(1,iorb,kp)*sg1*sg2
              hv(impi) = hv(impi) + htmp*vin(j)
           endif
           if( (ib(iorb)==0) .AND. (ib(ms+Ns)==1) )then
              call c(ms+Ns,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp = dmft_bath%u(1,iorb,kp)*sg1*sg2
              hv(impi) = hv(impi) + htmp*vin(j)
           endif
           ! IMP DW <--> BATH UP
           if( (ib(iorb+Ns)==1) .AND. (ib(ms)==0) )then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp = dmft_bath%u(Nspin,iorb,kp)*sg1*sg2
              hv(impi) = hv(impi) + htmp*vin(j)
           endif
           if( (ib(iorb+Ns)==0) .AND. (ib(ms)==1) )then
              call c(ms,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp = dmft_bath%u(Nspin,iorb,kp)*sg1*sg2
              hv(impi) = hv(impi) + htmp*vin(j)
           endif
        enddo
     enddo
  endif
