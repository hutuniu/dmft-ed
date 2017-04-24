  !Diagonal Elements, i.e. local part
  do iorb=1,Norb
     htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
     htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
  enddo


  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        if(iorb==jorb)cycle
        !UP
        Jcondition = &
             (impHloc(1,1,iorb,jorb)/=zero) .AND. &
             (ib(jorb)==1)                  .AND. &
             (ib(iorb)==0)
        if (Jcondition) then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j = binary_search(H%map,k2)
           htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
           !
           hv(impi) = hv(impi) + htmp*vin(j)
           !hv(j) = hv(j) + htmp*vin(impi)!??
           !
        endif
        !DW
        Jcondition = &
             (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
             (ib(jorb+Ns)==1)                  .AND. &
             (ib(iorb+Ns)==0)
        if (Jcondition) then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j = binary_search(H%map,k2)
           htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
           !
           hv(impi) = hv(impi) + htmp*vin(j)
           !hv(j) = hv(j) + htmp*vin(impi)!??
           !
        endif
        !
     enddo
  enddo
  !
  !
  !1. spin-flip part (only for the nonSU2 channel!)
  if(ed_mode=="nonsu2")then
     do iorb=1,Norb
        do jorb=1,Norb
           do ispin=1,Nspin
              !do jspin=1,Nspin
              jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
              alfa = iorb + (ispin-1)*Ns
              beta = jorb + (jspin-1)*Ns
              Jcondition=&
                   (impHloc(ispin,jspin,iorb,jorb)/=zero) .AND. &
                   (ib(beta)==1) .AND. (ib(alfa)==0)
              if(Jcondition)then
                 call c(beta,m,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = impHloc(ispin,jspin,iorb,jorb)*sg1*sg2
                 !
                 hv(impi) = hv(impi) + htmp*vin(j)
                 !hv(j) = hv(j) + htmp*vin(impi)!??
                 !
              endif
              !
              !enddo
           enddo
        enddo
     enddo
  endif

