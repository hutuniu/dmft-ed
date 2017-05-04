  !Diagonal Elements, i.e. local part
  htmp = zero
  htmp = htmp - xmu*(sum(nup)+sum(ndw))
  !
  do iorb=1,Norb
     htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
     htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
  enddo
  !
  hv(impi) = hv(impi) + htmp*vin(i)
  !
  
  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        !this loop considers only the orbital off-diagonal terms
        !because iorb=jorb can not have simultaneously
        !occupation 0 and 1, as required by this if Jcondition:
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
           !
        endif
     enddo
  enddo
  !
  !2. spin-flip part (only for the nonSU2 channel!)
  if(ed_mode=="nonsu2")then
     do ispin=1,Nspin
        jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
     !
     do iorb=1,Norb
        do jorb=1,Norb           
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
                 !
              endif
              !
           enddo
        enddo
     enddo
  endif

