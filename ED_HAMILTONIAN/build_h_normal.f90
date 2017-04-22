  states: do i=first_state,last_state
     m = H%map(i)
     impi = i-ishift
     ib = bdecomp(m,2*Ns)
     !
     htmp=zero
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     !
     !
     !CHEMICAL POTENTIAL
     htmp = htmp - xmu*(sum(nup)+sum(ndw))
     !
     !
     !IMPURITY  HAMILTONIAN (diagonal elements: local part)
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
     enddo
     !
     !
     !LOCAL INTERACTION
     !
     !density-density interaction: same orbital, opposite spins:
     ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Norb
        htmp = htmp + Uloc(iorb)*nup(iorb)*ndw(iorb)
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Norb
           htmp = htmp - 0.5d0*Uloc(iorb)*(nup(iorb)+ndw(iorb)) + 0.25d0*uloc(iorb)
        enddo
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp=htmp-0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*Ust
                 htmp=htmp-0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*(Ust-Jh)
              enddo
           enddo
        endif
     endif
     !
     !
     call sp_insert_element(spH0,htmp,impi,i)
     !
     !
     !
     !
     !IMPURITY  HAMILTONIAN (off-diagonal elements: non-local part)
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
              call sp_insert_element(spH0,htmp,impi,j)
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
              call sp_insert_element(spH0,htmp,impi,j)
              !
           endif
           !
        enddo
     enddo
     !
     !
     !
     ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
     !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
     !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
     if(Norb>1.AND.Jhflag)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (ib(jorb)==1).AND.&
                   (ib(iorb+Ns)==1).AND.&
                   (ib(jorb+Ns)==0).AND.&
                   (ib(iorb)==0))
              if(Jcondition)then
                 call c(jorb,m,k1,sg1)
                 call c(iorb+Ns,k1,k2,sg2)
                 call cdg(jorb+Ns,k2,k3,sg3)
                 call cdg(iorb,k3,k4,sg4)
                 j=binary_search(H%map,k4)
                 htmp = one*Jh*sg1*sg2*sg3*sg4
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
     if(Norb>1.AND.Jhflag)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (ib(jorb)==1).AND.&
                   (ib(jorb+Ns)==1).AND.&
                   (ib(iorb+Ns)==0).AND.&
                   (ib(iorb)==0))
              if(Jcondition)then
                 call c(jorb,m,k1,sg1)
                 call c(jorb+Ns,k1,k2,sg2)
                 call cdg(iorb+Ns,k2,k3,sg3)
                 call cdg(iorb,k3,k4,sg4)
                 j=binary_search(H%map,k4)
                 htmp = one*Jh*sg1*sg2*sg3*sg4
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
           enddo
        enddo
     endif
     !
     !
     !
     !
     !
     !
     !BATH HAMILTONIAN
     if(bath_type/="replica") then
        !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        !
        htmp=zero
        !
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              alfa=getBathStride(iorb,kp)
              htmp =htmp + &
                   dmft_bath%e(1,iorb,kp)*ib(alfa) +& !UP
                   dmft_bath%e(Nspin,iorb,kp)*ib(alfa+Ns) !DW
           enddo
        enddo
        !
        call sp_insert_element(spH0,htmp,impi,i)
        !
     else
        !diagonal elements
        !
        htmp=zero
        !
        do kp=1,Nbath
           do iorb=1,Norb
              alfa = getBathStride(iorb,kp) !iorb + kp*Norb
              htmp = htmp + &
                   dmft_bath%h(1,1,iorb,iorb,kp)*ib(alfa) + & !UP
                   dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ib(alfa+Ns) !DW
           enddo
        enddo
        !
        call sp_insert_element(spH0,htmp,impi,i)
        !
        !off-diagonal elements
        !
        htmp=zero
        !
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 if(iorb==jorb)cycle
                 !UP
                 alfa = getBathStride(iorb,kp)
                 beta = getBathStride(jorb,kp)
                 Jcondition = &
                      (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) .AND. &
                      (ib(beta)==1)                         .AND. &
                      (ib(alfa)==0)
                 if (Jcondition)then
                    call c(beta,m,k1,sg1)
                    call cdg(alfa,k1,k2,sg2)
                    j = binary_search(H%map,k2)
                    htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                    !
                    call sp_insert_element(spH0,htmp,impi,j)
                    !
                 endif
                 !DW
                 alfa = getBathStride(iorb,kp)
                 beta = getBathStride(jorb,kp)
                 Jcondition = &
                      (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                      (ib(beta+Ns)==1)                         .AND. &
                      (ib(alfa+Ns)==0)
                 if (Jcondition)then
                    call c(beta+Ns,m,k1,sg1)
                    call cdg(alfa+Ns,k1,k2,sg2)
                    j = binary_search(H%map,k2)
                    htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                    !
                    call sp_insert_element(spH0,htmp,impi,j)
                    !
                 endif
                 !
              enddo
           enddo
        enddo
        !
     endif
     !
     !
     !
     !
     !
     !
     !
     !IMPURITY- BATH HYBRIDIZATION
     !DIAGONAL
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
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(ms)==1) )then
              call c(ms,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           !IMP DW <--> BATH DW
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==1) .AND. (ib(ms+Ns)==0) )then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(ms+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==0) .AND. (ib(ms+Ns)==1) )then
              call c(ms+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
        enddo
     enddo

  enddo states

