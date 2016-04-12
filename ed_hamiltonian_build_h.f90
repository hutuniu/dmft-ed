  states: do i=first_state,last_state
     m=Hmap(i)
     impi = i-ishift
     call bdecomp(m,ib)
     htmp=0.d0
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     !
     !
     !##############################################
     !#                                            #
     !#          LOCAL HAMILTONIAN PART            #
     !#                                            #
     !##############################################
     !
     !
     !          ==> LOCAL ENERGIES <==
     !
     htmp = htmp - xmu*(sum(nup)+sum(ndw))
     htmp = htmp + dot_product(eloc(1,:),nup) + dot_product(eloc(Nspin,:),ndw)
     !
     !
     !     ==> IMPURITY INTERACTIONS <==
     !
     ! DENSITY-DENSITY INTERACTIONS: same orbital, opposite spins
     !    \sum_\a \pm U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Norb
        htmp = htmp + Uloc(iorb)*nup(iorb)*ndw(iorb)
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins
        !Eust=        U'   *sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        !    =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        !Eund = \sum_{i<j}    Und     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        !    "="\sum_{i<j} (Ust-Jh)   *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        !    "="\sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
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
     call sp_insert_element(spH0,htmp,impi,i)
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
                 j=binary_search(Hmap,k4)
                 htmp = Jh*sg1*sg2*sg3*sg4
                 call sp_insert_element(spH0,htmp,impi,j)
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
                 j=binary_search(Hmap,k4)
                 htmp = Jh*sg1*sg2*sg3*sg4
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
           enddo
        enddo
     endif
     !
     !
     !   ==> IMPURITY LOCAL HYBRIDIZATION <==
     !
     do iorb=1,Norb
        do jorb=1,Norb
           ! HYBRIDIZATION TERMS I: same or different orbitals, same spins.
           ! SPIN UP
           if((impHloc(1,1,iorb,jorb)/=0d0).AND.(ib(iorb)==0).AND.(ib(jorb)==1))then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           ! SPIN DW
           if((impHloc(Nspin,Nspin,iorb,jorb)/=0d0).AND.(ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           ! HYBRIDIZATION TERMS II: same or different orbitals, opposite spins.
           ! UP-DW
           if(impHloc(1,Nspin,iorb,jorb)/=0d0)then
              if( (ib(iorb)==0).AND.(ib(jorb+Ns)==1) )then
                 call c(jorb+Ns,m,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp = impHloc(1,Nspin,iorb,jorb)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
           endif
           ! DW-UP
           if(impHloc(Nspin,1,iorb,jorb)/=0d0)then
              if( (ib(iorb+Ns)==0).AND.(ib(jorb)==1) )then
                 call c(jorb,m,k1,sg1)
                 call cdg(iorb+Ns,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp = impHloc(Nspin,1,iorb,jorb)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
           endif
        enddo
     enddo
     !
     !
     !##############################################
     !#                                            #
     !#              BATH HAMILTONIAN              #
     !#                                            #
     !##############################################
     !
     !
     if(bath_type/="replica") then
        ! DIAGONAL BATH HAMILTONIAN: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
        htmp=0.0d0
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ib(ms),8) + dmft_bath%e(Nspin,iorb,kp)*real(ib(ms+Ns),8)
           enddo
        enddo
        call sp_insert_element(spH0,htmp,impi,i)
        !
     else
        ! CLUSTER REPLICA-HAMILTONIAN - no inter-cluster couplings
        do kp=1,Nbath
           !
           do iorb=1,Norb
              do jorb=1,Norb
                 do ispin=1,Nspin
                    do jspin=1,Nspin
                       !
                       if(dmft_bath%h(ispin,jspin,iorb,jorb,kp)/=0d0) then
                          !
                          ! STRIDE:  {Norb,up}[Norb,up(k=1)][Norb,up(k=2)]{Norb,dw}[Norb,dw(k=1)][Norb,dw(k=2)]
                          ! [ c+_(iorb,ispin)c_(jorb,jspin) ]_k + h.c. = c+_alpha c_beta + c+_beta c_alpha
                          !  
                          alfa = iorb + kp*Norb + (ispin-1)*Ns
                          beta = jorb + kp*Norb + (jspin-1)*Ns
                          !
                          if ((ib(beta)==1) .AND. (ib(alfa)==0)) then
                             call c(beta,m,k1,sg1)
                             call cdg(alfa,k1,k2,sg2)
                             j = binary_search(Hmap,k2)
                             htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                             call sp_insert_element(spH0,htmp,impi,j)
                          endif
                          !
                          if( (ib(alfa)==1) .AND. (ib(beta)==0) )then
                             call c(alfa,m,k1,sg1)
                             call cdg(beta,k1,k2,sg2)
                             j = binary_search(Hmap,k2)
                             htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                             call sp_insert_element(spH0,htmp,impi,j)
                          endif
                       endif
                       !
                    enddo
                 enddo
              enddo
           enddo
           !
        enddo
       !
     endif
     !
     !
     !##############################################
     !#                                            #
     !#        IMPURITY- BATH HYBRIDIZATION        #
     !#                                            #
     !##############################################
     !
     !
     !        ==> DIAGONAL-HYBRIDIZATION <==
     !
     do iorb=1,Norb
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           !
           ! IMP UP <--> BATH UP
           if( (dmft_bath%v(1,iorb,kp)/=0d0) .AND. (ib(iorb)==1) .AND. (ib(ms)==0) )then
              call c(iorb,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j = binary_search(Hmap,k2)
              htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           if( (dmft_bath%v(1,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(ms)==1) )then
              call c(ms,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           !IMP DW <--> BATH DW
           if( (dmft_bath%v(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==1) .AND. (ib(ms+Ns)==0) )then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(ms+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
           if( (dmft_bath%v(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb+Ns)==0) .AND. (ib(ms+Ns)==1) )then
              call c(ms+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
              call sp_insert_element(spH0,htmp,impi,j)
           endif
        enddo
     enddo
     !
     !      ==> OFF-DIAGONAL-HYBRIDIZATION <== 
     !
     if((ed_mode=="nonsu2").and.(bath_type/="replica"))then
        do iorb=1,Norb
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              !
              ! IMP UP <--> BATH DW
              if( (ib(iorb)==1) .AND. (ib(ms+Ns)==0) )then
                 call c(iorb,m,k1,sg1)
                 call cdg(ms+Ns,k1,k2,sg2)
                 j = binary_search(Hmap,k2)
                 htmp = dmft_bath%u(1,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              if( (ib(iorb)==0) .AND. (ib(ms+Ns)==1) )then
                 call c(ms+Ns,m,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp = dmft_bath%u(1,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              ! IMP DW <--> BATH UP
              if( (ib(iorb+Ns)==1) .AND. (ib(ms)==0) )then
                 call c(iorb+Ns,m,k1,sg1)
                 call cdg(ms,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp = dmft_bath%u(Nspin,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              if( (ib(iorb+Ns)==0) .AND. (ib(ms)==1) )then
                 call c(ms,m,k1,sg1)
                 call cdg(iorb+Ns,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp = dmft_bath%u(Nspin,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
           enddo
        enddo
     endif
     !
     !
     !##############################################
     !#                                            #
     !#    ANOMALOUS PAIR-CREATION/DESTRUCTION     #
     !#                                            #
     !##############################################
     !
     !
     if(ed_mode=="superc")then
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              !\Delta_l c_{\up,ms} c_{\dw,ms}
              if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==1) .AND. (ib(ms+Ns)==1) )then
                 call c(ms,m,k1,sg1)
                 call c(ms+Ns,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
              !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
              if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==0) .AND. (ib(ms+Ns)==0) )then
                 call cdg(ms+Ns,m,k1,sg1)
                 call cdg(ms,k1,k2,sg2)
                 j=binary_search(Hmap,k2)
                 htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2 !
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
           enddo
        enddo
     endif

    
  enddo states

