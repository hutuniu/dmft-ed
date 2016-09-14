  !
  !
  call build_sector(isector,Hup,Hdw)
  !
  if(spH0%status)call sp_delete_matrix(spH0)
  if(spH0up%status)call sp_delete_matrix(spH0up)
  if(spH0dw%status)call sp_delete_matrix(spH0dw)
  !
  !UP:
  dimup   = getDimUp(isector)
  call sp_init_matrix(spH0up,dimup)
  !DW
  dimdw   = getDimDw(isector)
  call sp_init_matrix(spH0dw,dimdw)
  !FULL:
  dim  = dimup*dimdw
  call sp_init_matrix(spH0,dim)
  ishift=0
  !
  !
  !##############################################
  !#                                            #
  !#          LOCAL HAMILTONIAN PART            #
  !#                                            #
  !##############################################  
  do idw=1,dimdw
     mdw = Hdw%map(idw)
     ndw = bdecomp(mdw,Ns)
     !
     do iup=1,dimdw
        mup = Hup%map(iup)
        nup = bdecomp(mup,Ns)
        !
        i    = iup + (idw-1)*dimup
        impi = i - ishift
        !
        htmp=0d0
        !
        ! ==> CHEMICAL POTENTIAL 
        htmp = htmp - xmu*sum(nup(1:Norb) + ndw(1:Norb))
        ! ==> IMPURITY LOCAL ENERGY
        do iorb=1,Norb
           htmp = htmp &
                + impHloc(1,1,iorb,iorb)*nup(iorb) &
                + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb) &
                - xmu*(nup(iorb)+ndw(iorb))
        enddo
        ! ==> BATH HAMILTONIAN
        do iorb=1,size(dmft_bath%e,2) !irrred=Norb OR hybrid=1
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              htmp = htmp & 
                   + dmft_bath%e(1,iorb,kp)*nup(ms) &
                   + dmft_bath%e(Nspin,iorb,kp)*ndw(ms)
           enddo
        enddo
        !
        !     ==> IMPURITY INTERACTIONS (no SE & PH) <==
        ! 1. density-density interaction: same orbital, opposite spins
        ! \sum_i U_i*[ n_{i,up}*n_{i,dw} ]
        do iorb=1,Norb
           htmp = htmp + Uloc(iorb)*nup(iorb)*ndw(iorb)
        enddo
        if(Norb>1)then
           !2. density-density interaction: different orbitals, opposite spins
           !  U' sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ], w/ U ==  (U-2*Jh)
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
              enddo
           enddo
           !3. density-density interaction: different orbitals, parallel spins
           ! U''  \sum_{i<j}  [ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ], w/ U''== (U'-Jh) == (U-3*Jh)
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
        call sp_insert_element(spH0,htmp,impi,i)
        !
     enddo
  enddo





  !##############################################
  !#                                            #
  !#  SPIN-EXCHANGE (S-E) & PAIR-HOPPING (P-H)  #
  !#                                            #
  !##############################################
  !
  if(Jhflag.AND.Norb>1)then
     do idw=1,dimdw
        mdw = Hdw%map(idw)
        ndw = bdecomp(mdw,Ns)
        !
        do iup=1,dimup
           mup = Hup%map(iup)
           nup = bdecomp(mup,Ns)
           !
           i    = iup + (idw-1)*dimup
           impi = i - ishift
           !
           !    S-E: J c^+_iorb_up c_jorb_up * c^+_jorb_dw c_iorb_dw
           do iorb=1,Norb
              do jorb=1,Norb
                 if(iorb==jorb)cycle
                 Jcondition=(&
                      (nup(jorb)==1).AND.&
                      (ndw(iorb)==1).AND.&
                      (ndw(jorb)==0).AND.&
                      (nup(iorb)==0))
                 if(Jcondition)then
                    call c(iorb,mdw,k1,sg1) !DW
                    call cdg(jorb,k1,k2,sg2) !DW
                    jdw=binary_search(Hdw%map,k2)
                    call c(jorb,mup,k1,sg3) !UP
                    call cdg(iorb,k1,k2,sg4)    !UP
                    jup=binary_search(Hup%map,k2)
                    htmp = Jh*sg1*sg2*sg3*sg4
                    j = jup + (jdw-1)*dimup
                    call sp_insert_element(spH0,htmp,impi,j)
                 endif
              enddo
           enddo
           !
           !    P-H: J c^+_iorb_up. c_jorb_up * c^+_iorb_dw. c_jorb_dw     (i.ne.j) 
           do iorb=1,Norb
              do jorb=1,Norb
                 if(iorb==jorb)cycle
                 Jcondition=(&
                      (nup(jorb)==1).AND.&
                      (ndw(jorb)==1).AND.&
                      (ndw(iorb)==0).AND.&
                      (nup(iorb)==0))
                 if(Jcondition)then
                    call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                    call cdg(iorb,k1,k2,sg2)       !c^+_iorb_dw
                    jdw = binary_search(Hdw%map,k2)
                    call c(jorb,mup,k1,sg1)       !c_jorb_up
                    call cdg(iorb,k1,k2,sg4)       !c^+_iorb_up
                    jup = binary_search(Hup%map,k2)
                    htmp = Jh*sg1*sg2*sg3*sg4
                    j = jup + (jdw-1)*dimup
                    call sp_insert_element(spH0,htmp,impi,j)
                 endif
              enddo
           enddo
           !
        enddo
     enddo
  endif




  !##############################################
  !#                                            #
  !#          NON-LOCAL HAMILTONIAN             #
  !#                                            #
  !##############################################
  !
  !HYBRIDIZATION SPIN UP:
  do iup=1,dimup
     mup = Hup%map(iup)
     nup = bdecomp(mup,Ns)
     !
     !==> INTRA-IMPURITY 
     do iorb=1,Norb
        do jorb=1,Norb
           if((impHloc(1,1,iorb,jorb)/=0d0).AND.(nup(iorb)==0).AND.(nup(jorb)==1))then
              call c(jorb,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup=binary_search(Hup%map,k2)
              htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,impi_up,jup)
              !
           endif
           !
        enddo
     enddo
     !
     !==> IMPURITY-BATH
     do iorb=1,Norb
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==1) .AND. (nup(ms)==0) )then
              call c(iorb,mup,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              jup = binary_search(Hup%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,impi_up,jup)
              !
           endif
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==0) .AND. (nup(ms)==1) )then
              call c(ms,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup=binary_search(Hup%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,impi_up,jup)
              !
           endif
        enddo
     enddo
     !
  enddo
  !
  !
  !HYBRIDIZATION SPIN DW:
  do idw=1,dimdw
     mdw = Hdw%map(idw)
     ndw = bdecomp(mdw,Ns)
     !
     !==> INTRA-IMPURITY 
     do iorb=1,Norb
        do jorb=1,Norb
           if((impHloc(Nspin,Nspin,iorb,jorb)/=0d0).AND.(ndw(iorb)==0).AND.(ndw(jorb)==1))then
              call c(jorb,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw=binary_search(Hdw%map,k2)
              htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,impi_dw,jdw)
              !
           endif
        enddo
     enddo
     !
     !==> IMPURITY-BATH 
     do iorb=1,Norb
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==1) .AND. (ndw(ms)==0) )then
              call c(iorb,mdw,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              jdw=binary_search(Hdw%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,impi_dw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==0) .AND. (ndw(ms)==1) )then
              call c(ms,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw=binary_search(Hdw%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,impi_dw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo
