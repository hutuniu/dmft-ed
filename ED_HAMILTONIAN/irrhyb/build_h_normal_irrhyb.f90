  !
  !
  call build_sector(isector,Hup,Hdw)
  !
  ! if(spH0up%status)call sp_delete_matrix(spH0up)
  ! if(spH0dw%status)call sp_delete_matrix(spH0dw)
  if(spH0%status)call sp_delete_matrix(spH0)
  !
  ! !UP:
  dimup   = getDimUp(isector)
  ! mpiQ_up = dimup/ED_MPI_SIZE
  ! mpiR_up = 0
  ! if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR_up=mod(dimup,ED_MPI_SIZE)
  ! call sp_init_matrix(spH0up,mpiQ_up + mpiR_up)
  ! ishift_up      = ED_MPI_ID*mpiQ_up
  ! first_state_up = ED_MPI_ID*mpiQ_up + 1
  ! last_state_up  = (ED_MPI_ID+1)*mpiQ_up + mpiR_up
  ! !DW
  dimdw   = getDimDw(isector)
  ! mpiQ_dw = dimdw/ED_MPI_SIZE
  ! mpiR_dw = 0
  ! if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR_dw=mod(dimdw,ED_MPI_SIZE)
  ! call sp_init_matrix(spH0dw,mpiQ_dw + mpiR_dw)
  ! ishift_dw      = ED_MPI_ID*mpiQ_dw
  ! first_state_dw = ED_MPI_ID*mpiQ_dw + 1
  ! last_state_dw  = (ED_MPI_ID+1)*mpiQ_dw + mpiR_dw
  !FULL:
  dim  = dimup*dimdw
  ! mpiQ = dim/ED_MPI_SIZE
  ! mpiR = 0
  ! if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
  call sp_init_matrix(spH0,dim)!mpiQ + mpiR)
  ! ishift      = ED_MPI_ID*mpiQ
  ! first_state = ED_MPI_ID*mpiQ + 1
  ! last_state  = (ED_MPI_ID+1)*mpiQ + mpiR
  !
  !
  !##############################################
  !#                                            #
  !#          LOCAL HAMILTONIAN PART            #
  !#                                            #
  !##############################################  
  do idw=1,dimdw!first_state_dw, last_state_dw
     mdw = Hdw%map(idw)
     ndw = bdecomp(mdw,Ns)
     !impi_dw = idw-ishift_dw
     !
     do iup=1,dimup!first_state_up, last_state_up
        mup = Hup%map(iup)
        nup = bdecomp(mup,Ns)
        !impi_up = iup-ishift_up
        !
        !
        i = iup + (idw-1)*dimup
        impi = i                !-ishift
        !
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
  !#            NON-LOCAL HAMILTONIAN           #
  !#                                            #
  !##############################################
  !


  do idw=1,dimdw!first_state_dw, last_state_dw
     mdw = Hdw%map(idw)
     ndw = bdecomp(mdw,Ns)
     !
     do iup=1,dimup!first_state_up, last_state_up
        mup = Hup%map(iup)
        nup = bdecomp(mup,Ns)
        !
        m = mup + mdw*2**Ns
        i = iup + (idw-1)*dimup
        impi = i                !-ishift
        !
        !==> INTRA-IMPURITY HYBRIDIZATION
        do iorb=1,Norb
           do jorb=1,Norb
              !SPIN UP
              if((impHloc(1,1,iorb,jorb)/=0d0).AND.(nup(iorb)==0).AND.(nup(jorb)==1))then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup=binary_search(Hup%map,k2)
                 j = jup + (idw-1)*dimup
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
              !
              !SPIN DW
              if((impHloc(Nspin,Nspin,iorb,jorb)/=0d0).AND.(ndw(iorb)==0).AND.(ndw(jorb)==1))then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw=binary_search(Hdw%map,k2)
                 j  = iup + (jdw-1)*dimup
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
           enddo
        enddo
        !
        !==> IMPURITY-BATH HYBRIDIZATION
        do iorb=1,Norb
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              !
              !SPIN UP
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==1) .AND. (nup(ms)==0) )then
                 call c(iorb,mup,k1,sg1)
                 call cdg(ms,k1,k2,sg2)
                 jup = binary_search(Hup%map,k2)
                 j   = jup + (idw-1)*dimup
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==0) .AND. (nup(ms)==1) )then
                 call c(ms,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup=binary_search(Hup%map,k2)
                 j   = jup + (idw-1)*dimup
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              !
              !SPIN DW
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==1) .AND. (ndw(ms)==0) )then
                 call c(iorb,mdw,k1,sg1)
                 call cdg(ms,k1,k2,sg2)
                 jdw=binary_search(Hdw%map,k2)
                 j  = iup + (jdw-1)*dimup
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ndw(iorb)==0) .AND. (ndw(ms)==1) )then
                 call c(ms,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jdw=binary_search(Hdw%map,k2)
                 j  = iup + (jdw-1)*dimup
                 htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 call sp_insert_element(spH0,htmp,impi,j)
              endif
           enddo
        enddo
        !
     enddo
  enddo
