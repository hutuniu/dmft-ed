!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN
  USE CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private

  !Get sparse sector Hamiltonian
  public                       :: ed_buildH_d,ed_buildH_c

contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_buildH_d(isector,Hmat)
    real(8),dimension(:,:),optional    :: Hmat
    real(8),dimension(:,:),allocatable :: Hredux
    integer                            :: isector
    integer,dimension(:),allocatable   :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Nlevels)         :: ivec
    integer                            :: mpiQ,mpiR                
    integer                            :: dim,iup,idw
    integer                            :: i,hi,j,m,ms,iorb,jorb,ispin
    integer                            :: kp,k1,k2,k3,k4
    real(8)                            :: sg1,sg2,sg3,sg4
    real(8)                            :: htmp
    real(8),dimension(Norb)            :: nup,ndw
    real(8),dimension(Nspin,Norb)      :: eloc
    logical                            :: Jcondition
    integer                            :: first_state,last_state
    !
    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
    !
    first_state= 1
    last_state = dim
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    first_state= ED_MPI_ID*mpiQ+1
    last_state = (ED_MPI_ID+1)*mpiQ+mpiR
#else
    call sp_init_matrix(spH0,dim)
#endif
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=dreal(Hloc(ispin,ispin,iorb,iorb))
       enddo
    enddo
    !
    do i=first_state,last_state
#ifdef _MPI
       hi=i-ED_MPI_ID*mpiQ
#else
       hi=i
#endif
       m=Hmap(i)
       call bdecomp(m,ivec)
       htmp=0.d0
       do iorb=1,Norb
          nup(iorb)=real(ivec(iorb),8)
          ndw(iorb)=real(ivec(iorb+Ns),8)
       enddo
       !
       !LOCAL HAMILTONIAN PART:
       !local energies -mu*\sum_\a (n_{\a,up}+n_{\a,dw}) + \sum_{\a\s}Eloc_{\a\s} n_{\a\s}
       htmp = -xmu*(sum(nup)+sum(ndw))  + dot_product(eloc(1,:),nup)  + dot_product(eloc(Nspin,:),ndw)
       !
       !energy of the bath=\sum_\a\sum_k \e_{\a,k} n_{\a,k}
       do iorb=1,size(dmft_bath%e,2)
          do kp=1,Nbath
             ms=getBathStride(iorb,kp)
             htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ivec(ms),8) + dmft_bath%e(Nspin,iorb,kp)*real(ivec(ms+Ns),8)
          enddo
       enddo
       !
       !Density-density interaction: same orbital, opposite spins
       !\sum_\a \pm U_\a*(n_{\a,up}*n_{\a,dw})
       if(ed_mode=="superc") then
          htmp = htmp - dot_product(uloc,nup*ndw) 
       else
          htmp = htmp + dot_product(uloc,nup*ndw) 
          if(Norb>1)then
             !density-density interaction: different orbitals, opposite spins
             ! U'\sum_{i,j}[ n_{i,up}*n_{j,dw} +  n_{j,up}*n_{i,dw} ]
             do iorb=1,Norb         
                do jorb=iorb+1,Norb
                   htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
                enddo
             enddo
             !density-density interaction: different orbitals, parallel spins
             ! (U'-Jh)'sum_{i,j} [ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
                enddo
             enddo
          endif
       end if
       !
       !Hartree-shift: mu=0 corresponds to half-filling condition.
       !sum up the contributions to hartree terms from expanding (n_up-1/2)(n_dw-1/2)
       if(hfmode)then
          if(ed_mode=="superc")then
             htmp=htmp + 0.5d0*dot_product(uloc,nup+ndw) - 0.25d0*sum(uloc)
          else
             htmp=htmp - 0.5d0*dot_product(uloc,nup+ndw) + 0.25d0*sum(uloc)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      htmp=htmp-0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*Ust
                      htmp=htmp-0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*(Ust-Jh)
                   enddo
                enddo
             endif
          endif
       endif
       !
       call sp_insert_element(spH0,htmp,hi,i)
       !
       !NON-LOCAL HAMILTONIAN
       !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
       do iorb=1,Norb
          do jorb=1,Norb
             !SPIN UP
             if((ivec(iorb)==0).AND.(ivec(jorb)==1))then
                call c(jorb,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(1,1,iorb,jorb)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !SPIN DW
             if((ivec(iorb+Ns)==0).AND.(ivec(jorb+Ns)==1))then
                call c(jorb+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
          enddo
       enddo
       !==> HYBRIDIZATION TERMS II: same or different orbitals, opposite spins.
       !this is active only in the non-SU(2) invariant channel
       !\sum_{\a,\b} H_{\a\b;\s\s`} c^+_{\a,\s}c_{\b,\s`} + c^+_{\a,\s`}c_{\b,\s}
       if(ed_mode=="nonsu2")then
          do iorb=1,Norb
             do jorb=1,Norb
                !UP-DW
                if((ivec(iorb)==0).AND.(ivec(jorb+Ns)==1))then
                   call c(jorb+Ns,m,k1,sg1)
                   call cdg(iorb,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp = Hloc(1,Nspin,iorb,jorb)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                !DW-UP
                if((ivec(iorb+Ns)==0).AND.(ivec(jorb)==1))then
                   call c(jorb,m,k1,sg1)
                   call cdg(iorb+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp = Hloc(Nspin,1,iorb,jorb)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !
       !
       !==> INTERACTION TERMS
       if(Norb>1.AND.Jhflag)then
          !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
          !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
          !it shoud rather be (not ordered product) but this changes sign in the code,
          !so it is more natura to NORMAL order the products of operators:
          !S-E: J c^+_iorb_up c_iorb_dw   c^+_jorb_dw    c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb}  c_{iorb+Ns} c^+_{jorb+Ns}  c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ivec(jorb)==1).AND.&
                     (ivec(iorb+Ns)==1).AND.&
                     (ivec(jorb+Ns)==0).AND.&
                     (ivec(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
          !
          !PAIR-HOPPING (P-H) TERMS
          !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
          !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ivec(jorb)==1).AND.&
                     (ivec(jorb+Ns)==1).AND.&
                     (ivec(iorb+Ns)==0).AND.&
                     (ivec(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !
       !IMPURITY - BATH TERMS
       !==> Hybridizations between imp. and bath, same spins
       do iorb=1,Norb
          do kp=1,Nbath
             ms=getBathStride(iorb,kp)
             !SPIN UP
             if(ivec(iorb) == 1 .AND. ivec(ms) == 0)then
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j = binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             if(ivec(iorb) == 0 .AND. ivec(ms) == 1)then
                call c(ms,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !SPIN DOWN
             if(ivec(iorb+Ns) == 1 .AND. ivec(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             if(ivec(iorb+Ns) == 0 .AND. ivec(ms+Ns) == 1)then
                call c(ms+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
          enddo
       enddo
       !
       !==> Hybridizations between imp. and bath,  spin flip. 
       if(ed_mode=="nonsu2")then
          do iorb=1,Norb
             do kp=1,Nbath
                ms=getBathStride(iorb,kp)
                !IMP UP <--> BATH DW
                if(ivec(iorb) == 1 .AND. ivec(ms+Ns) == 0)then
                   call c(iorb,m,k1,sg1)
                   call cdg(ms+Ns,k1,k2,sg2)
                   j = binary_search(Hmap,k2)
                   htmp = dmft_bath%w(1,Nspin,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                if(ivec(iorb) == 0 .AND. ivec(ms+Ns) == 1)then
                   call c(ms+Ns,m,k1,sg1)
                   call cdg(iorb,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp = dmft_bath%w(Nspin,1,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                !IMP DW <--> BATH UP
                if(ivec(iorb+Ns) == 1 .AND. ivec(ms) == 0)then
                   call c(iorb+Ns,m,k1,sg1)
                   call cdg(ms,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%w(Nspin,1,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                if(ivec(iorb+Ns) == 0 .AND. ivec(ms) == 1)then
                   call c(ms,m,k1,sg1)
                   call cdg(iorb+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%w(1,Nspin,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !
       !
       !PAIR-CREATION & DESTRUCTION (SUPERCONDUCTING CHANNEL)
       if(ed_mode=="superc")then
          do iorb=1,size(dmft_bath%e,2)
             do kp=1,Nbath
                ms=getBathStride(iorb,kp)
                !\Delta_l c_{\up,ms} c_{\dw,ms}
                if(ivec(ms)==1 .AND. ivec(ms+Ns)==1)then
                   call c(ms,m,k1,sg1)
                   call c(ms+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                !\Delta_l cdg_{\dw,ms} cdg_{\up,ms}
                if(ivec(ms)==0 .AND. ivec(ms+Ns)==0)then
                   call cdg(ms+Ns,m,k1,sg1)
                   call cdg(ms,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2 !
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !
    enddo
    !
    deallocate(Hmap)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "ED_HAMILTONIAN/ed_buildH_d: wrong dimensions in Hmat"
#ifdef _MPI
       allocate(Hredux(dim,dim))
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine ed_buildH_d












  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine ed_buildH_c(isector,Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(:),allocatable      :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Nlevels)            :: ivec
    integer                               :: mpiQ,mpiR                
    integer                               :: dim,iup,idw
    integer                               :: hi,i,j,m,ms,iorb,jorb,ispin
    integer                               :: kp,k1,k2,k3,k4
    real(8)                               :: sg1,sg2,sg3,sg4
    complex(8)                            :: htmp
    real(8),dimension(Norb)               :: nup,ndw
    complex(8),dimension(Nspin,Norb)      :: eloc
    logical                               :: Jcondition
    integer                               :: first_state,last_state
    !
    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
    !
    first_state= 1
    last_state = dim
    !
    if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
    mpiQ = dim/ED_MPI_SIZE
    mpiR = 0
    if(ED_MPI_ID==(ED_MPI_SIZE-1))mpiR=mod(dim,ED_MPI_SIZE)
    call sp_init_matrix(spH0,mpiQ+mpiR)
    first_state= ED_MPI_ID*mpiQ+1
    last_state = (ED_MPI_ID+1)*mpiQ+mpiR
#else
    call sp_init_matrix(spH0,dim)
#endif
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=Hloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    do i=first_state,last_state
#ifdef _MPI
       hi=i-ED_MPI_ID*mpiQ
#else
       hi=i
#endif
       m=Hmap(i)
       call bdecomp(m,ivec)
       htmp=0.d0
       do iorb=1,Norb
          nup(iorb)=real(ivec(iorb),8)
          ndw(iorb)=real(ivec(iorb+Ns),8)
       enddo
       !
       !LOCAL HAMILTONIAN PART:
       !local energies
       htmp = -xmu*(sum(nup)+sum(ndw))  + &
            dot_product(eloc(1,:),nup)  + &
            dot_product(eloc(Nspin,:),ndw)
       !Density-density interaction: same orbital, opposite spins
       htmp = htmp + dot_product(uloc,nup*ndw)!=\sum=i U_i*(n_u*n_d)_i
       if(hfmode)htmp=htmp - 0.5d0*dot_product(uloc,nup+ndw) + 0.25d0*sum(uloc)
       if(Norb>1)then
          !density-density interaction: different orbitals, opposite spins
          do iorb=1,Norb         ! n_up_i*n_dn_j
             do jorb=iorb+1,Norb ! n_up_j*n_dn_i
                htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
             enddo
          enddo
          !density-density interaction: different orbitals, parallel spins
          !Jhund effect: U``=U`-J smallest of the interactions
          do iorb=1,Norb         ! n_up_i*n_up_j
             do jorb=iorb+1,Norb ! n_dn_i*n_dn_j
                htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
             enddo
          enddo
       endif
       !
       !Hbath: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
       do iorb=1,size(dmft_bath%e,2)
          do kp=1,Nbath
             ! ms=Norb+(iorb-1)*Nbath + kp
             ! if(bath_type=='hybrid')ms=Norb+kp
             ms=getBathStride(iorb,kp)
             htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ivec(ms),8) + &
                  dmft_bath%e(Nspin,iorb,kp)*real(ivec(ms+Ns),8)
          enddo
       enddo
       !
       !
       call sp_insert_element(spH0,htmp,hi,i)
       !
       if(Norb>1.AND.Jhflag)then
          !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
          !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
          !it shoud rather be (not ordered product) but this changes sign in the code,
          !so it is more natura to NORMAL order the products of operators:
          !S-E: J c^+_iorb_up c_iorb_dw   c^+_jorb_dw    c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb}  c_{iorb+Ns} c^+_{jorb+Ns}  c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ivec(jorb)==1).AND.&
                     (ivec(iorb+Ns)==1).AND.&
                     (ivec(jorb+Ns)==0).AND.&
                     (ivec(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
          !
          !PAIR-HOPPING (P-H) TERMS
          !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
          !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ivec(jorb)==1).AND.&
                     (ivec(jorb+Ns)==1).AND.&
                     (ivec(iorb+Ns)==0).AND.&
                     (ivec(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !

       !LOCAL HYBRIDIZATION
       do iorb=1,Norb
          do jorb=1,Norb
             !SPIN UP
             if((ivec(iorb)==0).AND.(ivec(jorb)==1))then
                call c(jorb,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(1,1,iorb,jorb)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !SPIN DW
             if((ivec(iorb+Ns)==0).AND.(ivec(jorb+Ns)==1))then
                call c(jorb+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
          enddo
       enddo

       !IMP-BATH HYBRIDIZATION
       do iorb=1,Norb
          do kp=1,Nbath
             ms=getBathStride(iorb,kp)
             if(ivec(iorb) == 1 .AND. ivec(ms) == 0)then
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j = binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !
             if(ivec(iorb) == 0 .AND. ivec(ms) == 1)then
                call c(ms,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !
             if(ivec(iorb+Ns) == 1 .AND. ivec(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
             !
             if(ivec(iorb+Ns) == 0 .AND. ivec(ms+Ns) == 1)then
                call c(ms+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                call sp_insert_element(spH0,htmp,hi,j)
             endif
          enddo
       enddo


       if(ed_supercond)then
          !Anomalous pair-creation/destruction
          do iorb=1,size(dmft_bath%e,2)
             do kp=1,Nbath
                ms=getBathStride(iorb,kp)
                !\Delta_l c_{\up,ms} c_{\dw,ms}
                if(ivec(ms)==1 .AND. ivec(ms+Ns)==1)then
                   call c(ms,m,k1,sg1)
                   call c(ms+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
                !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
                if(ivec(ms)==0 .AND. ivec(ms+Ns)==0)then
                   call cdg(ms,m,k1,sg1)
                   call cdg(ms+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   htmp=dmft_bath%d(1,iorb,kp)*sg1*sg2 !
                   call sp_insert_element(spH0,htmp,hi,j)
                endif
             enddo
          enddo
       endif
       !
    enddo
    !
    deallocate(Hmap)
    !
    if(present(Hmat))then
       if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "ED_HAMILTONIAN/ed_buildH_d: wrong dimensions in Hmat"
#ifdef _MPI
       allocate(Hredux(dim,dim))
       call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
       call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,ED_MPI_ERR)
#else
       call sp_dump_matrix(spH0,Hmat)
#endif
    endif
    !
  end subroutine ed_buildH_c


end MODULE ED_HAMILTONIAN
