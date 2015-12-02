!########################################################################
!purpose  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_ENERGY
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_AUX_FUNX
  USE ED_MATVEC
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_TIMER
  USE SF_LINALG, only: inv,eye,zeye
  USE SF_MISC,  only:assert_shape
  implicit none
  private


  interface ed_kinetic_energy
     module procedure kinetic_energy_impurity_normal_main
     module procedure kinetic_energy_impurity_superc_main
     module procedure kinetic_energy_impurity_normal_1B
     module procedure kinetic_energy_impurity_normal_MB
     module procedure kinetic_energy_impurity_superc_1B
     module procedure kinetic_energy_impurity_superc_MB
  end interface ed_kinetic_energy

  interface ed_kinetic_energy_lattice
     module procedure kinetic_energy_lattice_normal_main
     module procedure kinetic_energy_lattice_superc_main
     module procedure kinetic_energy_lattice_normal_1
     module procedure kinetic_energy_lattice_normal_2
     module procedure kinetic_energy_lattice_normal_1B
     module procedure kinetic_energy_lattice_superc_1
     module procedure kinetic_energy_lattice_superc_2
     module procedure kinetic_energy_lattice_superc_1B
  end interface ed_kinetic_energy_lattice

  !PUBLIC in DMFT
  public :: ed_kinetic_energy
  public :: ed_kinetic_energy_lattice

  !INTERNAL in DMFT
  public :: local_energy_impurity

  real(8),dimension(:),allocatable        :: wm

contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity()
    integer,dimension(Nlevels)                        :: ib
    integer                                           :: i,j
    integer                                           :: izero
    integer                                           :: isector
    integer                                           :: idim
    integer                                           :: iorb,jorb,ispin
    integer                                           :: numstates
    integer                                           :: m,k1,k2,k3,k4
    real(8)                                           :: sg1,sg2,sg3,sg4
    real(8)                                           :: Egs,gs_weight
    real(8)                                           :: Ei
    real(8)                                           :: peso
    real(8)                                           :: norm
    real(8),dimension(Norb)                           :: nup,ndw
    real(8),dimension(Nspin,Norb)                     :: eloc
    real(8),dimension(:),pointer                      :: gsvec
    complex(8),dimension(:),pointer                   :: gscvec
    integer,allocatable,dimension(:)                  :: Hmap
    logical                                           :: Jcondition
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    numstates=state_list%size
    do izero=1,numstates
       isector = es_return_sector(state_list,izero)
       Ei      = es_return_energy(state_list,izero)
       idim    = getdim(isector)
       !
       if(ed_type=='d')then
          gsvec  => es_return_vector(state_list,izero)
          norm=sqrt(dot_product(gsvec,gsvec))
       elseif(ed_type=='c')then
          gscvec  => es_return_cvector(state_list,izero)
          norm=sqrt(dot_product(gscvec,gscvec))
       endif
       if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       allocate(Hmap(idim))
       call build_sector(isector,Hmap)
       !
       do i=1,idim
          m=Hmap(i)
          call bdecomp(m,ib)
          !
          if(ed_type=='d')then
             gs_weight=peso*gsvec(i)**2
          elseif(ed_type=='c')then
             gs_weight=peso*abs(gscvec(i))**2
          endif
          !
          !Get operators:
          do iorb=1,Norb
             nup(iorb)= dble(ib(iorb))
             ndw(iorb)= dble(ib(iorb+Ns))
          enddo
          !
          !start evaluating the Tr(H_loc) to estimate potential energy
          !
          !LOCAL ENERGY
          ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup)*gs_weight + dot_product(eloc(Nspin,:),ndw)*gs_weight
          !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
          do iorb=1,Norb
             do jorb=1,Norb
                !SPIN UP
                if((ib(iorb)==0).AND.(ib(jorb)==1))then
                   call c(jorb,m,k1,sg1)
                   call cdg(iorb,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*gs_weight
                endif
                !SPIN DW
                if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                   call c(jorb+Ns,m,k1,sg1)
                   call cdg(iorb+Ns,k1,k2,sg2)
                   j=binary_search(Hmap,k2)
                   ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*gs_weight
                endif
             enddo
          enddo
          !==> HYBRIDIZATION TERMS II: same or different orbitals, opposite spins.
          if(ed_mode=="nonsu2")then
             do iorb=1,Norb
                do jorb=1,Norb
                   !UP-DW
                   if((impHloc(1,Nspin,iorb,jorb)/=zero).AND.(ib(iorb)==0).AND.(ib(jorb+Ns)==1))then
                      call c(jorb+Ns,m,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      j=binary_search(Hmap,k2)
                      ed_Eknot = ed_Eknot + impHloc(1,Nspin,iorb,jorb)*sg1*sg2*gs_weight
                   endif
                   !DW-UP
                   if((impHloc(Nspin,1,iorb,jorb)/=zero).AND.(ib(iorb+Ns)==0).AND.(ib(jorb)==1))then
                      call c(jorb,m,k1,sg1)
                      call cdg(iorb+Ns,k1,k2,sg2)
                      j=binary_search(Hmap,k2)
                      ed_Eknot = ed_Eknot + impHloc(Nspin,1,iorb,jorb)*sg1*sg2*gs_weight
                   endif
                enddo
             enddo
          endif
          !
          !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
          !Euloc=\sum=i U_i*(n_u*n_d)_i
          !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
          do iorb=1,Norb
             ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*gs_weight
          enddo
          !
          !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
          !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
          !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                enddo
             enddo
          endif
          !
          !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
          !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                enddo
             enddo
          endif
          !
          !SPIN-EXCHANGE (S-E) TERMS
          !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
          if(Norb>1.AND.Jhflag)then
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=((iorb/=jorb).AND.&
                        (ib(jorb)==1)      .AND.&
                        (ib(iorb+Ns)==1)   .AND.&
                        (ib(jorb+Ns)==0)   .AND.&
                        (ib(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,m,k1,sg1)
                      call c(iorb+Ns,k1,k2,sg2)
                      call cdg(jorb+Ns,k2,k3,sg3)
                      call cdg(iorb,k3,k4,sg4)
                      j=binary_search(Hmap,k4)
                      ed_Epot = ed_Epot + Jh*sg1*sg2*sg3*sg4*gs_weight
                      ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*gs_weight
                   endif
                enddo
             enddo
          endif
          !
          !PAIR-HOPPING (P-H) TERMS
          !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
          !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
          if(Norb>1.AND.Jhflag)then
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=((iorb/=jorb).AND.&
                        (ib(jorb)==1)      .AND.&
                        (ib(jorb+Ns)==1)   .AND.&
                        (ib(iorb+Ns)==0)   .AND.&
                        (ib(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,m,k1,sg1)
                      call c(jorb+Ns,k1,k2,sg2)
                      call cdg(iorb+Ns,k2,k3,sg3)
                      call cdg(iorb,k3,k4,sg4)
                      j=binary_search(Hmap,k4)
                      ed_Epot = ed_Epot + Jh*sg1*sg2*sg3*sg4*gs_weight
                      ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*gs_weight
                   endif
                enddo
             enddo
          endif
          !
          !HARTREE-TERMS CONTRIBUTION:
          if(hfmode)then
             !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
             do iorb=1,Norb
                ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*uloc(iorb)*gs_weight
             enddo
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*Ust*gs_weight
                      ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*(Ust-Jh)*gs_weight
                   enddo
                enddo
             endif
          endif
       enddo
       if(associated(gsvec))nullify(gsvec)
       if(associated(gscvec))nullify(gscvec)
       deallocate(Hmap)
    enddo
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ED_MPI_ID==0)then
       if(ed_verbose<0)then
          write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
          write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
          write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
          write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
          write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
          write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
          write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
          write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
       endif
    endif
    !
    !
    call write_energy_info()
    call write_energy()
  end subroutine local_energy_impurity








  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the lattice model, given 
  ! the Hamiltonian matrix H(k) and the DMFT self-energy Sigma.
  ! The main routine accept self-energy as:
  ! - Sigma: [Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [L]
  ! - Sigma: [Nspin][Nspin][Norb][Norb][L]
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_impurity_normal_main(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)                 :: Hk ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:)                 :: Sigma ![Nspin][Nspin][Norb][Norb][Liw]
    real(8),dimension(size(Hk,3))               :: Wtk   ![Lk]
    integer                                     :: Lk,Nso,Liw
    integer                                     :: i,ik
    !
    real(8),dimension(size(Hk,1),size(Hk,1))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Gk,Tk
    real(8)                                     :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                     :: H0,Hl,ed_Ekin,ed_Eloc
    real(8)                                     :: Eout(2)
    !
    Nso = size(Hk,1)
    Lk  = size(Hk,3)
    Liw = size(Sigma,3)
    call assert_shape(Hk,[Nso,Nso,Lk],"kinetic_energy_impurity_normal_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"kinetic_energy_impurity_normal_main","Sigma")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)call print_hloc(Hloc)
    !
    if(ED_MPI_ID==0)call start_timer()
    H0=0d0
    Hl=0d0
    do ik=1,Lk
       Ak = Hk(:,:,ik) - Hloc(:,:)
       Bk =-Hk(:,:,ik) - Sigma_HF(:,:)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*eye(Nso) - Sigma(:,:,i) - Hk(:,:,ik) 
          select case(Nso)
          case default
             call inv(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = eye(Nso)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak  ,Gk - Tk)
          Dk = matmul(Hloc,Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,Nso)
          Hl = Hl + Wtk(ik)*trace_matrix(Dk,Nso)
       enddo
       if(ED_MPI_ID==0)call eta(ik,Lk)
    enddo
    if(ED_MPI_ID==0)call stop_timer()
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    Hl=Hl/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nso)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    Lail0=spin_degeneracy*Lail0
    Lail1=spin_degeneracy*Lail1*beta
    ed_Ekin=H0+Tail0+Tail1
    ed_Eloc=Hl+Lail0+Lail1
    Eout = [ed_Ekin,ed_Eloc]
    deallocate(wm)
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_impurity_normal_main

  function kinetic_energy_impurity_normal_1B(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(:)               :: Hk
    real(8),dimension(size(Hk))           :: Wtk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8),dimension(2)                  :: Eout
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    Eout = kinetic_energy_impurity_normal_main(Hk_,Wtk,Sigma_)
  end function kinetic_energy_impurity_normal_1B

  function kinetic_energy_impurity_normal_MB(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)             :: Hk
    real(8),dimension(size(Hk,3))           :: Wtk
    complex(8),dimension(:,:,:,:,:)         :: Sigma ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable :: Sigma_
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,io,jo
    real(8),dimension(2)                    :: Eout
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"kinetic_energy_impurity_normal_MB","Sigma")
    allocate(Sigma_(Nso,Nso,Lmats))
    do i=1,Lmats
       Sigma_(:,:,i) = nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
    enddo
    Eout = kinetic_energy_impurity_normal_main(Hk,Wtk,Sigma_)
  end function kinetic_energy_impurity_normal_MB






  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the lattice case Superconducting, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energy Sigma.
  ! The main routine accept self-energy as
  ! - Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
  ! - Sigma: [L]
  ! - Sigma: [Nspin][Nspin][Norb][Norb][L]
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_impurity_superc_main(Hk,Wtk,Sigma,SigmaA) result(Eout)
    integer                                         :: Lk,Nso,Liw
    integer                                         :: i,ik,iorb,jorb,inambu,jnambu,n,m
    complex(8),dimension(:,:,:)                     :: Hk
    complex(8),dimension(:,:,:)                     :: Sigma,SigmaA
    real(8),dimension(size(Hk,3))                   :: Wtk
    !
    real(8),dimension(size(Hk,1),size(Hk,1))        :: Sigma_HF,SigmaA_HF
    complex(8),dimension(size(Hk,1),size(Hk,1))     :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1))     :: Gk,Tk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,1)) :: Gk_Nambu
    complex(8),dimension(2,2)                       :: Gk_Nambu_ij
    !
    real(8)                                         :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                         :: H0,Hl,ed_Ekin,ed_Eloc
    real(8)                                         :: Eout(2)
    !
    Nso = size(Hk,1)
    Lk  = size(Hk,3)
    Liw = size(Sigma,3)
    call assert_shape(Hk,[Nso,Nso,Lk],"kinetic_energy_impurity_superc_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"kinetic_energy_impurity_superc_main","Sigma")
    call assert_shape(SigmaA,[Nso,Nso,Liw],"kinetic_energy_impurity_superc_main","Self")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !    
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)call print_hloc(Hloc)
    !
    H0=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)-Hloc
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       do i=1,Liw
          Gk=zero          
          do iorb=1,Nso
             do jorb=1,Nso
                Gk_Nambu_ij=zero
                Gk_Nambu_ij(1,1) =  -Hk(iorb,jorb,ik) - Sigma(iorb,jorb,i)
                Gk_Nambu_ij(1,2) =                    - SigmaA(iorb,jorb,i)
                Gk_Nambu_ij(2,1) =                    - SigmaA(iorb,jorb,i)
                Gk_Nambu_ij(2,2) =   Hk(iorb,jorb,ik) + conjg(Sigma(iorb,jorb,i))!-conjg(Gk_Nambu_ij(1,1))
                if(iorb==jorb) then
                   Gk_Nambu_ij(1,1) = Gk_Nambu_ij(1,1) + xi*wm(i) + xmu
                   Gk_Nambu_ij(2,2) = Gk_Nambu_ij(1,1) + xi*wm(i) - xmu
                end if
                do inambu=1,2
                   do jnambu=1,2
                      m=(inambu-1)*Nso + iorb
                      n=(jnambu-1)*Nso + jorb
                      Gk_nambu(m,n)=Gk_nambu_ij(inambu,jnambu)
                   enddo
                enddo
             enddo
          enddo
          call inv(Gk_Nambu)
          inambu=1
          jnambu=1
          do iorb=1,Nso
             do jorb=1,Nso
                m=(inambu-1)*Nso + iorb
                n=(jnambu-1)*Nso + jorb
                Gk(iorb,jorb) =  Gk_Nambu(m,n)
             enddo
          enddo
          Tk = eye(Nso)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak,Gk - Tk)
          Dk = matmul(Hloc,Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,Nso)
          Hl = Hl + Wtk(ik)*trace_matrix(Dk,Nso)
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    Hl=Hl/beta*2.d0*spin_degeneracy          
    !
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc,Nso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nso)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    Lail0=spin_degeneracy*Lail0
    Lail1=spin_degeneracy*Lail1*beta
    ed_Ekin=H0+Tail0+Tail1
    ed_Eloc=Hl+Lail0+Lail1
    deallocate(wm)
    Eout = [ed_Ekin,ed_Eloc]
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_impurity_superc_main

  function kinetic_energy_impurity_superc_1B(Hk,Wtk,Sigma,Self) result(Eout)
    real(8),dimension(:)                  :: Hk
    real(8),dimension(size(Hk))           :: Wtk
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(size(Sigma))     :: Self
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Sigma)) :: Self_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8)                               :: Eout(2)
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    Eout = kinetic_energy_impurity_superc_main(Hk_,Wtk,Sigma_,Self_)
  end function kinetic_energy_impurity_superc_1B

  function kinetic_energy_impurity_superc_MB(Hk,Wtk,Sigma,Self) result(Eout)
    complex(8),dimension(:,:,:)             :: Hk
    real(8),dimension(size(Hk,3))           :: Wtk
    complex(8),dimension(:,:,:,:,:)         :: Sigma
    complex(8),dimension(:,:,:,:,:)         :: Self
    complex(8),dimension(:,:,:),allocatable :: Sigma_
    complex(8),dimension(:,:,:),allocatable :: Self_
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,io,jo
    real(8)                                 :: Eout(2)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"kinetic_energy_impurity_superc_MB","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"kinetic_energy_impurity_superc_MB","Self")    
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    do i=1,Lmats
       Sigma_(:,:,i) = nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
       Self_(:,:,i)  = nn2so_reshape(Self(:,:,:,:,i),Nspin,Norb)
    enddo
    Eout = kinetic_energy_impurity_superc_main(Hk,Wtk,Sigma_,Self_)
  end function kinetic_energy_impurity_superc_MB

















  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the general lattice case, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energy Sigma.
  ! The main routine accept:
  ! - Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
  ! - Sigma: [Nlat][L]
  ! - Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_lattice_normal_main(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                                   :: wtk       ! [Nk]
    complex(8),dimension(:,:,:)                                     :: Sigma     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    !aux
    integer                                                         :: Lk,Nlso,Liw
    integer                                                         :: ik
    integer                                                         :: i,iorb,ilat,ispin,io,is
    integer                                                         :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,1))                        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Ak,Bk,Ck,Dk,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Gk
    real(8)                                                         :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                                         :: H0,H0k,H0ktmp,Hl,Hlk,Hlktmp
    real(8)                                                         :: ed_Ekin_lattice,ed_Eloc_lattice
    real(8)                                                         :: Eout(2)
    !
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Liw  = size(Sigma,3)
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_normal_main","Hk")
    call assert_shape(Sigma,[Nlso,Nlso,Liw],"kinetic_energy_lattice_normal_main","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
    ! ...we should make single-site routines uniform with this procedure...
    Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
    Hloc=0.d0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hloc(is,js)=Hloc_tmp(is,js) 
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)call print_hloc(Hloc)
    !
    !Get HF part of the self-energy
    Sigma_HF(:,:) = dreal(Sigma(:,:,Liw))
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    ed_Eloc_lattice = 0d0
    H0              = 0d0
    Hl              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    =  Hk(:,:,ik) - Hloc(:,:)
       Bk    = -Hk(:,:,ik) - Sigma_HF(:,:)
       H0ktmp= 0d0
       Hlktmp= 0d0
       H0k   = 0d0
       Hlk   = 0d0
       do i=1+mpiID,Lmats,mpiSIZE
          Gk = (xi*wm(i)+xmu)*eye(Nlso) - Sigma(:,:,i) - Hk(:,:,ik)
          call inv(Gk(:,:))
          Tk = zeye(Nlso)/(xi*wm(i)) - Bk/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck,Nlso)
          Hlktmp = Hlktmp + Wtk(ik)*trace_matrix(Dk,Nlso)
       enddo
#ifdef _MPI_INEQ
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       call MPI_ALLREDUCE(Hlktmp,Hlk,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
       H0k=H0ktmp
       Hlk=Hlktmp
#endif
       H0 = H0 + H0k
       Hl = Hl + Hlk
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    Hl = Hl/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    =  Hk(:,:,ik) - Hloc(:,:)
       Bk    = -Hk(:,:,ik) - Sigma_HF(:,:)
       Ck= matmul(Ak,Bk)
       Dk= matmul(Hloc,Bk)
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nlso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nlso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nlso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nlso)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    Lail0=spin_degeneracy*Lail0
    Lail1=spin_degeneracy*Lail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Eloc_lattice=Hl+Lail0+Lail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
    ed_Eloc_lattice=ed_Eloc_lattice/dble(Nlat)
    Eout = [ed_Ekin_lattice,ed_Eloc_lattice]
    deallocate(wm)
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_lattice_normal_main

  function kinetic_energy_lattice_normal_1(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:)                             :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,4)) :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                                   :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
    integer                                                   :: Nlat,Nlso,Nso,Lk,Liw
    real(8)                                                   :: ed_Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Nso  = size(Sigma,2)
    Lk   = size(Hk,3)
    Liw  = size(Sigma,4)
    Nlso = Nlat*Nso
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_normal_1","Hk")
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"kinetic_energy_lattice_normal_1","Sigma")
    Sigma_ = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Sigma_(is,js,:) = Sigma(ilat,io,jo,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_normal_main(Hk,Wtk,Sigma_)
  end function kinetic_energy_lattice_normal_1

  function kinetic_energy_lattice_normal_2(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)                         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,6)) :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                                   :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
    integer                                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                                   :: ed_Ekin_lattice(2)
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_normal_2","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_normal_2","Sigma")
    Sigma_ = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Sigma_(is,js,:) = Sigma(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_normal_main(Hk,Wtk,Sigma_)
  end function kinetic_energy_lattice_normal_2

  function kinetic_energy_lattice_normal_1B(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                             :: wtk    ! [Nk]
    complex(8),dimension(:,:)                                 :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,2)) :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb]][L]
    integer                                                   :: ilat
    integer                                                   :: Nlat,Nlso,Nso,Lk
    real(8)                                                   :: ed_Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Nlso = Nlat*1*1
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_normal_1B","Hk")
    Sigma_ = zero
    do ilat=1,Nlat
       Sigma_(ilat,ilat,:) = Sigma(ilat,:)
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_normal_main(Hk,Wtk,Sigma_)
  end function kinetic_energy_lattice_normal_1B










  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the general lattice case, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energy Sigma.
  ! The main routine accept:
  ! - Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
  ! - Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [Nlat][Nspin][Nspin][Norb][Norb][L]
  ! - Sigma: [Nlat][L]
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_lattice_superc_main(Hk,Wtk,Sigma,SigmaA) result(Eout)
    complex(8),dimension(:,:,:)                     :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                   :: wtk    ! [Nk]
    complex(8),dimension(:,:,:)                     :: Sigma  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:)                     :: SigmaA ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]  
    !aux
    integer                                         :: Lk,Nlso,Liw
    integer                                         :: ik
    integer                                         :: i,iorb,ilat,ispin,io,is
    integer                                         :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Ak,Bk,Ck,Dk,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Gk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2)) :: Gknambu
    real(8)                                         :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                         :: H0,H0k,H0ktmp,Hl,Hlk,Hlktmp
    real(8)                                         :: ed_Ekin_lattice,ed_Eloc_lattice
    real(8)                                         :: Eout(2)
    !Get generalized Lattice-Spin-Orbital index
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Liw  = size(Sigma,3)
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_superc_main","Hk")
    call assert_shape(Sigma,[Nlso,Nlso,Liw],"kinetic_energy_lattice_superc_main","Sigma")
    call assert_shape(SigmaA,[Nlso,Nlso,Liw],"kinetic_energy_lattice_superc_main","Self")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
    ! ...we should make single-site routines uniform with this procedure...
    Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
    Hloc=0.d0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hloc(is,js)=Hloc_tmp(is,js) 
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)call print_hloc(Hloc)
    !
    !Get HF part of the self-energy
    Sigma_HF(1:Nlso,1:Nlso) = dreal(Sigma(1:Nlso,1:Nlso,Liw))
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    ed_Eloc_lattice = 0d0
    H0              = 0d0
    Hl              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       Hlktmp= 0d0      
       Hlk   = 0d0
       !
       do i=1+mpiID,Lmats,mpiSIZE
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) -       Sigma(:,:,i)  - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =                            -       SigmaA(:,:,i)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =                            -       SigmaA(:,:,i)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + conjg(Sigma(:,:,i)) + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          !extract the 11-block component:
          Gk(:,:) = Gknambu(1:Nlso,1:Nlso)
          Tk = zeye(Nlso)/(xi*wm(i)) - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck,Nlso)
          Hlktmp = Hlktmp + Wtk(ik)*trace_matrix(Dk,Nlso)
       enddo
#ifdef _MPI_INEQ
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       call MPI_ALLREDUCE(Hlktmp,Hlk,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
       H0k = H0ktmp
       Hlk = Hlktmp
#endif
       H0 = H0 + H0k
       Hl = Hl + Hlk
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    Hl = Hl/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nlso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nlso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nlso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nlso)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    Lail0=spin_degeneracy*Lail0
    Lail1=spin_degeneracy*Lail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Eloc_lattice=Hl+Lail0+Lail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
    ed_Eloc_lattice=ed_Eloc_lattice/dble(Nlat)
    Eout = [ed_Ekin_lattice,ed_Eloc_lattice]
    deallocate(wm)
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_lattice_superc_main

  function kinetic_energy_lattice_superc_1(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:)                             :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:,:)                             :: SigmaA ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,4)) :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,4)) :: SigmaA_! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                                   :: i,iorb,ilat,ispin,io,is
    integer                                                   :: j,jorb,jlat,jspin,jo,js
    integer                                                   :: Nlat,Nlso,Nso,Lk,Liw
    real(8)                                                   :: ed_Ekin_lattice(2)
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Nso  = size(Sigma,2)
    Lk   = size(Hk,3)
    Liw  = size(Sigma,4)
    Nlso = Nlat*Nso
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_superc_1","Hk")
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"kinetic_energy_lattice_superc_1","Sigma")
    call assert_shape(SigmaA,[Nlat,Nso,Nso,Liw],"kinetic_energy_lattice_superc_1","Self")
    Sigma_=zero
    SigmaA_=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb !spin-orbit stride
                   jo = jorb + (jspin-1)*Norb !spin-orbit stride
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_(is,js,:) = Sigma(ilat,io,jo,:)
                   SigmaA_(is,js,:)= SigmaA(ilat,io,jo,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_superc_main(Hk,Wtk,Sigma_,SigmaA_)
  end function kinetic_energy_lattice_superc_1

  function kinetic_energy_lattice_superc_2(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)                         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)                         :: SigmaA ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,6)) :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,6)) :: SigmaA_! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                                   :: i,iorb,ilat,ispin,io,is
    integer                                                   :: j,jorb,jlat,jspin,jo,js
    integer                                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                                   :: ed_Ekin_lattice(2)
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_superc_2","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Sigma")
    call assert_shape(SigmaA,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Self")
    Sigma_=zero
    SigmaA_=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_(is,js,:) = Sigma(ilat,ispin,jspin,iorb,jorb,:)
                   SigmaA_(is,js,:)= SigmaA(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_superc_main(Hk,Wtk,Sigma_,SigmaA_)
  end function kinetic_energy_lattice_superc_2

  function kinetic_energy_lattice_superc_1B(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                                 :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                               :: wtk    ! [Nk]
    complex(8),dimension(:,:)                                   :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Sigma,1),size(Sigma,2))           :: SigmaA ! [Nlat][L]  
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,2))   :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(size(Hk,1),size(Hk,1),size(Sigma,2))   :: SigmaA_! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                                     :: i,iorb,ilat,ispin,io,is
    integer                                                     :: j,jorb,jlat,jspin,jo,js
    integer                                                     :: Nlat,Nlso,Nso,Lk,Liw
    real(8)                                                     :: ed_Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Nlso = Nlat*1*1
    call assert_shape(Hk,[Nlso,Nlso,Lk],"kinetic_energy_lattice_superc_1B","Hk")
    Sigma_=zero
    SigmaA_=zero
    do ilat=1,Nlat
       Sigma_(ilat,ilat,:)  = Sigma(ilat,:)
       SigmaA_(ilat,ilat,:) = SigmaA(ilat,:)
    enddo
    ed_Ekin_lattice = kinetic_energy_lattice_superc_main(Hk,Wtk,Sigma_,SigmaA_)
  end function kinetic_energy_lattice_superc_1B









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8)                    :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>"
    close(unit)
  end subroutine write_energy_info

  subroutine write_kinetic_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="kinetic_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",reg(txtfy(1))//"<K>",reg(txtfy(2))//"<Eloc>"
    close(unit)
  end subroutine write_kinetic_info



  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy

  subroutine write_kinetic(Ekin)
    real(8) :: Ekin(2)
    integer :: unit
    unit = free_unit()
    open(unit,file="kinetic_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")Ekin(1),Ekin(2)
    close(unit)
  end subroutine write_kinetic


end MODULE ED_ENERGY
