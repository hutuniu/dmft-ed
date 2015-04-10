!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_ENERGY
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_TIMER
  USE SF_LINALG, only: matrix_inverse
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_MATVEC
  implicit none
  private


  interface ed_kinetic_energy
     module procedure kinetic_energy_impurity_normal
     module procedure kinetic_energy_impurity_normal_1B
     module procedure kinetic_energy_impurity_normal_MB
     module procedure kinetic_energy_impurity_superc
     module procedure kinetic_energy_impurity_superc_1B
     module procedure kinetic_energy_impurity_superc_MB
  end interface ed_kinetic_energy

  public  :: ed_kinetic_energy     !PUBLIC in DMFT
  public  :: local_energy_impurity !INTERNAL in DMFT

contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity()
    integer,dimension(Ntot)                           :: ib
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
          !
          !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
          !Euloc=\sum=i U_i*(n_u*n_d)_i
          !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
          do iorb=1,Norb
             ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*gs_weight
          enddo
          ! if(.not.ed_supercond) then
          !    ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
          ! else
          !    ed_Epot = ed_Epot - uloc(1)*(nup(1)-0.5d0)*(ndw(1)-0.5d0)*gs_weight
          ! end if
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
    !
    !
    call write_energy_info()
    call write_energy()
  end subroutine local_energy_impurity






  !-------------------------------------------------------------------------------------------
  !PURPOSE:  comment
  !-------------------------------------------------------------------------------------------
  subroutine kinetic_energy_impurity_normal_1B(Hk,Wtk,Sigma)
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(:)               :: Hk
    real(8),dimension(size(Hk))           :: Wtk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call kinetic_energy_impurity_normal(Hk_,Wtk,Sigma_)
  end subroutine kinetic_energy_impurity_normal_1B
  !
  subroutine kinetic_energy_impurity_normal_MB(Hk,Wtk,Sigma)
    complex(8),dimension(:,:,:)                               :: Hk
    real(8),dimension(size(Hk,3))                             :: Wtk
    complex(8),dimension(:,:,:,:,:)                           :: Sigma
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(Sigma,5)) :: Sigma_
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Sigma_(io,jo,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call kinetic_energy_impurity_normal(Hk,Wtk,Sigma_)
  end subroutine kinetic_energy_impurity_normal_MB
  !
  subroutine kinetic_energy_impurity_normal(Hk,Wtk,Sigma)
    integer                                  :: Lk,Nso,Liw
    integer                                  :: i,j,ik,iorb
    complex(8),dimension(:,:,:)              :: Hk
    complex(8),dimension(:,:,:)              :: Sigma
    real(8),dimension(size(Hk,3))            :: Wtk
    !
    real(8),dimension(:,:),allocatable       :: Sigma_HF
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk,Ck,Zk,Dk,Hloc
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
    real(8)                                  :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                  :: H0,Hl,ed_Ekin,ed_Eloc
    !
    Nso = size(Hk,1)
    Lk = size(Hk,3)
    Liw= size(Sigma,3)
    if(Nso/=size(Hk,2))stop "kinetic_energy_impurity_normal error: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(Nso/=size(Sigma,1).OR.Nso/=size(Sigma,2))stop "kinetic_energy_impurity_normal error: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    !
    allocate(wm(Liw))
    allocate(Sigma_HF(Nso,Nso))
    allocate(Ak(Nso,Nso),Bk(Nso,Nso),Ck(Nso,Nso),Dk(Nso,Nso),Hloc(Nso,Nso),Zk(Nso,Nso),Zeta(Nso,Nso),Gk(Nso,Nso),Tk(Nso,Nso))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)then
       call print_hloc(Hloc,LOGfile)
       ! write(LOGfile,"(A)")"Hloc"
       ! do i=1,Nso
       !    write(LOGfile,"(90F21.12)")(dreal(Hloc(i,j)),j=1,Nso)
       ! enddo
       ! write(LOGfile,*)""
       ! do i=1,Nso
       !    write(LOGfile,"(90F21.12)")(dimag(Hloc(i,j)),j=1,Nso)
       ! enddo
    endif

    !
    if(ED_MPI_ID==0)call start_timer()
    H0=0d0
    Hl=0d0
    Zk=0d0 ; forall(i=1:Nso)Zk(i,i)=1d0
    do ik=1,Lk
       Ak = Hk(:,:,ik) - Hloc(:,:)
       Bk =-Hk(:,:,ik) - Sigma_HF(:,:)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik) - Sigma(:,:,i)
          select case(Nso)
          case default
             call matrix_inverse(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
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
    ! Tail0=0d0
    ! Tail1=0d0
    ! do ik=1,Lk
    !    Ak= Hk(:,:,ik)
    !    Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
    !    Ck= matmul(Ak,Bk)
    !    Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nso)
    !    Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nso)
    ! enddo
    ! Tail0=spin_degeneracy*Tail0
    ! Tail1=spin_degeneracy*Tail1*beta
    ! ed_Ekin=H0+Tail0+Tail1
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
    deallocate(wm,Sigma_HF,Ak,Bk,Ck,Dk,Hloc,Zk,Zeta,Gk,Tk)
    call write_kinetic_info()
    call write_kinetic([ed_Ekin,ed_Eloc])
  end subroutine kinetic_energy_impurity_normal









  !-------------------------------------------------------------------------------------------
  !PURPOSE:  comment
  !-------------------------------------------------------------------------------------------
  subroutine kinetic_energy_impurity_superc_1B(Hk,Wtk,Sigma,Self)
    real(8),dimension(:)                     :: Hk
    real(8),dimension(size(Hk))              :: Wtk
    complex(8),dimension(:)                  :: Sigma,Self
    complex(8),dimension(1,1,size(Sigma))    :: Sigma_
    complex(8),dimension(1,1,size(Self))     :: Self_
    complex(8),dimension(1,1,size(Hk))       :: Hk_
    real(8),dimension(size(Hk))              :: Wtk_
    if(size(Sigma)/=size(Self)) stop "ed_kinetic_energy_sc: Normal and Anomalous self-energies have different size!"
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    call kinetic_energy_impurity_superc(Hk_,Wtk,Sigma_,Self_)
  end subroutine kinetic_energy_impurity_superc_1B

  subroutine kinetic_energy_impurity_superc_MB(Hk,Wtk,Sigma,Self)
    complex(8),dimension(:,:,:)                               :: Hk
    real(8),dimension(size(Hk,3))                             :: Wtk
    complex(8),dimension(:,:,:,:,:)                           :: Sigma
    complex(8),dimension(:,:,:,:,:)                           :: Self
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(Sigma,5)) :: Sigma_
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(Self,5))  :: Self_
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Sigma_(io,jo,:) = Sigma(ispin,jspin,iorb,jorb,:)
                Self_(io,jo,:)  = Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call kinetic_energy_impurity_superc(Hk,Wtk,Sigma_,Self_)
  end subroutine kinetic_energy_impurity_superc_MB

  subroutine kinetic_energy_impurity_superc(Hk,Wtk,Sigma,SigmaA)
    integer                                  :: Lk,Nso,Liw
    integer                                  :: i,j,ik,iorb,jorb,inambu,jnambu,n,m
    complex(8),dimension(:,:,:)              :: Hk
    complex(8),dimension(:,:,:)              :: Sigma,SigmaA
    real(8),dimension(size(Hk,3))            :: Wtk
    !
    real(8),dimension(:,:),allocatable       :: Sigma_HF,SigmaA_HF
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk,Ck,Dk,Zk,Hloc
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk,Gk_Nambu
    complex(8),dimension(2,2)                :: Gk_Nambu_ij
    !
    real(8)                                  :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                  :: H0,Hl,ed_Ekin,ed_Eloc
    !
    Nso = size(Hk,1)
    Lk = size(Hk,3)
    Liw= size(Sigma,3)
    if(Nso/=size(Hk,2))stop "kinetic_energy_impurity_superc error: size(Hk,1)!=size(Hk,2) [Nsorb_total]"
    if(Nso/=size(Sigma,1).OR.Nso/=size(Sigma,2))stop "kinetic_energy_impurity_superc error: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(Nso/=size(SigmaA,1).OR.Nso/=size(SigmaA,2))stop "kinetic_energy_impurity_superc error: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(Lk/=size(Wtk))stop "kinetic_energy_impurity_superc error: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Sigma_HF(Nso,Nso),SigmaA_HF(Nso,Nso))
    allocate(Ak(Nso,Nso),Bk(Nso,Nso),Ck(Nso,Nso),Dk(Nso,Nso),Zk(Nso,Nso),Hloc(Nso,Nso),Zeta(Nso,Nso),Gk(Nso,Nso),Tk(Nso,Nso))
    allocate(Gk_Nambu(2*Nso,2*Nso))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    if(ED_MPI_ID==0)then
       call print_hloc(Hloc,LOGfile)
       ! write(LOGfile,"(A)")"Hloc"
       ! do i=1,Nso
       !    write(LOGfile,"(90F21.12)")(dreal(Hloc(i,j)),j=1,Nso)
       ! enddo
       ! write(LOGfile,*)""
       ! do i=1,Nso
       !    write(LOGfile,"(90F21.12)")(dimag(Hloc(i,j)),j=1,Nso)
       ! enddo
    endif

    !
    H0=0d0
    Zk=0d0 ; forall(i=1:Nso)Zk(i,i)=1d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)-Hloc
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       do i=1,Liw
          Gk=zero          
          do iorb=1,Nso
             do jorb=1,Nso
                Gk_Nambu_ij=zero
                Gk_Nambu_ij(1,1) =  -Hk(iorb,jorb,ik)-Sigma(iorb,jorb,i)
                Gk_Nambu_ij(1,2) = -SigmaA(iorb,jorb,i)
                Gk_Nambu_ij(2,1) = -SigmaA(iorb,jorb,i)
                Gk_Nambu_ij(2,2) = -conjg(Gk_Nambu_ij(1,1))
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
          !
          call matrix_inverse(Gk_Nambu)
          !
          inambu=1
          jnambu=1
          do iorb=1,Nso
             do jorb=1,Nso
                m=(inambu-1)*Nso + iorb
                n=(jnambu-1)*Nso + jorb
                Gk(iorb,jorb) =  Gk_Nambu(m,n)
             enddo
          enddo
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
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
    ! Tail0=0d0
    ! Tail1=0d0
    ! do ik=1,Lk
    !    Ak= Hk(:,:,ik)
    !    Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
    !    Ck= matmul(Ak,Bk)
    !    do iorb=1,Nso
    !       Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak(iorb,iorb)!trace_matrix(Ak,Nso)
    !       Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck(iorb,iorb)!trace_matrix(Ck,Nso)
    !    end do
    ! enddo
    ! Tail0=spin_degeneracy*Tail0
    ! Tail1=spin_degeneracy*Tail1*beta
    ! ed_Ekin=H0+Tail0+Tail1
    ! deallocate(wm,Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
    ! call write_kinetic_info(ed_Ekin)
    ! call write_kinetic(ed_Ekin)
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
    deallocate(wm,Sigma_HF,Ak,Bk,Ck,Dk,Hloc,Zk,Zeta,Gk,Tk)
    call write_kinetic_info()
    call write_kinetic([ed_Ekin,ed_Eloc])
  end subroutine kinetic_energy_impurity_superc









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
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
    integer :: unit,iorb,jorb,ispin
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
    real(8) :: Ekin
    integer :: unit,iorb,jorb,ispin
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
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy

  subroutine write_kinetic(Ekin)
    real(8) :: Ekin(2)
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="kinetic_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")Ekin(1),Ekin(2)
    close(unit)
  end subroutine write_kinetic


end MODULE ED_ENERGY
