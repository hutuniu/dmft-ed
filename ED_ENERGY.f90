!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_ENERGY
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG, only: matrix_inverse
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_MATVEC
  implicit none
  private
  !
  interface kinetic_energy_impurity
     module procedure kinetic_energy_impurity_normal,kinetic_energy_impurity_superc
  end interface kinetic_energy_impurity
  !
  public  :: local_energy_impurity
  public  :: kinetic_energy_impurity
  !
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
          ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
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
             ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
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





  subroutine kinetic_energy_impurity_normal(Hk,Wtk,Sigma)
    integer                                  :: Lk,No,Liw
    integer                                  :: i,ik,iorb
    complex(8),dimension(:,:,:)              :: Hk
    complex(8),dimension(:,:,:)              :: Sigma
    real(8),dimension(:)                     :: Wtk
    !
    real(8),dimension(:,:),allocatable       :: Sigma_HF
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: H0
    !
    No = size(Hk,1)
    Lk = size(Hk,3)
    Liw= size(Sigma,3)
    if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(No/=size(Sigma,1).OR.No/=size(Sigma,2))stop "get_kinetic_energy: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Sigma_HF(No,No))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    H0=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    do ik=1,Lk
       ! Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik) - Sigma(:,:,i)
          select case(No)
          case default
             call matrix_inverse(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Hk(:,:,ik),Gk - Tk)!matmul(Ak,Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       Ck= matmul(Ak,Bk)
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin=H0+Tail0+Tail1
    deallocate(wm,Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)

    call write_energy_info(ed_Ekin)
    call write_energy(ed_Ekin)

  contains

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

  end subroutine kinetic_energy_impurity_normal







  subroutine kinetic_energy_impurity_superc(Hk,Wtk,Sigma,SigmaA)
    integer                                  :: Lk,No,Liw
    integer                                  :: i,ik,iorb,jorb,inambu,jnambu,n,m
    complex(8),dimension(:,:,:)              :: Hk
    complex(8),dimension(:,:,:)              :: Sigma,SigmaA
    real(8),dimension(:)                     :: Wtk
    !
    real(8),dimension(:,:),allocatable       :: Sigma_HF,SigmaA_HF
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk,Gk_Nambu
    complex(8),dimension(2,2)                :: Gk_Nambu_ij
    !
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: H0
    !
    No = size(Hk,1)
    Lk = size(Hk,3)
    Liw= size(Sigma,3)
    if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(No/=size(Sigma,1).OR.No/=size(Sigma,2))stop "get_kinetic_energy: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(No/=size(SigmaA,1).OR.No/=size(SigmaA,2))stop "get_kinetic_energy: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Sigma_HF(No,No),SigmaA_HF(No,No))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
    allocate(Gk_Nambu(2*No,2*No))
    !
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    H0=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       do i=1,Liw
          Gk=zero          
          do iorb=1,No
             do jorb=1,No
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
                      m=(inambu-1)*No + iorb
                      n=(jnambu-1)*No + jorb
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
          do iorb=1,No
             do jorb=1,No
                m=(inambu-1)*No + iorb
                n=(jnambu-1)*No + jorb
                Gk(iorb,jorb) =  Gk_Nambu(m,n)
             enddo
          enddo
          !
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak,Gk - Tk)
          !
          do iorb=1,Norb
             H0 = H0 + Wtk(ik)*Ck(iorb,iorb)
          enddo
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       Ck= matmul(Ak,Bk)
       do iorb=1,No
          Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak(iorb,iorb)!trace_matrix(Ak,No)
          Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck(iorb,iorb)!trace_matrix(Ck,No)
       end do
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin=H0+Tail0+Tail1
    deallocate(wm,Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)

    call write_energy_info(ed_Ekin)
    call write_energy(ed_Ekin)


  end subroutine kinetic_energy_impurity_superc



  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_energy_info(Ekin)
    real(8),optional :: Ekin
    integer          :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="energy_info.ed")
    if(present(Ekin))then
       write(unit,"(A1,90(A14,1X))")"#",&
            reg(txtfy(1))//"<K>",&
            reg(txtfy(2))//"<Hi>",&
            reg(txtfy(3))//"<H>",&
            reg(txtfy(4))//"<V>=<Hi-Ehf>",&
            reg(txtfy(5))//"<E0>",&
            reg(txtfy(6))//"<Ehf>",&
            reg(txtfy(7))//"<Dst>",&
            reg(txtfy(8))//"<Dnd>",&
            reg(txtfy(9))//"<Dse>",&
            reg(txtfy(10))//"<Dph>"
    else
       write(unit,"(A1,90(A14,1X))")"#",&
            reg(txtfy(1))//"<Hi>",&
            reg(txtfy(2))//"<V>=<Hi-Ehf>",&
            reg(txtfy(3))//"<E0>",&
            reg(txtfy(4))//"<Ehf>",&
            reg(txtfy(5))//"<Dst>",&
            reg(txtfy(6))//"<Dnd>",&
            reg(txtfy(7))//"<Dse>",&
            reg(txtfy(8))//"<Dph>"
    endif
    close(unit)
  end subroutine write_energy_info



  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_energy(Ekin)
    real(8),optional :: Ekin
    integer          :: unit
    integer          :: iorb,jorb,ispin
    if(present(Ekin))then
       unit = free_unit()
       open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")ed_Ekin,ed_Epot,ed_Ekin+ed_Epot+ed_Eknot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
       close(unit)
    else
       unit = free_unit()
       open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
       close(unit)
    endif
  end subroutine write_energy

end MODULE ED_ENERGY
