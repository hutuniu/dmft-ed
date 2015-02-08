module ED_WRAP_ENERGY
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_ARRAYS,    only:arange
  USE SF_IOTOOLS,   only:reg,sread,free_unit
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym,matrix_diagonalize,matrix_inverse_gj
  USE SF_TIMER
  !Impurity solver interface
  USE DMFT_ED
  !
  implicit none
  private


  interface kinetic_energy_lattice_normal
     module procedure &
          kinetic_energy_lattice_normal_case1,&
          kinetic_energy_lattice_normal_case2,&
          kinetic_energy_lattice_normal_case3
  end interface kinetic_energy_lattice_normal
  interface kinetic_energy_lattice_superc
     module procedure &
          kinetic_energy_lattice_superc_case1,&
          kinetic_energy_lattice_superc_case2,&
          kinetic_energy_lattice_superc_case3
  end interface kinetic_energy_lattice_superc
  public :: kinetic_energy_lattice_normal
  public :: kinetic_energy_lattice_superc

  !OBSOLETE:
  interface kinetic_energy_lattice
     module procedure kinetic_energy_lattice_normal,kinetic_energy_lattice_superc
  end interface kinetic_energy_lattice
  public :: kinetic_energy_lattice


  real(8),dimension(:),allocatable        :: wm

contains



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the general lattice case, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energy Sigma.
  !    
  ! we distinguish three different cases based on the shape of the Sigma array
  ! we assume Hk has always the following form:
  ! Hk: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  ! case I  : Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
  ! case II : Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  ! case III: Sigma: [Nlat][Nspin][Nspin][Norb][Norb][L]
  ! in no cases we will need Hloc, because this is included in Hk. If not the user must 
  ! take care of this condition
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_lattice_normal_case1(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                 :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))               :: wtk    ! [Nk]
    complex(8),dimension(:,:,:)                 :: Sigma  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    !aux
    complex(8),allocatable                      :: zeta_mats(:,:,:,:) ![Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                     :: Lk,Nlso,Liw
    integer                                     :: ik
    integer                                     :: i,iorb,ilat,ispin,io,is
    integer                                     :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Bk
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Gk
    real(8)                                     :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                     :: H0,H0k,H0ktmp
    real(8)                                     :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,3)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    allocate(zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Liw))
    zeta_mats=zero
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_normal error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_normal error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_normal error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlso) stop "kinetic_energy_lattice_normal error: size[Sigma,1] != Nlat*Nspin*Norb"
    if(size(Sigma,2)/=size(Sigma,1)) stop "kinetic_energy_lattice_normal error: size[Sigma,1] != size[Sigma,2] "
    !
    !Get HF part of the self-energy
    Sigma_HF(1:Nlso,1:Nlso) = dreal(Sigma(1:Nlso,1:Nlso,Liw))
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Sigma(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       do i=1+mpiID,Lmats,mpiSIZE
          Gk(:,:) = blocks_to_matrix(zeta_mats(:,:,:,i)) - Hk(:,:,ik)
          ! forall(j=1:Nlso)Gk(j,j) = xi*wm(i)+xmu
          ! Gk(:,:) = Gk(:,:) - Hk(:,:,ik) - Sigma(:,:,i)
          call matrix_inverse(Gk(:,:))
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_normal_case1

  function kinetic_energy_lattice_normal_case2(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                 :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))               :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:)               :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    !aux
    complex(8),allocatable                      :: zeta_mats(:,:,:,:) ![Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                     :: Lk,Nso,Nlso,Liw
    integer                                     :: ik
    integer                                     :: i,iorb,ilat,ispin,io,is
    integer                                     :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Gk
    real(8)                                     :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                     :: H0,H0k,H0ktmp
    real(8)                                     :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    Nso  = Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,4)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_normal error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_normal error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_normal error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlat) stop "kinetic_energy_lattice_normal error: size[Sigma,1] != Nlat"
    if(size(Sigma,2)/=Nso) stop "kinetic_energy_lattice_normal error: size[Sigma,2] != Nspin*Norb"
    if(size(Sigma,3)/=size(Sigma,2)) stop "kinetic_energy_lattice_normal error: size[Sigma,3] != size[Sigma,2] "
    !
    !Get HF part of the self-energy
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb !spin-orbit stride
                   jo = jorb + (jspin-1)*Norb !spin-orbit stride
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_HF(is,js) = dreal(Sigma(ilat,io,jo,Liw))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Sigma(ilat,io,jo,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       do i=1+mpiID,Lmats,mpiSIZE
          Gk(:,:) = blocks_to_matrix(zeta_mats(:,:,:,i)) - Hk(:,:,ik)
          ! forall(j=1:Nlso)Gk(j,j) = xi*wm(i)+xmu
          ! Gk(:,:) = Gk(:,:) - Hk(:,:,ik) !- Sigma(:,:,i)
          ! !<DEBUG  comment
          ! do ilat=1,Nlat
          !    do ispin=1,Nspin
          !       do jspin=1,Nspin
          !          do iorb=1,Norb
          !             do jorb=1,Norb
          !                io = iorb + (ispin-1)*Norb !spin-orbit stride
          !                jo = jorb + (jspin-1)*Norb !spin-orbit stride
          !                is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
          !                js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
          !                Gk(is,js) = Gk(is,js) - Sigma(ilat,io,jo,i)
          !             enddo
          !          enddo
          !       enddo
          !    enddo
          ! enddo
          ! !>DEBUG
          call matrix_inverse(Gk(:,:))
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2 !this can be merged into the prev do loops
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik) - Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_normal_case2

  function kinetic_energy_lattice_normal_case3(Hk,Wtk,Sigma) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                 :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))               :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)           :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    !aux
    complex(8),allocatable                      :: zeta_mats(:,:,:,:) ![Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                     :: Lk,Nso,Nlso,Liw
    integer                                     :: ik
    integer                                     :: i,iorb,ilat,ispin,io,is
    integer                                     :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Gk
    real(8)                                     :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                     :: H0,H0k,H0ktmp
    real(8)                                     :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    Nso  = Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,6)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_normal error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_normal error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_normal error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlat) stop "kinetic_energy_lattice_normal error: size[Sigma,1] != Nlat"
    if(size(Sigma,2)/=Nspin) stop "kinetic_energy_lattice_normal error: size[Sigma,2] != Nspin"
    if(size(Sigma,3)/=size(Sigma,2)) stop "kinetic_energy_lattice_normal error: size[Sigma,3] != size[Sigma,2] "
    if(size(Sigma,4)/=Norb) stop "kinetic_energy_lattice_normal error: size[Sigma,4] != Norb"
    if(size(Sigma,5)/=size(Sigma,4)) stop "kinetic_energy_lattice_normal error: size[Sigma,5] != size[Sigma,4] "
    !
    !Get HF part of the self-energy
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_HF(is,js) = dreal(Sigma(ilat,ispin,jspin,iorb,jorb,Liw))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Sigma(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       do i=1+mpiID,Lmats,mpiSIZE
          Gk(:,:) = blocks_to_matrix(zeta_mats(:,:,:,i)) - Hk(:,:,ik)
          ! forall(j=1:Nlso)Gk(j,j) = xi*wm(i)+xmu
          ! Gk(:,:) = Gk(:,:) - Hk(:,:,ik) !- Sigma(:,:,i)
          ! !<DEBUG  comment
          ! do ilat=1,Nlat
          !    do ispin=1,Nspin
          !       do jspin=1,Nspin
          !          do iorb=1,Norb
          !             do jorb=1,Norb
          !                is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
          !                js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
          !                Gk(is,js) = Gk(is,js) - Sigma(ilat,ispin,jspin,iorb,jorb,i)
          !             enddo
          !          enddo
          !       enddo
          !    enddo
          ! enddo
          ! !>DEBUG
          call matrix_inverse(Gk(:,:))
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik) - Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_normal_case3









  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the general lattice case, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energies Sigma (normal) and
  ! SigmaA (anomalous).
  !    
  ! we distinguish three different cases based on the shape of the Sigma array
  ! we assume Hk has always the following form:
  ! Hk: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  ! case I  : Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
  ! case II : Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  ! case III: Sigma: [Nlat][Nspin][Nspin][Norb][Norb][L]
  ! in no cases we will need Hloc, because this is included in Hk. If not the user must 
  ! take care of this condition
  !-------------------------------------------------------------------------------------------
  function kinetic_energy_lattice_superc_case1(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                     :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                   :: wtk    ! [Nk]
    complex(8),dimension(:,:,:)                     :: Sigma  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:)                     :: SigmaA ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]  
    !aux
    complex(8),allocatable                          :: zeta_mats(:,:,:,:,:,:) ![2][2][Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                         :: Lk,Nso,Nlso,Liw
    integer                                         :: ik
    integer                                         :: i,iorb,ilat,ispin,io,is
    integer                                         :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Gk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2)) :: Gknambu
    real(8)                                         :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                         :: H0,H0k,H0ktmp
    real(8)                                         :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    Nso  = Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,3)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    allocate(zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Liw))
    zeta_mats=zero
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_superc error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_superc error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_superc error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlso) stop "kinetic_energy_lattice_superc error: size[Sigma,1] != Nlat"
    if(size(Sigma,2)/=size(Sigma,1)) stop "kinetic_energy_lattice_superc error: size[Sigma,2] != size[Sigma,1] "
    if(size(SigmaA,1)/=Nlso) stop "kinetic_energy_lattice_superc error: size[SigmaA,1] != Nlat"
    if(size(SigmaA,2)/=size(SigmaA,1)) stop "kinetic_energy_lattice_superc error: size[SigmaA,2] != size[SigmaA,1] "
    !
    !Get HF part of the self-energy
    Sigma_HF(1:Nlso,1:Nlso) = dreal(Sigma(1:Nlso,1:Nlso,Liw))
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu !+ Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Sigma(is,js,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - SigmaA(is,js,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - SigmaA(is,js,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg(Sigma(is,js,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       !
       do i=1+mpiID,Lmats,mpiSIZE
          Gknambu(1:Nlso,1:Nlso)               = blocks_to_matrix(zeta_mats(1,1,:,:,:,i)) - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        = blocks_to_matrix(zeta_mats(1,2,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        = blocks_to_matrix(zeta_mats(2,1,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_mats(2,2,:,:,:,i)) + Hk(:,:,ik)
          call matrix_inverse(Gknambu(:,:))
          !extract the 11-block component:
          Gk(:,:) = Gknambu(1:Nlso,1:Nlso)
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2 !this can be merged into the prev do loops
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik) - Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_superc_case1

  function kinetic_energy_lattice_superc_case2(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                     :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                   :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:)                   :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:,:)                   :: SigmaA ! [Nlat][Nspin*Norb][Nspin*Norb][L]  
    !aux
    complex(8),allocatable                          :: zeta_mats(:,:,:,:,:,:)
    integer                                         :: Lk,Nso,Nlso,Liw
    integer                                         :: ik
    integer                                         :: i,iorb,ilat,ispin,io,is
    integer                                         :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Gk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2)) :: Gknambu
    real(8)                                         :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                         :: H0,H0k,H0ktmp
    real(8)                                         :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    Nso  = Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,4)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    allocate(zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Liw))
    zeta_mats=zero
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_superc error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_superc error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_superc error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlat) stop "kinetic_energy_lattice_superc error: size[Sigma,1] != Nlat"
    if(size(Sigma,2)/=Nso) stop "kinetic_energy_lattice_superc error: size[Sigma,2] != Nspin*Norb"
    if(size(Sigma,3)/=size(Sigma,2)) stop "kinetic_energy_lattice_superc error: size[Sigma,3] != size[Sigma,2] "
    if(size(SigmaA,1)/=Nlat) stop "kinetic_energy_lattice_superc error: size[SigmaA,1] != Nlat"
    if(size(SigmaA,2)/=Nso) stop "kinetic_energy_lattice_superc error: size[SigmaA,2] != Nspin*Norb"
    if(size(SigmaA,3)/=size(SigmaA,2)) stop "kinetic_energy_lattice_superc error: size[SigmaA,3] != size[SigmaA,2] "
    !
    !Get HF part of the self-energy
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb !spin-orbit stride
                   jo = jorb + (jspin-1)*Norb !spin-orbit stride
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_HF(is,js) = dreal(Sigma(ilat,io,jo,Liw))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu !+ Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Sigma(ilat,io,jo,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - SigmaA(ilat,io,jo,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - SigmaA(ilat,io,jo,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg(Sigma(ilat,io,jo,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       !
       do i=1+mpiID,Lmats,mpiSIZE
          Gknambu(1:Nlso,1:Nlso)               = blocks_to_matrix(zeta_mats(1,1,:,:,:,i)) - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        = blocks_to_matrix(zeta_mats(1,2,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        = blocks_to_matrix(zeta_mats(2,1,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_mats(2,2,:,:,:,i)) + Hk(:,:,ik)
          call matrix_inverse(Gknambu(:,:))
          !extract the 11-block component:
          Gk(:,:) = Gknambu(1:Nlso,1:Nlso)
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2 !this can be merged into the prev do loops
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik) - Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_superc_case2

  function kinetic_energy_lattice_superc_case3(Hk,Wtk,Sigma,SigmaA) result(ed_Ekin_lattice)
    complex(8),dimension(:,:,:)                     :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))                   :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)               :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)               :: SigmaA ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    !aux
    complex(8),allocatable                          :: zeta_mats(:,:,:,:,:,:)
    integer                                         :: Lk,Nso,Nlso,Liw
    integer                                         :: ik
    integer                                         :: i,iorb,ilat,ispin,io,is
    integer                                         :: j,jorb,jlat,jspin,jo,js
    real(8),dimension(size(Hk,1),size(Hk,2))        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Ck
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Gk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2)) :: Gknambu
    real(8)                                         :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                         :: H0,H0k,H0ktmp
    real(8)                                         :: ed_Ekin_lattice
    !Get generalized Lattice-Spin-Orbital index
    Nlso = Nlat*Nspin*Norb
    Nso  = Nspin*Norb
    !
    !Get number of k-points:
    Lk=size(Hk,3)
    !
    !Get number of matsubara freq.
    Liw=size(Sigma,6)
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    allocate(zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Liw))
    zeta_mats=zero
    !
    ! check Hk and Sigma dimensions 
    if(size(Hk,1)/=Nlso) stop "kinetic_energy_lattice_superc error: size[Hk,1] != Nlat*Nspin*Norb"
    if(size(Hk,2)/=size(Hk,1)) stop "kinetic_energy_lattice_superc error: size[Hk,1] != size[Hk,2]"
    if(size(Wtk)/=Lk) stop "kinetic_energy_lattice_superc error: size[Wtk] != size[Hk,3] = Lk"
    if(size(Sigma,1)/=Nlat) stop "kinetic_energy_lattice_superc error: size[Sigma,1] != Nlat"
    if(size(Sigma,2)/=Nspin) stop "kinetic_energy_lattice_superc error: size[Sigma,2] != Nspin"
    if(size(Sigma,3)/=size(Sigma,2)) stop "kinetic_energy_lattice_superc error: size[Sigma,3] != size[Sigma,2] "
    if(size(Sigma,4)/=Norb) stop "kinetic_energy_lattice_superc error: size[Sigma,4] != Norb"
    if(size(Sigma,5)/=size(Sigma,4)) stop "kinetic_energy_lattice_superc error: size[Sigma,5] != size[Sigma,4] "
    if(size(SigmaA,1)/=Nlat) stop "kinetic_energy_lattice_superc error: size[SigmaA,1] != Nlat"
    if(size(SigmaA,2)/=Nspin) stop "kinetic_energy_lattice_superc error: size[SigmaA,2] != Nspin"
    if(size(SigmaA,3)/=size(SigmaA,2)) stop "kinetic_energy_lattice_superc error: size[SigmaA,3] != size[SigmaA,2] "
    if(size(SigmaA,4)/=Norb) stop "kinetic_energy_lattice_superc error: size[SigmaA,4] != Norb"
    if(size(SigmaA,5)/=size(SigmaA,4)) stop "kinetic_energy_lattice_superc error: size[SigmaA,5] != size[SigmaA,4] "
    !
    !Get HF part of the self-energy
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Sigma_HF(is,js) = dreal(Sigma(ilat,ispin,jspin,iorb,jorb,Liw))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Get the k-INDEPENDENT part of Gk: zeta= iw+mu-Sigma(iw)
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu !+ Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Sigma(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - SigmaA(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - SigmaA(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg(Sigma(ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Start the timer:
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    ed_Ekin_lattice = 0d0
    H0              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk       
       Gk    = zero
       H0ktmp= 0d0      
       H0k   = 0d0
       !
       do i=1+mpiID,Lmats,mpiSIZE
          Gknambu(1:Nlso,1:Nlso)               = blocks_to_matrix(zeta_mats(1,1,:,:,:,i)) - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        = blocks_to_matrix(zeta_mats(1,2,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        = blocks_to_matrix(zeta_mats(2,1,:,:,:,i))
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_mats(2,2,:,:,:,i)) + Hk(:,:,ik)
          call matrix_inverse(Gknambu(:,:))
          !extract the 11-block component:
          Gk(:,:) = Gknambu(1:Nlso,1:Nlso)
          Tk = zero
          forall(j=1:Nlso)Tk(j,j)=one/(xi*wm(i))
          Tk = Tk - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2 !this can be merged into the prev do loops
          Ck = matmul(Hk(:,:,ik),Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
       enddo
       call MPI_ALLREDUCE(H0ktmp,H0k,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       H0 = H0 + H0k
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ck= matmul(Hk(:,:,ik),(-Hk(:,:,ik) - Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Hk(:,:,ik))
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=H0+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  end function kinetic_energy_lattice_superc_case3






























  !###################################################################################################
  !###################################################################################################
  !               POSSIBLY OBSOLETE ROUTINES NOW SUPERSEDED BY TOP ROUTINE HERE
  !               GETTING KINETIC ENERGY FOR ANY Norb*Nspin*Nlat*Nk NORMAL PHASE
  !                (these routines are left temporarily for back-compatibility)
  !                                     ( to be removed)
  !###################################################################################################
  !###################################################################################################

  function kinetic_energy_lattice_normal(Hk,wtk,Hloc,Sigma) result(ed_Ekin_lattice)
    !+- ONLY SPIN DIAGONAL QUANTITIES -+!
    complex(8),dimension(:,:,:)              :: Hk     ! [Norb*Nlat][Norb*Nlat][Nk]
    real(8),dimension(:)                     :: wtk    ! [Lk]
    complex(8),dimension(:,:,:)              :: Hloc   ! [Nlat][Norb][Norb]    
    complex(8),dimension(:,:,:,:)            :: Sigma  ! [Nlat][Norb][Norb][Lmats]
    integer                                  :: Lk,No
    integer                                  :: i,ik,iorb,ilat,io,jorb,jo
    real(8),dimension(:,:,:),allocatable     :: Sigma_HF
    real(8),dimension(:),allocatable         :: elocal_
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Tk
    complex(8),dimension(:,:),allocatable    :: Gk
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: ed_Ekin_lattice,sumMats
    real(8)                                  :: sumMatstmp,sumMatsk
    !
    ! checks !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    No = Norb*Nlat  ! generalized "orbital index"
    write(*,*) No,Norb,Nlat
    Lk=size(Hk,3)
    ! check Hk !
    if(size(Hk,1)/=No) stop "get_kinetic_energy: size Hk(:,*,*) /= No"
    if(size(Hk,2)/=size(Hk,1)) stop "get_kinetic_energy: size Hk(*,:,*) /= Hk(:,*,*) "
    ! check wt !
    if(size(wtk)/=Lk) stop "get_kinetic_energy: size wt(*) /= Lk"
    ! check Hloc !
    if(size(Hloc,1)/=Nlat) stop "get_kinetic_energy: size Hloc(:,*,*) /= Nlat"
    if(size(Hloc,2)/=size(Hloc,3)) stop "get_kinetic_energy: size Hloc(*,:,*) /= Hloc(*,*,:) "
    if(size(Hloc,2)/=Norb) stop "get_kinetic_energy: size Hloc(*,:,*) /= Norb"
    ! check Smats !
    if(size(Sigma,1)/=Nlat) stop "get_kinetic_energy: size Smats(:,*,*,*) /= Nlat"
    if(size(Sigma,2)/=size(Sigma,3)) stop "get_kinetic_energy: size Smats(*,:,*,*) /= Smats(*,*,:,*) "
    if(size(Sigma,2)/=Norb) stop "get_kinetic_energy: size Smats(*,:,*,*) /= Norb"
    if(size(Sigma,4)/=Lmats) stop "get_kinetic_energy: size Smats(*,*,*,:) /= Lmats"
    !
    allocate(Sigma_HF(Nlat,Norb,Norb))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Tk(No,No))
    allocate(Gk(No,No))
    !
    Sigma_HF = dreal(Sigma(:,:,:,Lmats))
    !
    ed_Ekin_lattice=0.d0
    sumMats=0.d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    do ik=1,Lk       
       Ak =  Hk(:,:,ik) !- xmu*Zk
       Bk = -Hk(:,:,ik) !+ xmu*Zk      
       do ilat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io = (ilat-1)*Norb + iorb
                jo = (ilat-1)*Norb + jorb
                Bk(io,jo) = Bk(io,jo) - Sigma_HF(ilat,iorb,jorb) 
                Bk(io,jo) = Bk(io,jo) - Hloc(ilat,iorb,jorb)
             end do
          end do
       end do
       !Gk(iw) computation
       Gk=zero
       sumMatstmp=0.d0      
       sumMatsk=0.d0      
       do i=1+mpiID,Lmats,mpiSIZE
          Gk(:,:) = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   io = (ilat-1)*Norb + iorb
                   jo = (ilat-1)*Norb + jorb
                   Gk(io,jo) = Gk(io,jo) - Sigma(ilat,iorb,jorb,i)
                   Gk(io,jo) = Gk(io,jo) - Hloc(ilat,iorb,jorb)
                end do
             end do
          end do
          !
          call matrix_inverse(Gk(:,:))
          !
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          do io=1,No
             sumMatstmp = sumMatstmp + Wtk(ik)*Ck(io,io)
          end do
       enddo
       call MPI_ALLREDUCE(sumMatstmp,sumMatsk,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       sumMats=sumMats+sumMatsk
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    sumMats=sumMats/beta*2.d0*spin_degeneracy
    !
    ! subtract tails
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik) !- xmu*Zk
       Bk=-Hk(:,:,ik) !+ xmu*Zk      
       do ilat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io = (ilat-1)*Norb + iorb
                jo = (ilat-1)*Norb + jorb
                ! Bk(io,jo) = Bk(io,jo) - Sigma_HF(iorb,jorb,ilat)
                ! if(iorb==jorb) Bk(io,jo) = Bk(io,jo) - elocal_(ilat) ![tmp...think to something more clever]
                Bk(io,jo) = Bk(io,jo) - Sigma_HF(ilat,iorb,jorb)
                Bk(io,jo) = Bk(io,jo) - Hloc(ilat,iorb,jorb)
             end do
          end do
       end do
       Ck= matmul(Ak,Bk)
       do io=1,No
          Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak(io,io)
          Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck(io,io)
       end do
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ed_Ekin_lattice=sumMats+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
    deallocate(Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  end function kinetic_energy_lattice_normal





  function kinetic_energy_lattice_superc(Hk,Wtk,Sigma,SigmaA)  result(ed_Ekin_lattice)
    integer                                  :: Lk,No,Liw,Norb
    integer                                  :: i,ik,iorb,jorb,ilat,io,jo
    integer                                  :: inambu,jnambu,n,m
    complex(8),dimension(:,:,:)              :: Hk              ! (Norb*Nlat,Norb*Nlat,Lk)
    complex(8),dimension(:,:,:,:)            :: Sigma,SigmaA    ! (Norb,Norb,Nlat,Lw)
    real(8),dimension(:)                     :: Wtk             ! (Lk)
    !
    real(8),dimension(:,:,:),allocatable     :: Sigma_HF
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Tk,Gk_Nambu
    complex(8),dimension(2,2)                :: Gk_Nambu_ij
    complex(8),dimension(:,:),allocatable    :: Gk
    !
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: ed_Ekin_lattice,sumMats
    real(8)                                  :: sumMatstmp,sumMatsk
    complex(8),dimension(:,:),allocatable    :: tmpGloc
    integer                                  :: is
    if(size(Sigma,1) /= size(SigmaA,1))stop "get_kinetic_energy: size1 Sigma / size sigmaA"
    if(size(Sigma,2) /= size(SigmaA,2))stop "get_kinetic_energy: size2 Sigma / size sigmaA"
    if(size(Sigma,3) /= size(SigmaA,3))stop "get_kinetic_energy: size3 Sigma / size sigmaA"
    if(size(Sigma,4) /= size(SigmaA,4))stop "get_kinetic_energy: size4 Sigma / size sigmaA"
    !
    Norb = size(Sigma,1)
    Nlat = size(Sigma,3)
    Liw  = size(Sigma,4)
    Lk = size(Hk,3)
    No = Norb*Nlat  ! generalized "orbital index"    
    if(No/=size(Hk,1).or.No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb*Nlat]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"    !
    !
    allocate(Sigma_HF(Norb,Norb,Nlat))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Tk(No,No))
    allocate(Gk_Nambu(2*No,2*No))
    allocate(Gk(No,No))
    !
    Sigma_HF = dreal(Sigma(:,:,:,Liw))
    !
    sumMats=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
    if(mpiID==0)call start_timer
    do ik=1,Lk     
       Ak =  Hk(:,:,ik) !- xmu*Zk
       Bk = -Hk(:,:,ik) !+ xmu*Zk      
       do ilat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io = (ilat-1)*Norb + iorb
                jo = (ilat-1)*Norb + jorb
                Bk(io,jo) = Bk(io,jo) - Sigma_HF(iorb,jorb,ilat)
             end do
          end do
       end do
       Gk=zero          
       sumMatstmp=0.d0      
       sumMatsk=0.d0      
       ! Compute Gk(iw)
       do i=1+mpiID,Liw,mpiSIZE
          Gk_Nambu=zero
          ! Build not interacting part of the Gk_Nambu matrix
          do io=1,No   
             do jo=1,No
                Gk_Nambu_ij=zero                
                Gk_Nambu_ij(1,1) =  - Hk(io,jo,ik) 
                Gk_Nambu_ij(1,2) =  zero
                Gk_Nambu_ij(2,1) =  zero
                Gk_Nambu_ij(2,2) =  + Hk(io,jo,ik) ! CONJG(H(-k))=H(K)  
                if(io==jo) then
                   Gk_Nambu_ij(1,1) = Gk_Nambu_ij(1,1) + xi*wm(i)+xmu
                   Gk_Nambu_ij(2,2) = Gk_Nambu_ij(2,2) + xi*wm(i)-xmu
                end if
                ! Fill the 2*Norb*Nlat x 2*Norb*Nlat  Gk_nambu matrix
                do inambu=1,2
                   do jnambu=1,2
                      m=get_index(inambu,io)
                      n=get_index(jnambu,jo)
                      Gk_nambu(m,n)=Gk_nambu_ij(inambu,jnambu)
                   enddo
                enddo
             enddo
          enddo
          ! Add the Self-energy contribution
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   io = (ilat-1)*Norb + iorb
                   jo = (ilat-1)*Norb + jorb
                   Gk_Nambu_ij=zero
                   Gk_Nambu_ij(1,1) = -Sigma(iorb,jorb,ilat,i)
                   Gk_Nambu_ij(1,2) = -SigmaA(iorb,jorb,ilat,i)
                   Gk_Nambu_ij(2,1) = -SigmaA(iorb,jorb,ilat,i)
                   Gk_Nambu_ij(2,2) =  conjg(Sigma(iorb,jorb,ilat,i))
                   ! Add contribution to the 2*Norb*Nlat Gk_nambu matrix
                   do inambu=1,2
                      do jnambu=1,2
                         m=get_index(inambu,io)
                         n=get_index(jnambu,jo)
                         Gk_nambu(m,n)=Gk_nambu(m,n)+Gk_nambu_ij(inambu,jnambu)
                      enddo
                   enddo
                end do
             end do
          end do
          ! Invert and select the (1,1) Nambu submatrix
          call matrix_inverse(Gk_Nambu)
          inambu=1
          jnambu=1
          do io=1,No
             do jo=1,No                
                m=get_index(inambu,io)
                n=get_index(jnambu,jo)
                Gk(io,jo) =  Gk_Nambu(m,n)
             enddo
          enddo
          !
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak,Gk - Tk)
          !
          do io=1,No
             sumMatstmp = sumMatstmp + Wtk(ik)*Ck(io,io)
          enddo
          !
       enddo
       call MPI_ALLREDUCE(sumMatstmp,sumMatsk,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
       sumMats = sumMats + sumMatsk
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    enddo
    if(mpiID==0)call stop_timer
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    sumMats=sumMats/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak =  Hk(:,:,ik) !- xmu*Zk
       Bk = -Hk(:,:,ik) !+ xmu*Zk        
       do ilat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io = (ilat-1)*Norb + iorb
                jo = (ilat-1)*Norb + jorb
                Bk(io,jo) = Bk(io,jo) - Sigma_HF(iorb,jorb,ilat)
             end do
          end do
       end do
       Ck= matmul(Ak,Bk)
       do io=1,No
          Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak(io,io)
          Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck(io,io)
       end do
    enddo
    Tail0 = spin_degeneracy*Tail0
    Tail1 = spin_degeneracy*Tail1*beta
    ed_Ekin_lattice = sumMats+Tail0+Tail1
    ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
    deallocate(Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  contains
    function get_index(inambu,io) result(m)
      integer :: inambu,io,m
      m=(inambu-1)*No + io
      ! This is the formulation for 2x2 matrix made of No x No matrices !
    end function get_index
  end function kinetic_energy_lattice_superc





  !-------------------------------------------------------------------------------------------
  !PURPOSE: evaluate the trace of a square matrix
  !-------------------------------------------------------------------------------------------
  function trace_matrix(M) result(tr)
    complex(8),dimension(:,:) :: M
    integer                   :: dim
    complex(8)                :: tr
    integer                   :: i
    tr=dcmplx(0d0,0d0)
    dim=size(M,1)
    if(size(M,2)/=dim)stop "trace_matrix error: wrong Matrix dimensions"
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix




end module ED_WRAP_ENERGY
