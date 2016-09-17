module ED_WRAP_ENERGY
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_ARRAYS,    only:arange
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym,zeye
  USE SF_IOTOOLS,   only:free_unit,reg,txtfy
  USE SF_TIMER
  !
  implicit none
  private


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
  public :: ed_kinetic_energy_lattice

  real(8),dimension(:),allocatable        :: wm



contains



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy for the general lattice case, given 
  ! the Hamiltonian matrix Hk and the DMFT self-energy Sigma.
  ! The main routine accept:
  ! - Sigma: [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
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
    real(8),dimension(size(Hk,1),size(Hk,2))                        :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Gk
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
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! Get the block diagonal part of the local Hamiltonian 
    ! This term gives rise to an additional E_loc = Hloc*n
    Hloc = zero
    Hloc(:,:) = sum(Hk(:,:,:),dim=3)/Lk
    !
    !Get HF part of the self-energy
    Sigma_HF(:,:) = dreal(Sigma(:,:,Liw))
    !
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
          call matrix_inverse(Gk(:,:))
          Tk = zeye(Nlso)/(xi*wm(i)) - Bk/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
          Hlktmp = Hlktmp + Wtk(ik)*trace_matrix(Dk)
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
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:))
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk)
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
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_lattice_normal_main







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
    complex(8),dimension(size(Hk,1),size(Hk,2))     :: Ak,Bk,Ck,Dk,Hloc
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
    if(allocated(wm))deallocate(wm)
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    ! Get the block diagonal part of the local Hamiltonian 
    ! This term gives rise to an additional E_loc = Hloc*n
    ! This latter was already evaluated in the previous version of the code
    ! by means of Eii = Epot + Eknot in the WRAP_ED module.
    ! Thus we had to remove the corresponding part of the Hamiltonian.
    ! Get the block diagonal part of the local Hamiltonian 
    ! This term gives rise to an additional E_loc = Hloc*n
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
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
          call matrix_inverse(Gknambu(:,:))
          !extract the 11-block component:
          Gk(:,:) = Gknambu(1:Nlso,1:Nlso)
          Tk = zeye(Nlso)/(xi*wm(i)) - (-Hk(:,:,ik) - Sigma_HF(:,:))/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          H0ktmp = H0ktmp + Wtk(ik)*trace_matrix(Ck)
          Hlktmp = Hlktmp + Wtk(ik)*trace_matrix(Dk)
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
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:))
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk)
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
    call write_kinetic_info()
    call write_kinetic(Eout)
  end function kinetic_energy_lattice_superc_main















  !+-----------------------------------------------------------------------------+!
  !PURPOSE: additional interfaces for different allocation of the required 
  ! functions.
  ! We distinguish different interface according to other shapes of the self-energy:
  ! - Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [Nlat][L]
  ! - Sigma: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !+-----------------------------------------------------------------------------+!

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
  !
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
  !
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



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
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
  subroutine write_kinetic(Ekin)
    real(8) :: Ekin(2)
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="kinetic_last.ed")
    write(unit,"(90F15.9)")Ekin(1),Ekin(2)
    close(unit)
  end subroutine write_kinetic










  ! !###################################################################################################
  ! !###################################################################################################
  ! !               POSSIBLY OBSOLETE ROUTINES NOW SUPERSEDED BY TOP ROUTINE HERE
  ! !               GETTING KINETIC ENERGY FOR ANY Norb*Nspin*Nlat*Nk NORMAL PHASE
  ! !                (these routines are left temporarily for back-compatibility)
  ! !                                     ( to be removed)
  ! !###################################################################################################
  ! !###################################################################################################
  ! function kinetic_energy_lattice_OLD(Hk,wtk,Hloc,Sigma) result(ed_Ekin_lattice)
  !   !+- ONLY SPIN DIAGONAL QUANTITIES -+!
  !   complex(8),dimension(:,:,:)              :: Hk     ! [Norb*Nlat][Norb*Nlat][Nk]
  !   real(8),dimension(:)                     :: wtk    ! [Lk]
  !   complex(8),dimension(:,:,:)              :: Hloc   ! [Nlat][Norb][Norb]    
  !   complex(8),dimension(:,:,:,:)            :: Sigma  ! [Nlat][Norb][Norb][Lmats]
  !   integer                                  :: Lk,No
  !   integer                                  :: i,ik,iorb,ilat,io,jorb,jo
  !   real(8),dimension(:,:,:),allocatable     :: Sigma_HF
  !   real(8),dimension(:),allocatable         :: elocal_
  !   complex(8),dimension(:,:),allocatable    :: Ak,Bk
  !   complex(8),dimension(:,:),allocatable    :: Ck,Zk
  !   complex(8),dimension(:,:),allocatable    :: Zeta,Tk
  !   complex(8),dimension(:,:),allocatable    :: Gk
  !   real(8)                                  :: Tail0,Tail1,spin_degeneracy
  !   !
  !   real(8)                                  :: ed_Ekin_lattice,sumMats
  !   real(8)                                  :: sumMatstmp,sumMatsk
  !   !
  !   ! checks !
  !   if(allocated(wm))deallocate(wm)
  !   allocate(wm(Lmats))
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   No = Norb*Nlat  ! generalized "orbital index"
  !   write(*,*) No,Norb,Nlat
  !   Lk=size(Hk,3)
  !   ! check Hk !
  !   if(size(Hk,1)/=No) stop "get_kinetic_energy: size Hk(:,*,*) /= No"
  !   if(size(Hk,2)/=size(Hk,1)) stop "get_kinetic_energy: size Hk(*,:,*) /= Hk(:,*,*) "
  !   ! check wt !
  !   if(size(wtk)/=Lk) stop "get_kinetic_energy: size wt(*) /= Lk"
  !   ! check Hloc !
  !   if(size(Hloc,1)/=Nlat) stop "get_kinetic_energy: size Hloc(:,*,*) /= Nlat"
  !   if(size(Hloc,2)/=size(Hloc,3)) stop "get_kinetic_energy: size Hloc(*,:,*) /= Hloc(*,*,:) "
  !   if(size(Hloc,2)/=Norb) stop "get_kinetic_energy: size Hloc(*,:,*) /= Norb"
  !   ! check Smats !
  !   if(size(Sigma,1)/=Nlat) stop "get_kinetic_energy: size Smats(:,*,*,*) /= Nlat"
  !   if(size(Sigma,2)/=size(Sigma,3)) stop "get_kinetic_energy: size Smats(*,:,*,*) /= Smats(*,*,:,*) "
  !   if(size(Sigma,2)/=Norb) stop "get_kinetic_energy: size Smats(*,:,*,*) /= Norb"
  !   if(size(Sigma,4)/=Lmats) stop "get_kinetic_energy: size Smats(*,*,*,:) /= Lmats"
  !   !
  !   allocate(Sigma_HF(Nlat,Norb,Norb))
  !   allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Tk(No,No))
  !   allocate(Gk(No,No))
  !   !
  !   Sigma_HF = dreal(Sigma(:,:,:,Lmats))
  !   !
  !   ed_Ekin_lattice=0.d0
  !   sumMats=0.d0
  !   Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
  !   if(mpiID==0) write(LOGfile,*) "Kinetic energy computation"
  !   if(mpiID==0)call start_timer
  !   do ik=1,Lk       
  !      Ak =  Hk(:,:,ik) !- xmu*Zk
  !      Bk = -Hk(:,:,ik) !+ xmu*Zk      
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               Bk(io,jo) = Bk(io,jo) - Sigma_HF(ilat,iorb,jorb) 
  !               Bk(io,jo) = Bk(io,jo) - Hloc(ilat,iorb,jorb)
  !            end do
  !         end do
  !      end do
  !      !Gk(iw) computation
  !      Gk=zero
  !      sumMatstmp=0.d0      
  !      sumMatsk=0.d0      
  !      do i=1+mpiID,Lmats,mpiSIZE
  !         Gk(:,:) = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
  !         do ilat=1,Nlat
  !            do iorb=1,Norb
  !               do jorb=1,Norb
  !                  io = (ilat-1)*Norb + iorb
  !                  jo = (ilat-1)*Norb + jorb
  !                  Gk(io,jo) = Gk(io,jo) - Sigma(ilat,iorb,jorb,i)
  !                  Gk(io,jo) = Gk(io,jo) - Hloc(ilat,iorb,jorb)
  !               end do
  !            end do
  !         end do
  !         !
  !         call matrix_inverse(Gk(:,:))
  !         !
  !         Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
  !         Ck = matmul(Ak,Gk(:,:) - Tk)
  !         do io=1,No
  !            sumMatstmp = sumMatstmp + Wtk(ik)*Ck(io,io)
  !         end do
  !      enddo
  !      call MPI_ALLREDUCE(sumMatstmp,sumMatsk,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !      sumMats=sumMats+sumMatsk
  !      if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !   enddo
  !   if(mpiID==0)call stop_timer
  !   spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  !   sumMats=sumMats/beta*2.d0*spin_degeneracy
  !   !
  !   ! subtract tails
  !   Tail0=0d0
  !   Tail1=0d0
  !   do ik=1,Lk
  !      Ak= Hk(:,:,ik) !- xmu*Zk
  !      Bk=-Hk(:,:,ik) !+ xmu*Zk      
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               ! Bk(io,jo) = Bk(io,jo) - Sigma_HF(iorb,jorb,ilat)
  !               ! if(iorb==jorb) Bk(io,jo) = Bk(io,jo) - elocal_(ilat) ![tmp...think to something more clever]
  !               Bk(io,jo) = Bk(io,jo) - Sigma_HF(ilat,iorb,jorb)
  !               Bk(io,jo) = Bk(io,jo) - Hloc(ilat,iorb,jorb)
  !            end do
  !         end do
  !      end do
  !      Ck= matmul(Ak,Bk)
  !      do io=1,No
  !         Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak(io,io)
  !         Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck(io,io)
  !      end do
  !   enddo
  !   Tail0=spin_degeneracy*Tail0
  !   Tail1=spin_degeneracy*Tail1*beta
  !   ed_Ekin_lattice=sumMats+Tail0+Tail1
  !   ed_Ekin_lattice=ed_Ekin_lattice/dble(Nlat)
  !   deallocate(Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  ! end function kinetic_energy_lattice_OLD


end module ED_WRAP_ENERGY
