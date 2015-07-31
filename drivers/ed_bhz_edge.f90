program ed_bhz_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none

  integer                                       :: iloop
  integer                                       :: Nlso
  integer                                       :: Nso
  integer                                       :: ilat
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath(:,:),Bath_Prev(:,:)
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc,S0
  !gamma matrices:
  complex(8),allocatable,dimension(:,:)         :: gamma1
  complex(8),allocatable,dimension(:,:)         :: gamma2
  complex(8),allocatable,dimension(:,:)         :: gamma5

  real(8),allocatable,dimension(:)              :: Wtk
  real(8),allocatable,dimension(:)              :: kxgrid
  integer                                       :: Nk,Ly,Nkpath
  real(8)                                       :: e0,mh,lambda,wmixing
  logical                                       :: spinsym
  character(len=16)                             :: finput
  character(len=32)                             :: hkfile

#ifdef _MPI_INEQ
  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(nkpath,"NKPATH",finput,default=250)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))


  !set the global number of lattice sites equal to the number of layers along the y-axis
  Nlat = Ly

  !set the local number of total spin-orbitals (4)
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso  = Nspin*Norb

  !set the total lattice-spin-orbit dimension:
  Nlso=Nlat*Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb))
  S0=zero


  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)
  if(mpiID==0)then
     do ilat=1,Nlat
        call write_Hloc(Hloc(ilat,:,:,:,:))
     enddo
  endif

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb), Bath_Prev(Nlat,Nb) )
  call ed_init_solver_lattice(bath)

  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   
     bath_prev=bath
     ! solve the impurities on each y-layer
     call ed_solve_lattice(bath,Hloc)
     ! retrieve the self-energies
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     call ed_get_sigma_real_lattice(Sreal,Nlat)
     S0 = Smats(:,:,:,:,:,1)
     ! compute the local gf:
     call ed_get_gloc_lattice(Hkr,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     ! compute the Weiss field
     call ed_get_weiss_lattice(Nlat,Gmats,Smats,Delta,Hloc)
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf_lattice(bath,Delta,Hloc,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf_lattice(bath,Delta,Hloc,ispin=2)
     else
        do ilat=1,Nlat
           call spin_symmetrize_bath(bath(ilat,:))
        enddo
     endif
     bath=wmixing*bath + (1.d0-wmixing)*bath_prev
     if(mpiID==0)converged = check_convergence(delta(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
  end do

  if(mpiID==0)call build_eigenbands()

#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif

contains



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build the BHZ Hamiltonian H(k_x,R_y) on the STRIPE along Y
  !+-----------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional          :: file
    !
    !SETUP THE GAMMA MATRICES:
    allocate(gamma1(Nso,Nso),gamma2(Nso,Nso),gamma5(Nso,Nso))
    gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x )
    gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y )
    gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z )
    !
    !SETUP THE H(kx,Ry):
    if(mpiID==0)then
       write(LOGfile,*)"Build H(k_x,r_y) for BHZ-stripe:"
       write(*,*)"# of k-points     :",Nk
       write(*,*)"# of layers       :",Nlat
       write(*,*)"# of SO-bands     :",Nso
    endif
    !
    allocate(kxgrid(Nk))
    allocate(Hkr(Nlso,Nlso,Nk))
    kxgrid = kgrid(Nk)
    Hkr    = build_Hkr_model(bhz_edge_model,Ly,Nso,kxgrid,[0d0],[0d0],pbc=.false.)
    call write_hk_w90("Hkrfile.data",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Ly,&
         Hk=Hkr,&
         kxgrid=kxgrid,kygrid=[0d0],kzgrid=[0d0])
    allocate(Wtk(Nk))
    Wtk = 1d0/Nk
    !
    !SETUP THE LOCAL PART Hloc(Ry)
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = extract_Hloc(Hkr,Nlat,Nspin,Norb)
    !
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    if(mpiID==0)call build_eigenbands()
    !
  end subroutine build_hkr



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve H_BHZ(k_x,R_y) along the 1d -pi:pi path in the BZ.
  !+-----------------------------------------------------------------------------+!  
  subroutine build_eigenbands()
    real(8),dimension(:,:),allocatable :: kpath
    integer                            :: Npts
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    Npts=3
    allocate(Kpath(Npts,1))
    kpath(1,:)=[-1]*pi
    kpath(2,:)=[ 0]*pi
    kpath(3,:)=[ 1]*pi
    call solve_HkR_along_BZpath(bhz_edge_model,Ly,Nso,kpath,Nkpath,"hkr_Eigenbands.nint",pbc=.false.)
  end subroutine build_eigenbands



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ-edge model hamiltonian
  !+-----------------------------------------------------------------------------+!
  !BHZ on a stripe geometry;
  function bhz_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    Hmat=h0_rk_bhz(kx,N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax) = Hmat + select_block(i,S0)
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 + i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax) = Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax) = TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax) = TmatH
       Hrk(Itmin:Itmax,1:N) = Tmat
    endif
  end function bhz_edge_model

  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = (mh-e0*cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = -0.5d0*e0*gamma5 + xi*0.5d0*lambda*gamma2
  end function T0_rk_bhz

end program ed_bhz_edge
