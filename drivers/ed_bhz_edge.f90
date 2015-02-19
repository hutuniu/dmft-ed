program ed_bhz_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none

  integer                                       :: iloop
  integer                                       :: Npts,Nlso
  integer                                       :: Nso
  integer                                       :: ilat
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb(2)
  real(8),allocatable                           :: Bath(:,:,:),Bath_(:,:,:)
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !hamiltonian input:
  complex(8),allocatable                        :: Hkr(:,:,:),bhzHloc(:,:),Hloc(:,:,:,:,:),S0(:,:,:,:,:)
  real(8),allocatable                           :: Wtk(:)
  real(8),allocatable                           :: kxgrid(:)
  integer                                       :: Nk,Ly,Nkpath
  real(8)                                       :: mh,lambda,wmixing
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
#ENDIF

  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=10)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)

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

  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = reshape_Hloc(bhzHloc,Nlat,Nspin,Norb)
  if(mpiID==0)then
     do ilat=1,Nlat
        call write_Hloc(Hloc(ilat,:,:,:,:))
     enddo
  endif

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb(1),Nb(2)))
  allocate(Bath_(Nlat,Nb(1),Nb(2)))
  call ed_init_solver_lattice(bath)


  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   
     bath_=bath
     ! solve the impurities on each y-layer
     call ed_solve_lattice(bath,Hloc)
     ! retrieve the self-energies
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     call ed_get_sigma_real_lattice(Sreal,Nlat)
     ! compute the local gf:
     call ed_get_gloc_lattice(Hkr,Wtk,Gmats,Greal,Smats,Sreal,iprint=0)
     ! compute the Weiss field
     call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc)
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf_lattice(bath,Delta,Hloc,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf_lattice(bath,Delta,Hloc,ispin=2)
     else
        do ilat=1,Nlat
           call spin_symmetrize_bath(bath(ilat,:,:))
        enddo
     endif
     bath=wmixing*bath + (1.d0-wmixing)*bath_
     if(mpiID==0)converged = check_convergence(delta(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
  end do


#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#ENDIF

contains



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build the BHZ-edge Hamiltonian H(k_x,R_y)
  !+-----------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional           :: file
    real(8),dimension(:,:),allocatable  :: kpath
    if(mpiID==0)write(LOGfile,*)"Build H(k_x,r_y) for BHZ-stripe:"
    if(mpiID==0)write(*,*)"# of k-points     :",Nk
    if(mpiID==0)write(*,*)"# of layers       :",Nlat
    if(mpiID==0)write(*,*)"# of SO-bands     :",Nso
    allocate(kxgrid(Nk))
    allocate(Hkr(Nlso,Nlso,Nk))
    kxgrid = kgrid(Nk)
    Hkr    = build_Hkr_model(hkr_model,Ly,Nso,kxgrid,[0d0],[0d0],pbc=.false.)
    call write_hk_w90("Hkrfile.data",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Ly,&
         Hk=Hkr,&
         kxgrid=kxgrid,kygrid=[0d0],kzgrid=[0d0])
    allocate(Wtk(Nk))
    Wtk = 1d0/Nk
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = extract_Hloc(Hkr,Nlat,Nspin,Norb)
    Npts=3
    allocate(Kpath(Npts,1))
    kpath(1,:)=[-1]*pi
    kpath(2,:)=[ 0]*pi
    kpath(3,:)=[ 1]*pi
    call solve_HkR_along_BZpath(hkr_model,Ly,Nso,kpath,Nkpath,"hkr_Eigenbands.nint",pbc=.false.)
  end subroutine build_hkr




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ and BHZ-edge model hamiltonians
  !+-----------------------------------------------------------------------------+!
  !BHZ:
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
  end function hk_model
  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky
    complex(8),dimension(2,2) :: hk
    hk = (mh-cos(kx)-cos(ky))*pauli_tau_z + lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y 
  end function hk_bhz2x2

  !BHZ on a stripe geometry;
  function hkr_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=h0_rk_bhz(kx,N)
    enddo
    do i=1,Nlat-1
       Idmin=1+(i-1)*N
       Idmax=i*N
       Itmin=i*N+1
       Itmax=(i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=t0_rk_bhz(N)
       Hrk(Itmin:Itmax,Idmin:Idmax)=conjg(transpose(t0_rk_bhz(N)))
    enddo
    if(pbc)then
       Itmin=(Nlat-1)*N+1
       Itmax= Nlat*N
       Hrk(1:N,Itmin:Itmax)=conjg(transpose(t0_rk_bhz(N)))!T0
       Hrk(Itmin:Itmax,1:N)=t0_rk_bhz(N)!conjg(transpose(T0))
    endif
  end function hkr_model

  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    complex(8),dimension(4,4)  :: gamma1,gamma5
    !gamma1=kronecker_product_pauli_matrices(pauli_sigma_x,pauli_sigma_z)
    !gamma5=kronecker_product_pauli_matrices(pauli_sigma_z,pauli_sigma_0)
    gamma1 = kron(pauli_z,pauli_x)
    gamma5 = kron(pauli_0,pauli_z)
    H = (mh-cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    complex(8),dimension(4,4)  :: gamma2,gamma5
    !gamma2=kronecker_product_pauli_matrices(-pauli_sigma_y,pauli_sigma_0)
    !gamma5=kronecker_product_pauli_matrices(pauli_sigma_z,pauli_sigma_0)
    gamma2 = kron(pauli_0,-pauli_y)
    gamma5 = kron(pauli_0,pauli_z)
    H = xi*0.5d0*lambda*gamma2 + 0.5d0*gamma5
  end function T0_rk_bhz



end program ed_bhz_edge
