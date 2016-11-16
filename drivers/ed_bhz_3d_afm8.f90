!       8UP-----7DW
!       /|      /|
!      / |  BG / |
!     /  |    /  |
!    /  5DW--/--6UP        
!  4DW--/--3UP  /
!   |  /    |  /
!   | / FG  | /
!   |/      |/
!  1UP-----2DW
program ed_bhz_afm_3d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none

  integer                                       :: i,j,ip,iloop,Nso,Nlso,ispin,iorb
  logical                                       :: converged
  integer                                       :: Nineq
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk]
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  real(8),allocatable,dimension(:)              :: Wtk
  !variables for the model:
  integer                                       :: Nktot,Nkx,Nkpath,unit
  real(8)                                       :: mh,lambda,wmixing
  character(len=16)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym,fullsym
  !Dirac matrices:
  complex(8),dimension(4,4)                     :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

#ifdef _MPI_INEQ
  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')
  call parse_input_variable(nkx,"NKX",finput,default=25)
  call parse_input_variable(Nkpath,"Nkpath",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"spinsym",finput,default=.false.)
  call parse_input_variable(fullsym,"fullsym",finput,default=.true.)
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call ed_read_input(trim(finput))

  Nlat=8                      !number of independent sites, 4 for AFM ordering
  Nineq=Nlat
  Nso=Nspin*Norb
  Nlso=Nlat*Nso!=16
  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=32)stop "Nlso != 32"

  if(fullsym)then
     Nineq=1
     if(mpiID==0)then
        write(*,*)"Using Nineq sites=",Nineq
        open(free_unit(unit),file="symmetries.used")
        write(unit,*)"Symmetries used are:"
        write(unit,*)"       8UP-----7DW"
        write(unit,*)"       /|      /| "
        write(unit,*)"      / |  BG / | "
        write(unit,*)"     /  |    /  | "
        write(unit,*)"    /  5DW--/--6UP"         
        write(unit,*)"  4DW--/--3UP  /  "
        write(unit,*)"   |  /    |  /   " 
        write(unit,*)"   | / FG  | /    " 
        write(unit,*)"   |/      |/     "
        write(unit,*)"  1UP-----2DW     "
        close(unit)
     endif
  endif
  if(spinsym)sb_field=0.d0

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  !Allocate Field:
  !
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver_lattice(Bath)
  do ip=1,Nineq
     call break_symmetry_bath(Bath(ip,:),sb_field,sb_field_sign(ip))
     Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  enddo

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")
     !
     ! Solve the impurity problems
     call ed_solve_lattice(Bath,Hloc_ineq,iprint=1)
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
     do ip=1,Nineq
        Smats(ip,:,:,:,:,:) = Smats_ineq(ip,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(ip,:,:,:,:,:)
     enddo
     if(fullsym)then
        call neel_symmetrize_function(Smats)
        call neel_symmetrize_function(Sreal)       
     endif
     !
     ! Compute the local gf:
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
        Greal_ineq(ip,:,:,:,:,:) = Greal(ip,:,:,:,:,:)
     enddo
     !
     ! Compute the Weiss field
     call ed_get_weiss_lattice(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=1)
     !
     ! Fit the new bath
     call ed_chi2_fitgf_lattice(Bath,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        do ip=1,Nineq
           call spin_symmetrize_bath(Bath(ip,:))
        enddo
     else
        call ed_chi2_fitgf_lattice(Bath,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     ! Check convergence
     if(mpiID==0)converged = check_convergence(Weiss_ineq(1,1,1,1,1,1:Lfit),dmft_error,nsuccess,nloop)
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     !
     ! End Dmft loop
     if(mpiID==0) call end_loop()
  enddo


#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif
  print*,"Bravo"




contains




  !--------------------------------------------------------------------!
  !PURPOSE: build the H(k) for the bhz-afm model.
  !--------------------------------------------------------------------!
  subroutine build_hk(file,bool)
    character(len=*)                   :: file
    integer                            :: Npts
    integer                            :: i,j,k,ik,iorb,jorb
    integer                            :: ix,iy,iz
    real(8)                            :: kx,ky,kz
    real(8),dimension(:),allocatable   :: kxgrid,wm,wr
    real(8),dimension(:,:),allocatable :: kpath
    real(8)                            :: n(Nlso)
    complex(8)                         :: w
    logical,optional                   :: bool
    complex(8)                         :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    !
    Nktot=Nkx*Nkx*Nkx
    allocate(Hk(Nlso,Nlso,Nktot))
    allocate(kxgrid(Nkx))
    if(mpiID==0)write(LOGfile,*)"Build H(k) AFM-BHZ 3d:"
    if(mpiID==0)write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    kxgrid = kgrid(Nkx)
    Hk = build_hk_model(hk_model,Nlso,kxgrid,kxgrid,kxgrid)
    call write_hk_w90(trim(file),Nlso,&
         Nd=Nso,&
         Np=0,   &
         Nineq=8,&
         hk=Hk,  &
         kxgrid=kxgrid,&
         kygrid=kxgrid,&
         kzgrid=kxgrid)
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    allocate(Wtk(Nktot))
    Wtk=1.d0/dble(Nktot)
    !
    !
    Npts=8
    allocate(kpath(Npts,3))
    kpath(1,:)=[0,1,0]!X
    kpath(2,:)=[0,0,0]!G
    kpath(3,:)=[1,1,0]!M
    kpath(4,:)=[1,1,1]!R
    kpath(5,:)=[0,0,1]!Z
    kpath(6,:)=[1,0,1]!A
    kpath(7,:)=[0,0,0]!G
    kpath(8,:)=[0,0,1]!Z
    kpath=kpath*pi
    call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
         colors_name=[character(len=10) :: &
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue',&
         'red','blue','red','blue'],&
         points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
         file="Eigenbands_afm.nint")
  end subroutine build_hk





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Return the model Hamiltonian H(k) for a specific k point.
  !
  ! Hk =  | Hk^2d   Hk^z  |
  !       | Hk^z^+  Hk^2d |
  ! 
  ! Hk^2d =  | m1                 tx + tx^+.e^-ikx    0                   ty^+ + ty.e^iky |
  !          | tx^+ + tx.e^-ikx   m2                  ty^+ + ty.e^-iky    0               |
  !          | 0                  ty + ty^+e^-iky     m3                  tx^+ + tx.e^ikx |
  !          | ty + ty^+e^-iky    0                   tx + tx^+.e^-ikx    m4              |
  !
  ! Hk^z  =  | tz^+ + tz.e^-ikz   0                   0                   0               |
  !          | 0                  tz^+ + tz.e^-ikz    0                   0               |
  !          | 0                  0                   tz^+ + tz.e^-ikz    0               |
  !          | 0                  0                   0                   tz^+ + tz.e^-ikz|

  !+-----------------------------------------------------------------------------+!
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)              :: kpoint
    integer                           :: N
    real(8)                           :: kx,ky,kz
    complex(8),dimension(N,N)         :: hk
    complex(8),dimension(4*Nso,4*Nso) :: h2d,hz,hzh
    complex(8),dimension(Nso,Nso)     :: M
    complex(8),dimension(Nso,Nso)     :: tx,ty,tz,thx,thy,thz
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = 2d0*kpoint(1)
    ky = 2d0*kpoint(2)
    kz = 2d0*kpoint(3)
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    tz = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma3
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    thz= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma3
    !
    h2d(1:4,1:4)    = M
    h2d(1:4,5:8)    = tx  + thx*exp(-xi*kx)
    h2d(1:4,9:12)   = zero
    h2d(1:4,13:16)  = thy + ty*exp(xi*ky)
    !
    h2d(5:8,1:4)    = thx + tx*exp(xi*kx)
    h2d(5:8,5:8)    = M
    h2d(5:8,9:12)   = thy + ty*exp(xi*ky) 
    h2d(5:8,13:16)  = zero
    !
    h2d(9:12,1:4)   = zero
    h2d(9:12,5:8)   = ty  + thy*exp(-xi*ky)
    h2d(9:12,9:12)  = M
    h2d(9:12,13:16) = thx + tx*exp(xi*kx)
    !
    h2d(13:16,1:4)  = ty  + thy*exp(-xi*ky)
    h2d(13:16,5:8)  = zero
    h2d(13:16,9:12) = tx  + thx*exp(-xi*kx)
    h2d(13:16,13:16)= M
    !
    Hz              = zero
    Hz(1:4,1:4)     = thz + tz*exp(-xi*kz)
    Hz(5:8,5:8)     = thz + tz*exp(-xi*kz)
    Hz(9:12,9:12)   = thz + tz*exp(-xi*kz)
    Hz(13:16,13:16) = thz + tz*exp(-xi*kz)
    !
    Hzh             = zero
    Hzh(1:4,1:4)    = tz + thz*exp(xi*kz)
    Hzh(5:8,5:8)    = tz + thz*exp(xi*kz)
    Hzh(9:12,9:12)  = tz + thz*exp(xi*kz)
    Hzh(13:16,13:16)= tz + thz*exp(xi*kz)
    !
    Hk(1:16,1:16)   = H2d
    Hk(1:16,17:32)  = Hz
    Hk(17:32,1:16)  = Hzh
    Hk(17:32,17:32) = H2d
  end function hk_model






  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Neel symmetrize a function of rank: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !  using relation among the 8 inequivalent sites in the unit cell.
  !+-----------------------------------------------------------------------------+!
  subroutine neel_symmetrize_function(Flat)
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Flat
    integer                                         :: ispin,L
    L=size(Flat,6)
    call assert_shape(Flat,[Nlat,Nspin,Nspin,Norb,Norb,L],"neel_symmetrize_function","Flat")
    do ispin=1,2
       Flat(2,ispin,ispin,:,:,:)=Flat(1,3-ispin,3-ispin,:,:,:)
    enddo
    Flat(3,:,:,:,:,:) = Flat(1,:,:,:,:,:)
    Flat(4,:,:,:,:,:) = Flat(2,:,:,:,:,:)
    Flat(5,:,:,:,:,:) = Flat(2,:,:,:,:,:)
    Flat(6,:,:,:,:,:) = Flat(1,:,:,:,:,:)
    Flat(7,:,:,:,:,:) = Flat(2,:,:,:,:,:)
    Flat(8,:,:,:,:,:) = Flat(1,:,:,:,:,:)
  end subroutine neel_symmetrize_function



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  return the direction of the symmetry breaking field  per site
  !+-----------------------------------------------------------------------------+!
  function sb_field_sign(ip) result(segno)
    integer :: ip
    integer :: i,j,ilat
    real(8) :: sgn(8),segno
    do j=0,1
       do i=1,4
          ilat=i + 4*j
          sgn(ilat)=-1d0**(i+j+1)
       enddo
    enddo
    segno = sgn(ip)
  end function sb_field_sign


end program ed_bhz_afm_3d



