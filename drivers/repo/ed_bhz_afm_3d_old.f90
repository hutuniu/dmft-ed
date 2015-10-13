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
  integer                   :: i,j,ip,iloop,Nso,Nlso,ispin,iorb
  logical                   :: converged
  integer                   :: Nineq
  !Bath:
  integer                   :: Nb
  real(8),allocatable       :: Bath(:,:),Bath_ineq(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:)    :: Weiss,Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)    :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)    :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)    :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)    :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable    :: Hk(:,:,:) ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk]
  complex(8),allocatable    :: bhzHloc(:,:),Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  real(8),allocatable       :: Wtk(:)
  !variables for the model:
  integer                   :: Nktot,Nkx,Nkpath,unit
  real(8)                   :: mh,lambda,wmixing
  character(len=16)         :: finput
  character(len=32)         :: hkfile
  logical                   :: waverage,spinsym,fullsym
  !Dirac matrices:
  complex(8),dimension(4,4) :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

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
  call parse_input_variable(waverage,"WAVERAGE",finput,default=.false.)
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
        write(unit,*)"(site=2,l,s)=(site=1,l,-s)"
        write(unit,*)"(site=3,l,s)=(site=1,l,s)"
        write(unit,*)"(site=4,l,s)=(site=1,l,-s)"
        write(unit,*)"(site=5,l,s)=(site=1,l,-s)"
        write(unit,*)"(site=6,l,s)=(site=1,l,s)"
        write(unit,*)"(site=7,l,s)=(site=1,l,-s)"
        write(unit,*)"(site=8,l,s)=(site=1,l,s)"
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
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
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
  allocate(Bath(Nlat,Nb))
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  do j=0,1
     do i=1,4
        ip=i + 4*j
        ed_file_suffix="_site"//reg(txtfy(ip))
        call ed_init_solver(bath(ip,:))
        call break_symmetry_bath(bath(ip,:),sb_field,(-1d0)**(i+j+1) )
        write(LOGfile,*)"Updated Hloc, site:"//reg(txtfy(ip))
     enddo
  enddo
  ! call ed_init_solver_lattice(Bath_ineq)
  ! do ip=1,Nineq
  !    call break_symmetry_bath(Bath_ineq(ip,:),sb_field,sb_field_sign(ip))
  !    Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  ! enddo

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     ! NEW
     ! call ed_solve_lattice(Bath_ineq,Hloc_ineq,iprint=1)
     ! call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     ! call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
     ! call ineq2lattice(Smats_ineq,Smats)
     ! call ineq2lattice(Sreal_ineq,Sreal)
     ! OLD
     do ip=1,Nineq
        write(LOGfile,*)"Solving site:",ip
        ed_file_suffix="_site"//reg(txtfy(ip))
        call set_hloc(Hloc(ip,:,:,:,:))
        call ed_solve(bath(ip,:))
        call ed_get_sigma_matsubara(Smats(ip,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(ip,:,:,:,:,:))
     enddo
     if(fullsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
        Smats(3,:,:,:,:,:) = Smats(1,:,:,:,:,:)
        Smats(4,:,:,:,:,:) = Smats(2,:,:,:,:,:)
        Smats(5,:,:,:,:,:) = Smats(2,:,:,:,:,:)
        Smats(6,:,:,:,:,:) = Smats(1,:,:,:,:,:)
        Smats(7,:,:,:,:,:) = Smats(2,:,:,:,:,:)
        Smats(8,:,:,:,:,:) = Smats(1,:,:,:,:,:)
        Sreal(3,:,:,:,:,:) = Sreal(1,:,:,:,:,:)
        Sreal(4,:,:,:,:,:) = Sreal(2,:,:,:,:,:)
        Sreal(5,:,:,:,:,:) = Sreal(2,:,:,:,:,:)
        Sreal(6,:,:,:,:,:) = Sreal(1,:,:,:,:,:)
        Sreal(7,:,:,:,:,:) = Sreal(2,:,:,:,:,:)
        Sreal(8,:,:,:,:,:) = Sreal(1,:,:,:,:,:)
     endif
     !
     ! compute the local gf:
     ! call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     ! call lattice2ineq(Gmats,Gmats_ineq)
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     !
     ! compute the Weiss field
     ! call ed_get_weiss_lattice(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=1)
     call ed_get_weiss_lattice(Gmats,Smats,Weiss,Hloc,iprint=1)
     !
     ! fit the new bath
     ! NEW
     ! call ed_chi2_fitgf_lattice(Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     ! if(spinsym)then
     !    do ip=1,Nineq
     !       call spin_symmetrize_bath(Bath_ineq(ip,:))
     !    enddo
     ! else
     !    call ed_chi2_fitgf_lattice(Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     ! endif
     ! OLD
     do ip=1,Nineq
        ed_file_suffix="_site"//reg(txtfy(ip))
        call set_hloc(Hloc(ip,:,:,:,:))
        call ed_chi2_fitgf(Weiss(ip,1,1,:,:,:),Bath(ip,:),ispin=1)
        if(spinsym)then
           call spin_symmetrize_bath(Bath(ip,:))
        else
           call ed_chi2_fitgf(Weiss(ip,2,2,:,:,:),Bath(ip,:),ispin=2)
        endif
     enddo
     !
     ! manipulate the baths according to the symmetry
     if(fullsym)then
        do ispin=1,2
           call copy_component_bath(bath(1,:),3-ispin,bath(2,:),ispin)
        enddo
        bath(3,:)=bath(1,:)
        bath(4,:)=bath(2,:)
        bath(5,:)=bath(2,:)
        bath(6,:)=bath(1,:)
        bath(7,:)=bath(2,:)
        bath(8,:)=bath(1,:)
     elseif(.not.fullsym)then
        if(waverage)then
           bath(1,:)=(bath(1,:)+bath(3,:))/2.d0
           bath(2,:)=(bath(2,:)+bath(4,:))/2.d0
        endif
        bath(3,:)=bath(1,:)
        bath(4,:)=bath(2,:)
        bath(5,:)=bath(2,:)
        bath(6,:)=bath(1,:)
        bath(7,:)=bath(2,:)
        bath(8,:)=bath(1,:)
     endif
     ! call ineq_symmetry_bath(bath_ineq)
     !
     ! Mixing:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     ! Convergence
     if(mpiID==0)converged = check_convergence(Weiss_ineq(1,1,1,1,1,1:Lfit),dmft_error,nsuccess,nloop)
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(mpiID==0) call end_loop()
  enddo


#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif
  print*,"Bravo"


contains




  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
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
    !Build the local GF:
    if(present(bool).AND.bool)then
       allocate(wm(Lmats),wr(Lreal))
       wm = pi/beta*real(2*arange(1,Lmats)-1,8)
       wr = linspace(wini,wfin,Lreal)
       Greal = zero
       Gmats = zero 
       do ik=1,Nktot
          do i=1,Lreal
             w = dcmplx(wr(i),eps)+xmu
             Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
          enddo
          do i=1,Lmats
             w = xi*wm(i)+xmu
             Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
          enddo
       enddo
       do iorb=1,Nlso
          call splot("Gloc_lso"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
          call splot("Gloc_lso"//reg(txtfy(iorb))//"_realw.nint",wr,&
               -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
          n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
       enddo
       !
       !plot observables
       if(mpiID==0)then
          open(10,file="observables.nint")
          write(10,"(33F20.12)")(n(iorb),iorb=1,Nlso),sum(n)/size(n)
          close(10)
          write(*,"(A,33F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)/size(n)
       endif
    endif
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
         colors_name=[character(len=10) :: 'red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue'],&
         points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
         file="Eigenbands_afm.nint")


  end subroutine build_hk





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
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









  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function inverse_gk(zeta,hk) result(gk)
    complex(8),dimension(Nlso,Nlso)   :: zeta,hk
    complex(8),dimension(Nlso,Nlso)   :: gk
    integer :: i
    gk = -Hk
    forall(i=1:Nlso)gk(i,i)=zeta(i,i) - hk(i,i)
    call matrix_inverse(gk)
  end function inverse_gk
  !
  function inverse_g0k(zeta,hk) result(g0k)
    complex(8)                      :: zeta
    complex(8),dimension(Nlso,Nlso) :: hk
    complex(8),dimension(Nlso,Nlso) :: g0k
    integer :: i
    g0k = -hk
    forall(i=1:Nlso)g0k(i,i)=zeta - hk(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Extend or Retrieve by Neel symmetry a function of the form 
  !  Fin : [Nineq][Nspin][Nspin][Norb][Norb][L]
  !  to another function/bath of the form
  !  Flat: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !  using relation among the 8 inequivalent sites in the unit cell.
  !
  !  Note: in this version Nineq=1, so we obtain all the components from Fin[1]
  !+-----------------------------------------------------------------------------+!
  subroutine ineq2lattice(Fineq,Flat)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Fineq
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Flat
    integer :: ispin,L
    L=size(Fineq,6)
    call assert_shape(Fineq,[Nineq,Nspin,Nspin,Norb,Norb,L],"ineq2lattice","Fineq")
    call assert_shape(Flat,[Nlat,Nspin,Nspin,Norb,Norb,L],"ineq2lattice","Flat")
    if(fullsym)then
       Flat(1,:,:,:,:,:) = Fineq(1,:,:,:,:,:)
       do ispin=1,2
          Flat(2,ispin,ispin,:,:,:)=Flat(1,3-ispin,3-ispin,:,:,:)
       enddo
       Flat(3,:,:,:,:,:) = Flat(1,:,:,:,:,:)
       Flat(4,:,:,:,:,:) = Flat(2,:,:,:,:,:)
       Flat(5,:,:,:,:,:) = Flat(2,:,:,:,:,:)
       Flat(6,:,:,:,:,:) = Flat(1,:,:,:,:,:)
       Flat(7,:,:,:,:,:) = Flat(2,:,:,:,:,:)
       Flat(8,:,:,:,:,:) = Flat(1,:,:,:,:,:)
    else
       if(Nineq/=Nlat)stop "ineq2lattice_function error: Nineq != Nlat "
       Flat=Fineq
    endif
  end subroutine ineq2lattice
  !
  subroutine lattice2ineq(Flat,Fineq)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Flat
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Fineq
    integer                                         :: ispin,L
    L=size(Fineq,6)
    call assert_shape(Fineq,[Nineq,Nspin,Nspin,Norb,Norb,L],"ineq2lattice","Fineq")
    call assert_shape(Flat,[Nlat,Nspin,Nspin,Norb,Norb,L],"ineq2lattice","Flat")
    if(fullsym)then
       Fineq(1,:,:,:,:,:) = Flat(1,:,:,:,:,:)
    else
       if(Nineq/=Nlat)stop "ineq2lattice_function error: Nineq != Nlat "
       Fineq = Flat
    endif
  end subroutine lattice2ineq



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Symmetrize the Bath of inequivalent sites.
  ! - if fullsym : do nothing (there is only 1 sites and is the one actually used 
  ! in the calculation)
  ! - else: 
  !   if waverage: get the average and symmetrize accordingly
  !
  !+-----------------------------------------------------------------------------+!
  subroutine ineq_symmetry_bath(bath)
    real(8),dimension(:,:),intent(inout) :: bath
    integer                              :: ispin,Nb
    Nb=size(Bath,2)
    call assert_shape(Bath,[Nineq,Nb],"ineq_symmetry_bath","Bath")
    if(waverage)then
       Bath(1,:)=(Bath(1,:)+Bath(3,:))/2.d0
       Bath(2,:)=(Bath(2,:)+Bath(4,:))/2.d0
       Bath(3,:)=Bath(1,:)
       Bath(4,:)=Bath(2,:)
       Bath(5,:)=Bath(2,:)
       Bath(6,:)=Bath(1,:)
       Bath(7,:)=Bath(2,:)
       Bath(8,:)=Bath(1,:)
    endif
  end subroutine ineq_symmetry_bath






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



