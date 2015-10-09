!  4DW----3UP
!   |      |
!   |      |
!  1UP----2DW
program ed_bhz_afm_2d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                   :: ip,iloop,Lk,Nso,Nlso,ispin,iorb
  logical                   :: converged
  integer                   :: Nineq
  !Bath:
  integer                   :: Nb
  real(8),allocatable       :: Bath(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable    :: Delta(:,:,:,:,:,:)
  complex(8),allocatable    :: Smats(:,:,:,:,:,:)
  complex(8),allocatable    :: Sreal(:,:,:,:,:,:)
  complex(8),allocatable    :: Gmats(:,:,:,:,:,:)
  complex(8),allocatable    :: Greal(:,:,:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable    :: Hk(:,:,:),bhzHloc(:,:),Hloc(:,:,:,:,:)
  real(8),allocatable       :: Wtk(:)
  !variables for the model:
  integer                   :: Nktot,Nkx,Nkpath,unit
  real(8)                   :: mh,lambda,wmixing
  character(len=16)         :: finput
  character(len=32)         :: hkfile
  logical                   :: waverage,spinsym,fullsym
  !Dirac matrices:
  complex(8),dimension(4,4) :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

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

  Nlat=4                      !number of independent sites, 4 for AFM ordering
  Nineq=Nlat
  Nso=Nspin*Norb
  Nlso=Nlat*Nso!=16
  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=16)stop "Nlso != 16"
  if(sitesym.AND.fullsym)stop "sitesym AND fullsym both .true.! chose one symmetry!"

  if(fullsym)then
     Nineq=1
     write(*,*)"Using Nineq sites=",Nineq
     open(free_unit(unit),file="symmetries.used")
     write(unit,*)"Symmetries used are:"
     write(unit,*)"(site=2,l,s)=(site=1,l,-s)"
     write(unit,*)"(site=3,l,s)=(site=1,l,s)"
     write(unit,*)"(site=4,l,s)=(site=1,l,-s)"
     close(unit)
  endif

  if(spinsym)sb_field=0.d0

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  !Allocate Weiss Field:
  allocate(Delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  do ip=1,Nlat
     ed_file_suffix="_site"//reg(txtfy(ip))
     call ed_init_solver(bath(ip,:))
     call break_symmetry_bath(bath(ip,:),sb_field,(-1d0)**(ip+1))
     write(LOGfile,*)"Updated Hloc, site:",ip
  enddo


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
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
        Sreal(3,:,:,:,:,:) = Sreal(1,:,:,:,:,:)
        Smats(4,:,:,:,:,:) = Smats(2,:,:,:,:,:)
        Sreal(4,:,:,:,:,:) = Sreal(2,:,:,:,:,:)
     endif
     !
     ! compute the local gf:
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     !
     ! compute the Weiss field
     call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc,iprint=1)
     !
     ! fit the new bath
     do ip=1,Nineq
        ed_file_suffix="_is"//reg(txtfy(ip))
        call set_hloc(Hloc(ip,:,:,:,:))
        call ed_chi2_fitgf(delta(ip,1,1,:,:,:),bath(ip,:),ispin=1)
        if(.not.spinsym)then
           call ed_chi2_fitgf(delta(ip,2,2,:,:,:),bath(ip,:),ispin=2)
        else
           call spin_symmetrize_bath(bath(ip,:))
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
     elseif(.not.fullsym)then
        if(waverage)then
           bath(1,:)=(bath(1,:)+bath(3,:))/2.d0
           bath(2,:)=(bath(2,:)+bath(4,:))/2.d0
        endif
        bath(3,:)=bath(1,:)
        bath(4,:)=bath(2,:)
     endif
     !
     ! Mixing:
     if(iloop>1)bath = wmixing*bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     ! Convergence
     converged = check_convergence(delta(1,1,1,1,1,1:Lfit),dmft_error,nsuccess,nloop)
     !
     call end_loop
  enddo
  print*,"Bravo"


contains




  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk(file)
    character(len=*)                        :: file
    integer                                 :: Npts
    integer                                 :: i,j,k,ik,iorb,jorb
    integer                                 :: ix,iy,iz
    real(8)                                 :: kx,ky,kz
    real(8),dimension(:),allocatable        :: kxgrid,wm,wr
    real(8),dimension(:,:),allocatable      :: kpath
    real(8)                                 :: n(Nlso)
    complex(8)                              :: w
    complex(8)                              :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    !
    Nktot=Nkx*Nkx
    allocate(Hk(Nlso,Nlso,Nktot))
    allocate(kxgrid(Nkx))
    write(LOGfile,*)"Build H(k) AFM-BHZ 2d:"
    write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    kxgrid = kgrid(Nkx)
    Hk = build_hk_model(hk_model,Nlso,kxgrid,kxgrid,[0d0])
    call write_hk_w90(trim(file),Nlso,&
         Nd=Nso,&
         Np=0,   &
         Nineq=4,&
         hk=Hk,  &
         kxgrid=kxgrid,&
         kygrid=kxgrid,&
         kzgrid=[0d0])
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    allocate(Wtk(Nktot))
    Wtk=1.d0/dble(Nktot)
    !
    !
    !Build the local GF:
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
    open(10,file="observables.nint")
    write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
    close(10)
    write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)/size(n)
    !
    !
    !solve along the standard path in the 2D BZ.
    Npts=4
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_Gamma
    kpath(2,:)=kpoint_X1
    kpath(3,:)=kpoint_M1
    kpath(4,:)=kpoint_Gamma
    call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
         colors_name=[character(len=10) :: 'red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue'],&
         points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
         file="Eigenbands_afm.nint")
  end subroutine build_hk





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  !+-----------------------------------------------------------------------------+!
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    complex(8),dimension(Nso,Nso) :: M
    complex(8),dimension(Nso,Nso) :: tx,ty,thx,thy
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = 2d0*kpoint(1)
    ky = 2d0*kpoint(2)
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    ! H4 =  | m1                 tx + tx^+.e^-ikx    0                   ty^+ + ty.e^iky |
    !       | tx^+ + tx.e^-ikx   m2                  ty^+ + ty.e^-iky    0               |
    !       | 0                  ty + ty^+e^-iky     m3                  tx^+ + tx.e^ikx |
    !       | ty + ty^+e^-iky    0                   tx + tx^+.e^-ikx    m4              |
    !
    hk(1:4,1:4)    = M
    hk(1:4,5:8)    = tx  + thx*exp(-xi*kx)
    hk(1:4,9:12)   = zero
    hk(1:4,13:16)  = thy + ty*exp(xi*ky)
    !
    hk(5:8,1:4)    = thx + tx*exp(xi*kx)
    hk(5:8,5:8)    = M
    hk(5:8,9:12)   = thy + ty*exp(xi*ky) 
    hk(5:8,13:16)  = zero
    !
    hk(9:12,1:4)   = zero
    hk(9:12,5:8)   = ty  + thy*exp(-xi*ky)
    hk(9:12,9:12)  = M
    hk(9:12,13:16) = thx + tx*exp(xi*kx)
    !
    hk(13:16,1:4)  = ty  + thy*exp(-xi*ky)
    hk(13:16,5:8)  = zero
    hk(13:16,9:12) = tx  + thx*exp(-xi*kx)
    hk(13:16,13:16)= M
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
  !PURPOSE:  Read Self-energies in Matsubara and Real freq.
  !+-----------------------------------------------------------------------------+!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nlat)stop "read_sigma: error in dim 1. Nlat"
    if(size(sigma,2)/=Nspin)stop "read_sigma: error in dim 2. Nspin"
    if(size(sigma,4)/=Norb)stop "read_sigma: error in dim 4. Norb"
    L=size(sigma,6);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 6. Lmats/Lreal"
    do ip=1,Nineq
       do iorb=1,Norb
          unit=free_unit()
          if(L==Lreal)then
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw_is"//reg(txtfy(ip))//".ed"
          elseif(L==Lmats)then
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw_is"//reg(txtfy(ip))//".ed"
          endif
          write(*,*)"read from file=","impSigma"//reg(suffix)
          open(unit,file="impSigma"//reg(suffix),status='old')
          do i=1,L
             read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
             forall(ispin=1:Nspin)sigma(ip,ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
          enddo
          close(unit)
       enddo
    enddo
    if(sitesym)then
       Sigma(4,:,:,:,:,:) = Sigma(1,:,:,:,:,:)
       Sigma(3,:,:,:,:,:) = Sigma(2,:,:,:,:,:)
    endif
    if(fullsym)then
       do ispin=1,2
          Sigma(2,ispin,ispin,:,:,:)=Sigma(1,3-ispin,3-ispin,:,:,:)
       enddo
       Sigma(4,:,:,:,:,:) = Sigma(1,:,:,:,:,:)
       Sigma(3,:,:,:,:,:) = Sigma(2,:,:,:,:,:)
    endif
  end subroutine read_sigma


end program ed_bhz_afm_2d



