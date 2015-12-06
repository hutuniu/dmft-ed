
program ed_STO
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb,unit
  real(8),allocatable    :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:),Ti3dt2g_Hloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin
  real(8)                :: soc,ivb,wmixing,sumdens
  logical                :: surface,Hk_test
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym
  !convergence function
  complex(8),allocatable :: delta_conv(:,:,:),delta_conv_avrg(:)
  !density matrix
  real(8),allocatable    :: dm_eig(:)
  complex(8),allocatable :: density_matrix(:,:),dm_rot(:,:)
  !
  real(8),dimension(2)   :: Eout
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,   "FINPUT",           default='inputED_STO.in')
  call parse_input_variable(hkfile, "HKFILE",finput,    default="hkfile.in")
  call parse_input_variable(nk,     "NK",finput,        default=10)
  call parse_input_variable(nkpath, "NKPATH",finput,    default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,   default=0.5d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,   default=.true.)
  call parse_input_variable(soc,    "SOC",finput,       default=0.25d0)
  call parse_input_variable(ivb,    "IVB",finput,       default=0.02d0)
  call parse_input_variable(surface,"SURFACE",finput,   default=.false.)
  call parse_input_variable(Hk_test,"HK_TEST",finput,   default=.false.)
  call ed_read_input(trim(finput))

  !if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(density_matrix(Nspin*Norb,Nspin*Norb))
  allocate(dm_eig(Nspin*Norb),dm_rot(Nspin*Norb,Nspin*Norb))
  allocate(delta_conv(Nso,Nso,Lmats))
  allocate(delta_conv_avrg(Lmats))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  call build_hk_path

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(bath)
  call set_hloc(reshape_A1_to_A2(Ti3dt2g_Hloc))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call build_hk_path
     call ed_get_sigma_real(Sreal)
     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=3)
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=reshape_A1_to_A2(Ti3dt2g_Hloc),iprint=3)
     !density matrix
     if(ED_MPI_ID==0)then
        call ed_get_density_matrix(density_matrix,2,dm_eig,dm_rot)
     endif
     !Fit the new bath, starting from the old bath + the supplied delta 
     Bath_=bath
     if (ed_mode=="normal") then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        call spin_symmetrize_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(delta,bath)
     endif

     !MIXING:
     if(iloop>1) Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     delta_conv=zero
     delta_conv_avrg=zero
     do i=1,Lmats
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    if((ispin.eq.jspin).and.(iorb.eq.jorb)) then
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       delta_conv(io,jo,i)=delta(ispin,jspin,iorb,jorb,i)
                       delta_conv_avrg(i)=delta_conv_avrg(i)+delta_conv(io,jo,i)
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     delta_conv_avrg=delta_conv_avrg/Nso

     if(ED_MPI_ID==0) converged = check_convergence(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(ED_MPI_ID==0) converged = check_convergence_global(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(ED_MPI_ID==0) converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !if(ED_MPI_ID==0)converged = check_convergence_global(delta_conv(:,:,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif

     sumdens=sum(ed_get_dens())
     write(*,*) "sumdens",sumdens,"xmu",xmu,"converged",converged
     if(nread/=0.d0)call search_chemical_potential(xmu,sumdens,converged)
     write(*,*) "sumdens",sumdens,"xmu",xmu,"converged",converged

     if(ED_MPI_ID==0)call end_loop
  enddo
#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif
contains



!_______________________________________________________________________
!                            HAMILTONIAN
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: H(k) file for main program and write G0_loc
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky,kz    
    integer                             :: io,jo
    integer                             :: iorb,jorb,ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    complex(8),dimension(6,6)           :: theta
    real(8)                             :: wm(Lmats),wr(Lreal),dw

    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for STO:"

    if(.not.surface) Lk=Nk**3
    if(surface)      Lk=Nk**2
    
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nso,Nso,Lk));allocate(wtk(Lk));allocate(kxgrid(Nk),kygrid(Nk),kzgrid(Nk))
    kxgrid=0.0d0;kygrid=0.0d0;kzgrid=0.0d0
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    if(.not.surface) kzgrid = kgrid(Nk)

    if(.not.surface) Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,kzgrid)
    if(surface)      Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,[0d0])

    wtk = 1.0d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
        if(.not.surface) call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=kzgrid)
        if(surface)      call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=[0d0] )
    endif

    allocate(Ti3dt2g_Hloc(Nso,Nso))
    Ti3dt2g_Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=0d0

    if(ED_MPI_ID==0) then
       call write_Hloc(Ti3dt2g_Hloc)
       write(*,*)
       write(*,*) "Sum over k of H(k) nella versione A1"
       write(*,*) "real"
       do i=1,Nso
          write(*,'(6F10.4)') (real(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(*,*) "complex"
       do i=1,Nso
          write(*,'(6F10.4)') (dimag(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(*,*)
    endif

    !Build the local GF in the spin-orbital Basis:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu , Hk(:,:,ik) )/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,Gmats(io,jo,:))
                call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    !Build the local GF in the J-Jz Basis: theta funziona per Gloc nell'ordine di Z o per H *solo al punto gamma*
    if(.not.Hk_test) then
    call build_rotation(theta)
    do i=1,Lmats
       Gmats(:,:,i)=matmul(transpose(conjg(theta)),matmul(Gmats(:,:,i),theta))
    enddo
    do i=1,Lreal
       Greal(:,:,i)=matmul(transpose(conjg(theta)),matmul(Greal(:,:,i),theta))
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed",wm,Gmats(io,jo,:))
                call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    endif

  end subroutine build_hk

  !---------------------------------------------------------------------
  !PURPOSE: GET STO HAMILTONIAN in the A1 shape
  !---------------------------------------------------------------------
  function hk_Ti3dt2g(kvec,N) result(hk)
    real(8),dimension(:)        :: kvec
    complex(8),dimension(N,N)   :: hk
    complex(8),dimension(2,2)   :: s_x,s_y,s_z
    complex(8),dimension(2,2)   :: t_inter
    real(8)                     :: kx,ky,kz
    integer                     :: N,ndx
    real(8),dimension(Norb,0:6) :: HoppingMatrix

    s_x=cmplx(0.0d0,0.0d0);s_y=cmplx(0.0d0,0.0d0);s_z=cmplx(0.0d0,0.0d0)
    s_x(1,2)=cmplx(1.0d0,0.0d0); s_x(2,1)=cmplx(1.0d0,0.0d0)
    s_y(1,2)=cmplx(0.0d0,-1.0d0);s_y(2,1)=cmplx(0.0d0,1.0d0)
    s_z(1,1)=cmplx(1.0d0,0.0d0); s_z(2,2)=cmplx(-1.0d0,0.0d0)
    kx=kvec(1);ky=kvec(2);kz=kvec(3)

    Hk=zero
    if(Hk_test) then
       do i=1,Norb
          ndx=2*i-1
          Hk(ndx:ndx+1,ndx:ndx+1) = band_cos_omo(kx,ky,kz)
       enddo
    else
       call get_hopping(HoppingMatrix)
       do i=1,Norb
          ndx=2*i-1
          Hk(ndx:ndx+1,ndx:ndx+1) = orbital_dispersion(kx,ky,kz,HoppingMatrix(i,:))
       enddo
       !upper triangle
       Hk(1:2,3:4)= +xi * s_z * soc/2.
       Hk(1:2,5:6)= -xi * s_y * soc/2. + ivb*2*xi*sin(kx)*eye(2)
       Hk(3:4,5:6)= +xi * s_x * soc/2. + ivb*2*xi*sin(ky)*eye(2)
       !lower triangle
       do i=1,Nspin*Norb
          do j=1,Nspin*Norb
             Hk(j,i)=conjg(Hk(i,j))
          enddo
       enddo
    endif

    !A1 shape: [Norb*Norb]*Nspin
    Hk = reshape_Z_to_A1(Hk)
  
  end function hk_Ti3dt2g

  function hk_Ti3dt2g_Hartree(kvec,N) result(hk)
    real(8),dimension(:)        :: kvec
    complex(8),dimension(N,N)   :: hk
    complex(8),dimension(2,2)   :: s_x,s_y,s_z
    complex(8),dimension(2,2)   :: t_inter
    real(8)                     :: kx,ky,kz
    integer                     :: N,ndx
    real(8),dimension(Norb,0:6) :: HoppingMatrix

    s_x=cmplx(0.0d0,0.0d0);s_y=cmplx(0.0d0,0.0d0);s_z=cmplx(0.0d0,0.0d0)
    s_x(1,2)=cmplx(1.0d0,0.0d0); s_x(2,1)=cmplx(1.0d0,0.0d0)
    s_y(1,2)=cmplx(0.0d0,-1.0d0);s_y(2,1)=cmplx(0.0d0,1.0d0)
    s_z(1,1)=cmplx(1.0d0,0.0d0); s_z(2,2)=cmplx(-1.0d0,0.0d0)
    kx=kvec(1);ky=kvec(2);kz=kvec(3)

    Hk=zero
    if(Hk_test) then
       do i=1,Norb
          ndx=2*i-1
          Hk(ndx:ndx+1,ndx:ndx+1) = band_cos_omo(kx,ky,kz)
       enddo
    else
       call get_hopping(HoppingMatrix)
       do i=1,Norb
          ndx=2*i-1
          Hk(ndx:ndx+1,ndx:ndx+1) = orbital_dispersion(kx,ky,kz,HoppingMatrix(i,:))
       enddo
       !upper triangle
       Hk(1:2,3:4)= +xi * s_z * soc/2.
       Hk(1:2,5:6)= -xi * s_y * soc/2. + ivb*2*xi*sin(kx)*eye(2)
       Hk(3:4,5:6)= +xi * s_x * soc/2. + ivb*2*xi*sin(ky)*eye(2)
       !lower triangle
       do i=1,Nspin*Norb
          do j=1,Nspin*Norb
             Hk(j,i)=conjg(Hk(i,j))
          enddo
       enddo
    endif

    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = ispin + (iorb-1)*Nspin
                jo = jspin + (jorb-1)*Nspin
                Hk(io,jo) = Hk(io,jo) + Smats(ispin,jspin,iorb,jorb,1)
             enddo
          enddo
       enddo
    enddo

    !A1 shape: [Norb*Norb]*Nspin
    Hk = reshape_Z_to_A1(Hk)
  
  end function hk_Ti3dt2g_Hartree

  !---------------------------------------------------------------------
  !PURPOSE: various 2x2 band BULK structures
  !---------------------------------------------------------------------
  function band_cos_omo(kx,ky,kz) result(hk)
    real(8)                         :: kx,ky,kz
    complex(8),dimension(2,2)       :: hk
    hk = zero
    hk(1,1) = -2.*(cos(kx)+cos(ky)+cos(kz))
    hk(2,2) = hk(1,1)
  end function band_cos_omo

  function orbital_dispersion(kx,ky,kz,t) result(hk)
    real(8)                           :: kx,ky,kz
    real(8),intent(in),dimension(0:6) :: t
    complex(8),dimension(2,2)         :: hk
    !perovskite dispersion
    hk = zero
    if(.not.surface)then
       hk(1,1) = t(0)                      & !onsite
                 -2.*t(1)*cos(kx)          & !t_100
                 -2.*t(2)*cos(ky)          & !t_010
                 -2.*t(3)*cos(kz)          & !t_001
                 -4.*t(4)*cos(ky)*cos(kz)  & !t_011
                 -4.*t(5)*cos(kx)*cos(kz)  & !t_101
                 -4.*t(6)*cos(kx)*cos(ky)    !t_110
       hk(2,2) = hk(1,1)
    else
       hk(1,1) = t(0)                      & !onsite
                 -2.*t(1)*cos(kx)          & !t_100
                 -2.*t(2)*cos(ky)          & !t_010
                 -1.*t(3)                  & !t_001
                 -2.*t(4)*cos(ky)          & !t_011
                 -2.*t(5)*cos(kx)          & !t_101
                 -4.*t(6)*cos(kx)*cos(ky)    !t_110
       hk(2,2) = hk(1,1)

    endif
  end function orbital_dispersion

  !---------------------------------------------------------------------
  !PURPOSE: Get the STO(bulk) Hamiltonian along path
  !---------------------------------------------------------------------
  subroutine build_hk_path()
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),allocatable :: kpath
    if(.not.surface)then
       if(ED_MPI_ID==0)then
          write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
          write(LOGfile,*)
       endif
       Npts = 7
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_M1
       kpath(2,:)=kpoint_R
       kpath(3,:)=kpoint_Gamma
       kpath(4,:)=kpoint_M1
       kpath(5,:)=kpoint_X1
       kpath(6,:)=kpoint_Gamma
       kpath(7,:)=kpoint_X1
       if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,&
            colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
            points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],&
            file="Eigenband_bulk.nint")
       if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_Ti3dt2g_Hartree,Nso,kpath,Lk,&
            colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
            points_name=[character(len=20) :: 'M', 'R', 'G', 'M', 'X', 'G', 'X'],&
            file="Eigenband_bulk_H.nint")


    else
       if(ED_MPI_ID==0)then
          if(ED_MPI_ID==0)write(LOGfile,*)"Build surface H(k) along the path M-X-G-X"
          write(LOGfile,*)
       endif
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_M1
       kpath(2,:)=kpoint_X1
       kpath(3,:)=kpoint_Gamma
       kpath(4,:)=kpoint_X1
       if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_Ti3dt2g,Nso,kpath,Lk,&
            colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
            points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],&
            file="Eigenband_surf.nint")
       if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_Ti3dt2g_Hartree,Nso,kpath,Lk,&
            colors_name=[character(len=20) :: 'red','green','blue','red','green','blue'],&
            points_name=[character(len=20) :: 'M', 'X', 'G', 'X'],&
            file="Eigenband_surf_H.nint")
    endif
  end subroutine build_hk_path

  !---------------------------------------------------------------------
  !PURPOSE: Build the hopping integrals in k-space for realistic bandstructure
  !---------------------------------------------------------------------
  subroutine get_hopping(T)
  real(8),dimension(Norb,0:6),intent(out)   ::  T
  real(8),dimension(Norb,0:6)               ::  T_bulk,T_VACSTO,T_LAOSTO
  real(8)                                   ::  Eo,t1,t2,t3
  real(8)                                   ::  t_010_yz,t_001_yz
  real(8)                                   ::  t_100_zx,t_001_zx
  real(8)                                   ::  t_100_xy,t_010_xy,t_001_xy

  !pristine lattice
  Eo = 3.31
  t1 = 0.276536
  t2 = 0.031329
  t3 = 0.076842

  !lattice distortion
  t_010_yz = 0.232 !se c'è solo l'abbassamento del Ti questo dovrebbe essere uguale a t1, magari c'è anche altro dovuto ad LAO
  t_001_yz = 0.475

  t_100_zx = 0.232
  t_001_zx = 0.475

  t_100_xy = 0.286
  t_010_xy = 0.286
  t_001_xy = 0.03

  !  BULK STO
  !orbital1 = YZ
  T_bulk(1,0) = Eo
  T_bulk(1,1) = t2
  T_bulk(1,2) = t1
  T_bulk(1,3) = t1
  T_bulk(1,4) = t3
  T_bulk(1,5) = 0.d0
  T_bulk(1,6) = 0.d0
  !orbital1 = ZX
  T_bulk(2,0) = Eo
  T_bulk(2,1) = t1
  T_bulk(2,2) = t2
  T_bulk(2,3) = t1
  T_bulk(2,4) = 0.d0
  T_bulk(2,5) = t3
  T_bulk(2,6) = 0.d0
  !orbital1 = XY
  T_bulk(3,0) = Eo
  T_bulk(3,1) = t1
  T_bulk(3,2) = t1
  T_bulk(3,3) = t2
  T_bulk(3,4) = 0.d0
  T_bulk(3,5) = 0.d0
  T_bulk(3,6) = t3

  ! VAC/STO
  T_VACSTO=T_bulk
  
  !  LAO/STO
  !orbital1 = YZ
  T_LAOSTO(1,0) = 1.087
  T_LAOSTO(1,1) = t2
  T_LAOSTO(1,2) = t_010_yz
  T_LAOSTO(1,3) = t_001_yz
  T_LAOSTO(1,4) = t3
  T_LAOSTO(1,5) = 0.d0
  T_LAOSTO(1,6) = 0.d0
  !orbital1 = ZX
  T_LAOSTO(2,0) = 1.087
  T_LAOSTO(2,1) = t_100_zx
  T_LAOSTO(2,2) = t2
  T_LAOSTO(2,3) = t_001_zx
  T_LAOSTO(2,4) = 0.d0
  T_LAOSTO(2,5) = t3
  T_LAOSTO(2,6) = 0.d0
  !orbital1 = XY
  T_LAOSTO(3,0) = 1.035
  T_LAOSTO(3,1) = t_100_xy
  T_LAOSTO(3,2) = t_010_xy
  T_LAOSTO(3,3) = t_001_xy
  T_LAOSTO(3,4) = 0.d0
  T_LAOSTO(3,5) = 0.d0
  T_LAOSTO(3,6) = t3

  if(.not.surface) then
     T=T_bulk
  else
     T=T_LAOSTO
  endif

  end subroutine get_hopping

  !---------------------------------------------------------------------
  !PURPOSE: Build the rotation "theta" that diagonalizes SOC
  !---------------------------------------------------------------------
  subroutine build_rotation(theta_)
    complex(8),dimension(6,6),intent(out) :: theta_
    theta_=cmplx(0.0d0,0.0d0)
    !J=1/2 jz=-1/2
    theta_(1,1)=-Xi
    theta_(3,1)=-1.0d0
    theta_(6,1)=+Xi
    theta_(:,1)=theta_(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    theta_(2,2)=-Xi
    theta_(4,2)=+1.0d0
    theta_(5,2)=-Xi
    theta_(:,2)=theta_(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    theta_(2,3)=-Xi
    theta_(4,3)=+1.0d0
    theta_(5,3)=+2.0d0*Xi
    theta_(:,3)=theta_(:,3)/sqrt(6.)
    !J=3/2 jz=-1/2
    theta_(1,4)=+Xi
    theta_(3,4)=-1.0d0
    theta_(:,4)=theta_(:,4)/sqrt(2.)
    !J=3/2 jz=+1/2
    theta_(2,5)=-Xi 
    theta_(4,5)=-1.0d0
    theta_(:,5)=theta_(:,5)/sqrt(2.)
    !J=3/2 jz=+3/2
    theta_(1,6)=+Xi
    theta_(3,6)=+1.0d0
    theta_(6,6)=+2.0d0*Xi
    theta_(:,6)=theta_(:,6)/sqrt(6.)
    theta_=reshape_Z_to_A1(theta_)
  end subroutine build_rotation


!_______________________________________________________________________
!                                    Gfs
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: G0_loc functions DA RIFARE ATTENZIONE CHE H(k) è nella forma A1
  !---------------------------------------------------------------------
  function inverse_g0k(iw,hk) result(g0k)
    implicit none
    complex(8)                                    :: iw
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: hk
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: g0k,g0k_tmp
    integer                                       :: i,ndx
    integer (kind=4), dimension(6)                :: ipiv
    integer (kind=1)                              :: ok
    integer (kind=4), parameter                   :: lwork=2000
    complex (kind=8), dimension(lwork)            :: work
    real    (kind=8), dimension(lwork)            :: rwork

    g0k=zero
    g0k_tmp=zero

    g0k=iw*eye(Nspin*Norb)-hk
    g0k_tmp=g0k

    call inv(g0k)
    call inversion_test(g0k,g0k_tmp,1.e-9)

  end function inverse_g0k



!_______________________________________________________________________
!                            reshape functions
!_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: reshape functions
  !  Z  = [Nspin,Nspin]*Norb
  !  A1 = [Norb*Norb]*Nspin
  !  A2 = [Nspin,Nspin,Norb,Norb]
  !---------------------------------------------------------------------
  function reshape_Z_to_A1(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !O-index
                   io1 = iorb + (ispin-1)*Norb
                   jo1 = jorb + (jspin-1)*Norb
                   !I-index
                   io2 = ispin + (iorb-1)*Nspin
                   jo2 = jspin + (jorb-1)*Nspin
                   !switch
                   g(io1,jo1)  = fg(io2,jo2)
                   !
                enddo
             enddo
          enddo
       enddo
  end function reshape_Z_to_A1
  function reshape_A1_to_Z(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !O-index
                   io1 = ispin + (iorb-1)*Nspin
                   jo1 = jspin + (jorb-1)*Nspin
                   !I-index
                   io2 = iorb + (ispin-1)*Norb
                   jo2 = jorb + (jspin-1)*Norb
                   !switch
                   g(io1,jo1)  = fg(io2,jo2)
                   !
                enddo
             enddo
          enddo
       enddo
  end function reshape_A1_to_Z
  function reshape_A1_to_A2(fg) result(g)
    complex(8),dimension(Nso,Nso)                   :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb)     :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin,io,jo
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   g(ispin,jspin,iorb,jorb)  = fg(io,jo)
                enddo
             enddo
          enddo
       enddo
  end function reshape_A1_to_A2

  !---------------------------------------------------------------------
  !PURPOSE: Inversion test
  !---------------------------------------------------------------------
  subroutine inversion_test(A,B,tol)
    implicit none
    complex (kind=8), intent(in)   ::   A(Nspin*Norb,Nspin*Norb)
    complex (kind=8), intent(in)   ::   B(Nspin*Norb,Nspin*Norb)
    real    (kind=4), intent(in)   ::   tol
    integer (kind=2)               ::   dime

    if (size(A).ne.size(B)) then
       write(*,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    if (abs(float(dime)-real(sum(matmul(A,B)))).gt.tol) write(*,'(A30)') "inversion test fail"
  end subroutine inversion_test


end program ed_STO

