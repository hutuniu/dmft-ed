
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
  complex(8),allocatable :: Hk(:,:,:)
  complex(8),allocatable :: Ti3dt2g_Hloc(:,:)
  complex(8),allocatable :: Ti3dt2g_Hloc_nn(:,:,:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin
  real(8)                :: soc,ivb,wmixing,sumdens
  logical                :: surface,Hk_test,rotateG0loc
  character(len=16)      :: finput
  character(len=32)      :: hkfile
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
  call parse_input_variable(soc,    "SOC",finput,       default=0.25d0)
  call parse_input_variable(ivb,    "IVB",finput,       default=0.02d0)
  call parse_input_variable(surface,"SURFACE",finput,   default=.false.)
  call parse_input_variable(Hk_test,"HK_TEST",finput,   default=.false.)
  call parse_input_variable(rotateG0loc,"ROTATEG0loc",finput, default=.false.)
  call ed_read_input(trim(finput))
  !
  Nso=Nspin*Norb
  !
  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Ti3dt2g_Hloc_nn(Nspin,Nspin,Norb,Norb))
  allocate(density_matrix(Nspin*Norb,Nspin*Norb))
  allocate(dm_eig(Nspin*Norb),dm_rot(Nspin*Norb,Nspin*Norb))
  allocate(delta_conv(Nso,Nso,Lmats))
  allocate(delta_conv_avrg(Lmats))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  call build_hk_path
  Ti3dt2g_Hloc_nn=reshape_A1_to_A2(Ti3dt2g_Hloc)

  !Setup solver impHloc non è ne' allocata ne' riempita
  Nb=get_bath_size()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))


  !old
  !1) alloco impHloc e la riempo se trovo il file
  call ed_init_solver(bath)
  !2) riempio impHloc
  call set_hloc(reshape_A1_to_A2(Ti3dt2g_Hloc))

  !new
  !1) alloco impHloc e la riempo se trovo il file, e la riempo con argomento opzionale, se non è presente tiene quella che ha letto dal file o a zero e posso cmq settarla fuori come prima
  call ed_init_solver(bath,Ti3dt2g_Hloc_nn)

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
     call mpi_barrier(MPI_COMM_WORLD,ED_MPI_ERR)
     if(ED_MPI_ID==0)call rotate_Gloc(Greal)
     if(ED_MPI_ID==1)call Quantum_operator()
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=reshape_A1_to_A2(Ti3dt2g_Hloc),iprint=3)
     !density matrix
     if(ED_MPI_ID==0)then
        call ed_get_density_matrix(density_matrix,2,dm_eig,dm_rot)
     endif
     !Fit the new bath, starting from the old bath + the supplied delta 
     Bath_=bath
     if (ed_mode=="normal") then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        call spin_symmetrize_bath(bath,save=.false.)
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
    complex(8),dimension(6,6)           :: theta,impHloc_rot
    real(8)                             :: wm(Lmats),wr(Lreal),dw

    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for STO:"
    !
    !Alloco quantità
    !
    if(surface) then
       Lk=Nk**2
    else
       Lk=Nk**3
    endif
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nso,Nso,Lk));allocate(wtk(Lk));allocate(kxgrid(Nk),kygrid(Nk),kzgrid(Nk))
    wtk = 1.0d0/Lk
    kxgrid=0.0d0;kygrid=0.0d0;kzgrid=0.0d0
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    if(.not.surface) kzgrid = kgrid(Nk)
    !
    !scrivo H(k)
    !
    if(surface) then
       Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,[0d0])
       if(ED_MPI_ID==0.AND.present(file)) call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=[0d0] )
    else
       Hk = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,kzgrid)
       if(ED_MPI_ID==0.AND.present(file)) call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=kzgrid)
    endif
    !
    !calcolo impHloc = Sum_k [ H(k) ]
    !
    allocate(Ti3dt2g_Hloc(Nso,Nso))
    Ti3dt2g_Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=0d0
    !
    !scrivo impHloc
    !
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
    !
    !Build the local GF in the spin-orbital Basis:
    !
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
    open(unit=106,file='sum_w_G0loc.dat',status='unknown',action='write',position='rewind')
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(Greal(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    !Build the local GF in the J-Jz Basis (theta is already in the A1 shape)
    !
    if((rotateG0loc).and.(.not.Hk_test)) then
       call build_rotation(theta,impHloc_rot)
       do i=1,Lreal
          Greal(:,:,i)=matmul(transpose(conjg(theta)),matmul(Greal(:,:,i),theta))
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot("G0locrot_a_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
       open(unit=106,file='sum_w_G0loc_rot_anlytc.dat',status='unknown',action='write',position='rewind')
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
       close(106)
       !
       do i=1,Lreal
          Greal(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(Greal(:,:,i),impHloc_rot))
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   call splot("G0locrot_H_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
       open(unit=106,file='sum_w_G0loc_rot_Hloc.dat',status='unknown',action='write',position='rewind')
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
       close(106)
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
  !PURPOSE: 2x2 band structures
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
            file="Eigenband_bulk_Hartree.nint")


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
            file="Eigenband_surf_Hartree.nint")
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
  subroutine build_rotation(theta_,impHloc_rot_)
    complex(8),dimension(6,6),intent(out)           :: theta_
    complex(8),dimension(6,6),optional,intent(out)  :: impHloc_rot_
    real(8),dimension(6)                            :: impHloc_eig
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

    if(present(impHloc_rot_)) then
       impHloc_rot_=zero
       impHloc_rot_=Ti3dt2g_Hloc
       call matrix_diagonalize(impHloc_rot_,impHloc_eig,'V','U')
    endif

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

  subroutine rotate_Gloc(Gsowr)
    implicit none
    complex(8),allocatable,intent(in)             ::   Gsowr(:,:,:,:,:)
    complex(8),allocatable                        ::   G_in(:,:,:),G_out(:,:,:)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   ::   theta,impHloc_rot
    integer                                       ::   io,jo
    integer                                       ::   ispin,jspin
    integer                                       ::   iorb,jorb
    real(8)                                       ::   wr(Lreal),dw
    !
    write(*,*) "A(w) rotation"
    !
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    call build_rotation(theta,impHloc_rot)
    if(allocated( G_in))deallocate( G_in);allocate( G_in(Nspin*Norb,Nspin*Norb,Lreal));G_in=zero
    if(allocated(G_out))deallocate(G_out);allocate(G_out(Nspin*Norb,Nspin*Norb,Lreal));G_out=zero
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                G_in(io,jo,:)=Gsowr(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_Gloc.dat',status='unknown',action='write',position='rewind')
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_in(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(G_in(:,:,i),impHloc_rot))
    enddo
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                call splot("Glocrot_H_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_Gloc_rot_Hloc.dat',status='unknown',action='write',position='rewind')
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_in(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
    G_out=zero
    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(theta)),matmul(G_in(:,:,i),theta))
    enddo
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                call splot("Glocrot_a_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    open(unit=106,file='sum_w_Gloc_rot_anlytc.dat',status='unknown',action='write',position='rewind')
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                write(106,*) io,jo,"---",ispin,jspin,iorb,jorb,sum(abs(G_in(io,jo,:)))
             enddo
          enddo
       enddo
    enddo
    close(106)
    !
  end subroutine rotate_Gloc

  !---------------------------------------------------------------------
  !PURPOSE: 
  !---------------------------------------------------------------------
  subroutine Quantum_operator()
    implicit none
    complex(8),allocatable             :: Gso(:,:,:,:,:)
    complex(8),allocatable             :: Stot(:,:,:),Ltot(:,:,:)
    complex(8)                         :: LdotS
    complex(8)                         :: Sx,Lx,Sy,Ly,Sz,Lz,jz
    integer                            :: ilat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    real(8)                            :: site_mag(Norb)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate( Gso(Nspin,Nspin,Norb,Norb,Lmats)); Gso=zero
    call ed_get_gimp_matsubara(Gso)
    !
    !##############################################################
    !
    !                              S
    !
    !##############################################################
    !
    write(*,*) "Computing total Spin operator per orbital"
    write(*,*) "Lmats used:",Lmats
    allocate(Stot(3,Norb,Norb));Stot=zero
    !
    do iorb=1,Norb
       do jorb=1,Norb
          !Sx
          Stot(1,iorb,jorb)=sum(    (Gso(1,2,iorb,jorb,:)+Gso(2,1,iorb,jorb,:) ))/beta
          !Sy
          Stot(2,iorb,jorb)=sum( xi*(Gso(2,1,iorb,jorb,:)-Gso(1,2,iorb,jorb,:) ))/beta
          !Sz
          Stot(3,iorb,jorb)=sum(    (Gso(1,1,iorb,jorb,:)-Gso(2,2,iorb,jorb,:) ))/beta
       enddo
    enddo
    !
    site_mag=0.d0
    site_mag=ed_get_mag()
    !
    open(unit=105,file='Stot.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(105,'(a100)') "#diagonal orbital"
    write(105,'(30a20)') "Re{Sx}_11","Re{Sx}_22","Re{Sx}_33" &
                        ,"Re{Sy}_11","Re{Sy}_22","Re{Sy}_33" &
                        ,"Re{Sz}_11","Re{Sz}_22","Re{Sz}_33" &
                        ,"Im{Sx}_11","Im{Sx}_22","Im{Sx}_33" &
                        ,"Im{Sy}_11","Im{Sy}_22","Im{Sy}_33" &
                        ,"Im{Sz}_11","Im{Sz}_22","Im{Sz}_33" &
                        ,"mag" 
    write(105,'(30F20.12)') real(Stot(1,1,1)), real(Stot(1,2,2)), real(Stot(1,3,3)) &
                          , real(Stot(2,1,1)), real(Stot(2,2,2)), real(Stot(2,3,3)) &
                          , real(Stot(3,1,1)), real(Stot(3,2,2)), real(Stot(3,3,3)) &
                          ,aimag(Stot(1,1,1)),aimag(Stot(1,2,2)),aimag(Stot(1,3,3)) &       
                          ,aimag(Stot(2,1,1)),aimag(Stot(2,2,2)),aimag(Stot(2,3,3)) &
                          ,aimag(Stot(3,1,1)),aimag(Stot(3,2,2)),aimag(Stot(3,3,3)) &
                          ,site_mag(1)/2.,site_mag(2)/2.,site_mag(3)/2.
    write(105,*)
    write(105,*)
    write(105,'(a100)') "#inter-orbital"
    write(105,'(30a20)') "#Sx(orb_1)","Sx(orb_2)","Sx(orb_3)","Sy(orb_1)","Sy(orb_2)","Sy(orb_3)","Sz(orb_1)","Sz(orb_2)","Sz(orb_3)"
    do iorb=1,Norb
       write(105,'(30F20.12)') (real(Stot(1,iorb,jorb)),jorb=1,Norb) &
                              ,(real(Stot(2,iorb,jorb)),jorb=1,Norb) &
                              ,(real(Stot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(105,*)
    do iorb=1,Norb
       write(105,'(30F20.12)') (aimag(Stot(1,iorb,jorb)),jorb=1,Norb) &
                              ,(aimag(Stot(2,iorb,jorb)),jorb=1,Norb) &
                              ,(aimag(Stot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    close(105)
    !
    !##############################################################
    !
    !                              L
    !
    !##############################################################
    !
    write(*,*) "Computing total Orbital operator per spin"
    write(*,*) "Lmats used:",Lmats
    allocate(Ltot(3,Nspin,Nspin));Ltot=zero
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          !Lx
          Ltot(1,ispin,jspin)=sum(  xi*(Gso(ispin,jspin,1,3,:)-Gso(ispin,jspin,3,1,:))  )/beta
          !Ly
          Ltot(2,ispin,jspin)=sum(  xi*(Gso(ispin,jspin,3,2,:)-Gso(ispin,jspin,2,3,:))  )/beta
          !Lz
          Ltot(3,ispin,jspin)=sum(  xi*(Gso(ispin,jspin,1,2,:)-Gso(ispin,jspin,2,1,:))  )/beta
       enddo
    enddo
    !
    open(unit=106,file='Ltot.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(106,'(a100)') "#diagonal spin"
    write(106,'(a8,30a20)') "Re{Lx}_1","Re{Lx}_2" &
                           ,"Re{Ly}_1","Re{Ly}_2" &
                           ,"Re{Lz}_1","Re{Lz}_2" &
                           ,"Im{Lx}_1","Im{Lx}_2" &
                           ,"Im{Ly}_1","Im{Ly}_2" &
                           ,"Im{Lz}_1","Im{Lz}_2"
    write(106,'(30F20.12)') real(Ltot(1,1,1)), real(Ltot(1,2,2)) &
                         ,  real(Ltot(2,1,1)), real(Ltot(2,2,2)) &
                         ,  real(Ltot(3,1,1)), real(Ltot(3,2,2)) &
                         , aimag(Ltot(1,1,1)),aimag(Ltot(1,2,2)) &       
                         , aimag(Ltot(2,1,1)),aimag(Ltot(2,2,2)) &
                         , aimag(Ltot(3,1,1)),aimag(Ltot(3,2,2))
    write(106,*)
    write(106,*)
    write(106,'(a100)') "#inter-spin"
    write(106,'(30a20)') "#Lx(spin_1)","Lx(spin_2)","Ly(spin_1)","Ly(spin_2)","Lz(spin_1)","Lz(spin_2)"
    do ispin=1,Nspin
       write(106,'(30F20.12)') (real(Ltot(1,ispin,jspin)),jspin=1,Nspin) &
                              ,(real(Ltot(2,ispin,jspin)),jspin=1,Nspin) &
                              ,(real(Ltot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(106,*)
    do ispin=1,Nspin
       write(106,'(30F20.12)') (aimag(Ltot(1,ispin,jspin)),jspin=1,Nspin) &
                              ,(aimag(Ltot(2,ispin,jspin)),jspin=1,Nspin) &
                              ,(aimag(Ltot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    close(106)
    !
    !##############################################################
    !
    !                              L.dot.S
    !
    !##############################################################
    !
    write(*,*) "Computing total L dot S operator per site"
    write(*,*) "Lmats used:",Lmats
    LdotS=zero
    !
    LdotS=sum(       +xi*Gso(1,1,1,2,:) &
                     +xi*Gso(1,2,1,3,:) &  
                     -xi*Gso(2,2,1,2,:) &  
                     +xi*Gso(2,1,1,3,:) &  
                     -xi*Gso(1,1,2,1,:) &  
                     -   Gso(1,2,2,3,:) &  
                     +xi*Gso(2,2,2,1,:) &  
                     +   Gso(2,1,2,3,:) &  
                     -xi*Gso(1,2,3,1,:) &  
                     +   Gso(1,2,3,2,:) &  
                     -xi*Gso(2,1,3,1,:) &  
                     -   Gso(2,1,3,2,:) &  
                      )/beta
    LdotS=LdotS/2.d0
    !
    open(unit=107,file='Jz.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(107,'(30a20)') "Re{Sz}_11","Re{Sz}_22","Re{Sz}_33","Im{Sz}_11","Im{Sz}_22","Im{Sz}_33","Re{Tr[Sz]}","Im{Tr[Sz]}" &
                                                ,"Re{Lz}_uu","Im{Lz}_uu","Re{Lz}_dd","Im{Lz}_dd","Re{Tr[Lz]}","Im{Tr[Lz]}" &
                                                ,"Re{jz}","Im{jz}","Re{L.S}","Im{L.S}"
    Sx=trace(Stot(1,:,:));Sy=trace(Stot(2,:,:));Sz=trace(Stot(3,:,:))
    Lx=trace(Ltot(1,:,:));Ly=trace(Ltot(2,:,:));Lz=trace(Ltot(3,:,:))
    jz=Sz+Lz
    write(107,'(30F20.12)') real(Stot(3,1,1)), real(Stot(3,2,2)), real(Stot(3,3,3)) &
                         , aimag(Stot(3,1,1)),aimag(Stot(3,2,2)),aimag(Stot(3,3,3)),real(Sz),aimag(Sz) &
                         ,  real(Ltot(3,1,1)), real(Ltot(3,2,2)) &
                         , aimag(Ltot(3,1,1)),aimag(Ltot(3,2,2)),real(Lz),aimag(Lz) &
                         ,  real(jz),aimag(jz),real(LdotS),aimag(LdotS)
    close(107)
    !
    deallocate(Ltot,Stot)
    !
  end subroutine Quantum_operator





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

