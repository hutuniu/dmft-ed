program ed_tddpam_lattice
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                :: iloop
  logical                :: converged
  integer                :: Npd
  !Bath:
  integer                :: Nb
  real(8),allocatable    :: Bath(:)
  complex(8),allocatable :: Delta(:,:,:),Delta_Old(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  complex(8),allocatable :: pamHloc(:,:)
  real(8),allocatable    :: Wtk(:)
  complex(8),allocatable :: Smats(:,:,:)
  !variables for the model:
  character(len=32)      :: hkfile,finput
  integer                :: Nx,Lk,ntype,nkpath
  real(8)                :: nobj,wmixing
  real(8)                :: alpha,tpp,ep0,tpd,v0,gzero,gzerop,gzerom,gmu
  logical                :: bool
  real(8),dimension(2)   :: Eout

  !parse additional variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ntype,"NTYPE",finput,default=0)
  call parse_input_variable(Nx,"NX",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(tpp,"TPP",finput,default=0.25d0)
  call parse_input_variable(alpha,"ALPHA",finput,default=0.d0)
  call parse_input_variable(tpd,"TPD",finput,default=0.d0)
  call parse_input_variable(v0,"V0",finput,default=0.d0)
  call parse_input_variable(ep0,"EP0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call ed_read_input(trim(finput))

  !Number of orbitals:
  Npd=2

  if(wmixing==0.d0)stop "error: wmixing=0 is not allowed."

  inquire(file="last_mu.restart",exist=bool)
  if(bool.AND.nread/=0.d0)then
     open(100,file="last_mu.restart")
     read(100,*)xmu
     close(100)
     print*,"UPDATE XMU:",xmu
  endif

  !this shift contain |ep0-ed0|
  gmu=xmu
  gzerop=0.5d0*(ep0 + dsqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - dsqrt(ep0**2 + 4.d0*tpd**2))
  gzero=0.d0
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))
  allocate(delta_old(Norb,Norb,Lmats))

  !Read/Build the H(k)
  call build_hk("hkfile.in")

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call ed_init_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)

     !Get the Weiss field/Delta function to be fitted (user defined)
     delta_old=delta
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     if(iloop>1)delta = wmixing*delta + (1.d0-wmixing)*delta_old
     call ed_chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop)
     if(nread/=0.d0)call search_chemical_potential(xmu,nobj,converged)
     call end_loop
  enddo

  allocate(Smats(Npd,Npd,Lmats))
  Smats=zero
  Smats(1,1,:)=impSmats(1,1,1,1,:)
  Eout = ed_kinetic_energy(Hk,Wtk,Smats)
  deallocate(Smats)

  if(nread/=0.d0)then
     open(100,file='last_mu.restart')
     write(100,*)xmu-gzero
     close(100)
  endif



contains



  !BUILD H(k) from model Hamiltonian
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    complex(8),dimension(Npd,Npd,Lmats) :: Gmats
    complex(8),dimension(Npd,Npd,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    real(8),allocatable                 :: kxgrid(:),kygrid(:)
    integer,allocatable                 :: ik2ix(:),ik2iy(:)
    integer                             :: Npts
    real(8),dimension(:,:),allocatable  :: kpath

    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k):"
    !
    !SOLVE ALONG THE 2D PATH:
    Npts = 4
    Lk=(Npts-1)*Nkpath
    if(ED_MPI_ID==0)write(*,*)"# of k-points PATH    :",Lk
    allocate(kpath(Npts,2))
    kpath(1,:)=kpoint_Gamma
    kpath(2,:)=kpoint_M1
    kpath(3,:)=kpoint_X1
    kpath(4,:)=kpoint_Gamma
    if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_model,Npd,kpath,Lk,&
         colors_name=[character(len=20) :: 'red','blue'],&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file="Eigenband.nint")
    !
    !BUILD THE H(k) HAMILTONIAN
    Lk = Nx*Nx
    if(ED_MPI_ID==0)write(*,*)"# of k-points  BZ   :",Lk
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    if(allocated(kxgrid))deallocate(kxgrid)
    if(allocated(kygrid))deallocate(kygrid)
    !
    allocate(Hk(Npd,Npd,Lk))
    allocate(wtk(Lk))
    allocate(kxgrid(Nx),kygrid(Nx))
    kxgrid = kgrid(Nx)
    kygrid = kgrid(Nx)
    Hk     = build_hk_model(hk_model,Npd,kxgrid,kygrid,[0d0])
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90(file,Npd,&
            Nd=1,&
            Np=1,   &
            Nineq=1,&
            hk=Hk,  &
            kxgrid=kxgrid,&
            kygrid=kygrid,&
            kzgrid=[0d0])
    endif
    !
    !GET LOCAL HAM:
    allocate(pamHloc(Npd,Npd))
    pamHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(pamHloc))<1.d-9)pamHloc=0d0
    if(ED_MPI_ID==0)call write_Hloc(pamHloc)
    !Build the local GF:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu,Hk(:,:,ik))/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
       enddo
    enddo
    do iorb=1,Npd
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,&
            -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
    enddo
  end subroutine build_hk


  !+----------------------------------------+


  !MODEL HAMILTONIAN
  function Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    real(8)                   :: kx,ky,epsik,vpsik
    complex(8),dimension(N,N) :: Hk
    if(N/=Npd)stop "hk_model error: N != N_pd"
    kx=kvec(1)
    ky=kvec(2)
    epsik = cos(kx)+cos(ky)
    vpsik = sin(kx)*sin(ky)
    Hk(1,1) =     - 2.d0*alpha*tpp*epsik
    Hk(2,2) = ep0 - 2.d0*tpp*epsik
    Hk(1,2) = tpd - 4.d0*v0*vpsik
    Hk(2,1) = tpd - 4.d0*v0*vpsik
  end function Hk_model


  !+----------------------------------------+



  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8)                              :: iw,zeta(2),fg(Npd,Npd)
    complex(8),dimension(:,:,:),allocatable :: gloc
    real(8)                                 :: wm(Lmats),wr(Lreal),npimp,ntotal
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta=zero
    allocate(gloc(Npd,Npd,Lmats))
    do i=1,Lmats
       zeta(1)= xi*wm(i)+xmu  - impSmats(1,1,1,1,i)
       zeta(2)= xi*wm(i)+xmu
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zeta,Hk(:,:,ik))*Wtk(ik)
       enddo
       gloc(:,:,i)  = fg
       if(cg_scheme=='weiss')then
          delta(1,1,i) = one/(one/fg(1,1) + impSmats(1,1,1,1,i))
       else
          delta(1,1,i) = zeta(1) - one/fg(1,1)
       endif
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(1,1,:))
    do iorb=1,Npd
       do jorb=iorb,Npd
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
       enddo
    enddo
    npimp=2*fft_get_density(gloc(2,2,:),beta)
    ntotal=ed_dens(1)+npimp
    write(*,"(A,F25.18)")"np  =",npimp
    write(*,"(A,F25.18)")"ntot=",ntotal
    call splot("nd_np_ntot_all.ed",ed_dens(1),npimp,ntotal,append=.true.)
    call splot("nd_np_ntot.ed",ed_dens(1),npimp,ntotal)
    deallocate(gloc)


    allocate(gloc(Npd,Npd,Lreal))
    do i=1,Lreal
       zeta(1) = dcmplx(wr(i),eps)+xmu - impSreal(1,1,1,1,i)
       zeta(2) = dcmplx(wr(i),eps)+xmu
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zeta,Hk(:,:,ik))*Wtk(ik)
       enddo
       gloc(:,:,i)  = fg
    enddo
    !Print:
    do iorb=1,Npd
       do jorb=iorb,Npd
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi,dreal(gloc(iorb,jorb,:)))
       enddo
    enddo
    deallocate(gloc)

    if(ntype==1)then
       nobj=ed_dens(1)
    elseif(ntype==2)then
       nobj=npimp
    else
       nobj=ntotal
    endif

  end subroutine get_delta


  !+----------------------------------------+

  !INVERT 2x2 MATRICES ROUTINES:
  function inverse_g0k(iw,hk) result(g0k)
    integer                   :: i,M
    complex(8),dimension(2,2)    :: hk
    complex(8)                :: iw
    complex(8),dimension(2,2) :: g0k
    complex(8)                :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k

  function inverse_gk(zeta,hk) result(gk)
    integer                   :: i,M
    complex(8),dimension(2,2)    :: hk
    complex(8),dimension(2)   :: zeta
    complex(8),dimension(2,2) :: gk
    complex(8)                :: delta,ppi,vmix
    gk=zero
    delta = zeta(1) - hk(1,1)
    ppi   = zeta(2) - hk(2,2)
    vmix  = -hk(1,2)
    gk(1,1) = one/(delta - abs(vmix)**2/ppi)
    gk(2,2) = one/(ppi - abs(vmix)**2/delta)
    gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    gk(2,1) = conjg(gk(1,2))
  end function inverse_gk

end program ed_tddpam_lattice



