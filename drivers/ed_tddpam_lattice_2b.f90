program ed_tddpam_lattice
  USE DMFT_ED
  USE CONSTANTS
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  USE ERROR
  USE ARRAYS
  USE MATRIX
  USE IOTOOLS
  USE PARSE_INPUT
  implicit none
  integer                :: iloop
  logical                :: converged
  integer                :: Npd
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  complex(8),allocatable :: Delta(:,:,:),Delta_Old(:,:,:)
  !Hamiltonian input:
  real(8),allocatable    :: Hk(:,:,:),Hloc(:,:)
  real(8),allocatable    :: Wtk(:)
  complex(8),allocatable :: Smats(:,:,:)
  !variables for the model:
  character(len=32)      :: hkfile,finput
  integer                :: Nx,Lk,ntype
  real(8)                :: nobj,wmixing
  real(8)                :: alpha,tpp,ep0,tpd,v0,gzero,gzerop,gzerom,gmu
  logical                :: bool

  !parse additional variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ntype,"NTYPE",finput,default=0)
  call parse_input_variable(Nx,"NX",finput,default=100)
  call parse_input_variable(tpp,"TPP",finput,default=0.25d0)
  call parse_input_variable(alpha,"ALPHA",finput,default=0.d0)
  call parse_input_variable(tpd,"TPD",finput,default=0.d0)
  call parse_input_variable(v0,"V0",finput,default=0.d0)
  call parse_input_variable(ep0,"EP0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call ed_read_input(trim(finput))

  !Number of orbitals:
  if(Norb/=2)stop "This driver requires Norb==2"
  Npd=Norb

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
  gzerop=0.5d0*(ep0 + sqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - sqrt(ep0**2 + 4.d0*tpd**2))
  gzero=0.d0
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))
  allocate(delta_old(Norb,Norb,Lmats))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

  !Read/Build the H(k)
  call build_hk()
  call set_Hloc(dcmplx(Hloc),1)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath)

     !Get the Weiss field/Delta function to be fitted (user defined)
     delta_old=delta
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     if(iloop>1)delta = wmixing*delta + (1.d0-wmixing)*delta_old
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop)
     if(nread/=0.d0)call search_chemical_potential(xmu,nobj,converged)
     call end_loop
  enddo

  allocate(Smats(Npd,Npd,Lmats))
  Smats=zero
  Smats=impSmats(1,1,:,:,:)
  call ed_kinetic_energy(Smats,Hk,wtk)
  deallocate(Smats)

  if(nread/=0.d0)then
     open(100,file='last_mu.restart')
     write(100,*)xmu-gzero
     close(100)
  endif


contains


  !MODEL HAMILTONIAN
  function Hk_model(kx,ky) result(Hk)
    real(8)                    :: kx,ky,epsik,vpsik
    real(8),dimension(Npd,Npd) :: Hk
    epsik = cos(kx)+cos(ky)
    vpsik = sin(kx)*sin(ky)
    Hk(1,1) =     - 2.d0*alpha*tpp*epsik
    Hk(2,2) = ep0 - 2.d0*tpp*epsik
    Hk(1,2) = tpd - 4.d0*v0*vpsik
    Hk(2,1) = tpd - 4.d0*v0*vpsik
  end function Hk_model


  !BUILD H(k) from model Hamiltonian
  subroutine build_hk()
    integer                :: ix,iy,i,j,ik,iorb,jorb,unit,units(2)
    real(8)                :: kx,ky
    real(8),allocatable    :: pEloc(:)
    real(8),allocatable    :: pHk(:,:)
    Lk=Nx*Nx
    allocate(Hk(Npd,Npd,Lk),pHk(Npd,Npd),peloc(Npd),Hloc(Npd,Npd))
    allocate(wtk(Lk))
    units=free_units(2)
    open(units(1),file="Eigenband_l1.ed")
    open(units(2),file="Eigenband_l2.ed")
    ik=0
    do ix=1,Nx
       kx=-pi + 2.d0*pi*dble(ix-1)/dble(nx)
       do iy=1,Nx
          ky=-pi + 2.d0*pi*dble(iy-1)/dble(Nx)
          ik=ik+1
          Hk(:,:,ik)=Hk_model(kx,ky)
       enddo
    enddo
    wtk = 1.d0/dble(Lk)
    Hloc = sum(Hk,dim=3)/dble(Lk)
    where(abs(Hloc)<1.d-9)Hloc=0.d0
    !
    ik=0
    do ix=1,Lk
       ik=ik+1
       kx = 0.d0 + pi*real(ix-1,8)/dble(Lk)
       ky = 0.d0
       pHk=Hk_model(kx,ky) 
       call matrix_diagonalize(pHk,peloc,'V')
       write(units(1),"(I,100F18.12)")ik,peloc(1),(pHk(i,1)**2,i=1,size(pHk(:,1)))
       write(units(2),"(I,100F18.12)")ik,peloc(2),(pHk(i,2)**2,i=1,size(pHk(:,2)))
    enddo
    !From X=(pi,0) to M=(pi,pi): 100 steps
    do iy=1,Lk
       ik=ik+1
       kx = pi
       ky = 0.d0 + pi*real(iy-1,8)/dble(Lk)
       pHk=Hk_model(kx,ky) 
       call matrix_diagonalize(pHk,peloc,'V')
       write(units(1),"(I,100F18.12)")ik,peloc(1),(pHk(i,1)**2,i=1,size(pHk(:,1)))
       write(units(2),"(I,100F18.12)")ik,peloc(2),(pHk(i,2)**2,i=1,size(pHk(:,2)))
    enddo
    !From M=(pi,pi) to \Gamma=(0,0): 100 steps
    do ix=1,Lk
       ik=ik+1
       iy=ix
       kx = pi - pi*real(ix-1,8)/dble(Lk)
       ky = pi - pi*real(iy-1,8)/dble(Lk)
       pHk=Hk_model(kx,ky) 
       call matrix_diagonalize(pHk,peloc,'V')
       write(units(1),"(I,100F18.12)")ik,peloc(1),(pHk(i,1)**2,i=1,size(pHk(:,1)))
       write(units(2),"(I,100F18.12)")ik,peloc(2),(pHk(i,2)**2,i=1,size(pHk(:,2)))
    enddo
  end subroutine build_hk


  !+----------------------------------------+



  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8)                              :: iw,zeta(Npd,Npd),fg(Npd,Npd),gweiss(Npd,Npd)
    complex(8),dimension(:,:,:),allocatable :: gloc
    real(8)                                 :: wm(Lmats),wr(Lreal),npimp,ntotal
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta=zero
    allocate(gloc(Npd,Npd,Lmats))
    do i=1,Lmats
       zeta=zero
       forall(iorb=1:Norb)zeta(iorb,iorb) = (xi*wm(i) + xmu)
       zeta(:,:) = zeta(:,:) - impSmats(1,1,:,:,i)
       !
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zeta,Hk(:,:,ik))*Wtk(ik)
       enddo
       gloc(:,:,i)  = fg
       call matrix_inverse(fg)
       if(cg_scheme=='weiss')then
          gweiss = fg + impSmats(1,1,:,:,i)
          call matrix_inverse(gweiss)
          delta(:,:,i) = gweiss
       else
          delta(:,:,i) = zeta(:,:) - Hloc(:,:) - fg(:,:)
       endif
    enddo
    !Print:
    do iorb=1,Npd
       do jorb=iorb,Npd
          call splot("Delta_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,Delta(iorb,jorb,:))
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
       enddo
    enddo
    deallocate(gloc)


    allocate(gloc(Npd,Npd,Lreal))
    do i=1,Lreal
       zeta=zero
       forall(iorb=1:Norb)zeta(iorb,iorb) = dcmplx(wr(i),eps)+xmu
       zeta(:,:) = zeta(:,:) - impSreal(1,1,:,:,i)
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
       nobj=ed_dens(1)
    else
       nobj=sum(ed_dens)
    endif

  end subroutine get_delta


  !+----------------------------------------+






  !INVERT 2x2 MATRICES ROUTINES:
  function inverse_g0k(iw,hk) result(g0k)
    integer                   :: i,M
    real(8),dimension(2,2)    :: hk
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
    real(8),dimension(2,2)    :: hk
    complex(8),dimension(2,2) :: zeta
    complex(8),dimension(2,2) :: gk,invg
    complex(8)                :: delta,ppi,vmix
    ! gk=zero
    ! delta = zeta(1,1) - hk(1,1)
    ! ppi   = zeta(2,2) - hk(2,2)
    ! vmix  = zeta(1,2) - hk(1,2)
    ! gk(1,1) = one/(delta - abs(vmix)**2/ppi)
    ! gk(2,2) = one/(ppi - abs(vmix)**2/delta)
    ! gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    ! gk(2,1) = conjg(gk(1,2))
    gk = zeta-hk
    call matrix_inverse(gk)
  end function inverse_gk


end program ed_tddpam_lattice



