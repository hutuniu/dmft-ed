!                    MODEL Hamiltonian is:
!
! |     h^{2x2}(k)             &          hso^{2x2}(k)        |
! |     [hso^{2x2}]*(k)        &         [h^{2x2}]*(-k)       |
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
! hso^{2x2}(k):=
! | xi*rh*(sin(kx)-xi*sin(ky))  &         \delta              |
! |         -\delta             &             0               |
program ed_bhz
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
  integer                :: Nb
  real(8),allocatable    :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                :: Nk,Nkpath
  real(8)                :: mh,lambda,wmixing,akrange,rh
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym,getak,getdeltaw,getpoles
  type(finter_type)      :: finter_func
  !
  real(8),dimension(2)   :: Eout
  real(8),allocatable :: dens(:)
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(getak,"GETAK",finput,default=.false.)
  call parse_input_variable(getdeltaw,"GETDELTAW",finput,default=.false.)
  call parse_input_variable(getpoles,"GETPOLES",finput,default=.false.)
  call parse_input_variable(akrange,"AKRANGE",finput,default=5.d0)
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(rh,"RH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  !
  call ed_read_input(trim(finput))

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(dens(Norb))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))


  !OPTIONAL CHANNELS:
  if(getak)then
     call get_Akw
     stop
  endif
  !
  if(getpoles)then
     call get_poles
     stop
  endif
  !
  if(getdeltaw)then
     call get_deltaw
     stop
  endif


  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(bath)
  call set_hloc(j2so(bhzHloc))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     dens= ed_get_dens()
     print*,dens

     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=j2so(bhzHloc),iprint=1)
     ! !Get the Weiss field/Delta function to be fitted (user defined)
     ! call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call ed_chi2_fitgf(delta(1,1,:,:,:),bath,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf(delta(2,2,:,:,:),bath,ispin=2)
     else
        call spin_symmetrize_bath(bath,save=.true.)
     endif

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     if(ED_MPI_ID==0)then
        converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
        if(nread/=0.d0)call search_chemical_potential(xmu,sum(dens),converged)
     endif
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0)call end_loop
  enddo

  Eout = ed_kinetic_energy(Hk,Wtk,Smats)
  print*,Eout

#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif
contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw,n0(Nso)
    call build_hk_GXMG()

    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**2
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    allocate(kxgrid(Nk),kygrid(Nk))
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    Hk     = build_hk_model(hk_bhz,Nso,kxgrid,kygrid,[0d0])
    wtk = 1d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90(file,Nso,&
            Nd=Norb,&
            Np=1,   &
            Nineq=1,&
            hk=Hk,  &
            kxgrid=kxgrid,&
            kygrid=kxgrid,&
            kzgrid=[0d0])
    endif
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0d0
    if(ED_MPI_ID==0)call write_Hloc(bhzHloc)
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
    do iorb=1,Nso
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,&
            -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
       n0(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
    enddo
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG THE Gamma-X-M-Gamma path
  !---------------------------------------------------------------------
  subroutine build_hk_GXMG(kpath_)
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    character(len=64)                      :: file
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    if(present(kpath_))then
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eig_path.nint"
    else
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_Gamma
       kpath(2,:)=kpoint_M1
       kpath(3,:)=kpoint_X1
       kpath(4,:)=kpoint_Gamma
       file="Eigenbands.nint"
    endif
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    Hk     = build_hk_model(hk_bhz,Nso,kpath,Nkpath)
    wtk = 1d0/Lk
    if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_bhz,Nso,kpath,Lk,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine build_hk_GXMG







  !---------------------------------------------------------------------
  !PURPOSE: GET DELTA FUNCTION
  !---------------------------------------------------------------------
  subroutine get_delta
    integer                                     :: i,j,ik,iorb,jorb,ispin,jspin,iso,unit
    complex(8),dimension(Nso,Nso)               :: zeta,fg,gdelta
    complex(8),dimension(:,:,:),allocatable     :: Smats
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc
    complex(8),dimension(:,:,:,:),allocatable   :: gk
    complex(8)                                  :: iw
    real(8)                                     :: wm(Lmats),wr(Lreal)
    character(len=20)                           :: suffix
    !
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    delta = zero
    !
    !MATSUBARA AXIS
    if(ED_MPI_ID==0)print*,"Get Gloc_iw:"
    allocate(gloc(Nspin,Nspin,Norb,Norb,Lmats))
    do i=1,Lmats
       iw = xi*wm(i)
       forall(iorb=1:Nso)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - so2j(impSmats(:,:,:,:,i),Nso)
       fg=zero
       if(lambda==0.d0)then
          do ik=1,Lk
             forall(iorb=1:Nso)&
                  fg(iorb,iorb) = fg(iorb,iorb) + Wtk(ik)/(zeta(iorb,iorb)-Hk(iorb,iorb,ik))
          enddo
          gloc(:,:,:,:,i) = j2so(fg)
          !Get Delta=\Delta or G_0
          gdelta=zero
          forall(iorb=1:Nso)fg(iorb,iorb)=one/fg(iorb,iorb)
          if(cg_scheme=='weiss')then
             gdelta = fg + so2j(impSmats(:,:,:,:,i),Nso)
             forall(iorb=1:Nso)gdelta(iorb,iorb)=one/gdelta(iorb,iorb)
          else
             forall(iorb=1:Nso)gdelta(iorb,iorb) = zeta(iorb,iorb) - bhzHloc(iorb,iorb) - fg(iorb,iorb)
          endif
       else
          do ik=1,Lk
             fg = fg + inverse_gk(zeta,Hk(:,:,ik))*Wtk(ik)
          enddo
          gloc(:,:,:,:,i) = j2so(fg)
          !Get Delta=\Delta or G_0
          call matrix_inverse(fg)
          if(cg_scheme=='weiss')then
             gdelta = fg + so2j(impSmats(:,:,:,:,i),Nso)
             call matrix_inverse(gdelta)
          else
             gdelta = zeta(:,:) - bhzHloc - fg(:,:)
          endif
       endif
       !
       delta(:,:,:,:,i) = j2so(gdelta(:,:))
       !
    enddo
    if(ED_MPI_ID==0)then
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call splot("Delta"//reg(suffix),wm,delta(ispin,ispin,iorb,iorb,:))
             call splot("Gloc"//reg(suffix),wm,gloc(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
    endif
    deallocate(gloc)
    !
    !REAL AXIS
    allocate(gloc(Nspin,Nspin,Norb,Norb,Lreal))
    if(ED_MPI_ID==0)print*,"Get Gloc_realw:"
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       forall(iorb=1:Nso)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - so2j(impSreal(:,:,:,:,i),Nso)
       fg=zero
       do ik=1,Lk         
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*wtk(ik)
       enddo
       gloc(:,:,:,:,i) = j2so(fg)
    enddo
    if(ED_MPI_ID==0)then
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call splot("Gloc"//reg(suffix),wr,-dimag(gloc(ispin,ispin,iorb,iorb,:))/pi,dreal(gloc(ispin,ispin,iorb,iorb,:)))
          enddo
       enddo
    endif
    deallocate(gloc)
    ! !Get Kinetic Energy too
    ! allocate(Smats(Nso,Nso,Lmats))
    ! do i=1,Lmats
    !    Smats(:,:,i)=so2j(impSmats(:,:,:,:,i),Nso)
    ! enddo
    ! call ed_kinetic_energy(Hk,wtk,Smats)
    ! deallocate(Smats)
  end subroutine get_delta





  !---------------------------------------------------------------------
  !PURPOSE: GET A(k,w)
  !---------------------------------------------------------------------
  subroutine get_Akw()
    integer,parameter                           :: Lw=250
    integer                                     :: i,j,ik=0
    integer                                     :: ix,iy
    integer                                     :: iorb,jorb
    integer                                     :: isporb,jsporb
    integer                                     :: ispin,jspin
    integer                                     :: iso,unit
    real(8)                                     :: foo
    real(8)                                     :: kx,ky    
    complex(8),dimension(Nso,Nso)               :: zeta,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal_,Sreal
    complex(8),dimension(:,:,:,:),allocatable   :: gk,gfoo,ReSmat
    complex(8),dimension(:,:,:),allocatable     :: Hktilde
    complex(8)                                  :: iw
    real(8)                                     :: wr_(Lreal),wr(Lw)
    real(8),dimension(:,:),allocatable          :: Kpath
    character(len=20)                           :: suffix
    if(ED_MPI_ID==0)then
       !
       wr_ = linspace(wini,wfin,Lreal)
       wr  = linspace(-akrange,akrange,Lw)
       !
       call build_hk_GXMG()
       ! !Uncomment these lines and comment the previous one to use a different path
       ! !attention: the path should be edited. now it's the same as GXMG
       ! allocate(kpath(4,3))
       ! kpath(1,:)=kpoint_Gamma
       ! kpath(2,:)=kpoint_M1
       ! kpath(3,:)=kpoint_X1
       ! kpath(4,:)=kpoint_Gamma
       ! call build_hk_GXMG(kpath)

       print*,"Get A(k,w):",Lk
       !
       allocate(Sreal_(Nspin,Nspin,Norb,Norb,Lreal))
       call read_sigma(Sreal_)
       allocate(Sreal(Nspin,Nspin,Norb,Norb,Lw))
       Sreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             call cubic_spline(wr_,Sreal_(ispin,ispin,iorb,iorb,:),wr,Sreal(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       allocate(gk(Lk,Nspin,Norb,Lw))
       allocate(gfoo(Nspin,Nspin,Norb,Norb))
       call start_progress(LOGfile)
       do i=1,Lw
          iw=dcmplx(wr(i),eps)
          zeta=zero
          forall(iso=1:Nso)zeta(iso,iso)=iw+xmu
          zeta(:,:) = zeta(:,:)-so2j(Sreal(:,:,:,:,i),Nso)
          do ik=1,Lk
             fgk = zeta - Hk(:,:,ik)!inverse_gk(zeta,Hk(:,:,ik))
             call inv(fgk)
             gfoo(:,:,:,:) = j2so(fgk(:,:))
             forall(ispin=1:Nspin,iorb=1:Norb)gk(ik,ispin,iorb,i) = gfoo(ispin,ispin,iorb,iorb)
          enddo
          call progress(i,Lw)
       enddo
       call stop_progress()
       !PRINT
       do ispin=1,1
          do iorb=1,Norb
             unit=free_unit()
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             print*,"printing: Ak"//reg(suffix)
             call splot3d("Ak"//reg(suffix),(/(dble(ik),ik=1,Lk)/),wr,-dimag(gk(:,ispin,iorb,:))/pi,ymin=-akrange,ymax=akrange,nosurface=.true.)
          enddo
       enddo
       deallocate(gk,gfoo)
    endif
    return
  end subroutine get_Akw





  !---------------------------------------------------------------------
  !PURPOSE: GET POLES ON THE REAL AXIS
  !---------------------------------------------------------------------
  subroutine get_poles
    USE IOFILE
    integer                                     :: i,j,ik,ix,iy
    integer                                     :: iorb,jorb
    integer                                     :: isporb,jsporb
    integer                                     :: ispin,jspin
    integer                                     :: iso,unit
    real(8),dimension(Nso)                      :: dzeta
    complex(8),dimension(Nso,Nso)               :: zeta,gfk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal,Smats
    complex(8),dimension(:,:,:,:),allocatable   :: gk,gfoo,ReSmat
    complex(8)                                  :: iw
    complex(8),dimension(:,:),allocatable       :: detGiw
    real(8)                                     :: wr(Lreal),wm(Lmats)
    real(8),dimension(Lreal)                    :: Den
    real(8),dimension(:),allocatable            :: Ipoles,Xcsign,Iweight
    real(8),dimension(:,:),allocatable          :: Ipoles3d
    real(8),dimension(:,:),allocatable          :: Mpoles,Mweight
    real(8),dimension(:,:,:),allocatable        :: Mpoles3d
    integer                                     :: Linterval
    integer                                     :: count,Ninterval,maxNinterval,int
    real(8)                                     :: sign,sign_old
    wr = linspace(wini,wfin,Lreal)
    wm = pi/beta*(2*arange(1,Lmats)-1)
    allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
    call read_sigma(Sreal)
    call read_sigma(Smats)
    !
    call build_hk_GXMG
    ! allocate(detGiw(Lk,Lmats))
    ! ! unit=free_unit()
    ! ! open(unit,file="detGk_iw.ed")
    ! do ik=1,Lk
    !    do i=1,Lmats
    !       forall(iorb=1:Nso)zeta(iorb,iorb)=xi*wm(i)+xmu
    !       zeta(:,:)    = zeta(:,:) - (so2j(Smats(:,:,:,:,i),Nso))
    !       detGiw(ik,i) = one/( (zeta(1,1) - Hk(1,1,ik))*(zeta(2,2) - Hk(2,2,ik)) - Hk(1,2,ik)*Hk(2,1,ik))
    !       write(unit,*)wm(i),dimag(detGiw(ik,i)),dreal(detGiw(ik,i))
    !    enddo
    !    write(unit,*)""
    ! enddo
    Linterval = 150 !Maximum number of allowed intervals to look for zeros&poles
    !
    allocate(Xcsign(0:Linterval))
    allocate(Ipoles(Lk),Iweight(Lk))
    allocate(Mpoles(Lk,Linterval),Mweight(Lk,Linterval))
    !
    !FINDING THE POLES:
    !assume \eps=0.d0 ==> the ImSigma(poles)=0 this condition should be automatically
    !verified at the pole from definition of the pole (the ImSigma at the pole is just
    !an artificial broadening of an otherwise delta function, whose position should be 
    !determined by ReSigma only.
    Ipoles=0.d0   
    Mpoles=0.d0
    write(LOGfile,*)"Solving for the poles..."
    maxNinterval=-1
    do ik=1,Lk
       !
       do i=1,Lreal
          ! forall(iorb=1:Nso)zeta(iorb,iorb)=wr(i)+xmu
          ! zeta(:,:) = (wr(i)+xmu)*eye(Nso) - dreal(so2j(Sreal(:,:,:,:,i),Nso))
          ! Den(i) = dreal((zeta(1,1) - Hk(1,1,ik))*(zeta(2,2) - Hk(2,2,ik))) - abs(Hk(1,2,ik)*Hk(2,1,ik))
          zeta = (wr(i)+xmu)*eye(Nso) - dreal(so2j(Sreal(:,:,:,:,i),Nso))
          gfk  = zeta - Hk(:,:,ik)
          Den(i) = det(Gfk)
       enddo
       Xcsign(0)=0.d0
       count=0
       sign_old=sgn(Den(Lreal/2+1))
       do i=Lreal/2+1,Lreal
          sign=sgn(Den(i))
          if(sign*sign_old<1)then
             count=count+1
             if(count>Linterval)stop"Allocate Xcsign to a larger array."
             Xcsign(count)=wr(i)
          endif
          sign_old=sign
       enddo
       Ninterval=count
       if(count>maxNinterval)maxNinterval=count
       call init_finter(finter_func,wr,Den,3)
       ! do int=1,Ninterval
       !    Mpoles(ik,int) = fzero_brentq(det_poles,Xcsign(int-1),Xcsign(int))
       !    Mweight(ik,int)= get_weight(hk(:,:,ik)-so2j(Smats(:,:,:,:,1),Nso))
       ! enddo
       ipoles(ik) = fzero_brentq(det_poles,-wr(Lreal),wr(Lreal))
       iweight(ik)= get_weight(hk(:,:,ik)-so2j(Smats(:,:,:,:,1),Nso))
       call delete_finter(finter_func)
    enddo
    call splot("BHZpoles.ed",(/(ik-1,ik=1,Lk)/),ipoles(:),iweight(:))
    unit=free_unit()
    ! open(unit,file="BHZpoles_all.ed")
    ! do int=1,maxNinterval
    !    if(any((Mpoles(:,int)/=0.d0)))then
    !       do ik=1,Lk
    !          if(Mpoles(ik,int)/=0.d0)write(unit,*)ik-1,Mpoles(ik,int),Mweight(ik,int)
    !       enddo
    !       write(unit,*)""
    !    endif
    ! enddo
    ! close(unit)
    !
    if(.false.)then
       call build_hk
       allocate(Ipoles3d(Nk,Nk))
       allocate(Mpoles3d(Nk,Nk,Linterval))
       write(LOGfile,*)"Solving for the 3d poles..."
       do ik=1,Lk
          ix=ik2ix(ik)
          iy=ik2iy(ik)
          do i=1,Lreal
             forall(iorb=1:Nso)zeta(iorb,iorb)=wr(i)+xmu
             zeta(:,:) = zeta(:,:) - dreal(so2j(Sreal(:,:,:,:,i),Nso))
             Den(i) = dreal((zeta(1,1) - Hk(1,1,ik))*(zeta(2,2) - Hk(2,2,ik))) - Hk(1,2,ik)*Hk(2,1,ik)
          enddo
          !
          Xcsign(0)=0.d0
          count=0
          sign_old=sgn(Den(Lreal/2+1))
          do i=Lreal/2+1,Lreal
             sign=sgn(Den(i))
             if(sign*sign_old<1)then
                count=count+1
                if(count>Linterval)stop"Allocate Xcsign to a larger array."
                Xcsign(count)=wr(i)
             endif
             sign_old=sign
          enddo
          Ninterval=count
          if(count>maxNinterval)maxNinterval=count
          call init_finter(finter_func,wr,Den,3)
          do int=1,Ninterval
             Mpoles3d(ix,iy,int) = fzero_brentq(det_poles,Xcsign(int-1),Xcsign(int))
          enddo
          ipoles3d(ix,iy) = fzero_brentq(det_poles,0.d0,wr(Lreal))
          call delete_finter(finter_func)
       enddo
       call splot3d("3dBHZpoles.ed",kxgrid,kygrid,ipoles3d)
       do int=1,maxNinterval
          call splot3d("3dBHZpoles_layer_"//reg(txtfy(int))//".ed",kxgrid,kygrid,Mpoles3d(:,:,int))
       enddo
    endif
  end subroutine get_poles

  function det_poles(w) result(det)
    real(8),intent(in) :: w
    real(8)            :: det
    det = finter(finter_func,w)
  end function det_poles

  function get_weight(hk) result(wt)
    complex(8),dimension(4,4) :: hk,foo
    real(8),dimension(4)      :: eigv
    real(8) :: wt
    foo = hk
    call matrix_diagonalize(foo,eigv)
    wt = sum(foo(:,1))
  end function Get_Weight






  !---------------------------------------------------------------------
  !PURPOSE: GET DELTA FUNCTION ON REAL AXIS
  !---------------------------------------------------------------------
  subroutine get_deltaw
    integer                                     :: i,j,ik,iorb,jorb,ispin,jspin,iso,unit
    complex(8),dimension(Nso,Nso)               :: zeta,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal,Smats
    complex(8),dimension(:,:,:,:),allocatable   :: gk,gfoo,ReSmat
    complex(8),dimension(:,:,:),allocatable     :: Hktilde
    complex(8)                                  :: iw
    real(8)                                     :: wm(Lmats),wr(Lreal)
    real(8),dimension(:,:),allocatable          :: Ktrim,Ev
    character(len=20)                           :: suffix
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    delta = zero
    !
    call build_hk()
    !
    print*,"Get Delta(w):"
    allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
    call read_sigma(Sreal)
    if(allocated(delta))deallocate(delta)
    allocate(delta(Nspin,Nspin,Norb,Norb,Lreal))
    call start_progress(LOGfile)
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       forall(iorb=1:Nso)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - so2j(Sreal(:,:,:,:,i),Nso)
       fg=zero
       do ik=1,Lk         
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*wtk(ik)
       enddo
       call matrix_inverse(fg)
       if(cg_scheme=='weiss')then
          gdelta = fg + so2j(Sreal(:,:,:,:,i),Nso)
          call matrix_inverse(gdelta)
       else
          gdelta = zeta(:,:) - bhzHloc - fg(:,:)
       endif
       delta(:,:,:,:,i) = j2so(gdelta(:,:))
       call progress(i,Lreal)
    enddo
    call stop_progress()
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Delta"//reg(suffix),wr,-dimag(delta(ispin,ispin,iorb,iorb,:))/pi,dreal(delta(ispin,ispin,iorb,iorb,:)))
       enddo
    enddo
    return
  end subroutine get_deltaw













  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    integer                   :: N
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
    ! Hk(1,4) = -delta ; Hk(4,1)=-delta
    ! Hk(2,3) =  delta ; Hk(3,2)= delta
    Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
    Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
  end function hk_bhz

  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky,epsik
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2

  function inverse_gk(zeta,hk) result(gk)
    complex(8)                  :: zita(2)
    complex(8),dimension(4,4)   :: zeta,hk
    complex(8),dimension(4,4)   :: gk
    integer                     :: i
    select case(lambda==0.d0)
    case default
       zita(1)=zeta(1,1);zita(2)=zeta(2,2)
       gk(1:2,1:2) = inverse_gk2x2(zita,hk(1:2,1:2))
       zita(1)=zeta(3,3);zita(2)=zeta(4,4)
       gk(3:4,3:4) = inverse_gk2x2(zita,hk(3:4,3:4))
    case (.true.)
       gk=zero
       forall(i=1:4)gk(i,i)=one/(zeta(i,i)-hk(i,i))
    end select
  end function inverse_gk
  !
  function inverse_gk2x2(zeta,hk) result(gk)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8),dimension(2)     :: zeta
    complex(8),dimension(2,2)   :: gk
    Complex(8)                  :: delta,ppi,vmix
    gk=zero
    delta = zeta(1) - hk(1,1)
    ppi   = zeta(2) - hk(2,2)
    vmix  = -hk(1,2)
    gk(1,1) = one/(delta - abs(vmix)**2/ppi)
    gk(2,2) = one/(ppi - abs(vmix)**2/delta)
    gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    gk(2,1) = -conjg(vmix)/(ppi*delta - abs(vmix)**2)
  end function inverse_gk2x2

  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
    g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
    ! else
    !    g0k = -hk
    !    forall(i=1:4)g0k(i,i) = iw + xmu + g0k(i,i)
    !    call matrix_inverse(g0k)
    ! endif
  end function inverse_g0k
  !
  function inverse_g0k2x2(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k2x2

  function Eigk(hk) result(eig)
    complex(8),dimension(4,4) :: hk
    real(8),dimension(4)      :: eig
    call matrix_diagonalize(hk,eig)
  end function Eigk

  function Eigk2x2(hk) result(eig)
    complex(8),dimension(2,2) :: hk
    real(8),dimension(2)      :: eig
    call matrix_diagonalize(hk,eig)
    eig(1)=hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig(2)=hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig = eig/2.d0
  end function Eigk2x2









  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do iorb=1,Norb
       unit=free_unit()
       if(L==Lreal)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw.ed"
       elseif(L==Lmats)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed"
       endif
       write(*,*)"read from file=","impSigma"//reg(suffix)
       open(unit,file="impSigma"//reg(suffix),status='old')
       do i=1,L
          read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
          forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
       enddo
       close(unit)
    enddo
  end subroutine read_sigma


  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop"error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop"error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg,Nso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nso,Nso)               :: g
    integer                                     :: Nso,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nso,Nso)               :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so


end program ed_bhz



