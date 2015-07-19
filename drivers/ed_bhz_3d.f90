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
  real(8),allocatable    :: Bath(:),Bath_Prev(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:),ik2iz(:)
  !variables for the model:
  integer                :: Nk,Nkpath
  real(8)                :: ez,mh,lambda
  real(8)                :: wmixing
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym
  !
  real(8),dimension(2)   :: Eout
  !
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  !
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  !
  call parse_input_variable(ez,"ez",finput,default=1d0)
  call parse_input_variable(mh,"MH",finput,default=0.25d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  !
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))
  !
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  !
  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  !
  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  !
  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb))
  allocate(Bath_Prev(Nb))
  call ed_init_solver(bath)
  call set_hloc(j2so(bhzHloc))
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call build_z2_indices( so2j(Smats(:,:,:,:,1),Nso) ) !so2j transforms a Nspin:Nspin:Norb:Norb into a Nso:Nso matrix

     !Get the Local GF and the Weiss field/Delta function to be fitted
     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=2)
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=j2so(bhzHloc),iprint=2)

     !Fit the new bath, starting from the old bath + the supplied delta
     call ed_chi2_fitgf(delta(1,1,:,:,:),bath,ispin=1)
     if(.not.spinsym)then
        call ed_chi2_fitgf(delta(2,2,:,:,:),bath,ispin=2)
     else
        call spin_symmetrize_bath(bath,save=.true.)
     endif

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
     Bath_Prev=Bath

     if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0)call end_loop
  enddo
  !
  Eout = ed_kinetic_energy(Hk,Wtk,Smats)
  !
#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif


contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN IN THE FULL BZ
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
    !
    !get H(k) and solve the non-interacting problem along a path in 3d:
    call build_hk_path()
    !
    !Get H(k) in the BZ:    
    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**3
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    if(allocated(kxgrid))deallocate(kxgrid)
    if(allocated(kygrid))deallocate(kygrid)
    if(allocated(kzgrid))deallocate(kzgrid)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    allocate(kxgrid(Nk),kygrid(Nk),kzgrid(Nk))
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    kzgrid = kgrid(Nk)
    Hk     = build_hk_model(hk_bhz,Nso,kxgrid,kygrid,kzgrid)
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90(file,Nso,&
            Nd=Norb,&
            Np=1,   &
            Nineq=1,&
            hk=Hk,  &
            kxgrid=kxgrid,&
            kygrid=kygrid,&
            kzgrid=kzgrid)
    endif
    !
    !Get the local part of H(k)
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0d0
    if(ED_MPI_ID==0)call write_Hloc(bhzHloc)
    !
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
       call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_realw.ed",wr,&
            -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
       n0(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
    enddo
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG A 3D PATH IN THE BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path()
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),allocatable :: kpath
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
    Npts = 9
    Lk=(Npts-1)*Nkpath
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    if(allocated(kxgrid))deallocate(kxgrid)
    if(allocated(kygrid))deallocate(kygrid)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    allocate(kxgrid(Lk))
    allocate(kygrid(Lk))
    allocate(kzgrid(Lk))
    allocate(kpath(Npts,3))
    kpath(1,:)=[0,0,0]!G
    kpath(2,:)=[1,0,0]!X
    kpath(3,:)=[1,1,0]!M
    kpath(4,:)=[0,0,0]!G
    kpath(5,:)=[1,1,1]!R
    kpath(6,:)=[1,0,1]!M
    kpath(7,:)=[0,0,1]!X
    kpath(8,:)=[1,1,1]!R
    kpath(9,:)=[0,0,0]!G
    kpath=kpath*pi
    kxgrid = kgrid_from_path(kpath,Npts,Nkpath,1)
    kygrid = kgrid_from_path(kpath,Npts,Nkpath,2)
    kzgrid = kgrid_from_path(kpath,Npts,Nkpath,3)
    Hk     = build_hk_model(hk_bhz,Nso,kxgrid,kygrid,kzgrid)
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0)  call solve_Hk_along_BZpath(hk_bhz,Nso,kpath,Lk,&
         colors_name=[character(len=20) :: 'red','blue','red','blue'],&
         points_name=[character(len=20) :: "G","X","M","G","R","M","X","R","G"],&
         file="Eigenband.nint")
  end subroutine build_hk_path





  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - cos(kx) - cos(ky) - ez*cos(kz))*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y
    Hk(3:4,3:4) = conjg( &
         (Mh-cos(-kx) - cos(-ky) - ez*cos(-kz))*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x
  end function hk_bhz








  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine build_z2_indices(sigma0)
    integer                                :: unit
    complex(8),dimension(Nso,Nso),optional :: sigma0(Nso,Nso)
    complex(8),dimension(Nso,Nso)          :: sigma0_
    !
    integer,dimension(4)                   :: z2
    !
    sigma0_=zero;if(present(sigma0))sigma0_=sigma0
    !Evaluate the Z2 index:
    !STRONG TI
    z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi,sigma0_)
    !WEAK TI
    !K=1: n_1=1, n_2,3=0,1
    z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    !K=2: n_2=1, n_1,2=0,1
    z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi,sigma0_)
    !k=3: n_3=1, n_1,2=0,1
    z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    unit=free_unit()
    open(unit,file="z2_invariant.ed")
    write(unit,*)z2
    close(unit)
  end subroutine build_z2_indices

  function z2_number(ktrims,sigma0) result(z2)
    real(8),dimension(:,:),intent(in)            :: ktrims
    complex(8),dimension(Nso,Nso)                :: sigma0
    complex(8),dimension(Nso,Nso,size(Ktrims,2)) :: Htrims
    complex(8),dimension(size(Ktrims,2))         :: Delta
    integer                                      :: z2
    integer                                      :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_bhz(Ktrims(:,itrim),Nso) + sigma0(:,:)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number

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

  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k = (iw + xmu)*zeye(4)-hk
    call inv(g0k)
  end function inverse_g0k



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



