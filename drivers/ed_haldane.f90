program ed_haldane
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso
  logical                                       :: converged

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: halHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk
  real(8),allocatable,dimension(:)              :: kxgrid
  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: b1,b2,b3

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,tsp,phi,delta,Mh,wmixing
  character(len=16)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0


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
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputHALDANE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputHALDANE.conf",default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH","inputHALDANE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputHALDANE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))

  if(Nspin/=1.OR.Norb/=1)stop "Wrong setup from input file: Nspin=Norb=1"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso

  !Following:
  !http://www-personal.umich.edu/~sunkai/teaching/Fall_2012/chapter3_part7.pdf
  !Lattice basis (a=1) is:
  !\nu_1 = [ sqrt3  , 0  ]
  !\nu=2 = [-sqrt3/2,3/2 ]
  !nearest neighbor: A-->B, B-->A
  a1 = [0d0           ,  1d0    ]
  a2 = [-sqrt(3d0)/2d0, -1d0/2d0]
  a3 = [ sqrt(3d0)/2d0, -1d0/2d0]

  !next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
  b1=[ sqrt(3d0), 0d0]
  b2=[-sqrt(3d0)/2d0, 3d0/2d0]
  b3=[-sqrt(3d0)/2d0,-3d0/2d0]



  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(halHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver_lattice(Bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve_lattice(bath,Hloc,iprint=1)

     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     call ed_get_sigma_real_lattice(Sreal,Nlat)

     ! compute the local gf:
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=1)

     ! compute the Weiss field (only the Nineq ones)
     call ed_get_weiss_lattice(Gmats,Smats,Weiss,Hloc,iprint=1)

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf_lattice(Bath,Weiss,Hloc,ispin=1)


     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     if(ED_MPI_ID==0)converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0)call end_loop
  enddo

#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif



contains




  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    real(8)                         :: kdota(3),kdotb(3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !(k.b_j)
    kdotb(1) = dot_product(kpoint,b1)
    kdotb(2) = dot_product(kpoint,b2)
    kdotb(3) = dot_product(kpoint,b3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdotb(:)) )
    hx =-ts*sum( cos(kdota(:)) )
    hy =-ts*sum( sin(kdota(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdotb(:)) ) + Mh 
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
  end function hk_haldane_model






  !---------------------------------------------------------------------
  !PURPOSE: Get Haldane Model Hamiltonian
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
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw,n0(Nlso)

    call build_hk_path()

    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk*Nk
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(Nlso,Nlso,Lk))
    allocate(wtk(Lk))
    allocate(kxgrid(Nk))
    kxgrid = kgrid(Nk)
    Hk     = TB_build_model(hk_haldane_model,Nlso,kxgrid,kxgrid,[0d0])
    wtk    = 1d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90("Hkrfile_Haldane.data",&
            No=Nlso,&
            Nd=Norb,&
            Np=0,&
            Nineq=1,&
            Hk=Hk,&
            kxgrid=kxgrid,kygrid=kxgrid,kzgrid=[0d0])
    endif
    !
    allocate(halHloc(Nlso,Nlso))
    halHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(halHloc))<1.d-9)halHloc=0d0
    if(ED_MPI_ID==0)call write_Hloc(halHloc)
    !
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
    do iorb=1,Nlso
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,&
            -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
       n0(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
    enddo
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: solve H(k) along a specified path in the BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path(kpath_)
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    character(len=64)                  :: file

    if(present(kpath_))then
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="bands_path"
    else
       if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) along the path G-K-K`-G:"
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_Gamma
       kpath(2,:)=[2d0/3d0/sqrt(3d0),0d0]*pi2
       kpath(3,:)=[1d0/3d0/sqrt(3d0),1/3d0]*pi2
       kpath(4,:)=kpoint_Gamma
       file="non_interacting_bands"
    endif
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(Nlso,Nlso,Lk))
    allocate(wtk(Lk))
    Hk  = TB_build_model(hk_haldane_model,Nlso,kpath,Nkpath)
    wtk = 1d0/Lk
    if(ED_MPI_ID==0) call TB_Solve_path(hk_haldane_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file=reg(file))
  end subroutine build_hk_path













  function inverse_g0k(zeta,hk) result(gk)
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
  end function inverse_g0k


end program ed_haldane



