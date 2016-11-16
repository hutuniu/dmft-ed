!     | 1Ao  ---  1Bo | = 1-2
!     | 2Ao  ---  2Bo | = 3-4
!     | 3Ao  ---  3Bo | = 5-6
!            ...  
!     | LyAo --- LyBo | = 2Ly-1,2Ly
!
!   o =  1-2_UP x 1-2_DW
!
program ed_bhz_2d_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI

  implicit none

  integer                                       :: iloop
  integer                                       :: Nineq
  integer                                       :: Ncell
  integer                                       :: Nso
  integer                                       :: ilat,iy,iorb,ispin,ineq,i
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath_ineq
  real(8),allocatable,dimension(:,:)            :: Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !Nineq:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq

  !Hmiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  real(8),allocatable,dimension(:)              :: Wtk
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq

  integer                                       :: Nk,Ly,Nkpath                           
  real(8)                                       :: e0,mh,lambda,wmixing
  logical                                       :: spinsym,tridiag,lrsym
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile

  integer                                       :: comm,rank
  logical                                       :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)



  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ_EDGE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput),comm)


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  !set the spin-orbit dimension:
  Nso   = Nspin*Norb

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Ncell= 2                      !two atoms in the unit cell
  Nineq= Ly                     !Number of inequivalent lattice sites (maybe Nineq<Ly)
  if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     if(master)print*,"Using L-R Symmetry. Solve",Nineq," of",Ly," sites."
     call sleep(2)
  endif


  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  !allocate to Ncell[2] times the number of layers[Ly]:1A-1B,...,NA-NB
  allocate(Hloc(Ly*Ncell,Nspin,Nspin,Norb,Norb));Hloc     =zero
  !allocate to Nineq: only consider 1 per cell and half of the layers: 1A,...,NA
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
  !spread the full local H0 to Ncell*Ly 4x4[Nso**2] hamiltonians:
  Hloc = lso2nnn_reshape(bhzHloc,Ncell*Ly,Nspin,Norb)
  !select only the inequivalent ones: 1,3,5,...,Ly-1
  do ineq=1,Nineq
     ilat = ineq2ilat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  !Allocate Full Functions:
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Weiss_ineq=zero
  !
  allocate(Smats(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lmats)) ;Smats=zero
  allocate(Sreal(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lreal)) ;Sreal=zero
  allocate(Gmats(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lmats)) ;Gmats=zero
  allocate(Greal(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lreal)) ;Greal=zero





  !======================> DMFT <========================+!
  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb) )
  allocate(Bath_prev(Nineq,Nb) )
  call ed_init_solver(Comm,Bath_ineq)
  do ineq=1,Nineq
     call break_symmetry_bath(Bath_ineq(ineq,:),sb_field,(-1d0)**(ineq+1))
  enddo

  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")   

     ! solve the impurities on each inequivalent y-layer
     call ed_solve(comm,Bath_ineq,Hloc_ineq,iprint=1)

     ! retrieve the self-energies
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     ! extend ineq. to full lattice + AFM symmetrization
     do ilat=1,Ly
        ineq = ilat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
     enddo
     !even,sigma ---> odd,bar_sigma
     do ilat=2,Ly,2
        do ispin=1,2
           Smats(ilat,ispin,ispin,:,:,:)=Smats(ilat-1,3-ispin,3-ispin,:,:,:)
        enddo
     enddo

     ! compute the local gf:
     call dmft_gloc_matsubara(Comm,Hkr,Wtk,Gmats,Smats,iprint=4) !tridiag option off for the time being
     do ineq=1,Nineq
        ilat = ineq2ilat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo


     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     else
        call dmft_delta(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     endif


     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     call Bcast_MPI(Comm,Bath_ineq)

     Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     if(master)converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)

     if(master)call end_loop
  enddo



  call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
  do ilat=1,Ly
     ineq = ilat2ineq(ilat)
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
  enddo
  do ilat=2,Ly,2
     do ispin=1,2
        Sreal(ilat,ispin,ispin,:,:,:)=Sreal(ilat-1,3-ispin,3-ispin,:,:,:)
     enddo
  enddo
  call dmft_gloc_realaxis(Comm,Hkr,Wtk,Greal,Sreal,iprint=4)


  call build_eigenbands()


  call finalize_MPI()



contains





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build the BHZ Hamiltonian H(k_x,R_y) on the STRIPE along Y
  !+-----------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional        :: file
    integer                          :: i,ik
    real(8),allocatable,dimension(:) :: kxgrid
    !
    !SETUP THE H(kx,Ry):
    if(allocated(Hkr))deallocate(Hkr)
    if(allocated(Wtk))deallocate(Wtk)
    write(LOGfile,*)"Build H(kx,R) for BHZ-AFM-stripe:"
    write(*,*)"# of kx-points     :",Nk
    write(*,*)"# of y-layers      :",Ly
    !
    allocate(Kxgrid(Nk))
    allocate(Hkr(Ncell*Ly*Nso,Ncell*Ly*Nso,Nk))
    allocate(Wtk(Nk))
    kxgrid = kgrid(Nk)
    Hkr = TB_build_model(bhz_afm2_edge_model,Ly,Ncell*Nso,kxgrid,[0d0],[0d0],pbc=.false.)
    Wtk = 1d0/Nk
    !
    allocate(bhzHloc(Ncell*Ly*Nso,Ncell*Ly*Nso))
    bhzHloc = sum(Hkr,dim=3)/Nk
    where(abs(bhzHloc)<1.d-6)bhzHloc=zero
    !
  end subroutine build_hkr




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve H_BHZ(k_x,R_y) along the 1d -pi:pi path in the BZ.
  !+-----------------------------------------------------------------------------+!  
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    integer                            :: Npts
    character(len=64)                  :: file
    if(.not.present(kpath_))then
       write(LOGfile,*)"Solve H(kx,y) along [-pi:pi]:"
       Npts=3
       allocate(Kpath(Npts,1))
       kpath(1,:)=[-1]*pi
       kpath(2,:)=[ 0]*pi
       kpath(3,:)=[ 1]*pi
       file="Eigenbands_Gamma_AFM.nint"
    else
       write(LOGfile,*)"Solve H(kx,y) along a given path:"
       Npts = size(kpath_,1)
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eigenbands_path.nint"
    endif
    call TB_solve_path(bhz_afm2_edge_model,Ly,Ncell*Nso,kpath,Nkpath,&
         colors_name=[gray88,gray88,gray88,gray88,gray88,gray88,gray88,gray88],&
         points_name=[character(len=10) :: "-pi","0","pi"],&
         file=reg(file),pbc=.false.)
  end subroutine build_eigenbands







  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ-edge model hamiltonian
  !+-----------------------------------------------------------------------------+!
  !BHZ on a stripe geometry;
  function bhz_afm2_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    complex(8),dimension(Nso,Nso)       :: M
    complex(8),dimension(Nso,Nso)       :: tx,ty,thx,thy
    complex(8),dimension(Nso,Nso)       :: gamma1,gamma2,gamma5
    !
    !SETUP THE GAMMA MATRICES:
    gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x )
    gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y )
    gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z )
    !
    if(N/=Ncell*Nso)stop "hk_model error: N != Ncell*Nso" 
    !
    kx=kpoint(1)
    !
    M  = Mh*Gamma5
    tx = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    Hmat(1:Nso,1:Nso)             = M
    Hmat(1:Nso,Nso+1:2*Nso)       = tx  + thx*exp(xi*2*kx)
    Hmat(Nso+1:2*Nso,1:Nso)       = thx + tx*exp(-xi*2*kx)
    Hmat(Nso+1:2*Nso,Nso+1:2*Nso) = M
    !
    Tmat(1:Nso,1:Nso)             = zero
    Tmat(1:Nso,Nso+1:2*Nso)       = thy*exp(xi*kx)
    Tmat(Nso+1:2*Nso,1:Nso)       = thy*exp(-xi*kx)
    Tmat(Nso+1:2*Nso,Nso+1:2*Nso) = zero
    !
    TmatH=conjg(transpose(Tmat))
    !
    Hrk=zero
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat !+ dreal(select_block(i,S0)) !< H(k) + Re(Sigma_iy(:Nso,:Nso;omega=0))
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
    !Hrk = matmul(Zmats,Hrk)
  end function bhz_afm2_edge_model








  function ineq2ilat(ineq) result(ilat)
    integer,intent(in) :: ineq
    integer            :: ilat
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
    ilat=2*ineq-1
  end function ineq2ilat


  function ilat2ineq(ilat) result(ineq)
    integer,intent(in) :: ilat
    integer            :: ineq
    ineq=(ilat-1)/2 + 1
    if( lrsym .AND. (ineq>Nineq) )ineq=Ly-ineq+1
  end function ilat2ineq




end program
