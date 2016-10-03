program ed_graphene
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: graphHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bk1,bk2,pointK,pointKp,bklen
  complex(8),dimension(4,4)                     :: GammaX,GammaY,GammaZ,Gamma0
  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputGRAPHENE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputGRAPHENE.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))


  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Norb/=1.OR.Nspin/=2)stop "Wrong setup from input file: Norb!=1 OR Nspin!=2 (This is AFM!)"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso

  Gamma0 = kron_pauli(pauli_0,pauli_0)
  GammaX = kron_pauli(pauli_0,pauli_x)
  GammaY = kron_pauli(pauli_0,pauli_y)
  GammaZ = kron_pauli(pauli_0,pauli_z)

  !FOLLOWING REV.MOD.PHYS.81.109(2009)
  !Lattice basis (a=1; a0=sqrt3*a) is:
  !\a_1 = a0 [ sqrt3/2 , 1/2 ]
  !\a_2 = a0 [ sqrt3/2 ,-1/2 ]
  !
  !
  !nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]
  !
  !next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
  a1=[ sqrt(3d0)/2d0, 1d0/2d0]
  a2=[ sqrt(3d0)/2d0,-1d0/2d0]
  a3=a2-a1

  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/3d0
  bk1=bklen*[ 1d0/2d0 ,  sqrt(3d0)/2d0 ]
  bk2=bklen*[ 1d0/2d0 , -sqrt(3d0)/2d0 ]

  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


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
  Hloc = lso2nnn_reshape(graphHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(Bath)
  ! break SU(2) symmetry for magnetic solutions
  call break_symmetry_bath(Bath(1,:),sb_field, 1.d0)
  call break_symmetry_bath(Bath(2,:),sb_field,-1.d0)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(Bath,Hloc,iprint=1)
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     Smats(2,1,1,:,:,:) = Smats(1,2,2,:,:,:)
     Smats(2,2,2,:,:,:) = Smats(1,1,1,:,:,:)

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=4)

     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc,iprint=4)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc,iprint=4)
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=2)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(1,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  call ed_get_sigma_real_lattice(Sreal,Nlat)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal,iprint=4)


contains




  !--------------------------------------------------------------------!
  !Graphene HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_graphene_afm_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                       :: h0,hx,hy,hz
    real(8)                       :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = 0d0!2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = 0d0!2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    hk = h0*Gamma0 + hx*GammaX + hy*GammaY + hz*GammaZ
  end function hk_graphene_afm_model






  !---------------------------------------------------------------------
  !PURPOSE: Get Graphene Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky  
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    integer                               :: unit
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal,fooSreal
    real(8),dimension(2)                  :: kvec
    real(8)                               :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)      :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable      :: KPath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k) Graphene:",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    allocate(kxgrid(Nk),kygrid(Nk))
    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          kxgrid(ix) = kvec(1)
          kygrid(iy) = kvec(2)
          Hk(:,:,ik) = hk_graphene_afm_model(kvec,Nlso)
       enddo
    enddo
    Wtk = 1d0/Lk


    if(present(file))then
       call write_hk_w90("Hkrfile_graphene.data",&
            No=Nlso,&
            Nd=Norb,&
            Np=0,&
            Nineq=1,&
            Hk=Hk,&
            kxgrid=kxgrid,kygrid=kygrid,kzgrid=[0d0])
    endif
    !
    allocate(graphHloc(Nlso,Nlso))
    graphHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(graphHloc))<1.d-4)graphHloc=0d0
    call write_Hloc(graphHloc)
    call write_Hloc(graphHloc,'Hloc.txt')
    !
    !


    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_path(hk_graphene_afm_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,orange1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.nint")



    !Build the local GF:
    Gmats=zero
    Greal=zero
    fooSmats =zero
    fooSreal =zero
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,fooSmats,iprint=1)
    call dmft_gloc_realaxis(Hk,Wtk,Greal,fooSreal,iprint=1)
    !
  end subroutine build_hk








end program ed_graphene



