!AFM in STRIPE GEOMETRY and 2 ATOMS IN THE UNIT CELL, SEE STRUCTURE IN A SINGLE RUN
program ed_bhz_2d_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop
  integer                                       :: Ly
  integer                                       :: Nineq
  integer                                       :: Ncell
  integer                                       :: Nlat
  integer                                       :: Nso
  integer                                       :: ilat,iy,iorb,ispin,ineq,i
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath_ineq
  real(8),allocatable,dimension(:,:)            :: Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !Nineq:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_ineq
  !
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  !Hmiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  real(8),allocatable,dimension(:)              :: Wtk,sbpattern
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  !
  integer                                       :: Nk,Nkpath                           
  real(8)                                       :: e0,mh,lambda,wmixing
  logical                                       :: spinsym,lysym,neelsym,rebuild_sigma
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  !
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
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(lysym,"LYSYM",finput,default=.true.,comment="Enforce symmetry on half of the stripe")
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.,comment="Enforce spin symmetry (PM solution)")
  call parse_input_variable(neelsym,"NEELSYM",finput,default=.false.,comment="Enforce Neel symmetry on ineq. atoms in the same cell")
  call parse_input_variable(rebuild_sigma,"REBUILD_SIGMA",finput,default=.false.)
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
  Nso   = Nspin*Norb

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Ncell= 2                      !two atoms in the unit cell
  Nlat = Ncell*Ly               !total number of lattice sites (with unit cell multiplicity)
  if(lysym .AND. mod(Ly,2)==0 )stop "Wrong setup: (lysym=T) & Ly%2=0. Cant reduce Ly"
  if(neelsym)then
     Nineq= Ly                     !Number of inequivalent lattice sites
     if(lysym)Nineq = (Ly-1)/2+1
  else
     Nineq=Nlat
     if(lysym)Nineq = Nlat/2+1
  endif
  print*,"Ly    =",Ly
  print*,"Nlat  =",Nlat
  print*,"Nineq =",Nineq
  print*,""


  if(master)then
     call print_structure(Ly)
     call print_structure(Ly,"Structure_Lattice.dat")
  endif
  call sleep(2)


  !Allocate Full Functions:
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Weiss_ineq=zero
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
  !
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ;Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ;Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ;Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)) ;Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)) ;Greal=zero
  !
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0
  allocate(Zmats(Nlat*Nso,Nlat*Nso));Zmats=eye(Nlat*Nso)


  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  !allocate to Ncell[2] times the number of layers[Ly]:1A-1B,...,NA-NB
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc     =zero
  !allocate to Nineq: only consider 1 per cell and half of the layers: 1A,...,NA
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
  !spread the full local H0 to Ncell*Ly 4x4[Nso**2] hamiltonians:
  Hloc = lso2nnn_reshape(bhzHloc,Ncell*Ly,Nspin,Norb)
  !select only the inequivalent ones: 1,3,5,...,Ly-1
  do ineq=1,Nineq
     ilat = ineq2ilat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo





  !======================> DMFT <========================+!
  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb) )
  allocate(Bath_prev(Nineq,Nb) )
  call ed_init_solver(Comm,Bath_ineq)
  call MPI_Barrier(Comm,ineq)
  write(*,"(A)")"Breaking Symmetry pattern:"
  allocate(sbpattern(Nineq))
  if(neelsym)then
     sbpattern=1d0
  else
     do ineq=1,Nineq
        sbpattern(ineq)=(-1d0)**(mod(ineq,2)+1)
     enddo
  endif
  do ineq=1,Nineq
     call break_symmetry_bath(Bath_ineq(ineq,:),sb_field,sbpattern(ineq))
  enddo
  call MPI_Barrier(Comm,ineq)

  if(master)then
     open(10,file="Symmetry_breaking_pattern.dat")
     do ineq=1,Nineq
        write(*,"(A5,I3,2X,A5,I3,A4,F9.3)")"Ineq=",ineq,"Ilat=",ineq2ilat(ineq),"SB=",sbpattern(ineq)
        write(10,"(A5,I3,2X,A5,I3,A4,F9.3)")"Ineq=",ineq,"Ilat=",ineq2ilat(ineq),"SB=",sbpattern(ineq)
     enddo
     close(10)
  endif




  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")   
     !
     !
     !
     ! solve the impurities on each inequivalent y-layer
     call ed_solve(comm,Bath_ineq,Hloc_ineq,iprint=1)
     !
     !
     !
     ! retrieve the self-energies
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     do ilat=1,Nlat
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ilat2ineq(ilat),:,:,:,:,:)
     enddo
     if(neelsym)then
        do ilat=2,Nlat,2
           do ispin=1,2
              Smats(ilat,ispin,ispin,:,:,:)=Smats(ilat-1,3-ispin,3-ispin,:,:,:)
           enddo
        enddo
     endif
     S0 = Smats(:,:,:,:,:,1)![Nlat][Nspin][Nspin][Norb][Norb]
     do ilat=1,Nlat
        Zfoo(ilat,:,:) = select_block(ilat,S0)
        do iorb=1,Nso
           i = iorb + (ilat-1)*Nso
           Zmats(i,i)  = 1.d0/( 1.d0 + abs( dimag(Zfoo(ilat,iorb,iorb))/(pi/beta) ))
        enddo
     enddo
     !
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Comm,Hkr,Wtk,Gmats,Smats,iprint=4) !tridiag option off
     do ineq=1,Nineq
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ineq2ilat(ineq),:,:,:,:,:)
     enddo
     !
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     else
        call dmft_delta(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=4)
     endif
     !
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     !
     !
     if(iloop>1)Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     !
     if(master)converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)
     !
     if(master)call end_loop
  enddo




  call ed_get_sigma_real(Sreal_ineq,Nineq)
  do ilat=1,Nlat
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ilat2ineq(ilat),:,:,:,:,:)
  enddo
  do ilat=2,Nlat,2
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
    integer                          :: i,ik,ix
    real(8),allocatable,dimension(:) :: kxgrid
    real(8),dimension(2)             :: bk1,bk2,kvec
    real(8)                          :: kx
    !
    bk1 = pi*[1,-1]
    bk2 = 2*pi*[0,1]
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
    ! kxgrid = kgrid(Nk)
    ! Hkr = TB_build_model(bhz_afm2_edge_model,Ly,Ncell*Nso,kxgrid,[0d0],[0d0],pbc=.false.)
    do ix=1,Nk
       kx=dble(ix-1)/Nk
       kvec = kx*bk1
       Hkr(:,:,ix) = bhz_afm2_edge_model(kvec,Ly,Ncell*Nso,.false.)
    enddo
    allocate(Wtk(Nk))
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
    real(8),dimension(:,:),optional            :: kpath_
    real(8),dimension(:,:),allocatable         :: kpath
    integer                                    :: Npts
    character(len=64)                          :: file
    type(rgb_color),dimension(:,:),allocatable :: colors
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
    !
    allocate(colors(Ly,Ncell*Nso))
    colors = gray88
    colors(1,:) = [red1,gray88,blue1,gray88,red1,gray88,blue1,gray88]
    colors(Ly,:) =[blue1,gray88,red1,gray88,blue1,gray88,red1,gray88]
    !
    call TB_solve_model(bhz_afm2_edge_model,Ly,Ncell*Nso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "-pi","0","pi"],&
         file=reg(file),pbc=.false.)
  end subroutine build_eigenbands









  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ-edge model hamiltonian
  !+-----------------------------------------------------------------------------+!
  function bhz_afm2_edge_model(kpoint,Ly,N,pbc) result(Hrk)
    real(8),dimension(:)            :: kpoint
    real(8)                         :: kx
    integer                         :: Ly,N
    complex(8),dimension(N,N)       :: Hmat,Tmat,TmatH
    real(8),dimension(N,N)          :: Sh
    complex(8),dimension(Ly*N,Ly*N) :: Hrk
    integer                         :: i,Idmin,Idmax,Itmin,Itmax,indx1,indx2
    logical                         :: pbc
    complex(8),dimension(Nso,Nso)   :: M
    complex(8),dimension(Nso,Nso)   :: tx,ty,thx,thy
    complex(8),dimension(Nso,Nso)   :: gamma1,gamma2,gamma5
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
    do i=1,Ly
       indx1 = 1 + (i-1)*Ncell
       indx2 = 2 + (i-1)*Ncell
       Sh    = zero
       Sh(1:Nso,1:Nso)             = dreal(select_block(indx1,S0))
       Sh(Nso+1:2*Nso,Nso+1:2*Nso) = dreal(select_block(indx2,S0))
       Idmin = 1+(i-1)*N
       Idmax =       i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat + Sh
    enddo
    do i=1,Ly-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Ly-1)*N
       Itmax=0+Ly*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
    Hrk = matmul(Zmats,Hrk)
  end function bhz_afm2_edge_model








  function ineq2ilat(ineq) result(ilat)
    integer,intent(in) :: ineq
    integer            :: ilat
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
    ilat=ineq
    if(neelsym)ilat=2*ineq-1
  end function ineq2ilat


  function ilat2ineq(ilat) result(ineq)
    integer,intent(in) :: ilat
    integer            :: ineq
    if(neelsym)then
       ineq=(ilat-1)/2 + 1
       if( lysym .AND. (ineq>Nineq) )ineq=Ly-ineq+1
    else
       ineq=ilat
       if(lysym .AND. ineq>Nineq)then
          ineq=Ly-(ineq-Nineq)-1
          if(mod(ineq,2)==0)ineq=ineq+2
       endif
    endif
  end function ilat2ineq


  subroutine print_structure(Nlayers,file)
    integer                               :: Nlayers
    character(len=*),optional             :: file
    integer                               :: unit,iy,ineq,ilat
    character(len=2),dimension(2) :: str_spin
    integer,dimension(Nlayers,2)          :: spin
    !
    unit=6;if(present(file))open(free_unit(unit),file=trim(file))
    !
    write(unit,"(A)")"Structure:"
    str_spin(1)="up"
    str_spin(2)="dw"
    do iy=1,Nlayers
       if(mod(iy,2)/=0)then
          write(unit,"(A1,I3,1X,A2,A7,I3,1X,A2,A7)")"|",1+(iy-1)*2,str_spin(1)," ---- ",2+(iy-1)*2,str_spin(2),"      |"
       else
          write(unit,"(A7,I3,1X,A2,A7,I3,1X,A2,A1)")"      |",1+(iy-1)*2,str_spin(1)," ---- ",2+(iy-1)*2,str_spin(2),"|"
       endif
    enddo
    write(unit,"(A)")""
    write(unit,"(A)")"Ineq2Ilat:"
    do ineq=1,Nineq
       write(unit,"(A5,I3,1X,A5,I3)")"Ineq=",ineq,"Ilat=",ineq2ilat(ineq)
    enddo
    write(unit,"(A)")""
    write(unit,"(A)")"Ilat2Ineq:"
    do ilat=1,Nlat
       write(unit,"(A5,I3,1X,A5,I3)")"Ilat=",ilat,"Ineq=",ilat2ineq(ilat)
    enddo
    write(unit,"(A)")""
    if(present(file))close(unit)
  end subroutine print_structure



  function select_block(ip,Matrix) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block




end program ed_bhz_2d_edge
