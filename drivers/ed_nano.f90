program ed_nano
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                         :: iloop
  logical                                         :: converged
  integer                                         :: ilat,ineq,ispin,iorb
  !bath:
  integer                                         :: Nb
  real(8),allocatable                             :: Bath_prev(:,:),Bath_ineq(:,:)
  !local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal,Sreal_ineq ![Nlat*(Nspin*Norb)**2*Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal,Greal_ineq
  real(8), allocatable,dimension(:)               :: dens,dens_ineq
  real(8), allocatable,dimension(:)               :: docc,docc_ineq
  !hamiltonian input:
  complex(8),allocatable                          :: Hij(:,:,:) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk==1]
  complex(8),allocatable                          :: nanoHloc(:,:),Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  integer                                         :: Nk,Nlso,Nineq
  integer,dimension(:),allocatable                :: lat2ineq,ineq2lat
  integer,dimension(:),allocatable                :: sb_field_sign
  !
  real(8)                                         :: wmixing,Eout(2)
  !input files:
  character(len=32)                               :: finput
  character(len=32)                               :: nfile,hijfile
  !
  logical                                         :: phsym
  logical                                         :: leads
  logical                                         :: kinetic,trans,jrkky,chi0ij
  !non-local Green's function:
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gijmats,Gijreal
  !hybridization function to environment
  complex(8),dimension(:,:,:),allocatable         :: Hyb_mats,Hyb_real ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats/Lreal]
  integer :: mpiID,mpiERR,mpiSIZE

#ifdef _MPI
  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='inputED_NANO.conf')
  call parse_input_variable(nfile,"NFILE",finput,default="nano.in")
  call parse_input_variable(hijfile,"HIJFILE",finput,default="hij.in")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(phsym,"phsym",finput,default=.false.)

  ! parse environment & transport flags
  call parse_input_variable(leads,"leads",finput,default=.false.)
  call parse_input_variable(trans,"trans",finput,default=.false.)
  call parse_input_variable(jrkky,"jrkky",finput,default=.false.)
  call parse_input_variable(chi0ij,"chi0ij",finput,default=.false.)
  call parse_input_variable(kinetic,"kinetic",finput,default=.false.)

  ! read input
  call ed_read_input(trim(finput),MPI_COMM_WORLD)

  ! set input structure hamiltonian
  call build_Hij([nfile,hijfile])

  ! allocate weiss field:
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate self-energy
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate Green's function
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate Hloc
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))
  ! allocate density
  allocate(dens(Nlat))
  allocate(dens_ineq(Nineq))
  ! allocate double occupations
  allocate(docc(Nlat))
  allocate(docc_ineq(Nineq))

  !Hloc = reshape_Hloc(nanoHloc,Nlat,Nspin,Norb)
  Hloc = lso2nnn_reshape(nanoHloc,Nlat,Nspin,Norb)

  ! allocate hybridization matrix
  if(leads) call set_hyb()



  ! postprocessing options

  ! evaluates the kinetic energy
  if(kinetic)then
     ! read converged self-energy
     call read_sigma(Smats_ineq,Sreal_ineq)
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo
     !

     !THIS SETS THE COMMUNICATOR FOR THE MODULES NOT DIRECTLY INCLUDED IN ED (i.e.KineticEnergy,GetGloc,Weiss,etc)
     call ed_set_MPI(MPI_COMM_WORLD)

     ! computes the kinetic energy
     Eout = ed_kinetic_energy(Hij,[1d0],Smats)
     print*,Eout

     call add_ctrl_var(Norb,"Norb")
     call add_ctrl_var(Nspin,"Nspin")
     call add_ctrl_var(beta,"beta")
     call add_ctrl_var(xmu,"xmu")
     call add_ctrl_var(wini,"wini")
     call add_ctrl_var(wfin,"wfin")
     Eout = dmft_kinetic_energy(MPI_COMM_WORLD,Hij,[1d0],Smats(:,1,1,1,1,:))
     print*,Eout

     stop
  endif


  ! computes conductance on the real-axis
  if(trans)then
     ! read converged self-energy
     call read_sigma(Smats_ineq,Sreal_ineq)
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo
     !
     ! allocates and extracts non-local Green's function
     allocate(Gijmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
     allocate(Gijreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
     if(leads)then
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=0,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
     else
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=1)
     endif
     !
     deallocate(Gijmats)
     !
     ! extract the linear response (zero-bias) transmission function
     ! i.e. the conductance in units of the quantum G0 [e^2/h]
     call ed_get_conductance(Gijreal)
     !
     deallocate(Gijreal)
     stop
  endif



  ! compute effective non-local exchange
  if(jrkky)then
     ! read converged self-energy
     call read_sigma(Smats_ineq,Sreal_ineq)
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo
     ! allocates and extracts non-local Green's function
     allocate(Gijmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
     allocate(Gijreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
     if(leads)then
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=0,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=0,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
     else
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=0)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=0)
     endif
     deallocate(Gijmats,Smats)
     ! compute effective exchange
     call ed_get_jeff(Gijreal,Sreal)
     deallocate(Gijreal,Sreal)
     stop
  endif



  ! compute effective non-local exchange
  if(chi0ij)then
     ! read converged self-energy
     call read_sigma(Smats_ineq,Sreal_ineq)
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo
     ! allocates and extracts non-local Green's function
     allocate(Gijmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
     allocate(Gijreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
     if(leads)then
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=0,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=0,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
     else
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=0)
        call ed_get_gij_lattice(Hij,[1d0],Gijmats,Gijreal,Smats,Sreal,iprint=0)
     endif
     deallocate(Gijmats,Smats)
     ! compute bare static non-local susceptibility
     call ed_get_chi0ij(Gijreal)
     deallocate(Gijreal)
     stop
  endif


  !###################################################################################################


  ! setup solver
  Nb=get_bath_dimension()

  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(MPI_COMM_WORLD,Bath_ineq)

  do ineq=1,Nineq
     ilat = ineq2lat(ineq)
     ! break SU(2) symmetry for magnetic solutions
     if(Nspin>1) call break_symmetry_bath(Bath_ineq(ineq,:),sb_field,dble(sb_field_sign(ineq)))
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  !THIS SETS THE COMMUNICATOR FOR THE MODULES NOT DIRECTLY INCLUDED IN ED (i.e.KineticEnergy,GetGloc,Weiss,etc)
  call ed_set_MPI(MPI_COMM_WORLD)

  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   
     bath_prev=bath_ineq

     ! solve impurities on each inequivalent site:
     call ed_solve(MPI_COMM_WORLD,bath_ineq,Hloc_ineq,iprint=0)

     ! retrieve self-energies and occupations(Nineq,Norb=1)
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
     dens_ineq = ed_get_dens_lattice(Nineq,1)
     docc_ineq = ed_get_docc_lattice(Nineq,1)

     !  
     ! spread self-energies and occupation to all lattice sites
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        dens(ilat) = dens_ineq(ineq)
        docc(ilat) = docc_ineq(ineq)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo

     ! compute the local gf:
     if(leads)then
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1,Gamma_mats=Hyb_mats,Gamma_real=Hyb_real)
     else
        call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1)
     endif
     do ineq=1,Nineq
        ilat = ineq2lat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo


     ! compute the Weiss field
     call ed_get_weiss_lattice(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,iprint=0)

     ! fit baths and mix result with old baths
     do ispin=1,Nspin
        call ed_chi2_fitgf(MPI_COMM_WORLD,bath_ineq,Weiss_ineq,Hloc_ineq,ispin)
     enddo

     if(phsym)then
        do ineq=1,Nineq
           call ph_symmetrize_bath(bath_ineq(ineq,:),save=.true.)
        enddo
     endif
     Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     if(mpiID==0)then
        converged = check_convergence(Weiss_ineq(1,1,1,1,1,:),dmft_error,nsuccess,nloop)
        ! alternative convergency criteria
        !converged = check_convergence_local(docc_ineq,dmft_error,nsuccess,nloop)
        if(NREAD/=0.d0) call search_chemical_potential(xmu,sum(dens)/Nlat,converged)
     endif
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
     if(mpiID==0) call end_loop()
  end do


  call save_sigma(Smats_ineq,Sreal_ineq)

  Eout = ed_kinetic_energy(Hij,[1d0],Smats)
  print*,Eout

  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  Eout = dmft_kinetic_energy(MPI_COMM_WORLD,Hij,[1d0],Smats(:,1,1,1,1,:))
  print*,Eout


#ifdef _MPI
  call MPI_FINALIZE(mpiERR)
#endif


contains



  !----------------------------------------------------------------------------------------!
  ! purpose: build real-space Hamiltonian for a nanostructure of size [Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine build_Hij(file)
    character(len=*)     :: file(2)
    integer              :: ilat,jlat,iorb,jorb,is,js,ispin,ie
    integer              :: i,isite,iineq,iineq0,isign
    integer              :: EOF
    character, parameter :: tab = achar ( 9 )
    integer              :: unit,ineq_count
    integer              :: Ns,Ne,Nb,Nk         ! #atoms, #inequivalent, #bands
    real(8)              :: ret,imt
    logical              :: blank_at_right
    character(len=1)     :: next,prev
    character(len=6)     :: site,sign
    if(mpiID==0)write(LOGfile,*)"Build H(R_i,R_j) for a NANO object:"
    ! readin generic input
    ! allocate & fill inequivalent list
    unit = free_unit()
    open(unit,file=trim(file(1)),status='old')
    read(unit,*)Ns,Ne,Nb
    !Checks:
    if(Nb/=Norb)stop "build_Hij error: Nb read from file != Norb in input.conf"
    Nk   = 1
    Nb   = Norb
    Nlat = Ns
    Nineq= Ne
    Nlso = Nlat*Nspin*Norb
    allocate(lat2ineq(Nlat),ineq2lat(Nineq))
    read(unit,"(A1)",advance='no',IOSTAT=EOF)next
    site  = next
    isite = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          site=trim(site)//next
       else
          read(site,"(I6)")isite
          site=""
          i=i+1
          if(i>Nlat)stop "build_Hij error: lattice index > Nlat read from file"
          lat2ineq(i)=isite+1
       endif
       if(EOF<0)exit
    enddo
    if(i<Nlat)stop "build_Hij error: lattice index < Nlat read from file"
    if(mpiID==0)write(*,*)"# of sites      :",Nlat
    if(mpiID==0)write(*,*)"# of ineq sites :",Nineq
    if(mpiID==0)write(*,*)"# of bands      :",Norb
    !
    ineq_count=1
    iineq=lat2ineq(Nlat)
    do i=Nlat,2,-1
       iineq0=lat2ineq(i-1)!iineq
       iineq =lat2ineq(i)
       if(iineq/=iineq0)then
          ineq2lat(iineq)=i
          ineq_count=ineq_count+1
       endif
       !if(ineq_count==Nineq)exit
    enddo
    iineq=lat2ineq(1)
    ineq2lat(1)=iineq
    !close(unit) ! do not close unit if readin info below
    !
    ! allocate & fill sign list of symmetry-breaking field
    allocate(sb_field_sign(Nineq))
    sign  = next
    isign = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          sign=trim(sign)//next
       else
          read(sign,"(I6)")isign
          sign=""
          i=i+1
          if(i>Nineq)stop "build_Hij error: lattice index > Nineq read from file"
          sb_field_sign(i)=isign
       endif
       if(EOF<0)exit
    enddo
    close(unit)
    !
    ! allocate and initialize H(r_i,r_j)
    allocate(Hij(Nlso,Nlso,Nk))
    Hij = zero 
    unit = free_unit()
    open(unit,file=trim(file(2)),status='old')
    do !while(EOF>=0)
       read(unit,*,IOSTAT=EOF)ilat,iorb,jlat,jorb,ret,imt
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       if(EOF<0)exit
       do ispin=1,Nspin
          is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
          js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
          Hij(is,js,1)=dcmplx(ret,imt) 
          Hij(js,is,1)=dcmplx(ret,imt) ! symmetrize hopping
       enddo
    enddo
    close(unit)
    !
    call write_hk_w90("Hij_nano.data",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Nineq,&
         Hk=Hij,&
         kxgrid=[0d0],kygrid=[0d0],kzgrid=[0d0])
    !
    allocate(nanoHloc(Nlso,Nlso))
    nanoHloc = extract_Hloc(Hij,Nlat,Nspin,Norb)
    !
    !save lat2ineq,ineq2lat arrays
    unit=free_unit()
    open(unit,file="lat2ineq.ed")
    do ilat=1,Nlat
       write(unit,*)ilat,lat2ineq(ilat)
    enddo
    close(unit)
    unit=free_unit()
    open(unit,file="ineq2lat.ed")
    do i=1,Nineq
       write(unit,*)i,ineq2lat(i)
    enddo
    close(unit)
  end subroutine build_Hij


  !----------------------------------------------------------------------------------------!
  ! purpose: save the local self-energy on disk
  !----------------------------------------------------------------------------------------!
  subroutine save_sigma(Smats,Sreal)
    complex(8),intent(inout)         :: Smats(:,:,:,:,:,:)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

    if(size(Smats,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Smats,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Smats,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Smats,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wm(Lmats))
    allocate(wr(Lreal))

    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    if(mpiID==0)then
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call store_data("LSigma"//trim(suffix),Smats(:,ispin,ispin,iorb,iorb,:),wm)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call store_data("LSigma"//trim(suffix),Sreal(:,ispin,ispin,iorb,iorb,:),wr)
          enddo
       enddo
    endif

  end subroutine save_sigma



  !----------------------------------------------------------------------------------------!
  ! purpose: read the local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma(Smats,Sreal)
    complex(8),intent(inout)         :: Smats(:,:,:,:,:,:)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

    if(size(Smats,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Smats,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Smats,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Smats,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wm(Lmats))
    allocate(wr(Lreal))

    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    if(mpiID==0)then
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call read_data("LSigma"//trim(suffix),Smats(:,ispin,ispin,iorb,iorb,:),wm)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call read_data("LSigma"//trim(suffix),Sreal(:,ispin,ispin,iorb,iorb,:),wr)
          enddo
       enddo
    endif

  end subroutine read_sigma


  !     !----------------------------------------------------------------------------------------!
  !     ! purpose: evaluate the Normal Green's functions G_ij for a given Hamiltonian matrix and
  !     ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  !     ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !     !----------------------------------------------------------------------------------------!
  !     subroutine ed_get_gij_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,Gamma_mats,Gamma_real)
  !       complex(8),dimension(:,:,:)      :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  !       real(8)                          :: Wtk(size(Hk,3)) ![Nk]
  !       complex(8),intent(inout)         :: Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
  !       complex(8),intent(inout)         :: Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
  !       complex(8),intent(inout)         :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
  !       complex(8),intent(inout)         :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
  !       complex(8)                       :: zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Lmats)
  !       complex(8)                       :: zeta_real(Nlat,Nspin*Norb,Nspin*Norb,Lreal)
  !       complex(8)                       :: Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
  !       complex(8)                       :: Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
  !       complex(8),optional              :: Gamma_mats(size(Hk,1),size(Hk,2),Lmats)
  !       complex(8),optional              :: Gamma_real(size(Hk,1),size(Hk,2),Lreal)
  !       integer                          :: i,j,ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !       real(8),dimension(:),allocatable :: wm,wr
  !       !
  !       Lk=size(Hk,3)
  !       Nlso=Nlat*Norb*Nspin
  !       if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "ed_get_gij_lattice error: wrong dimensions of Hk"
  !       !
  !       allocate(wm(Lmats))
  !       allocate(wr(Lreal))
  !       wm = pi/beta*(2*arange(1,Lmats)-1)
  !       wr = linspace(wini,wfin,Lreal)
  !       !
  !   !
  !       if(mpiID==0)write(*,*)"Get local GF (id=0):"
  !       !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
  !       !then we use a newly created function *blocks_to_matrix* to spread the blocks into
  !       !a matrix of rank 2 dimensions Nlso*Nlso
  !       !
  !       zeta_mats=zero
  !       zeta_real=zero
  !       do ilat=1,Nlat
  !          do ispin=1,Nspin
  !             do iorb=1,Norb
  !                io = iorb + (ispin-1)*Norb
  !                js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
  !                zeta_mats(ilat,io,io,:) = xi*wm(:)       + xmu
  !                zeta_real(ilat,io,io,:) = wr(:) + xi*eps + xmu
  !             enddo
  !          enddo
  !          do ispin=1,Nspin
  !             do jspin=1,Nspin
  !                do iorb=1,Norb
  !                   do jorb=1,Norb
  !                      io = iorb + (ispin-1)*Norb
  !                      jo = jorb + (jspin-1)*Norb
  !                      zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Smats(ilat,ispin,jspin,iorb,jorb,:)
  !                      zeta_real(ilat,io,jo,:) = zeta_real(ilat,io,jo,:) - Sreal(ilat,ispin,jspin,iorb,jorb,:)
  !                   enddo
  !                enddo
  !             enddo
  !          enddo
  !       enddo
  !       !
  !       !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
  !       if(mpiID==0)call start_timer
  !       Gmats=zero
  !       Greal=zero
  !       do ik=1,Lk
  !          if(present(Gamma_mats))then
  !             call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),Gkmats,Gembed=Gamma_mats)
  !          else
  !             call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),Gkmats)
  !          endif
  !          if(present(Gamma_real))then
  !             call add_to_gloc_normal(zeta_real,Hk(:,:,ik),Gkreal,Gembed=Gamma_real)
  !          else
  !             call add_to_gloc_normal(zeta_real,Hk(:,:,ik),Gkreal)
  !          endif
  !          !call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),Gkmats)      
  !          !call add_to_gloc_normal(zeta_real,Hk(:,:,ik),Gkreal)
  !          Gmats = Gmats + Gkmats*Wtk(ik)
  !          Greal = Greal + Gkreal*Wtk(ik)
  !          if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !       end do
  !       if(mpiID==0)call stop_timer
  !       if(mpiID==0)then
  !          write(LOGfile,*)"write spin-orbital diagonal elements:"
  !          ispin=1
  !          iorb=1
  !          call store_data("Gij_l1_s1_iw.ed",Gmats(:,:,ispin,ispin,iorb,iorb,:),wm)
  !          call store_data("Gij_l1_s1_realw.ed",Greal(:,:,ispin,ispin,iorb,iorb,:),wr)
  !       endif
  !     end subroutine ed_get_gij_lattice

  !     subroutine add_to_gloc_normal(zeta_site,Hk,Gkout,Gembed)
  !       complex(8)               :: zeta_site(:,:,:,:)              ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
  !       complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
  !       real(8)                  :: Wtk                    
  !       !output:
  !       complex(8),intent(inout) :: Gkout(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
  !       complex(8)               :: Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
  !       !
  !       complex(8)               :: Gmatrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
  !       complex(8),optional      :: Gembed(Nlat*Nspin*Norb,Nlat*Nspin*Norb,size(zeta_site,4)) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
  !       integer                  :: i,j,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
  !       if(size(zeta_site,1)/=Nlat)stop "get_gloc_kpoint error: zeta_site wrong size 1 = Nlat"
  !       if(size(zeta_site,2)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 2 = Nspin*Norb"
  !       if(size(zeta_site,3)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 3 = Nspin*Norb"
  !       Lfreq = size(zeta_site,4)
  !       Gktmp=zero
  !       do i=1+mpiID,Lfreq,mpiSIZE
  !          Gmatrix  = blocks_to_matrix(zeta_site(:,:,:,i)) - Hk !(z+mu)I_ij - Sigma_ij.I_ij - H(k)_ij
  !          if(present(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
  !          call matrix_inverse_sym(Gmatrix)
  !          !store the diagonal blocks directly into the tmp output 
  !          do ilat=1,Nlat
  !             do jlat=1,Nlat
  !                do ispin=1,Nspin
  !                   do jspin=1,Nspin
  !                      do iorb=1,Norb
  !                         do jorb=1,Norb
  !                            io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
  !                            jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
  !                            Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
  !                         enddo
  !                      enddo
  !                   enddo
  !                enddo
  !             enddo
  !          enddo
  !       enddo
  !       Gkout=zero
  ! #ifdef _MPI
  !       call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  ! #else
  !       Gkout = Gktmp
  ! #endif
  !     end subroutine add_to_gloc_normal


  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the conductance (without vertex corrections) for a nanostructure 
  ! on the real axis, given the non-local Green's function and the L/R hybridization matrix, 
  ! of size [Nlat*Nspin*Norb**2*Lreal]
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_conductance(Gret)
    complex(8),intent(inout)              :: Gret(:,:,:,:,:,:,:)  ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    ! auxiliary variables for matmul        
    complex(8),dimension(:,:),allocatable :: GR,HR,GA,HL,Re,Le,Te ![Nlat*Norb]**2
    complex(8),dimension(:,:),allocatable :: transe               ![Nspin][Lreal]
    !
    real(8),dimension(:),allocatable      :: jcurrs               ![Nspin]
    !
    integer,dimension(:),allocatable      :: rmask,lmask          ![Nlat]
    !
    real(8),dimension(:),allocatable      :: wr
    integer                               :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,is,js,i,Nlso,Nlo
    integer                               :: unit,lfile
    character(len=30)                     :: suffix
    !
    Nlso = Nlat*Nspin*Norb
    Nlo  = Nlat*Norb
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)

    ! allocate variables for matrix-matrix multiplication
    allocate(GR(Nlo,Nlo));GR=zero
    allocate(HR(Nlo,Nlo));HR=zero
    allocate(GA(Nlo,Nlo));GA=zero
    allocate(HL(Nlo,Nlo));HL=zero
    allocate(Re(Nlo,Nlo));Re=zero
    allocate(Le(Nlo,Nlo));Le=zero
    allocate(Te(Nlo,Nlo));Te=zero



    ! set masks
    allocate(lmask(Nlat),rmask(Nlat))
    lmask(:)=0
    rmask(:)=0
    if(mpiID==0)then
       lfile = file_length("lmask.in")
       unit = free_unit()
       open(unit,file='lmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat
          ilat=ilat+1
          lmask(ilat)=1
          write(6,*) ilat,lmask(ilat)
       enddo
       lfile = file_length("rmask.in")
       unit = free_unit()
       open(unit,file='rmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat
          ilat=ilat+1
          rmask(ilat)=1
          write(6,*) ilat,rmask(ilat)
       enddo
    endif

    ! allocate spin-resolved transmission coefficient
    allocate(transe(Nspin,Lreal))

    do ispin=1,Nspin
       do i=1,Lreal
          ! fill auxiliary matrix [Nlso]**2
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb  + (ilat-1)*Norb
                      jo = jorb  + (jlat-1)*Norb
                      is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
                      js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
                      !
                      ! retarded Green's function
                      GR(io,jo)=Gret(ilat,jlat,ispin,ispin,iorb,jorb,i)
                      ! set \Gamma matrix for L/R according to masks: {ilat,jlat} \in L-subset OR R-subset
                      HR(io,jo)=zero
                      if( (rmask(ilat)==1) .AND. (rmask(jlat)==1) )HR(io,jo) = cmplx(dimag(Hyb_real(is,js,i)),0d0)
                      HL(io,jo)=zero
                      if( (lmask(ilat)==1) .AND. (lmask(jlat)==1) )HL(io,jo) = cmplx(dimag(Hyb_real(is,js,i)),0d0)
                   enddo
                enddo
             enddo
          enddo
          ! advanced Green's function
          GA=conjg(transpose(GR))
          ! get transmission function as T(ispin,i)=Tr[Gadvc*Hybl*Gret*Hybr]
          Re = matmul(GR,HR)
          Le = matmul(GA,HL)
          Te = matmul(Le,Re)
          transe(ispin,i) = trace_matrix(Te,Nlo)
       enddo
       suffix="_s"//reg(txtfy(ispin))//"_realw.ed"
       call store_data("Te"//trim(suffix),transe(ispin,:),wr)
    enddo


    ! allocate spin-resolved current (transmission coefficient integrated)
    allocate(jcurrs(Nspin));jcurrs=0.d0

    ! evaluate spin-resolved current: 
    ! actually this formula is wrong, because in the zero-bias limit the current should be zero 
    ! what matters is the integral over the eenrgy window included 
    ! between the chemical potentials of the L/R leads: i.e., the formula should be 
    ! J = \int_{-\infty}^{\infty} de T(e)*(f_L(e)-f_R(e))
    do ispin=1,Nspin
       do i=1,Lreal
          jcurrs(ispin) = jcurrs(ispin) + transe(ispin,i)*fermi(wr(i),beta)
       enddo
    enddo
    !unit = free_unit()
    !open(unit,file="Jcurrs.ed")
    !do ispin=1,Nspin
    !   write(unit,'(i3,1f16.9)')ispin,Jcurrs(ispin)
    !enddo
    !close(unit)



    deallocate(GR,HR,GA,HL,rmask,lmask,Re,Le,Te,jcurrs) 

  end subroutine ed_get_conductance


  !----------------------------------------------------------------------------------------!
  ! purpose: define the hybridization matrix of size [Nlat][Nlat][Nspin][Norb][Norb][Lreal] 
  ! reading the parameters from an input file
  !----------------------------------------------------------------------------------------!
  subroutine set_hyb()
    integer                                 :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,i,Nlso
    integer                                 :: k,kmax
    integer                                 :: unit,l,lfile
    ! leads
    integer                                 :: ikind,ilead,Nlead
    real(8)                                 :: D,mu,V,epsk
    complex(8)                              :: ksum
    complex(8),dimension(:,:,:),allocatable :: lead_real,lead_mats ![Nlead][Nspin][Lreal/Lmats]
    real(8),dimension(:),allocatable        :: wr,wm
    character(50)                           :: suffix
    !
    Nlso = Nlat*Nspin*Norb
    !
    kmax=10000
    !
    allocate(wm(Lmats),wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)


    ! initialize embedding hybridization function
    allocate(Hyb_mats(Nlso,Nlso,Lmats))
    allocate(Hyb_real(Nlso,Nlso,Lreal))
    Hyb_mats=zero
    Hyb_real=zero

    ! determine Nleads & allocate lead matrix
    if(mpiID==0)then
       lfile = file_length("lead.in")
       unit = free_unit()
       open(unit,file='lead.in',status='old')
       read(unit,*)Nlead
       allocate(lead_real(Nlead,Nspin,Lreal))
       allocate(lead_mats(Nlead,Nspin,Lmats))
       lead_real(:,:,:)=zero
       ! lead file setup lead by kind, half-bandwitdh (D) and chemical potential (mu)
       do l=1,lfile-1 ! because Nlead was read separately above
          read(unit,*) ilead, ispin, D, mu, ikind
          ilead=ilead+1
          ispin=ispin+1
          if(ilead>Nlead)stop "set_hyb error: in input file 'lead.in' ilead > Nlead"
          if(ispin>Nspin)stop "set_hyb error: non-spin degenerate leads for Nspin=1 calculation"
          suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))
          !
          ! set the lead's Green's function, depending on ikind
          if(ikind==0)then
             ! flat DOS (analytic)
             write(*,*) "flat DOS (analytic)"
             lead_real(ilead,ispin,:)=dcmplx( log(abs((D+wr(:)+mu)/(D-wr(:)-mu))) , -pi*heaviside(D-abs(wr(:)+mu)) )/(2d0*D)
          elseif(ikind==1)then
             ! flat DOS (k-sum)
             write(*,*) "flat DOS (k-sum)"
             do i=1,Lreal
                ksum=zero
                do k=1,kmax
                   epsk = -D + 2*D/kmax*(k-1)
                   ksum = ksum + 1d0/( wr(i)+xi*0.01d0+mu - epsk)
                enddo
                lead_real(ilead,ispin,i)=ksum/kmax
             enddo
          elseif(ikind==2)then
             ! broad-band limit
             write(*,*) "broad-band limit (analytic)" 
             lead_real(ilead,ispin,:)=dcmplx(0d0,-1d0) ! not very elegant...
          elseif(ikind==3)then
             ! semicircular DOS (k-sum) 
             write(*,*) "semicircular DOS (k-sum)"
             do i=1,Lreal
                ksum=zero
                do k=1,kmax
                   epsk = -D + 2*D/kmax*(k-1)
                   ksum = ksum + (4d0/(pi*kmax))*sqrt(1d0-(epsk/D)**2)/( wr(i)+xi*0.01d0+mu - epsk)
                enddo
                lead_real(ilead,ispin,i)=ksum
             enddo
          elseif(ikind==4)then
             ! readin hk DOS
             write(*,*) "readin hk DOS to be implemented and benchmarked w/ w2dynamics"
             stop
          else
             write(*,*) "set_hyb error: in input file 'lead.in' invalid ikind"
             stop
          endif
          ! store lead(s) on disk
          suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call store_data("lead"//trim(suffix),lead_real(ilead,ispin,:),wr)
          call get_matsubara_gf_from_dos(wr,lead_real(ilead,ispin,:),lead_mats(ilead,ispin,:),beta)
          suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_iw.ed"
          call store_data("lead"//trim(suffix),lead_mats(ilead,ispin,:),wm)
       enddo
       close(unit)
       !
       ! hybridization file determine lead-site connections 
       lfile = file_length("vij.in")
       unit = free_unit()
       open(unit,file='vij.in',status='old')
       do i=1,lfile
          read(unit,*) ilat, iorb, jlat, jorb, ilead, V
          ilat=ilat+1
          iorb=iorb+1
          jlat=jlat+1
          jorb=jorb+1
          ilead=ilead+1
          if((iorb>Norb).or.(jorb>Norb))stop "set_hyb error: in input file 'vij.in' i/jorb > Norb"
          if((ilat>Nlat).or.(jlat>Nlat))stop "set_hyb error: in input file 'vij.in' i/jlat > Nlat"
          if(ilead>Nlead)stop "set_hyb error: in input file 'vij.in' ilead > Nlead"
          do ispin=1,Nspin
             io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
             jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
             Hyb_real(io,jo,:)=Hyb_real(io,jo,:)+lead_real(ilead,ispin,:)*V**2
             Hyb_mats(io,jo,:)=Hyb_mats(io,jo,:)+lead_mats(ilead,ispin,:)*V**2
             io = iorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
             jo = jorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
             Hyb_real(io,jo,:)=Hyb_real(io,jo,:)+lead_real(ilead,ispin,:)*V**2
             Hyb_mats(io,jo,:)=Hyb_mats(io,jo,:)+lead_mats(ilead,ispin,:)*V**2
             io = jorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
             jo = iorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
             Hyb_real(io,jo,:)=Hyb_real(io,jo,:)+lead_real(ilead,ispin,:)*V**2
             Hyb_mats(io,jo,:)=Hyb_mats(io,jo,:)+lead_mats(ilead,ispin,:)*V**2
             io = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
             jo = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
             Hyb_real(io,jo,:)=Hyb_real(io,jo,:)+lead_real(ilead,ispin,:)*V**2
             Hyb_mats(io,jo,:)=Hyb_mats(io,jo,:)+lead_mats(ilead,ispin,:)*V**2
             suffix="_i"//reg(txtfy(ilat))//"_j"//reg(txtfy(jlat))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call store_data("Hyb"//trim(suffix),Hyb_real(io,jo,:),wr)
          enddo
       enddo
       close(unit)
    endif
    deallocate(lead_real,lead_mats,wr,wm)
    !
    call MPI_Bcast(hyb_real,size(hyb_real),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_Bcast(hyb_mats,size(hyb_mats),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    !
  end subroutine set_hyb


  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the effective exchange as in Katsnelson PRB 61, 8906 (2000), eq. (21)
  ! given the non-local Green's function, the local (auxiliary) self-energy S_i = (S_iup-S_ido)/2 
  ! and the fermi distribution on the real axis. 
  ! Jeff_ij = 1/pi Im \int_{-infty}^{infty} S_i(w) G_ijup(w) S_j(w) G_ijdo(w) f(w) dw
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_jeff(Gret,Sret)
    complex(8),intent(inout)                  :: Gret(:,:,:,:,:,:,:) ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout)                  :: Sret(:,:,:,:,:,:)   ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:),allocatable :: Saux(:,:,:,:)       ![Nlat][Norb][Norb][Lreal]
    complex(8)                                :: kernel
    real(8),dimension(:,:),allocatable        :: jeff(:,:)           ![Nlat][Nlat]
    real(8),dimension(:),allocatable          :: wr
    integer                                   :: ilat,jlat,iorb,jorb,i
    !
    !I/O
    integer                                   :: unit
    character(len=30)                         :: suffix
    !

    ! check inouts dimensions
    call assert_shape(Gret,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_jeff","Gret")
    call assert_shape(Sret,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_jeff","Sret")

    allocate(Saux(Nlat,Norb,Norb,Lreal))
    Saux(:,:,:,:)=zero
    !
    allocate(jeff(Nlat,Nlat))
    jeff(:,:)=zero
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    write(*,*) "computing effective non-local exchange"
    !
    ! sanity checks
    if(Nspin/=2)stop "ed_get_jeff error: Nspin /= 2"
    if(Norb>1)stop "ed_get_jeff error: Norb > 1 (mutli-orbital case no timplmented yet)"
    !
    ! define auxiliary local spin-less self-energy
    do ilat=1,Nlat
       Saux(ilat,1,1,:) = (Sret(ilat,1,1,1,1,:)-Sret(ilat,2,2,1,1,:))/2.d0
       !Saux(ilat,1,1,:) = one
    enddo
    !unit = free_unit()
    !open(unit,file="Saux.ed")
    !do ilat=1,Nlat
    !   do i=1,Lreal
    !      write(unit,'(i5,7f16.9)')ilat,wr(i),Saux(ilat,1,1,i),Sret(ilat,1,1,1,1,i),Sret(ilat,2,2,1,1,i)
    !   enddo
    !enddo
    !close(unit)
    !
    ! compute effective exchange
    do ilat=1,Nlat
       do jlat=1,Nlat
          ! perform integral over frequency
          kernel=0.d0
          do i=1,Lreal
             ! jeff kernel: non-local Green's function and fermi function convolution
             !              in the multi-orbital case: trace over the orbitals required
             kernel = kernel + Saux(ilat,1,1,i)*Gret(ilat,jlat,1,1,1,1,i)*Saux(jlat,1,1,i)*Gret(jlat,ilat,2,2,1,1,i)*fermi(wr(i),beta)
          enddo
          jeff(ilat,jlat) = 1.d0*dimag(kernel)/pi
       enddo
    enddo
    !
    ! write effective exchange on disk
    unit = free_unit()
    open(unit,file="jeff_ij.ed")
    do ilat=1,Nlat
       do jlat=1,Nlat
          write(unit,*)ilat,jlat,jeff(ilat,jlat)
       enddo
    enddo
    close(unit)

    deallocate(Saux,jeff,wr) 

  end subroutine ed_get_jeff


  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the non-local bare static spin susceptibility
  ! given the non-local Green's function and the fermi distribution on the real axis. 
  ! chi0_ij = 1/pi Im \int_{-infty}^{infty} G_ij(w) G_ji(w) f(w) dw
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_chi0ij(Gret)
    complex(8),intent(inout)                  :: Gret(:,:,:,:,:,:,:) ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8)                                :: kernel
    real(8),dimension(:,:),allocatable        :: jeff(:,:)           ![Nlat][Nlat]
    real(8),dimension(:),allocatable          :: wr
    integer                                   :: ilat,jlat,iorb,jorb,i
    !
    !I/O
    integer                                   :: unit
    character(len=30)                         :: suffix
    !

    ! check inouts dimensions
    call assert_shape(Gret,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_jeff","Gret")

    allocate(jeff(Nlat,Nlat))
    jeff(:,:)=zero
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    write(*,*) "computing bare static non-local susceptibility"
    !
    ! sanity checks
    if(Nspin/=1)stop "ed_get_chi0ij error: Nspin /= 1"
    if(Norb>1)stop "ed_get_chi0ij error: Norb > 1 (mutli-orbital case no timplmented yet)"
    !
    ! compute bare static non-local susceptibility
    do ilat=1,Nlat
       do jlat=1,Nlat
          ! perform integral over frequency
          kernel=0.d0
          do i=1,Lreal
             ! jeff kernel: non-local Green's function and fermi function convolution
             !              in the multi-orbital case: trace over the orbitals required
             kernel = kernel + Gret(ilat,jlat,1,1,1,1,i)*Gret(jlat,ilat,1,1,1,1,i)*fermi(wr(i),beta)
          enddo
          jeff(ilat,jlat) = 1.d0*dimag(kernel)/pi
       enddo
    enddo
    !
    ! write bare static non-local susceptibility on disk
    unit = free_unit()
    open(unit,file="jeff_ij.ed")
    do ilat=1,Nlat
       do jlat=1,Nlat
          write(unit,*)ilat,jlat,jeff(ilat,jlat)
       enddo
    enddo
    close(unit)

    deallocate(jeff,wr) 

  end subroutine ed_get_chi0ij



end program ed_nano
