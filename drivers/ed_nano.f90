program ed_nano
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none

  integer                                       :: iloop
  logical                                       :: converged
  integer                                       :: ilat,ineq
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_prev(:,:),Bath_ineq(:,:)
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  real(8), allocatable,dimension(:)             :: dens,dens_ineq
  real(8), allocatable,dimension(:)             :: docc,docc_ineq
  !hamiltonian input:
  complex(8),allocatable                        :: Hij(:,:,:)  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk==1]
  complex(8),allocatable                        :: nanoHloc(:,:),Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  integer                                       :: Nk,Nlso,Nineq
  integer,dimension(:),allocatable              :: lat2ineq,ineq2lat
  integer,dimension(:),allocatable              :: sb_field_sign
  !
  real(8)                                       :: wmixing,Eout(2)
  !input files
  character(len=32)                             :: finput
  character(len=32)                             :: nfile,hijfile
  logical                                       :: phsym

#ifdef _MPI_INEQ
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

  call ed_read_input(trim(finput))

  call build_Hij([nfile,hijfile])

  !Allocate Weiss Field:
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))
  !
  allocate(dens(Nlat))
  allocate(dens_ineq(Nineq))
  !
  allocate(docc(Nlat))
  allocate(docc_ineq(Nineq))

  Hloc = reshape_Hloc(nanoHloc,Nlat,Nspin,Norb)


  !Setup solver
  Nb=get_bath_size()

  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver_lattice(Bath_ineq)

  do ineq=1,Nineq
     ilat = ineq2lat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   
     bath_prev=bath_ineq

     ! solve the impurities on each inequivalent site:
     call ed_solve_lattice(bath_ineq,Hloc_ineq)
     ! retrieve the self-energies and spread them to all lattice sites
     dens_ineq = ed_get_dens_lattice(Nineq,1)
     docc_ineq = ed_get_docc_lattice(Nineq,1)
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        dens(ilat) = dens_ineq(ineq)
        docc(ilat) = docc_ineq(ineq)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo

     ! compute the local gf:
     call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1)
     do ineq=1,Nineq
        ilat = ineq2lat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo
     ! compute the Weiss field
     call ed_get_weiss_lattice(Nineq,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)

     ! fit baths and mix result with old baths
     call ed_chi2_fitgf_lattice(bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     !call ed_chi2_fitgf_lattice(bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     if(phsym)then
        do ineq=1,Nineq
           call ph_symmetrize_bath(bath_ineq(ineq,:),save=.true.)
        enddo
     endif
     Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     if(mpiID==0)then
       converged = check_convergence(Weiss_ineq(1,1,1,1,1,:),dmft_error,nsuccess,nloop)

       if(NREAD/=0.d0) call search_chemical_potential(xmu,sum(dens)/Nlat,converged)
    endif
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(mpiID==0) call end_loop()
  end do

  Eout = ed_kinetic_energy_lattice(Hij,[1d0],Smats)

  call ed_get_gloc_lattice(Hij,[1d0],Gmats,Greal,Smats,Sreal,iprint=1)

#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif


contains

  subroutine build_Hij(file)
    character(len=*)     :: file(2)
    integer              :: ilat,jlat,iorb,jorb,is,js,ispin,ie
    integer              :: i,isite,iineq,iineq0,isign
    integer              :: EOF
    character, parameter :: tab = achar ( 9 )
    integer              :: unit,ineq_count
    integer              :: Ns,Ne,Nb         ! #atoms, #inequivalent, #bands
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
          Hij(js,is,1)=dcmplx(ret,imt) ! symmetrize hoppping
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




end program ed_nano
