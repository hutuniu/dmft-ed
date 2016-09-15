module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_IO,only: ed_print_impSigma
  USE ED_SETUP
  USE ED_BATH
  USE ED_MATVEC
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_ENERGY
  USE ED_DIAG
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,store_data,txtfy,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  implicit none
  private


  interface ed_init_solver
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
  end interface ed_init_solver


  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
  end interface ed_solve

  interface ed_rebuild_sigma
     module procedure :: ed_rebuild_sigma_single
     module procedure :: ed_rebuild_sigma_lattice
  end interface ed_rebuild_sigma


  public :: ed_init_solver
  public :: ed_solve
  public :: ed_rebuild_sigma
  !

  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix

contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath,Hloc,hwband,Hunit,MpiComm)
    !
    integer,optional                     :: MpiComm
    !
    real(8),dimension(:),intent(inout)   :: bath
    complex(8),optional,intent(in)       :: Hloc(Nspin,Nspin,Norb,Norb)
    real(8),optional,intent(in)          :: hwband
    character(len=*),optional,intent(in) :: Hunit
    real(8)                              :: hwband_
    character(len=64)                    :: Hunit_
    logical                              :: check 
    logical,save                         :: isetup=.true.
    integer                              :: i
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    if(ed_verbose<2.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    hwband_=2.d0;if(present(hwband))hwband_=hwband
    Hunit_='inputHLOC.in';if(present(Hunit))Hunit_=Hunit
    if(present(Hloc))then
       check = check_bath_dimension(bath,Hloc)
    else
       check = check_bath_dimension(bath)
    endif
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    bath = 0d0
    if(isetup)call init_ed_structure(Hunit_)
    if(present(Hloc))call set_hloc(Hloc)
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call init_dmft_bath(dmft_bath,hwband_)
    call get_dmft_bath(dmft_bath,bath)
    !
    !DEBUG>
    !call write_dmft_bath(dmft_bath,LOGfile)
    ! if(ED_MPI_ID==0)write(LOGfile,'(20(F12.6,1X))') bath
    !stop
    !>DEBUG
    !
    if(isetup)then
       select case(ed_mode)
       case default
          call setup_pointers_normal
       case ("superc")
          call setup_pointers_superc
       case ("nonsu2")
          call setup_pointers_nonsu2
       end select
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
  end subroutine ed_init_solver_single

  subroutine ed_init_solver_lattice(bath,Hloc,hwband,Hunit,MpiComm)
    !
    integer,optional                     :: MpiComm
    !
    real(8),dimension(:,:)               :: bath ![Nlat][:]
    real(8),optional,intent(in)          :: hwband
    complex(8),optional,intent(in)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    character(len=*),optional,intent(in) :: Hunit
    real(8)                              :: hwband_
    character(len=64)                    :: Hunit_
    integer                              :: ilat,Nineq,Nsect
    logical                              :: check_dim
    character(len=5)                     :: tmp_suffix
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    !
    hwband_=2.d0;if(present(hwband))hwband_=hwband
    Hunit_='inputHLOC.in';if(present(Hunit))Hunit_=Hunit
    !
    Nineq = size(bath,1)
    if(Nineq > Nlat)stop "init_lattice_bath error: size[bath,1] > Nlat"
    call setup_ed_dimensions() ! < Nsectors this is called first to allocate the following arrays
    if(allocated(neigen_sectorii))deallocate(neigen_sectorii)
    if(allocated(neigen_totalii))deallocate(neigen_totalii)
    allocate(neigen_sectorii(Nineq,Nsectors))
    allocate(neigen_totalii(Nineq))
    do ilat=1,Nineq
       if(present(Hloc))then
          check_dim = check_bath_dimension(bath(ilat,:),Hloc(ilat,:,:,:,:))
       else
          check_dim = check_bath_dimension(bath(ilat,:))
       endif
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       if(present(Hloc))then
          call ed_init_solver(bath(ilat,:),Hloc(ilat,:,:,:,:),hwband_,Hunit_)
       else
          call ed_init_solver(bath(ilat,:),hwband=hwband_,Hunit=Hunit_)
       endif
       neigen_sectorii(ilat,:) = neigen_sector(:)
       neigen_totalii(ilat)    = lanc_nstates_total
    end do
#ifdef _MPI
    call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
#endif
    ed_file_suffix=""
  end subroutine ed_init_solver_lattice










  !+------------------------------------------------------------------+
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+------------------------------------------------------------------+
  subroutine ed_solve_single(bath_,MpiComm)
    !
    integer,optional                :: MpiComm
    !
    real(8),dimension(:),intent(in) :: bath_
    logical                         :: check
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    check = check_bath_dimension(bath_,impHloc)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    if(ED_MPI_ID==0)then
       if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
       call save_dmft_bath(dmft_bath,used=.true.)
    endif
    !ASSOCIATE THE GLOBAL PROCEDURES TO THE SPECIFIC ONES
    !AS DICTATED BY THE INPUT 
    if(associated(ed_buildh_d))nullify(ed_buildh_d)
    if(associated(ed_buildh_c))nullify(ed_buildh_c)
    select case(ed_mode)
    case ('normal')
       ed_buildh_d=>build_H_normal_d
       ed_buildh_c=>build_H_normal_c
    case default
       ed_buildh_d=>build_H_all_d
       ed_buildh_c=>build_H_all_c
    end select
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity         !find target states by digonalization of Hamiltonian
    call observables_impurity         !obtain impurity observables as thermal averages.  
    call buildgf_impurity             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity !build the local susceptibilities (spin [todo charge])
    call local_energy_impurity        !obtain the local energy of the effective impurity problem.
    call deallocate_dmft_bath(dmft_bath)   
    call es_delete_espace(state_list) 
  end subroutine ed_solve_single

  subroutine ed_solve_lattice(bath,Hloc,iprint,MpiComm,Uloc_ii,Ust_ii,Jh_ii)
    !
    integer,optional :: MpiComm
    !
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer          :: iprint
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Fmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Freal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)          :: nii_tmp(size(bath,1),Norb)
    real(8)          :: dii_tmp(size(bath,1),Norb)
    real(8)          :: mii_tmp(size(bath,1),Norb)
    real(8)          :: pii_tmp(size(bath,1),Norb)
    real(8)          :: eii_tmp(size(bath,1),4)
    real(8)          :: ddii_tmp(size(bath,1),4)
    !
    integer          :: neigen_sectortmp(size(bath,1),Nsectors)
    integer          :: neigen_totaltmp(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    integer          :: MPI_COLOR
    integer          :: MPI_COLOR_RANK
    integer          :: MPI_COLOR_SIZE
    integer          :: MPI_COLOR_COMM
    logical          :: MPI_COLOR_MASTER
    !
    integer          :: MPI_MASTERS_COLOR
    integer          :: MPI_MASTERS_RANK
    integer          :: MPI_MASTERS_SIZE
    integer          :: MPI_MASTERS_COMM
    !
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    print*,Nsites,ED_MPI_ID
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    if(allocated(Fmatsii))deallocate(Fmatsii)
    if(allocated(Frealii))deallocate(Frealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Frealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    if(size(neigen_sectorii,1)<Nsites)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nsites"
    if(size(neigen_totalii)<Nsites)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nsites"
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+ED_MPI_ID,Nsites,ED_MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero ; Smats_tmp  = zero
    Srealii  = zero ; Sreal_tmp  = zero
    SAmatsii = zero ; SAmats_tmp = zero
    SArealii = zero ; SAreal_tmp = zero
    Gmatsii  = zero ; Gmats_tmp  = zero
    Grealii  = zero ; Greal_tmp  = zero
    Fmatsii  = zero ; Fmats_tmp  = zero
    Frealii  = zero ; Freal_tmp  = zero
    nii      = 0d0  ; nii_tmp    = 0d0
    dii      = 0d0  ; dii_tmp    = 0d0
    mii      = 0d0  ; mii_tmp    = 0d0
    pii      = 0d0  ; pii_tmp    = 0d0
    eii      = 0d0  ; eii_tmp    = 0d0
    ddii     = 0d0  ; ddii_tmp   = 0d0
    !
    if(ED_MPI_MASTER)call start_timer
    if(.not.ED_MPI_MASTER)LOGfile = 800+ED_MPI_ID
    !
#ifdef _MPI
    if(MPI_Colors==0.OR.MPI_Colors>ED_MPI_SIZE)MPI_Colors=ED_MPI_SIZE
    if(ED_MPI_SIZE<2)MPI_Colors=1
    if(ED_MPI_MASTER)then
       write(LOGfile,*)"MPI_SIZE      =",ED_MPI_SIZE
       write(LOGfile,*)"MPI_COLORS    =",MPI_Colors
       write(LOGfile,*)"MPI_COLOR_SIZE=",ED_MPI_SIZE/MPI_Colors
    endif
    !
    MPI_COLOR = MPI_UNDEFINED
    MPI_COLOR = mod(ED_MPI_ID,MPI_Colors)
    !
    !Split the user provided communicator into MPI_Colors groups.
    !Each group (or color) communicate via the MPI communicator MPI_Color_Comm
    call MPI_Comm_split(ED_MPI_COMM, MPI_Color, ED_MPI_ID, MPI_Color_Comm, ED_MPI_ERR)
    MPI_COLOR_SIZE=get_Size_MPI(MPI_COLOR_COMM)
    MPI_COLOR_RANK=get_Rank_MPI(MPI_COLOR_COMM)
    MPI_COLOR_MASTER=get_Master_MPI(MPI_COLOR_COMM)
    do i=0,ED_MPI_SIZE-1
       if(ED_MPI_ID==i)write(*,*)ED_MPI_ID," is now ",MPI_COLOR_RANK," in color group: ",MPI_Color
       call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
    enddo
    call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
    if(ED_MPI_MASTER)write(*,*)""
#else
    MPI_Colors=1
    MPI_Color=0
    MPI_COLOR_SIZE=1
    MPI_COLOR_RANK=0
    MPI_COLOR_MASTER=.true.
#endif

    ! call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
    ! do i=0,ED_MPI_SIZE-1
    !    if(ED_MPI_ID==i)write(*,*)ED_MPI_ID,MPI_Color,MPI_Colors
    !    call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
    ! enddo
    ! call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)
    ! write(*,*)""

    do ilat=1+MPI_Color,Nsites,MPI_Colors
       if(MPI_COLOR_MASTER)write(*,*)"Solving site:"//reg(txtfy(ilat,Npad=4))//" by group: "//txtfy(mpi_color)
       write(*,*)reg(txtfy(ED_MPI_ID))//" solves site: "//reg(txtfy(ilat,Npad=4))//" as "//reg(txtfy(MPI_COLOR_RANK))//" in group: "//txtfy(mpi_color)

       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       !
       !Solve the impurity problem for the ilat-th site
       neigen_sector(:)   = neigen_sectorii(ilat,:)
       lanc_nstates_total = neigen_totalii(ilat)
#ifdef _MPI       
       call ed_solve_single(bath(ilat,:), MpiComm=MPI_COLOR_COMM)
       print*,ED_MPI_ID,MPI_COLOR_RANK,ed_dens(1)
#else
       call ed_solve_single(bath(ilat,:))
#endif
       neigen_sectortmp(ilat,:)   = neigen_sector(:)
       neigen_totaltmp(ilat)      = lanc_nstates_total
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Fmats_tmp(ilat,:,:,:,:,:)  = impFmats(:,:,:,:,:)
       Freal_tmp(ilat,:,:,:,:,:)  = impFreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)       = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)       = ed_docc(1:Norb)
       mii_tmp(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       pii_tmp(ilat,1:Norb)       = ed_phisc(1:Norb)
       eii_tmp(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii_tmp(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]

    enddo

    call MPI_Barrier(ED_MPI_COMM,ED_MPI_ERR)

    if(ED_MPI_MASTER)call stop_timer
    ed_file_suffix=""

#ifdef _MPI
    !Now we need to collect the results: This is done in three steps (I could not realize
    !a better strategy so far):
    !1. We split the original communicators in  masters/non-masters, i.e. the processes
    !    with color_rank 0 or >0. By construction the overall master (ED_MPI_ID=0)
    !    belongs to this group. 
    !2. We perform an AllReduce operation among these two groups. This will merge the
    !    copies on each master into a single array. The copies on the slaves are identical
    !    to that of their master so we can disregard them.
    !3. We let the overall master (ED_MPI_ID=0) Bcast the result of the reduction to
    !    all the other nodes.
    neigen_sectorii=0
    neigen_totalii =0
    !
    !
    !1. split ED_MPI_COMM into MASTER/NON-MASTERS
    if(ED_MPI_MASTER)write(LOGfile,"(A)")"Creating Masters/non-Masters group:"
    MPI_Masters_Color=1
    if(MPI_COLOR_MASTER)MPI_Masters_Color=0
    call MPI_Comm_Split(ED_MPI_COMM, MPI_Masters_Color, ED_MPI_ID, MPI_MASTERS_COMM, ED_MPI_ERR)
    do j=0,ED_MPI_SIZE-1
       MPI_MASTERS_RANK=get_Rank_MPI(MPI_MASTERS_COMM)
       if(MPI_COLOR_MASTER)then
          IF(ED_MPI_ID==j)write(*,"(A)")reg(str(ED_MPI_ID))//"(~ master in group"//reg(str(MPI_Color))//&
               ") is now "//reg(str(MPI_MASTERS_RANK))//" in masters group: "//reg(str(MPI_Masters_Color))
       else
          IF(ED_MPI_ID==j)write(*,"(A)")reg(str(ED_MPI_ID))//"(~ slave in group"//reg(str(MPI_Color))//&
               ") is now "//reg(str(MPI_MASTERS_RANK))//" in NON-masters group: "//reg(str(MPI_Masters_Color))
       endif
       call MPI_Barrier(ED_MPI_COMM, ED_MPI_ERR)
    enddo
    call MPI_Barrier(ED_MPI_COMM, ED_MPI_ERR)
    if(ED_MPI_MASTER)write(LOGfile,"(A)")""
    !
    !2. Reduce in the masters/non-masters group
    call MPI_AllReduce(neigen_sectortmp, neigen_sectorii, Nsites*Nsectors, MPI_INTEGER, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(neigen_totaltmp, neigen_totalii, Nsites, MPI_INTEGER, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Smats_tmp, Smatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Sreal_tmp, Srealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(SAmats_tmp, SAmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(SAreal_tmp, SArealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Gmats_tmp, Gmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Greal_tmp, Grealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Fmats_tmp, Fmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(Freal_tmp, Frealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(nii_tmp, nii, Nsites*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(dii_tmp, dii, Nsites*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(mii_tmp, mii, Nsites*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(pii_tmp, pii, Nsites*Norb, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(eii_tmp, eii, Nsites*4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    call MPI_AllReduce(ddii_tmp, ddii, Nsites*4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MASTERS_COMM, ED_MPI_ERR)
    !
    !3. have the overall master bcast the results
    call MPI_Bcast(neigen_sectorii, Nsites*Nsectors, MPI_INTEGER, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(neigen_totalii, Nsites, MPI_INTEGER, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Smatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Srealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(SAmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(SArealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Gmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Grealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Fmatsii, Nsites*Nspin*Nspin*Norb*Norb*Lmats, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(Frealii, Nsites*Nspin*Nspin*Norb*Norb*Lreal, MPI_DOUBLE_COMPLEX, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(nii, Nsites*Norb, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(dii, Nsites*Norb, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(mii, Nsites*Norb, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(pii, Nsites*Norb, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(eii, Nsites*4, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Bcast(ddii, Nsites*4, MPI_DOUBLE_PRECISION, 0, ED_MPI_COMM, ED_MPI_ERR)
    !
    !
    ! call MPI_ALLREDUCE(neigen_sectortmp,neigen_sectorii,Nsites*Nsectors,MPI_INTEGER,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(neigen_totaltmp,neigen_totalii,Nsites,MPI_INTEGER,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Fmats_tmp,Fmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(Freal_tmp,Frealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    ! call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    !
    !
    if(present(MpiComm))then
       ED_MPI_COMM=MpiComm
    else
       ED_MPI_COMM=MPI_COMM_WORLD
    endif
    call MPI_Barrier(ED_MPI_COMM, ED_MPI_ERR)
    call MPI_Comm_Free(MPI_MASTERS_COMM,ED_MPI_ERR)
    call MPI_Comm_Free(MPI_COLOR_COMM,ED_MPI_ERR)
#else
    neigen_sectorii=neigen_sectortmp
    neigen_totalii =neigen_totaltmp
    Smatsii  =  Smats_tmp
    Srealii  =  Sreal_tmp
    SAmatsii = SAmats_tmp
    SArealii = SAreal_tmp
    Gmatsii  = Gmats_tmp
    Grealii  = Greal_tmp
    Fmatsii  = Fmats_tmp
    Frealii  = Freal_tmp
    nii      = nii_tmp
    dii      = dii_tmp
    mii      = mii_tmp
    pii      = pii_tmp
    eii      = eii_tmp
    ddii     = ddii_tmp
#endif
    call ed_print_impSigma(iprint)
  end subroutine ed_solve_lattice











  !+------------------------------------------------------------------+
  !PURPOSE: Rebuild the Self-energy (normal mode only) 
  !+------------------------------------------------------------------+
  subroutine ed_rebuild_sigma_single(bath_,MpiComm)
    !
    integer,optional                :: MpiComm
    !
    real(8),dimension(:),intent(in) :: bath_
    logical                         :: check
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    if(ed_mode/="normal")stop "WARNING: ed_rebuild_sigma works only in normal mode"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    if(ED_MPI_MASTER)then
       if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
       call save_dmft_bath(dmft_bath,used=.true.)
    endif
    call rebuildgf_impurity             !build the one-particle impurity Green's functions & Self-energy
    call deallocate_dmft_bath(dmft_bath)   
  end subroutine ed_rebuild_sigma_single


  subroutine ed_rebuild_sigma_lattice(bath,Hloc,iprint,MpiComm)
    !
    integer,optional :: MpiComm
    !
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer          :: iprint
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    ! 
    integer          :: ilat,iorb,jorb,ispin,jspin,i
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
#ifdef _MPI
    if(present(MpiComm))then
       ED_MPI_COMM = MpiComm
    else
       ED_MPI_COMM = MPI_COMM_WORLD
    end if
    ED_MPI_ID     = get_Rank_MPI(ED_MPI_COMM)
    ED_MPI_SIZE   = get_Size_MPI(ED_MPI_COMM)
    ED_MPI_MASTER = get_Master_MPI(ED_MPI_COMM)
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Check the dimensions of the bath are ok:
    if(ED_MPI_ID==0)then
       write(LOGfile,*)"Rebuilding Sigma: have you moved .used bath files to .restart ones?! "
       call sleep(3)
    end if
    do ilat=1+ED_MPI_ID,Nsites,ED_MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero ; Smats_tmp  = zero
    Srealii  = zero ; Sreal_tmp  = zero
    SAmatsii = zero ; SAmats_tmp = zero
    SArealii = zero ; SAreal_tmp = zero
    !
    if(ED_MPI_MASTER)call start_timer
    if(.not.ED_MPI_MASTER)LOGfile = 800+ED_MPI_ID
    do ilat=1+ED_MPI_ID,Nsites,ED_MPI_SIZE
       if(ED_MPI_MASTER)write(LOGfile,*)"Solving site:"//reg(txtfy(ilat,Npad=4))
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       !
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       ! 
       !Rebuild for the ilat-th site
       call ed_rebuild_sigma(bath(ilat,:))
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
    enddo
    if(ED_MPI_MASTER)call stop_timer
    ed_file_suffix=""
#ifdef _MPI
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Smatsii  =  Smats_tmp
    Srealii  =  Sreal_tmp
    SAmatsii = SAmats_tmp
    SArealii = SAreal_tmp
#endif
    call ed_print_impSigma(iprint)

  end subroutine ed_rebuild_sigma_lattice










end module ED_MAIN

