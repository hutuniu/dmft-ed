module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH_DMFT
  USE ED_BATH_USER
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_ENERGY
  USE ED_DIAG
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: reg,store_data,txtfy,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  implicit none
  private

  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_1
     module procedure ed_get_sigma_matsubara_2
     module procedure ed_get_sigma_matsubara_3
  end interface ed_get_sigma_matsubara

  interface ed_get_self_matsubara
     module procedure ed_get_self_matsubara_1
     module procedure ed_get_self_matsubara_2
     module procedure ed_get_self_matsubara_3
  end interface ed_get_self_matsubara

  interface ed_get_sigma_real
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
     module procedure ed_get_sigma_real_3
  end interface ed_get_sigma_real

  interface ed_get_self_real
     module procedure ed_get_self_real_1
     module procedure ed_get_self_real_2
     module procedure ed_get_self_real_3
  end interface ed_get_self_real

  interface ed_get_sigma_matsubara_lattice
     module procedure ed_get_sigma_matsubara_lattice_1
     module procedure ed_get_sigma_matsubara_lattice_2
     module procedure ed_get_sigma_matsubara_lattice_3
     module procedure ed_get_sigma_matsubara_lattice_11
     module procedure ed_get_sigma_matsubara_lattice_21
     module procedure ed_get_sigma_matsubara_lattice_31
  end interface ed_get_sigma_matsubara_lattice

  interface ed_get_self_matsubara_lattice
     module procedure ed_get_self_matsubara_lattice_1
     module procedure ed_get_self_matsubara_lattice_2
     module procedure ed_get_self_matsubara_lattice_3
     module procedure ed_get_self_matsubara_lattice_11
     module procedure ed_get_self_matsubara_lattice_21
     module procedure ed_get_self_matsubara_lattice_31
  end interface ed_get_self_matsubara_lattice

  interface ed_get_sigma_real_lattice
     module procedure ed_get_sigma_real_lattice_1
     module procedure ed_get_sigma_real_lattice_2
     module procedure ed_get_sigma_real_lattice_3
     module procedure ed_get_sigma_real_lattice_11
     module procedure ed_get_sigma_real_lattice_21
     module procedure ed_get_sigma_real_lattice_31
  end interface ed_get_sigma_real_lattice

  interface ed_get_self_real_lattice
     module procedure ed_get_self_real_lattice_1
     module procedure ed_get_self_real_lattice_2
     module procedure ed_get_self_real_lattice_3
     module procedure ed_get_self_real_lattice_11
     module procedure ed_get_self_real_lattice_21
     module procedure ed_get_self_real_lattice_31
  end interface ed_get_self_real_lattice



  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
     module procedure ed_get_gimp_matsubara_3
  end interface ed_get_gimp_matsubara

  interface ed_get_fimp_matsubara
     module procedure ed_get_fimp_matsubara_1
     module procedure ed_get_fimp_matsubara_2
     module procedure ed_get_fimp_matsubara_3
  end interface ed_get_fimp_matsubara

  interface ed_get_gimp_real
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
     module procedure ed_get_gimp_real_3
  end interface ed_get_gimp_real

  interface ed_get_fimp_real
     module procedure ed_get_fimp_real_1
     module procedure ed_get_fimp_real_2
     module procedure ed_get_fimp_real_3
  end interface ed_get_fimp_real

  interface ed_get_gimp_matsubara_lattice
     module procedure ed_get_gimp_matsubara_lattice_1
     module procedure ed_get_gimp_matsubara_lattice_2
     module procedure ed_get_gimp_matsubara_lattice_3
     module procedure ed_get_gimp_matsubara_lattice_11
     module procedure ed_get_gimp_matsubara_lattice_21
     module procedure ed_get_gimp_matsubara_lattice_31
  end interface ed_get_gimp_matsubara_lattice

  interface ed_get_fimp_matsubara_lattice
     module procedure ed_get_fimp_matsubara_lattice_1
     module procedure ed_get_fimp_matsubara_lattice_2
     module procedure ed_get_fimp_matsubara_lattice_3
     module procedure ed_get_fimp_matsubara_lattice_11
     module procedure ed_get_fimp_matsubara_lattice_21
     module procedure ed_get_fimp_matsubara_lattice_31
  end interface ed_get_fimp_matsubara_lattice

  interface ed_get_gimp_real_lattice
     module procedure ed_get_gimp_real_lattice_1
     module procedure ed_get_gimp_real_lattice_2
     module procedure ed_get_gimp_real_lattice_3
     module procedure ed_get_gimp_real_lattice_11
     module procedure ed_get_gimp_real_lattice_21
     module procedure ed_get_gimp_real_lattice_31
  end interface ed_get_gimp_real_lattice

  interface ed_get_fimp_real_lattice
     module procedure ed_get_fimp_real_lattice_1
     module procedure ed_get_fimp_real_lattice_2
     module procedure ed_get_fimp_real_lattice_3
     module procedure ed_get_fimp_real_lattice_11
     module procedure ed_get_fimp_real_lattice_21
     module procedure ed_get_fimp_real_lattice_31
  end interface ed_get_fimp_real_lattice



  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
  end interface ed_get_docc

  interface ed_get_phisc
     module procedure ed_get_phisc_1
     module procedure ed_get_phisc_2
  end interface ed_get_phisc

  interface ed_get_dens_lattice
     module procedure ed_get_dens_lattice_1
     module procedure ed_get_dens_lattice_2
  end interface ed_get_dens_lattice

  interface ed_get_mag_lattice
     module procedure ed_get_mag_lattice_1
     module procedure ed_get_mag_lattice_2
  end interface ed_get_mag_lattice

  interface ed_get_docc_lattice
     module procedure ed_get_docc_lattice_1
     module procedure ed_get_docc_lattice_2
  end interface ed_get_docc_lattice

  interface ed_get_phisc_lattice
     module procedure ed_get_phisc_lattice_1
     module procedure ed_get_phisc_lattice_2
  end interface ed_get_phisc_lattice

  public :: ed_init_solver
  public :: ed_init_solver_lattice
  public :: ed_solve
  public :: ed_solve_lattice
  !
  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_real
  public :: ed_get_self_real
  !
  public :: ed_get_sigma_matsubara_lattice
  public :: ed_get_self_matsubara_lattice
  public :: ed_get_sigma_real_lattice
  public :: ed_get_self_real_lattice
  !
  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_real
  public :: ed_get_fimp_real
  !
  public :: ed_get_gimp_matsubara_lattice
  public :: ed_get_fimp_matsubara_lattice
  public :: ed_get_gimp_real_lattice
  public :: ed_get_fimp_real_lattice
  !
  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc
  !
  public :: ed_get_dens_lattice
  public :: ed_get_mag_lattice
  public :: ed_get_docc_lattice
  public :: ed_get_phisc_lattice
  !
  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot
  !
  public :: ed_get_eimp_lattice
  public :: ed_get_epot_lattice
  public :: ed_get_eint_lattice
  public :: ed_get_ehartree_lattice
  public :: ed_get_eknot_lattice
  !
  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph
  !
  public :: ed_get_doubles_lattice
  public :: ed_get_dust_lattice
  public :: ed_get_dund_lattice
  public :: ed_get_dse_lattice
  public :: ed_get_dph_lattice
  !
  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators



  real(8),dimension(:,:),allocatable,save            :: nii,dii,mii,pii,ddii,eii ![Nlat][Norb/4]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: SAmatsii,SArealii        ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Fmatsii,Frealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  integer,allocatable,dimension(:,:)                 :: neigen_sectorii          ![Nlat][Nsectors]
  integer,allocatable,dimension(:)                   :: neigen_totalii           ![Nlat]
  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix

contains

  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver(bath_,himp_,hwband,Hunit)
    real(8),dimension(:),intent(inout)   :: bath_
    real(8),optional,intent(in)          :: hwband
    complex(8),allocatable,optional,intent(in)     :: himp_(:,:,:,:)
    real(8)                              :: hwband_
    complex(8)                           :: himp(Nspin,Nspin,Norb,Norb)
    character(len=*),optional,intent(in) :: Hunit
    character(len=64)                    :: Hunit_
    logical                              :: check 
    logical,save                         :: isetup=.true.
    integer :: i

    hwband_=2.d0;if(present(hwband))hwband_=hwband
    Hunit_='inputHLOC.in';if(present(Hunit))Hunit_=Hunit
    if(ed_verbose<2.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    himp=zero
    if(present(himp_))then
       check = check_size_bath(bath_,himp_)
       himp=himp_
    else
       check = check_size_bath(bath_)
    endif
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    bath_ = 0.d0

    !qui alloco le gf, impHloc=0 e provo a leggerla da file , setto le dimensioni dei blocchi
    if(isetup)call init_ed_structure(Hunit_)
    call set_hloc(himp)

    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call init_dmft_bath(dmft_bath,hwband_)
    call get_dmft_bath(dmft_bath,bath_)
    !
    !call write_dmft_bath(dmft_bath,LOGfile)
    !write(*,*) bath_
    !
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
  end subroutine ed_init_solver

  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:)             :: bath ![Nlat][:]
    integer                            :: ilat,Nineq,Nsect
    logical                            :: check_dim
    character(len=5)                   :: tmp_suffix
    Nineq = size(bath,1)
    if(Nineq > Nlat)stop "init_lattice_bath error: size[bath,1] > Nlat"
    call setup_ed_dimensions() ! < Nsectors
    if(allocated(neigen_sectorii))deallocate(neigen_sectorii)
    if(allocated(neigen_totalii))deallocate(neigen_totalii)
    allocate(neigen_sectorii(Nineq,Nsectors))
    allocate(neigen_totalii(Nineq))
    do ilat=1,Nineq
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       call ed_init_solver(bath(ilat,:))
       neigen_sectorii(ilat,:) = neigen_sector(:)
       neigen_totalii(ilat)    = lanc_nstates_total
    end do
#ifdef _MPI_INEQ
    call MPI_Barrier(MPI_COMM_WORLD,mpiERR)
#endif
  end subroutine ed_init_solver_lattice







  !+------------------------------------------------------------------+
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !-------------------------------------------------------------------------------------------
  subroutine ed_solve(bath_)
    real(8),dimension(:),intent(in) :: bath_
    logical                         :: check
    check = check_bath_dimension(bath_,impHloc)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    if(ED_MPI_ID==0)then
       if(ed_verbose<2)call write_dmft_bath(dmft_bath,LOGfile)
       call save_dmft_bath(dmft_bath,used=.true.)
    endif
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity         !find target states by digonalization of Hamiltonian
    call buildgf_impurity             !build the one-particle impurity Green's functions
    if(chiflag)call buildchi_impurity !build the local susceptibilities (spin [todo charge])
    call observables_impurity         !obtain impurity observables as thermal averages.  
    call local_energy_impurity        !obtain the local energy of the effective impurity problem.
    !
    call deallocate_dmft_bath(dmft_bath)   
    call es_delete_espace(state_list) 
  end subroutine ed_solve

  subroutine ed_solve_lattice(bath,Hloc,iprint,Uloc_ii,Ust_ii,Jh_ii)
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
    integer          :: ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    ! Check dimensions !
    Nsites=size(bath,1)
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
    if(size(neigen_sectorii,1)/=Nsites)stop "ed_solve_lattice error: size(neigen_sectorii,1)!=Nsites"
    if(size(neigen_totalii)/=Nsites)stop "ed_solve_lattice error: size(neigen_totalii,1)!=Nsites"
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+mpiID,Nsites,mpiSIZE
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
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       ! 
       !Solve the impurity problem for the ilat-th site\
       neigen_sector(:)   = neigen_sectorii(ilat,:)
       lanc_nstates_total = neigen_totalii(ilat)
       call ed_solve(bath(ilat,:))
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
    if(mpiID==0)call stop_timer
#ifdef _MPI_INEQ
    neigen_sectorii=0
    neigen_totalii =0
    call MPI_ALLREDUCE(neigen_sectortmp,neigen_sectorii,Nsites*Nsectors,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(neigen_totaltmp,neigen_totalii,Nsites,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Fmats_tmp,Fmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Freal_tmp,Frealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
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
    if(ed_verbose>4)then
       if(mpiID==0)then
          if(allocated(wm))deallocate(wm)
          if(allocated(wr))deallocate(wr)
          allocate(wm(Lmats))
          allocate(wr(Lreal))
          wm = pi/beta*(2*arange(1,Lmats)-1)
          wr = linspace(wini,wfin,Lreal)
          select case(iprint)
          case (0)
             write(LOGfile,*)"Sigma not written on file."
          case(1)                  !print only diagonal elements
             write(LOGfile,*)"write spin-orbital diagonal elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,iorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,iorb,:),wr)
                   if(ed_mode=="superc")then
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,iorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,iorb,:),wr)
                   endif
                enddo
             enddo
          case(2)                  !print spin-diagonal, all orbitals 
             write(LOGfile,*)"write spin diagonal and all orbitals elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,jorb,:),wr)
                      if(ed_mode=="superc")then
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                         call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                         call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,jorb,:),wr)
                      endif
                   enddo
                enddo
             enddo
          case default                  !print all off-diagonals
             write(LOGfile,*)"write all elements:"
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                         call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                         call store_data("LSigma"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         if(ed_mode=="superc")then
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                            call store_data("LSelf"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                            call store_data("LSelf"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          end select
       endif
    endif
  end subroutine ed_solve_lattice










  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_sigma_matsubara_1(Smats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_1
  !
  subroutine ed_get_sigma_matsubara_2(Smats)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Smats(io,jo,:) = impSmats(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_matsubara_2
  !
  subroutine ed_get_sigma_matsubara_3(Smats,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lmats),intent(inout) :: Smats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Smats(:) = impSmats(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_matsubara_3

  !ANOMALous, MATSUBARA SELF-ENERGY
  subroutine ed_get_self_matsubara_1(SAmats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
    SAmats(:,:,:,:,:) = impSAmats(:,:,:,:,:)
  end subroutine ed_get_self_matsubara_1
  !
  subroutine ed_get_self_matsubara_2(SAmats)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAmats(io,jo,:) = impSAmats(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_matsubara_2
  !
  subroutine ed_get_self_matsubara_3(SAmats,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lmats),intent(inout) :: SAmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    SAmats(:) = impSAmats(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_matsubara_3


  !NORMAL, REAL SELF-ENERGY
  subroutine ed_get_sigma_real_1(Sreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
  end subroutine ed_get_sigma_real_1
  !
  subroutine ed_get_sigma_real_2(Sreal)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Sreal(io,jo,:) = impSreal(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_real_2
  !
  subroutine ed_get_sigma_real_3(Sreal,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lreal),intent(inout) :: Sreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Sreal(:) = impSreal(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_real_3

  !ANOMALous, REAL SELF-ENERGY
  subroutine ed_get_self_real_1(SAreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
    SAreal(:,:,:,:,:) = impSAreal(:,:,:,:,:)
  end subroutine ed_get_self_real_1
  !
  subroutine ed_get_self_real_2(SAreal)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAreal(io,jo,:) = impSAreal(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_real_2
  !
  subroutine ed_get_self_real_3(SAreal,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lreal),intent(inout) :: SAreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    SAreal(:) = impSAreal(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_real_3


  !LATTICE:
  !NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_sigma_matsubara_lattice_1(Smats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(1:Nsites,:,:,:,:,:) = Smatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_lattice_1
  !
  subroutine ed_get_sigma_matsubara_lattice_2(Smats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Smats(ilat,io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_matsubara_lattice_2
  !
  subroutine ed_get_sigma_matsubara_lattice_3(Smats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Smats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Smats(1:Nsites,:) = Smatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_matsubara_lattice_3
  !
  subroutine ed_get_sigma_matsubara_lattice_11(Smats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(:,:,:,:,:) = Smatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_lattice_11
  !
  subroutine ed_get_sigma_matsubara_lattice_21(Smats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Smats(io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_matsubara_lattice_21
  !
  subroutine ed_get_sigma_matsubara_lattice_31(Smats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Smats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Smats(:) = Smatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_matsubara_lattice_31



  !ANOMALous, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_self_matsubara_lattice_1(SAmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
    SAmats(1:Nsites,:,:,:,:,:) = SAmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_self_matsubara_lattice_1
  !
  subroutine ed_get_self_matsubara_lattice_2(SAmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   SAmats(ilat,io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_matsubara_lattice_2
  !
  subroutine ed_get_self_matsubara_lattice_3(SAmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: SAmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAmats(1:Nsites,:) = SAmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_matsubara_lattice_3
  !
  subroutine ed_get_self_matsubara_lattice_11(SAmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
    SAmats(:,:,:,:,:) = SAmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_self_matsubara_lattice_11
  !
  subroutine ed_get_self_matsubara_lattice_21(SAmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAmats(io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_matsubara_lattice_21
  !
  subroutine ed_get_self_matsubara_lattice_31(SAmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: SAmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAmats(:) = SAmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_matsubara_lattice_31



  !NORMAL, REAL SELF-ENEGRGY
  subroutine ed_get_sigma_real_lattice_1(Sreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(1:Nsites,:,:,:,:,:) = Srealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_sigma_real_lattice_1
  !
  subroutine ed_get_sigma_real_lattice_2(Sreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Sreal(ilat,io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_real_lattice_2
  !
  subroutine ed_get_sigma_real_lattice_3(Sreal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Sreal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Sreal(1:Nsites,:) = Srealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_real_lattice_3
  !
  subroutine ed_get_sigma_real_lattice_11(Sreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(:,:,:,:,:) = Srealii(ilat,:,:,:,:,:)
  end subroutine ed_get_sigma_real_lattice_11
  !
  subroutine ed_get_sigma_real_lattice_21(Sreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Sreal(io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_real_lattice_21
  !
  subroutine ed_get_sigma_real_lattice_31(Sreal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Sreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Sreal(:) = Srealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_real_lattice_31



  !ANOMALous, REAL SELF-ENERGY
  subroutine ed_get_self_real_lattice_1(SAreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
    SAreal(1:Nsites,:,:,:,:,:) = SArealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_self_real_lattice_1
  !
  subroutine ed_get_self_real_lattice_2(SAreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   SAreal(ilat,io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_real_lattice_2
  !
  subroutine ed_get_self_real_lattice_3(SAreal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: SAreal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAreal(1:Nsites,:) = SArealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_real_lattice_3
  !
  subroutine ed_get_self_real_lattice_11(SAreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
    SAreal(:,:,:,:,:) = SArealii(ilat,:,:,:,:,:)
  end subroutine ed_get_self_real_lattice_11
  !
  subroutine ed_get_self_real_lattice_21(SAreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAreal(io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_real_lattice_21
  !
  subroutine ed_get_self_real_lattice_31(SAreal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: SAreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAreal(:) = SArealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_real_lattice_31








  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine ed_get_gimp_matsubara_1(Gmats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_1
  !
  subroutine ed_get_gimp_matsubara_2(Gmats)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Gmats(io,jo,:) = impGmats(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_matsubara_2
  !
  subroutine ed_get_gimp_matsubara_3(Gmats,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lmats),intent(inout) :: Gmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Gmats(:) = impGmats(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_matsubara_3


  !ANOMALous, MATSUBARA GREEN'S FUNCTION
  subroutine ed_get_fimp_matsubara_1(Fmats)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
    Fmats(:,:,:,:,:) = impFmats(:,:,:,:,:)
  end subroutine ed_get_fimp_matsubara_1
  !
  subroutine ed_get_fimp_matsubara_2(Fmats)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Fmats(io,jo,:) = impFmats(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_matsubara_2
  !
  subroutine ed_get_fimp_matsubara_3(Fmats,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lmats),intent(inout) :: Fmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Fmats(:) = impFmats(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_matsubara_3



  !NORMAL, REAL GREEN'S FUNCTION
  subroutine ed_get_gimp_real_1(Greal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
  end subroutine ed_get_gimp_real_1
  !
  subroutine ed_get_gimp_real_2(Greal)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Greal(io,jo,:) = impGreal(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_real_2
  !
  subroutine ed_get_gimp_real_3(Greal,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lreal),intent(inout) :: Greal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Greal(:) = impGreal(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_real_3

  !ANOMALous, REAL GREEN'S FUNCTION
  subroutine ed_get_fimp_real_1(Freal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
    Freal(:,:,:,:,:) = impFreal(:,:,:,:,:)
  end subroutine ed_get_fimp_real_1
  !
  subroutine ed_get_fimp_real_2(Freal)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
    integer  :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Freal(io,jo,:) = impFreal(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_real_2
  !
  subroutine ed_get_fimp_real_3(Freal,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lreal),intent(inout) :: Freal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1 ; if(present(iorb))iorb_=iorb
    jorb_=1 ; if(present(jorb))jorb_=jorb
    Freal(:) = impFreal(ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_real_3


  !LATTICE
  !NORMAL, MATSUBARA GREEN'S FUNCTION
  subroutine ed_get_gimp_matsubara_lattice_1(Gmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(1:Nsites,:,:,:,:,:) = Gmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_lattice_1
  !
  subroutine ed_get_gimp_matsubara_lattice_2(Gmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gmats(ilat,io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_matsubara_lattice_2
  !
  subroutine ed_get_gimp_matsubara_lattice_3(Gmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Gmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Gmats(1:Nsites,:) = Gmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_matsubara_lattice_3
  !
  subroutine ed_get_gimp_matsubara_lattice_11(Gmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(:,:,:,:,:) = Gmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_lattice_11
  !
  subroutine ed_get_gimp_matsubara_lattice_21(Gmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Gmats(io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_matsubara_lattice_21
  !
  subroutine ed_get_gimp_matsubara_lattice_31(Gmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Gmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Gmats(:) = Gmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_matsubara_lattice_31



  !ANOMALous, MATSUBARA GREEN'S FUNCTION
  subroutine ed_get_fimp_matsubara_lattice_1(Fmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
    Fmats(1:Nsites,:,:,:,:,:) = Fmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_fimp_matsubara_lattice_1
  !
  subroutine ed_get_fimp_matsubara_lattice_2(Fmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Fmats(ilat,io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_matsubara_lattice_2
  !
  subroutine ed_get_fimp_matsubara_lattice_3(Fmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Fmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Fmats(1:Nsites,:) = Fmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_matsubara_lattice_3
  !
  subroutine ed_get_fimp_matsubara_lattice_11(Fmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
    Fmats(:,:,:,:,:) = Fmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_fimp_matsubara_lattice_11
  !
  subroutine ed_get_fimp_matsubara_lattice_21(Fmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Fmats(io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_matsubara_lattice_21
  !
  subroutine ed_get_fimp_matsubara_lattice_31(Fmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Fmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Fmats(:) = Fmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_matsubara_lattice_31



  !NORMAL, REAL GREEN'S FUNCTION
  subroutine ed_get_gimp_real_lattice_1(Greal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(1:Nsites,:,:,:,:,:) = Grealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_gimp_real_lattice_1
  !
  subroutine ed_get_gimp_real_lattice_2(Greal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Greal(ilat,io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_real_lattice_2
  !
  subroutine ed_get_gimp_real_lattice_3(Greal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Greal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Greal(1:Nsites,:) = Grealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_real_lattice_3
  !
  subroutine ed_get_gimp_real_lattice_11(Greal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(:,:,:,:,:) = Grealii(ilat,:,:,:,:,:)
  end subroutine ed_get_gimp_real_lattice_11
  !
  subroutine ed_get_gimp_real_lattice_21(Greal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Greal(io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_real_lattice_21
  !
  subroutine ed_get_gimp_real_lattice_31(Greal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Greal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Greal(:) = Grealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_real_lattice_31



  !ANOMALous, REAL GREEN'S FUNCTION
  subroutine ed_get_fimp_real_lattice_1(Freal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
    Freal(1:Nsites,:,:,:,:,:) = Frealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_fimp_real_lattice_1
  !
  subroutine ed_get_fimp_real_lattice_2(Freal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Freal(ilat,io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_real_lattice_2
  !
  subroutine ed_get_fimp_real_lattice_3(Freal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Freal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Freal(1:Nsites,:) = Frealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_real_lattice_3
  !
  subroutine ed_get_fimp_real_lattice_11(Freal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
    Freal(:,:,:,:,:) = Frealii(ilat,:,:,:,:,:)
  end subroutine ed_get_fimp_real_lattice_11
  !
  subroutine ed_get_fimp_real_lattice_21(Freal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Freal(io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_real_lattice_21
  !
  subroutine ed_get_fimp_real_lattice_31(Freal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Freal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Freal(:) = Frealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_real_lattice_31









  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  function ed_get_dens_1() result(dens)
    real(8),dimension(Norb) :: dens
    dens = ed_dens
  end function ed_get_dens_1
  function ed_get_dens_2(iorb) result(dens)
    real(8)   :: dens
    integer   :: iorb
    if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
    dens = ed_dens(iorb)
  end function ed_get_dens_2

  function ed_get_dens_lattice_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0    
    if(allocated(nii))then
       if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
       yii=nii
    end if
  end function ed_get_dens_lattice_1
  function ed_get_dens_lattice_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
    yii=0d0
    if(allocated(nii))then
       if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
       yii=nii(:,iorb)
    endif
  end function ed_get_dens_lattice_2





  function ed_get_mag_1() result(mag)
    real(8),dimension(Norb) :: mag
    mag = (ed_dens_up-ed_dens_dw)
  end function ed_get_mag_1
  function ed_get_mag_2(iorb) result(mag)
    real(8)   :: mag
    integer   :: iorb
    if(iorb>Norb)stop "ed_get_dens_up error: orbital index > N_orbital"
    mag = ed_dens_up(iorb)-ed_dens_dw(iorb)
  end function ed_get_mag_2

  function ed_get_mag_lattice_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(mii))then
       if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
       yii=mii
    endif
  end function ed_get_mag_lattice_1
  function ed_get_mag_lattice_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_mag error: orbital index > N_orbital"
    yii=0d0
    if(allocated(mii))then
       if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
       yii=mii(:,iorb)
    endif
  end function ed_get_mag_lattice_2





  function ed_get_docc_1() result(docc)
    real(8),dimension(Norb) :: docc
    docc = ed_docc
  end function ed_get_docc_1
  function ed_get_docc_2(iorb) result(docc)
    real(8)   :: docc
    integer   :: iorb
    if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
    docc = ed_docc(iorb)
  end function ed_get_docc_2

  function ed_get_docc_lattice_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(dii))then
       if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
       yii=dii
    endif
  end function ed_get_docc_lattice_1
  function ed_get_docc_lattice_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
    yii=0d0
    if(allocated(dii))then
       if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
       yii=dii(:,iorb)
    endif
  end function ed_get_docc_lattice_2




  function ed_get_phisc_1() result(phisc)
    real(8),dimension(Norb) :: phisc
    phisc = ed_phisc
  end function ed_get_phisc_1
  function ed_get_phisc_2(iorb) result(phisc)
    real(8)   :: phisc
    integer   :: iorb
    if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
    phisc = ed_phisc(iorb)
  end function ed_get_phisc_2

  function ed_get_phisc_lattice_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(pii))then
       if(Nlat>size(pii,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"   
       yii=pii
    endif
  end function ed_get_phisc_lattice_1
  function ed_get_phisc_lattice_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
    yii=0d0
    if(allocated(pii))then
       if(Nlat>size(pii,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"
       yii=pii(:,iorb)
    endif
  end function ed_get_phisc_lattice_2





  function ed_get_eimp() result(eimp)
    real(8),dimension(4) :: eimp
    eimp = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
  end function ed_get_eimp
  function ed_get_epot() result(eimp)
    real(8) :: eimp
    eimp = ed_Epot
  end function ed_get_epot
  function ed_get_eint() result(eimp)
    real(8) :: eimp
    eimp = ed_Eint
  end function ed_get_eint
  function ed_get_ehartree() result(eimp)
    real(8) :: eimp
    eimp = ed_Ehartree
  end function ed_get_ehartree
  function ed_get_eknot() result(eimp)
    real(8) :: eimp
    eimp = ed_Eknot
  end function ed_get_eknot

  function ed_get_eimp_lattice(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,4)    :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_eimp error: required N_sites > evaluated N_sites"
       yii=eii
    endif
  end function ed_get_eimp_lattice
  function ed_get_epot_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_epot error: required N_sites > evaluated N_sites"
       yii=eii(:,1)
    endif
  end function ed_get_epot_lattice
  function ed_get_eint_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_eint error: required N_sites > evaluated N_sites"
       yii=eii(:,2)
    endif
  end function ed_get_eint_lattice
  function ed_get_ehartree_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
       yii=eii(:,3)
    endif
  end function ed_get_ehartree_lattice
  function ed_get_eknot_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_knot error: required N_sites > evaluated N_sites"
       yii=eii(:,4)
    endif
  end function ed_get_eknot_lattice





  function ed_get_doubles() result(docc)
    real(8),dimension(4) :: docc
    docc = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
  end function ed_get_doubles
  function ed_get_dust() result(docc)
    real(8) :: docc
    docc = ed_Dust
  end function ed_get_dust
  function ed_get_dund() result(docc)
    real(8) :: docc
    docc = ed_Dund
  end function ed_get_dund
  function ed_get_dse() result(docc)
    real(8) :: docc
    docc = ed_Dse
  end function ed_get_dse
  function ed_get_dph() result(docc)
    real(8) :: docc
    docc = ed_Dph
  end function ed_get_dph

  function ed_get_doubles_lattice(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,4)    :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
       yii=ddii(:,:)
    endif
  end function ed_get_doubles_lattice
  function ed_get_dust_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dust error: required N_sites > evaluated N_sites"
       yii=ddii(:,1)
    endif
  end function ed_get_dust_lattice
  function ed_get_dund_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dund error: required N_sites > evaluated N_sites"
       yii=ddii(:,2)
    endif
  end function ed_get_dund_lattice
  function ed_get_dse_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dse error: required N_sites > evaluated N_sites"
       yii=ddii(:,3)
    endif
  end function ed_get_dse_lattice
  function ed_get_dph_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dph error: required N_sites > evaluated N_sites"
       yii=ddii(:,4)
    endif
  end function ed_get_dph_lattice

  !DEBUG>>
  subroutine ed_get_density_matrix(dm_,iprint,dm_eig,dm_rot)
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    integer,intent(in)                           :: iprint
    real(8),allocatable,intent(out),optional     :: dm_eig(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot(:,:)
    !internal
    integer                                      :: unit  
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    !
    if (((.not.present(dm_eig)).or.(.not.present(dm_rot))).and.(iprint/=0)) then
       write(*,*) "iprint/=0 but matrices not allocated"
       stop
    elseif (.not.allocated(imp_density_matrix)) then
       write(*,*) "for some reason imp_density_matrix is not allocated"
       stop
    endif
    !
    dm_ = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                dm_(io,jo) =  imp_density_matrix(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
    !
    unit = free_unit()
    open(unit,file="imp_density_matrix.dat",action="write",position="rewind",status='unknown')
    if(iprint==0) then
       write(unit,"(A10)")"# Re{rho}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       write(unit,"(A10)")"# Im{rho}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
       enddo
    endif
    if(iprint>=1)then
       Tr=zero;Tr=trace(dm_)
       write(*,'(A25,6F15.7)') 'test #1: Tr[rho]:    ',real(Tr),aimag(Tr)
       dm_rot=zero;dm_rot=matmul(dm_,dm_)
       Tr=zero;Tr=trace(dm_rot)
       write(*,'(A25,6F15.7)') 'test #2: Tr[rho^2]:  ',real(Tr),aimag(Tr)
       dm_eig=0.0d0;dm_rot=zero;dm_rot=dm_
       call matrix_diagonalize(dm_rot,dm_eig,'V','U')
       Tr=zero;Tr=sum(dm_eig)
       write(*,'(A25,6F15.7)') 'test #3: Tr[rot(rho)]:',Tr
    endif
    if(iprint>=2)then
       dm_eig=0.0d0;dm_rot=zero;dm_rot=dm_
       call matrix_diagonalize(dm_rot,dm_eig,'V','U')
       write(unit,"(A10)")
       write(unit,"(A10)")"# rho_tilda"
       write(unit,'(10F22.12)') dm_eig
       write(unit,"(A10)")
       write(unit,"(A100)")"# Re{theta}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       write(unit,"(A100)")"# Im{theta}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot(io,jo)),jo=1,Nspin*Norb)
       enddo
    endif
    close(unit)
  end subroutine ed_get_density_matrix


  subroutine ed_get_quantum_SOC_operators(S_,L_,j_)
    !passed
    complex(8),allocatable,optional,intent(out)  ::  S_(:,:,:)
    complex(8),allocatable,optional,intent(out)  ::  L_(:,:,:)
    complex(8),allocatable,optional,intent(out)  ::  j_(:)
    integer                                      ::  unit_
    integer                                      ::  iorb,ispin,jorb,jspin
    real(8)                                      ::  Lxsq,Lysq,Lzsq,Lsq
    real(8)                                      ::  Sxsq,Sysq,Szsq,Ssq
    real(8)                                      ::  jxsq,jysq,jzsq,jsq
    if(Norb/=3)stop"SOC_operators implemented for 3 orbitals"
    if(present(S_).and.((size(S_,dim=1)/=3).or.(size(S_,dim=2)/=3).or.(size(S_,dim=3)/=3)))stop"wrong S size (3,3,3)"
    if(present(L_).and.((size(L_,dim=1)/=3).or.(size(L_,dim=2)/=2).or.(size(L_,dim=3)/=2)))stop"wrong L size (3,2,2)"
    if(present(j_).and.(size(j_)/=3))stop"wrong j size (3)"
    !
    if(present(S_))        S_=impStot
    if(present(L_))        L_=impLtot
    if(present(j_))        j_=impj_aplha
    !
    Sxsq = (trace(impStot(1,:,:)))*conjg(trace(impStot(1,:,:)))
    Sysq = (trace(impStot(2,:,:)))*conjg(trace(impStot(2,:,:)))
    Szsq = (trace(impStot(3,:,:)))*conjg(trace(impStot(3,:,:)))
    Ssq  = Sxsq + Sysq + Szsq
    !
    Lxsq = (trace(impLtot(1,:,:)))*conjg(trace(impStot(1,:,:)))
    Lysq = (trace(impLtot(2,:,:)))*conjg(trace(impLtot(2,:,:)))
    Lzsq = (trace(impLtot(3,:,:)))*conjg(trace(impLtot(3,:,:)))
    Lsq  = Lxsq + Lysq + Lzsq
    !
    jxsq = (impj_aplha(1))*conjg(impj_aplha(1))
    jysq = (impj_aplha(2))*conjg(impj_aplha(2))
    jzsq = (impj_aplha(3))*conjg(impj_aplha(3))
    jsq  = jxsq + jysq + jzsq
    !
    unit_ = free_unit()
    open(unit=unit_,file='S_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}",&
                                "7-|Sx|^2","8-|Sy|^2","9-|Sz|^2","10-|S|^2"
    write(unit_,'(30(F20.12,1X))') real(trace(impStot(1,:,:))),aimag(trace(impStot(1,:,:))),&
                                   real(trace(impStot(2,:,:))),aimag(trace(impStot(2,:,:))),&
                                   real(trace(impStot(3,:,:))),aimag(trace(impStot(3,:,:))),&
                                   Sxsq,Sysq,Szsq,Ssq

    write(unit_,*)
    write(unit_,'(30(a20,1X))') "#Sx(orb_1)","Sx(orb_2)","Sx(orb_3)","Sy(orb_1)","Sy(orb_2)","Sy(orb_3)","Sz(orb_1)","Sz(orb_2)","Sz(orb_3)"
    write(unit_,*)
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(impStot(1,iorb,jorb)),jorb=1,Norb) &
                                     ,(real(impStot(2,iorb,jorb)),jorb=1,Norb) &
                                     ,(real(impStot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (aimag(impStot(1,iorb,jorb)),jorb=1,Norb) &
                                     ,(aimag(impStot(2,iorb,jorb)),jorb=1,Norb) &
                                     ,(aimag(impStot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    close(unit_)
    !
    unit_ = free_unit()
    open(unit=unit_,file='L_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}",&
                               "7-|Lx|^2","8-|Ly|^2","9-|Lz|^2","10-|L|^2"
    write(unit_,'(30(F20.12,1X))') real(trace(impLtot(1,:,:))),aimag(trace(impLtot(1,:,:))),&
                                   real(trace(impLtot(2,:,:))),aimag(trace(impLtot(2,:,:))),&
                                   real(trace(impLtot(3,:,:))),aimag(trace(impLtot(3,:,:))),&
                                   Lxsq,Lysq,Lzsq,Lsq
    write(unit_,*)
    write(unit_,'(30(a20,1X))') "#Lx(spin_1)","Lx(spin_2)","Ly(spin_1)","Ly(spin_2)","Lz(spin_1)","Lz(spin_2)"
    write(unit_,*)
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(impLtot(1,ispin,jspin)),jspin=1,Nspin) &
                                     ,(real(impLtot(2,ispin,jspin)),jspin=1,Nspin) &
                                     ,(real(impLtot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (aimag(impLtot(1,ispin,jspin)),jspin=1,Nspin) &
                                     ,(aimag(impLtot(2,ispin,jspin)),jspin=1,Nspin) &
                                     ,(aimag(impLtot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    close(unit_)
    !
    unit_ = free_unit()
    open(unit=unit_,file='Jz_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "#1-Re{jx}","2-Im{jx}","3-Re{jy}","4-Im{jy}","5-Re{jz}","6-Im{jz}",&
                                 "7-|jx|^2","8-|jy|^2","9-|jz|^2","10-|j|^2","11-Re{L.S}","12-Im{L.S}"
    write(unit_,'(30(F20.12,1X))') real(impj_aplha(1)),aimag(impj_aplha(1)),real(impj_aplha(2)),aimag(impj_aplha(2)),real(impj_aplha(3)),aimag(impj_aplha(3)) &
                                  ,jxsq,jysq,jzsq,jsq,real(impLdotS),aimag(impLdotS)
    close(unit_)
    !
    unit_ = free_unit()
    open(unit=unit_,file='operators_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "#1-Re{Tr[Lz]}","2-|Lz|^2","3-|L|^2","4-Re{Tr[Sz]}","5-|Sz|^2","6-|S|^2",&
                                 "7-Re{jz}","7-|jz|^2","8-|j|^2","10-Re{L.S}"
    write(unit_,'(30(F20.12,1X))') real(trace(impLtot(3,:,:))),Lzsq,Lsq,real(trace(impStot(3,:,:))),Szsq,Ssq,&
                                   real(impj_aplha(3)),jzsq,jsq,real(impLdotS)
    close(unit_)
  end subroutine ed_get_quantum_SOC_operators
  !<<DEBUG

end module ED_MAIN

