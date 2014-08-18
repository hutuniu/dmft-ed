module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_ENERGY
  USE ED_DIAG
  implicit none
  private

  interface ed_kinetic_energy
     module procedure &
          ed_kinetic_energy_d1,ed_kinetic_energy_d1_,&
          ed_kinetic_energy_dm,ed_kinetic_energy_dm_,&
          ed_kinetic_energy_c1,ed_kinetic_energy_c1_,&
          ed_kinetic_energy_cm,ed_kinetic_energy_cm_
  end interface ed_kinetic_energy

  public :: init_ed_solver
  public :: ed_solver
  public :: ed_kinetic_energy

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine init_ed_solver(bath_,hwband,Hunit)
    real(8),dimension(:,:),intent(inout) :: bath_
    real(8),optional,intent(in)          :: hwband
    real(8)                              :: hwband_
    character(len=*),optional,intent(in) :: Hunit
    character(len=64)                    :: Hunit_
    logical                              :: check 
    logical,save                         :: isetup=.true.
    hwband_=2.d0;if(present(hwband))hwband_=hwband
    Hunit_='inputHLOC.in';if(present(Hunit))Hunit_=Hunit
    if(ed_verbose<2.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    if(isetup)call init_ed_structure(Hunit_)
    bath_ = 0.d0
    check = check_bath_dimension(bath_)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_bath(dmft_bath)
    call init_bath_ed(dmft_bath,hwband_)
    call copy_bath(dmft_bath,bath_)
    if(isetup)then
       if(.not.ed_supercond)then
          call setup_pointers
       else
          call setup_pointers_sc
       endif
    endif
    call deallocate_bath(dmft_bath)
    isetup=.false.
  end subroutine init_ed_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solver(bath_)
    real(8),dimension(:,:),intent(in) :: bath_
    integer                           :: unit
    logical                           :: check
    check = check_bath_dimension(bath_)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_bath(dmft_bath)
    call set_bath(bath_,dmft_bath)
    if(ED_MPI_ID==0)then
       if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
       call save_bath(dmft_bath,file=trim(Hfile)//trim(ed_file_suffix))
    endif
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity         !find target states by digonalization of Hamiltonian
    call buildgf_impurity             !build the one-particle impurity Green's functions
    if(chiflag)call buildchi_impurity !build the local susceptibilities (spin [todo charge])
    call observables_impurity         !obtain impurity observables as thermal averages.  
    call local_energy_impurity        !obtain the local energy of the effective impurity problem.
    !
    call deallocate_bath(dmft_bath)   
    call es_delete_espace(state_list) 
  end subroutine ed_solver




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_kinetic_energy_d1(Sigma,Hk)
    complex(8),dimension(:)               :: Sigma
    real(8),dimension(:)                  :: Hk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8),dimension(size(Hk))           :: Wtk_
    Wtk_          = 1d0/dble(size(Hk))
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call kinetic_energy_impurity(Hk_,Wtk_,Sigma_)
  end subroutine ed_kinetic_energy_d1
  !
  subroutine ed_kinetic_energy_dm(Sigma,Hk)
    complex(8),dimension(:,:,:)   :: Sigma
    real(8),dimension(:,:,:)      :: Hk
    real(8),dimension(size(Hk,3)) :: Wtk_
    Wtk_          = 1d0/dble(size(Hk))
    call kinetic_energy_impurity(dcmplx(1d0,0d0)*Hk,Wtk_,Sigma)
  end subroutine ed_kinetic_energy_dm
  !
  subroutine ed_kinetic_energy_d1_(Sigma,Hk,Wtk)
    complex(8),dimension(:)               :: Sigma
    real(8),dimension(:)                  :: Hk
    real(8),dimension(:)                  :: Wtk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call kinetic_energy_impurity(Hk_,Wtk,Sigma_)
  end subroutine ed_kinetic_energy_d1_
  !
  subroutine ed_kinetic_energy_dm_(Sigma,Hk,Wtk)
    complex(8),dimension(:,:,:) :: Sigma
    real(8),dimension(:,:,:)    :: Hk
    real(8),dimension(:)        :: Wtk
    call kinetic_energy_impurity(dcmplx(1d0,0d0)*Hk,Wtk,Sigma)
  end subroutine ed_kinetic_energy_dm_
  !
  subroutine ed_kinetic_energy_c1(Sigma,Hk)
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(:)               :: Hk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8),dimension(size(Hk))           :: Wtk_
    Wtk_          = 1d0/dble(size(Hk))
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call kinetic_energy_impurity(Hk_,Wtk_,Sigma_)
  end subroutine ed_kinetic_energy_c1
  !
  subroutine ed_kinetic_energy_cm(Sigma,Hk)
    complex(8),dimension(:,:,:)   :: Sigma
    complex(8),dimension(:,:,:)   :: Hk
    real(8),dimension(size(Hk,3)) :: Wtk_
    Wtk_          = 1d0/dble(size(Hk))
    call kinetic_energy_impurity(Hk,Wtk_,Sigma)
  end subroutine ed_kinetic_energy_cm
  !
  subroutine ed_kinetic_energy_c1_(Sigma,Hk,Wtk)
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(:)               :: Hk
    real(8),dimension(:)                  :: Wtk
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call kinetic_energy_impurity(Hk_,Wtk,Sigma_)
  end subroutine ed_kinetic_energy_c1_
  !
  subroutine ed_kinetic_energy_cm_(Sigma,Hk,Wtk)
    complex(8),dimension(:,:,:) :: Sigma
    complex(8),dimension(:,:,:) :: Hk
    real(8),dimension(:)        :: Wtk
    call kinetic_energy_impurity(Hk,Wtk,Sigma)
  end subroutine ed_kinetic_energy_cm_

end module ED_MAIN



! MODULE DMFT_ED
!   USE ED_INPUT_VARS
!   USE ED_VARS_GLOBAL, only: Hloc,&
!        impSmats,impSAmats,impSreal,impSAreal,&
!        ed_dens,ed_docc,ed_phisc,ed_Ekin,ed_Epot,ed_Ehartree
!   USE ED_AUX_FUNX, only: print_Hloc,search_chemical_potential
!   USE ED_BATH, only: get_bath_size,spin_symmetrize_bath,ph_symmetrize_bath,break_symmetry_bath
!   USE ED_MAIN, only: init_ed_solver,ed_solver,ed_kinetic_energy
!   USE ED_CHI2FIT
! END MODULE DMFT_ED
