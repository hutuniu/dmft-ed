
module DMFT_ED
  USE IOTOOLS, only:free_unit,reg
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_CHI2FIT
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  implicit none

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
    if(ed_verbose<2)write(LOGfile,"(A)")"INIT SOLVER FOR "//reg(ed_file_suffix)
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
    if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
    unit=free_unit()
    open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".used")
    call write_bath(dmft_bath,unit)
    close(unit)
    call lanc_ed_diag
    call lanc_ed_getgf
    if(chiflag)call lanc_ed_getchi
    call ed_getobs
    call deallocate_bath(dmft_bath)
    call es_delete_espace(state_list)
  end subroutine ed_solver


end module DMFT_ED
