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



  public :: ed_init_solver
  public :: ed_solve
  !
  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_real
  public :: ed_get_self_real
  !
  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_real
  public :: ed_get_fimp_real
  !
  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc
  !
  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot
  !
  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_init_solver(bath_,hwband,Hunit)
    real(8),dimension(:),intent(inout)   :: bath_
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
       select case(ed_mode)
       case default
          call setup_pointers_normal
       case ("superc")
          call setup_pointers_superc
       case ("nonsu2")
          call setup_pointers_nonsu2
       end select
    endif
    call deallocate_bath(dmft_bath)
    isetup=.false.
  end subroutine ed_init_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solve(bath_)
    real(8),dimension(:),intent(in) :: bath_
    integer                         :: unit
    logical                         :: check
    check = check_bath_dimension(bath_)
    if(.not.check)stop "init_ed_solver: wrong bath dimensions"
    call allocate_bath(dmft_bath)
    call set_bath(bath_,dmft_bath)
    if(ED_MPI_ID==0)then
       if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
       !call save_bath(dmft_bath,file=trim(Hfile)//trim(ed_file_suffix))
       call save_bath(dmft_bath,used=.true.)
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
  end subroutine ed_solve






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





end module ED_MAIN

