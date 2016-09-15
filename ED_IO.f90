MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: reg,store_data,txtfy,free_unit,free_units
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


  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_real
  public :: ed_get_self_real

  public :: ed_get_sigma_matsubara_lattice
  public :: ed_get_self_matsubara_lattice
  public :: ed_get_sigma_real_lattice
  public :: ed_get_self_real_lattice

  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_real
  public :: ed_get_fimp_real

  public :: ed_get_gimp_matsubara_lattice
  public :: ed_get_fimp_matsubara_lattice
  public :: ed_get_gimp_real_lattice
  public :: ed_get_fimp_real_lattice

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc

  public :: ed_get_dens_lattice
  public :: ed_get_mag_lattice
  public :: ed_get_docc_lattice
  public :: ed_get_phisc_lattice

  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot

  public :: ed_get_eimp_lattice
  public :: ed_get_epot_lattice
  public :: ed_get_eint_lattice
  public :: ed_get_ehartree_lattice
  public :: ed_get_eknot_lattice

  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph
  public :: ed_get_doubles_lattice
  public :: ed_get_dust_lattice
  public :: ed_get_dund_lattice
  public :: ed_get_dse_lattice
  public :: ed_get_dph_lattice
  !
  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators


  interface ed_print_impSigma
     module procedure :: ed_print_impSigma_single
     module procedure :: ed_print_impSigma_lattice
  end interface ed_print_impSigma

  interface  ed_print_impG
     module procedure :: ed_print_impG_single
     module procedure :: ed_print_impG_lattice
  end interface ed_print_impG

  public :: ed_print_PolesWeights
  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impChi

  !FOR INTERNAL USE ONLY:
  public :: print_poles_weights_normal
  public :: print_poles_weights_superc
  public :: print_poles_weights_nonsu2
  !
  public :: print_impSigma_normal
  public :: print_impSigma_superc
  public :: print_impSigma_nonsu2
  !
  public :: print_impSigma_normal_lattice
  public :: print_impSigma_superc_lattice
  public :: print_impSigma_nonsu2_lattice
  !
  public :: print_impG_normal
  public :: print_impG_superc
  public :: print_impG_nonsu2
  !
  public :: print_impG_normal_lattice
  public :: print_impG_superc_lattice
  public :: print_impG_nonsu2_lattice
  !
  public :: print_impG0_normal
  public :: print_impG0_superc
  public :: print_impG0_nonsu2
  !
  public :: print_chi_spin
  public :: print_chi_dens
  public :: print_chi_dens_mix
  public :: print_chi_dens_tot
  public :: print_chi_pair



  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm


  character(len=64)                :: suffix


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case:
  ! - Poles&Weights
  ! - impSigma
  ! - impG
  ! - impG0
  ! NORMAL - SUPERConducting - NONSU2
  !+------------------------------------------------------------------+
  subroutine ed_print_PolesWeights
    select case(ed_mode)
    case ("normal")
       call print_poles_weights_normal
    case ("superc")
       call print_poles_weights_superc
    case ("nonsu2")
       call print_poles_weights_nonsu2
    case default
       stop "ed_print_PolesWeights error: ed_mode not in the list"
    end select
  end subroutine ed_print_PolesWeights


  subroutine ed_print_impSigma_single
    select case(ed_mode)
    case ("normal")
       call print_impSigma_normal
    case ("superc")
       call print_impSigma_superc
    case ("nonsu2")
       call print_impSigma_nonsu2
    case default
       stop "ed_print_impSigma error: ed_mode not in the list"
    end select
  end subroutine ed_print_impSigma_single
  !
  subroutine ed_print_impSigma_lattice(iprint)
    integer :: iprint
    select case(ed_mode)
    case ("normal")
       call print_impSigma_normal_lattice(iprint)
    case ("superc")
       call print_impSigma_superc_lattice(iprint)
    case ("nonsu2")
       call print_impSigma_nonsu2_lattice(iprint)
    case default
       stop " print_impSigma_lattice error: ed_mode not in the list"
    end select
  end subroutine ed_print_impSigma_lattice



  subroutine ed_print_impG_single
    select case(ed_mode)
    case ("normal")
       call print_impG_normal
    case ("superc")
       call print_impG_superc
    case ("nonsu2")
       call print_impG_nonsu2
    case default
       stop "ed_print_impG error: ed_mode not in the list"
    end select
  end subroutine ed_print_impG_single
  !
  subroutine ed_print_impG_lattice(iprint)
    integer :: iprint
    select case(ed_mode)
    case ("normal")
       call print_impG_normal_lattice(iprint)
    case ("superc")
       call print_impG_superc_lattice(iprint)
    case ("nonsu2")
       call print_impG_nonsu2_lattice(iprint)
    case default
       stop " print_impG_lattice error: ed_mode not in the list"
    end select
  end subroutine ed_print_impG_lattice



  subroutine ed_print_impG0
    select case(ed_mode)
    case ("normal")
       call print_impG0_normal
    case ("superc")
       call print_impG0_superc
    case ("nonsu2")
       call print_impG0_nonsu2
    case default
       stop "ed_print_impG0 error: ed_mode not in the list"
    end select
  end subroutine ed_print_impG0



  subroutine ed_print_impChi
    call print_chi_spin
    call print_chi_dens
    call print_chi_dens_mix
    call print_chi_dens_tot
    call print_chi_pair
  end subroutine ed_print_impChi




  include "ED_IO/print_PolesWeights.f90"
  include "ED_IO/print_impSigma.f90"
  include "ED_IO/print_impG.f90"
  include "ED_IO/print_impG0.f90"
  include "ED_IO/print_impChi.f90"














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
  subroutine ed_get_density_matrix(dm_,dm_eig_,dm_rot_)
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    real(8),allocatable,intent(out),optional     :: dm_eig_(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
    !internal
    integer                                      :: unit  
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    !
    if (.not.allocated(imp_density_matrix)) then
       write(*,*) "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(present(dm_eig_))then
       if(allocated(dm_eig_))deallocate(dm_eig_);allocate(dm_eig_(Nspin*Norb));dm_eig_ = 0.0d0
    endif
    if(present(dm_rot_))then
       if(allocated(dm_rot_))deallocate(dm_rot_);allocate(dm_rot_(Nspin*Norb,Nspin*Norb));dm_rot_ = zero
    endif
    if(allocated(dm_))deallocate(dm_);allocate(dm_(Nspin*Norb,Nspin*Norb));dm_ = zero
    !
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
    dm_eig_=0.0d0;dm_rot_=zero;dm_rot_=dm_
    call matrix_diagonalize(dm_rot_,dm_eig_,'V','U')
    !
    unit = free_unit()
    open(unit,file="imp_density_matrix.dat",action="write",position="rewind",status='unknown')
    write(unit,"(A10)")"# Re{rho}: [Norb*Norb]*Nspin"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    write(unit,"(A10)")"# Im{rho}: [Norb*Norb]*Nspin"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    write(unit,"(A10)")"# rho_tilda"
    write(unit,'(10F22.12)') dm_eig_
    write(unit,"(A10)")
    write(unit,"(A100)")"# Re{theta}: [Norb*Norb]*Nspin"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    write(unit,"(A100)")"# Im{theta}: [Norb*Norb]*Nspin"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nspin*Norb)
    enddo
    close(unit)
    !
    if(ed_verbose<1) then
       Tr=zero;Tr=trace(dm_)
       if(ED_MPI_ID==0)write(LOGfile,'(A25,6F15.7)') 'test #1: Tr[rho]:    ',real(Tr),aimag(Tr)
       Tr=zero;Tr=trace(matmul(dm_,dm_))
       if(ED_MPI_ID==0)write(LOGfile,'(A25,6F15.7)') 'test #2: Tr[rho^2]:  ',real(Tr),aimag(Tr)
       Tr=zero;Tr=sum(dm_eig_)
       if(ED_MPI_ID==0)write(LOGfile,'(A25,6F15.7)') 'test #3: Tr[rot(rho)]:',Tr
    endif
  end subroutine ed_get_density_matrix


  subroutine ed_get_quantum_SOC_operators(Simp,Limp,jimp)
    complex(8),optional,intent(inout) ::  Simp(:,:,:)
    complex(8),optional,intent(inout) ::  Limp(:,:,:)
    complex(8),optional,intent(inout) ::  Jimp(:)
    integer                           ::  unit_
    integer                           ::  iorb,ispin,jorb,jspin
    real(8)                           ::  Lxsq,Lysq,Lzsq,Lsq
    real(8)                           ::  Sxsq,Sysq,Szsq,Ssq
    real(8)                           ::  jxsq,jysq,jzsq,Jsq
    if(Norb/=3)stop"SOC_operators implemented for 3 orbitals"
    if(present(Simp).and.((size(Simp,dim=1)/=3).or.(size(Simp,dim=2)/=3).or.(size(Simp,dim=3)/=3)))stop"wrong S size (3,3,3)"
    if(present(Limp).and.((size(Limp,dim=1)/=3).or.(size(Limp,dim=2)/=2).or.(size(Limp,dim=3)/=2)))stop"wrong L size (3,2,2)"
    if(present(jimp).and.(size(jimp)/=3))stop"wrong j size (3)"
    !
    if(present(Simp))        Simp=impStot
    if(present(Limp))        Limp=impLtot
    if(present(jimp))        jimp=impj_aplha
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
    Jsq  = jxsq + jysq + jzsq
    !
    !   IMPURITY SPIN OPERATOR - (S)
    !
    unit_ = free_unit()
    open(unit=unit_,file='S_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Sx_(iorb,jorb)}, Im{Sx_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(impStot(1,iorb,jorb)),jorb=1,Norb),(aimag(impStot(1,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sy_(iorb,jorb)}, Im{Sy_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(impStot(2,iorb,jorb)),jorb=1,Norb),(aimag(impStot(2,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sz_(iorb,jorb)}, Im{Sz_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(impStot(3,iorb,jorb)),jorb=1,Norb),(aimag(impStot(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}","7-|Sx|^2","8-|Sy|^2","9-|Sz|^2","10-|S|^2"
    write(unit_,'(30(F20.12,1X))') real(trace(impStot(1,:,:))),aimag(trace(impStot(1,:,:))),&
         real(trace(impStot(2,:,:))),aimag(trace(impStot(2,:,:))),&
         real(trace(impStot(3,:,:))),aimag(trace(impStot(3,:,:))),&
         Sxsq,Sysq,Szsq,Ssq
    close(unit_)
    !
    !   IMPURITY ORBITAL ANGULAR MOMENTUM OPERATOR - (L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='L_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Lx_(ipin,jspin)}, Im{Lx_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(impLtot(1,ispin,jspin)),jspin=1,Nspin),(aimag(impLtot(1,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Ly_(ipin,jspin)}, Im{Ly_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(impLtot(2,ispin,jspin)),jspin=1,Nspin),(aimag(impLtot(2,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Lz_(ipin,jspin)}, Im{Lz_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(impLtot(3,ispin,jspin)),jspin=1,Nspin),(aimag(impLtot(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}","7-|Lx|^2","8-|Ly|^2","9-|Lz|^2","10-|L|^2"
    write(unit_,'(30(F20.12,1X))') real(trace(impLtot(1,:,:))),aimag(trace(impLtot(1,:,:))),&
         real(trace(impLtot(2,:,:))),aimag(trace(impLtot(2,:,:))),&
         real(trace(impLtot(3,:,:))),aimag(trace(impLtot(3,:,:))),&
         Lxsq,Lysq,Lzsq,Lsq
    close(unit_)
    !
    !   IMPURITY TOTAL ANGULAR MOMENTUM OPERATOR - (J = S + L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='J_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "#1-Re{jx}","2-Im{jx}","3-Re{jy}","4-Im{jy}","5-Re{jz}","6-Im{jz}","7-|jx|^2","8-|jy|^2","9-|jz|^2","10-|j|^2","11-Re{L.S}","12-Im{L.S}"
    write(unit_,'(30(F20.12,1X))') real(impj_aplha(1)),aimag(impj_aplha(1)),real(impj_aplha(2)),aimag(impj_aplha(2)),real(impj_aplha(3)),aimag(impj_aplha(3)) &
         ,jxsq,jysq,jzsq,Jsq,real(impLdotS),aimag(impLdotS)
    close(unit_)
    !
  end subroutine ed_get_quantum_SOC_operators



END MODULE ED_IO
