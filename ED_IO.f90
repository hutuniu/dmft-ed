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
  public :: ed_get_density_matrix


  public :: ed_print_PolesWeights
  public :: ed_print_impSigma
  public :: ed_print_impSigma_lattice
  public :: ed_print_impG
  public :: ed_print_impG_lattice
  public :: ed_print_impG0

  public :: print_poles_weights_normal
  public :: print_poles_weights_superc
  public :: print_poles_weights_nonsu2
  public :: print_impSigma_normal
  public :: print_impSigma_superc
  public :: print_impSigma_nonsu2
  public :: print_impG_normal
  public :: print_impG_superc
  public :: print_impG_nonsu2
  public :: print_impG0_normal
  public :: print_impG0_superc
  public :: print_impG0_nonsu2


  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix


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



  subroutine ed_print_impSigma
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
  end subroutine ed_print_impSigma
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



  subroutine ed_print_impG
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
  end subroutine ed_print_impG
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



  include "ED_IO/ed_io_print_PolesWeights.f90"
  include "ED_IO/ed_io_print_impSigma.f90"
  include "ED_IO/ed_io_print_impG.f90"
  include "ED_IO/ed_io_print_impG0.f90"









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
    complex(8),allocatable,intent(inout)           :: dm_(:,:)
    integer,intent(in)                             :: iprint
    real(8),allocatable,intent(inout),optional     :: dm_eig(:)
    complex(8),allocatable,intent(inout),optional  :: dm_rot(:,:)
    !internal
    integer                                        :: unit  
    integer                                        :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                     :: Tr

    if (ed_mode/="nonsu2") then
       write(*,*) "Not tested for ed_mode different from nonsu2"
       stop
    elseif (bath_type/="hybrid") then
       write(*,*) "Not tested for bath_type different from hybrid"
       stop
    elseif (ed_type/="c") then
       write(*,*) "Not tested for ed_type different from c"
       stop
    elseif (((.not.present(dm_eig)).or.(.not.present(dm_rot))).and.(iprint/=0)) then
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
    if(iprint==0) then
       unit = free_unit()
       open(unit,file="imp_density_matrix.ed",action="write",position="rewind",status='unknown')
       write(unit,"(A10)")"# Re{rho}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       write(unit,"(A10)")"# Im{rho}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
       enddo
       close(unit)
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
       unit = free_unit()
       open(unit,file="imp_density_matrix.ed",action="write",position="append",status='unknown')
       write(unit,"(A10)")
       write(unit,"(A10)")"# rho_tilda"
       write(unit,'(10F22.12)') dm_eig
       write(unit,"(A10)")
       write(unit,"(A10)")"# Re{theta}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       write(unit,"(A10)")"# Im{theta}: [Norb*Norb]*Nspin"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot(io,jo)),jo=1,Nspin*Norb)
       enddo
       close(unit)
    endif
  end subroutine ed_get_density_matrix
  !<<DEBUG


END MODULE ED_IO
