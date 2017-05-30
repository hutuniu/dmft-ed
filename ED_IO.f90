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
     module procedure ed_get_sigma_matsubara_lattice_1
     module procedure ed_get_sigma_matsubara_lattice_2
     module procedure ed_get_sigma_matsubara_lattice_3
     module procedure ed_get_sigma_matsubara_lattice_11
     module procedure ed_get_sigma_matsubara_lattice_21
     module procedure ed_get_sigma_matsubara_lattice_31
  end interface ed_get_sigma_matsubara

  interface ed_get_self_matsubara
     module procedure ed_get_self_matsubara_1
     module procedure ed_get_self_matsubara_2
     module procedure ed_get_self_matsubara_3
     module procedure ed_get_self_matsubara_lattice_1
     module procedure ed_get_self_matsubara_lattice_2
     module procedure ed_get_self_matsubara_lattice_3
     module procedure ed_get_self_matsubara_lattice_11
     module procedure ed_get_self_matsubara_lattice_21
     module procedure ed_get_self_matsubara_lattice_31
  end interface ed_get_self_matsubara

  interface ed_get_sigma_real
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
     module procedure ed_get_sigma_real_3
     module procedure ed_get_sigma_real_lattice_1
     module procedure ed_get_sigma_real_lattice_2
     module procedure ed_get_sigma_real_lattice_3
     module procedure ed_get_sigma_real_lattice_11
     module procedure ed_get_sigma_real_lattice_21
     module procedure ed_get_sigma_real_lattice_31
  end interface ed_get_sigma_real

  interface ed_get_self_real
     module procedure ed_get_self_real_1
     module procedure ed_get_self_real_2
     module procedure ed_get_self_real_3
     module procedure ed_get_self_real_lattice_1
     module procedure ed_get_self_real_lattice_2
     module procedure ed_get_self_real_lattice_3
     module procedure ed_get_self_real_lattice_11
     module procedure ed_get_self_real_lattice_21
     module procedure ed_get_self_real_lattice_31
  end interface ed_get_self_real




  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
     module procedure ed_get_gimp_matsubara_3
     module procedure ed_get_gimp_matsubara_lattice_1
     module procedure ed_get_gimp_matsubara_lattice_2
     module procedure ed_get_gimp_matsubara_lattice_3
     module procedure ed_get_gimp_matsubara_lattice_11
     module procedure ed_get_gimp_matsubara_lattice_21
     module procedure ed_get_gimp_matsubara_lattice_31
  end interface ed_get_gimp_matsubara

  interface ed_get_fimp_matsubara
     module procedure ed_get_fimp_matsubara_1
     module procedure ed_get_fimp_matsubara_2
     module procedure ed_get_fimp_matsubara_3
     module procedure ed_get_fimp_matsubara_lattice_1
     module procedure ed_get_fimp_matsubara_lattice_2
     module procedure ed_get_fimp_matsubara_lattice_3
     module procedure ed_get_fimp_matsubara_lattice_11
     module procedure ed_get_fimp_matsubara_lattice_21
     module procedure ed_get_fimp_matsubara_lattice_31
  end interface ed_get_fimp_matsubara

  interface ed_get_gimp_real
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
     module procedure ed_get_gimp_real_3
     module procedure ed_get_gimp_real_lattice_1
     module procedure ed_get_gimp_real_lattice_2
     module procedure ed_get_gimp_real_lattice_3
     module procedure ed_get_gimp_real_lattice_11
     module procedure ed_get_gimp_real_lattice_21
     module procedure ed_get_gimp_real_lattice_31
  end interface ed_get_gimp_real

  interface ed_get_fimp_real
     module procedure ed_get_fimp_real_1
     module procedure ed_get_fimp_real_2
     module procedure ed_get_fimp_real_3
     module procedure ed_get_fimp_real_lattice_1
     module procedure ed_get_fimp_real_lattice_2
     module procedure ed_get_fimp_real_lattice_3
     module procedure ed_get_fimp_real_lattice_11
     module procedure ed_get_fimp_real_lattice_21
     module procedure ed_get_fimp_real_lattice_31
  end interface ed_get_fimp_real


  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
     module procedure ed_get_dens_lattice_1
     module procedure ed_get_dens_lattice_2
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
     module procedure ed_get_mag_lattice_1
     module procedure ed_get_mag_lattice_2
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
     module procedure ed_get_docc_lattice_1
     module procedure ed_get_docc_lattice_2
  end interface ed_get_docc

  interface ed_get_phisc
     module procedure ed_get_phisc_1
     module procedure ed_get_phisc_2
     module procedure ed_get_phisc_lattice_1
     module procedure ed_get_phisc_lattice_2
  end interface ed_get_phisc

  interface ed_get_eimp
     module procedure :: ed_get_eimp_
     module procedure :: ed_get_eimp_lattice
  end interface ed_get_eimp

  interface ed_get_epot
     module procedure :: ed_get_epot_
     module procedure :: ed_get_epot_lattice
  end interface ed_get_epot

  interface ed_get_eint
     module procedure :: ed_get_eint_
     module procedure :: ed_get_eint_lattice
  end interface ed_get_eint

  interface ed_get_ehartree
     module procedure :: ed_get_ehartree_
     module procedure :: ed_get_ehartree_lattice
  end interface ed_get_ehartree

  interface ed_get_eknot
     module procedure :: ed_get_eknot_
     module procedure :: ed_get_eknot_lattice
  end interface ed_get_eknot

  interface ed_get_doubles
     module procedure :: ed_get_doubles_
     module procedure :: ed_get_doubles_lattice
  end interface ed_get_doubles

  interface ed_get_dust
     module procedure :: ed_get_dust_
     module procedure :: ed_get_dust_lattice
  end interface ed_get_dust

  interface ed_get_dund
     module procedure :: ed_get_dund_
     module procedure :: ed_get_dund_lattice
  end interface ed_get_dund

  interface ed_get_dse
     module procedure :: ed_get_dse_
     module procedure :: ed_get_dse_lattice
  end interface ed_get_dse

  interface ed_get_dph
     module procedure :: ed_get_dph_
     module procedure :: ed_get_dph_lattice
  end interface ed_get_dph



  ! interface ed_get_sigma_matsubara_lattice
  !    module procedure ed_get_sigma_matsubara_lattice_1
  !    module procedure ed_get_sigma_matsubara_lattice_2
  !    module procedure ed_get_sigma_matsubara_lattice_3
  !    module procedure ed_get_sigma_matsubara_lattice_11
  !    module procedure ed_get_sigma_matsubara_lattice_21
  !    module procedure ed_get_sigma_matsubara_lattice_31
  ! end interface ed_get_sigma_matsubara_lattice

  ! interface ed_get_self_matsubara_lattice
  !    module procedure ed_get_self_matsubara_lattice_1
  !    module procedure ed_get_self_matsubara_lattice_2
  !    module procedure ed_get_self_matsubara_lattice_3
  !    module procedure ed_get_self_matsubara_lattice_11
  !    module procedure ed_get_self_matsubara_lattice_21
  !    module procedure ed_get_self_matsubara_lattice_31
  ! end interface ed_get_self_matsubara_lattice

  ! interface ed_get_sigma_real_lattice
  !    module procedure ed_get_sigma_real_lattice_1
  !    module procedure ed_get_sigma_real_lattice_2
  !    module procedure ed_get_sigma_real_lattice_3
  !    module procedure ed_get_sigma_real_lattice_11
  !    module procedure ed_get_sigma_real_lattice_21
  !    module procedure ed_get_sigma_real_lattice_31
  ! end interface ed_get_sigma_real_lattice

  ! interface ed_get_self_real_lattice
  !    module procedure ed_get_self_real_lattice_1
  !    module procedure ed_get_self_real_lattice_2
  !    module procedure ed_get_self_real_lattice_3
  !    module procedure ed_get_self_real_lattice_11
  !    module procedure ed_get_self_real_lattice_21
  !    module procedure ed_get_self_real_lattice_31
  ! end interface ed_get_self_real_lattice


  ! interface ed_get_gimp_matsubara_lattice
  !    module procedure ed_get_gimp_matsubara_lattice_1
  !    module procedure ed_get_gimp_matsubara_lattice_2
  !    module procedure ed_get_gimp_matsubara_lattice_3
  !    module procedure ed_get_gimp_matsubara_lattice_11
  !    module procedure ed_get_gimp_matsubara_lattice_21
  !    module procedure ed_get_gimp_matsubara_lattice_31
  ! end interface ed_get_gimp_matsubara_lattice

  ! interface ed_get_fimp_matsubara_lattice
  !    module procedure ed_get_fimp_matsubara_lattice_1
  !    module procedure ed_get_fimp_matsubara_lattice_2
  !    module procedure ed_get_fimp_matsubara_lattice_3
  !    module procedure ed_get_fimp_matsubara_lattice_11
  !    module procedure ed_get_fimp_matsubara_lattice_21
  !    module procedure ed_get_fimp_matsubara_lattice_31
  ! end interface ed_get_fimp_matsubara_lattice

  ! interface ed_get_gimp_real_lattice
  !    module procedure ed_get_gimp_real_lattice_1
  !    module procedure ed_get_gimp_real_lattice_2
  !    module procedure ed_get_gimp_real_lattice_3
  !    module procedure ed_get_gimp_real_lattice_11
  !    module procedure ed_get_gimp_real_lattice_21
  !    module procedure ed_get_gimp_real_lattice_31
  ! end interface ed_get_gimp_real_lattice

  ! interface ed_get_fimp_real_lattice
  !    module procedure ed_get_fimp_real_lattice_1
  !    module procedure ed_get_fimp_real_lattice_2
  !    module procedure ed_get_fimp_real_lattice_3
  !    module procedure ed_get_fimp_real_lattice_11
  !    module procedure ed_get_fimp_real_lattice_21
  !    module procedure ed_get_fimp_real_lattice_31
  ! end interface ed_get_fimp_real_lattice

  ! interface ed_get_dens_lattice
  !    module procedure ed_get_dens_lattice_1
  !    module procedure ed_get_dens_lattice_2
  ! end interface ed_get_dens_lattice

  ! interface ed_get_mag_lattice
  !    module procedure ed_get_mag_lattice_1
  !    module procedure ed_get_mag_lattice_2
  ! end interface ed_get_mag_lattice

  ! interface ed_get_docc_lattice
  !    module procedure ed_get_docc_lattice_1
  !    module procedure ed_get_docc_lattice_2
  ! end interface ed_get_docc_lattice

  ! interface ed_get_phisc_lattice
  !    module procedure ed_get_phisc_lattice_1
  !    module procedure ed_get_phisc_lattice_2
  ! end interface ed_get_phisc_lattice



  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_real
  public :: ed_get_self_real

  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_real
  public :: ed_get_fimp_real

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc

  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot

  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph

  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators




  ! public :: ed_get_eimp_lattice
  ! public :: ed_get_epot_lattice
  ! public :: ed_get_eint_lattice
  ! public :: ed_get_ehartree_lattice
  ! public :: ed_get_eknot_lattice

  ! public :: ed_get_doubles_lattice
  ! public :: ed_get_dust_lattice
  ! public :: ed_get_dund_lattice
  ! public :: ed_get_dse_lattice
  ! public :: ed_get_dph_lattice




  ! public :: ed_get_sigma_matsubara_lattice
  ! public :: ed_get_self_matsubara_lattice
  ! public :: ed_get_sigma_real_lattice
  ! public :: ed_get_self_real_lattice

  ! public :: ed_get_gimp_matsubara_lattice
  ! public :: ed_get_fimp_matsubara_lattice
  ! public :: ed_get_gimp_real_lattice
  ! public :: ed_get_fimp_real_lattice

  ! public :: ed_get_dens_lattice
  ! public :: ed_get_mag_lattice
  ! public :: ed_get_docc_lattice
  ! public :: ed_get_phisc_lattice




  interface ed_print_impSigma
     module procedure :: ed_print_impSigma_single
     module procedure :: ed_print_impSigma_lattice
  end interface ed_print_impSigma

  interface  ed_print_impG
     module procedure :: ed_print_impG_single
     module procedure :: ed_print_impG_lattice
  end interface ed_print_impG

  ! public :: ed_print_PolesWeights
  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impChi

  !FOR INTERNAL USE ONLY:
  ! public :: print_poles_weights_normal
  ! public :: print_poles_weights_superc
  ! public :: print_poles_weights_nonsu2
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
  ! subroutine ed_print_PolesWeights
  !   select case(ed_mode)
  !   case ("normal")
  !      call print_poles_weights_normal
  !   case ("superc")
  !      call print_poles_weights_superc
  !   case ("nonsu2")
  !      call print_poles_weights_nonsu2
  !   case default
  !      stop "ed_print_PolesWeights error: ed_mode not in the list"
  !   end select
  ! end subroutine ed_print_PolesWeights


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




  ! include "ED_IO/print_PolesWeights.f90"
  include "ED_IO/print_impSigma.f90"
  include "ED_IO/print_impG.f90"
  include "ED_IO/print_impG0.f90"
  include "ED_IO/print_impChi.f90"





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_sigma_matsubara.f90"
  include "ED_IO/get_self_matsubara.f90"
  include "ED_IO/get_sigma_realaxis.f90"
  include "ED_IO/get_self_realaxis.f90"








  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_gimp_matsubara.f90"
  include "ED_IO/get_fimp_matsubara.f90"
  include "ED_IO/get_gimp_realaxis.f90"
  include "ED_IO/get_fimp_realaxis.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  include "ED_IO/get_dens.f90"
  include "ED_IO/get_mag.f90"
  include "ED_IO/get_docc.f90"
  include "ED_IO/get_phisc.f90"
  include "ED_IO/get_eimp.f90"
  include "ED_IO/get_doubles.f90"





  subroutine ed_get_density_matrix(dm_,custom_rot,dm_eig_,dm_rot_)
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    complex(8),allocatable,intent(in)            :: custom_rot(:,:)
    real(8),allocatable,intent(out),optional     :: dm_eig_(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
    !internal
    integer                                      :: unit  
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_diag(:,:)
    !
    if (.not.allocated(imp_density_matrix)) then
       write(*,*) "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(present(dm_eig_))then
       if(allocated(dm_eig_))deallocate(dm_eig_)
       allocate(dm_eig_(Nspin*Norb));dm_eig_ = 0.0d0
    endif
    if(present(dm_rot_))then
       if(allocated(dm_rot_))deallocate(dm_rot_)
       allocate(dm_rot_(Nspin*Norb,Nspin*Norb));dm_rot_ = zero
    endif
    if(allocated(dm_))    deallocate(dm_);    allocate(dm_(Nspin*Norb,Nspin*Norb));    dm_ = zero
    if(allocated(dm_diag))deallocate(dm_diag);allocate(dm_diag(Nspin*Norb,Nspin*Norb));dm_diag = zero
    !
    if(bath_type=="replica")then
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
       if(present(dm_eig_).and.present(dm_rot_))then
          dm_rot_=dm_
          call eigh(dm_rot_,dm_eig_,'V','U')
       endif
       !
       ! rotate from the impurity basis to the J diagonal basis
       dm_diag=matmul(transpose(conjg(custom_rot)),matmul(dm_,custom_rot))
       !
    elseif(bath_type=="normal")then
       !
       do iorb=1,Norb
          do ispin=1,Nspin
             io = iorb + (ispin-1)*Norb
             dm_diag(io,io) =  imp_density_matrix(ispin,ispin,iorb,iorb)
          enddo
       enddo
       !
       dm_=matmul(custom_rot,matmul(dm_diag,transpose(conjg(custom_rot))))
       !
       if(present(dm_eig_).and.present(dm_rot_))then
          dm_rot_=dm_
          call eigh(dm_rot_,dm_eig_,'V','U')
       endif
       !
    endif
    !
    unit = free_unit()
    open(unit,file="imp_density_matrix.dat",action="write",position="rewind",status='unknown')
    !
    write(unit,"(A30)")"# Re{rho_asbs'}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A30)")"# Im{rho_asbs'}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    if(present(dm_eig_).and.present(dm_rot_))then
       write(unit,"(A30)")"# rho_tilda"
       write(unit,'(10F22.12)') dm_eig_
       write(unit,*)
       !
       write(unit,"(A30)")"# Re{theta}:"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,*)
       !
       write(unit,"(A30)")"# Im{theta}:"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,*)
    endif
    !
    write(unit,"(A30)")"# Re{rho_Jj}"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_diag(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A30)")"# Im{rho_Jj}"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_diag(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    write(unit,"(A30)")"# J basis densities"
    write(unit,"(90(F15.9,1X))") (real(dm_diag(io,io)),io=1,Nspin*Norb)
    !
    close(unit)
    !
  end subroutine ed_get_density_matrix


  subroutine ed_get_quantum_SOC_operators()
    complex(8)                        ::  Simp(3,Norb,Norb)
    complex(8)                        ::  Limp(3,Nspin,Nspin)
    complex(8)                        ::  Jimp(3)
    complex(8)                        ::  Jimp_sq(3)
    complex(8),allocatable            ::  rho_so(:,:),U(:,:),Udag(:,:)
    complex(8),allocatable            ::  rho_nn(:,:,:,:)
    complex(8)                        ::  Jsq
    complex(8)                        ::  LSimp
    complex(8)                        ::  LSbth(Nbath)
    integer                           ::  unit_
    integer                           ::  iorb,ispin,jorb,jspin,io,ibath
    !
    if(Norb/=3)stop "SOC_operators implemented for 3 orbitals"
    !
    if((bath_type=="replica").and.(.not.Jz_basis))then
       !the impurity dm is in the {a,Sz} basis
       Simp=impStot
       Limp=impLtot
       jimp=impj_aplha
       jimp_sq=impj_aplha_sq
       Jsq=sum(jimp_sq)
       LSimp=impLdotS
       do ibath=1,Nbath
          LSbth(ibath)=bthLdotS(ibath)
       enddo
    elseif((bath_type=="replica".and.Jz_basis).or.(bath_type=="replica"))then
       allocate(U(Nspin*Norb,Nspin*Norb));U=zero
       allocate(Udag(Nspin*Norb,Nspin*Norb));Udag=zero
       allocate(rho_so(Nspin*Norb,Nspin*Norb));rho_so=zero
       allocate(rho_nn(Nspin,Nspin,Norb,Norb));rho_nn=zero
       !
       rho_so=nn2so_reshape(imp_density_matrix,Nspin,Norb)
       !
       if(bath_type=="replica".and.Jz_basis)then
          !
          !the impurity dm is in the {Lz,Sz} basis
          !rotation {Lz,Sz}-->{a,Sz}
          U=transpose(conjg(orbital_Lz_rotation_NorbNspin()))
          Udag=transpose(conjg(U))
          !
       elseif(bath_type=="normal")then
          !
          !the impurity dm is in the {J} basis
          !rotation {J}-->{a,Sz}
          U=transpose(conjg(atomic_SOC_rotation()))
          Udag=transpose(conjg(U))
          !
       endif
       !
       rho_so=matmul(Udag,matmul(rho_so,U))
       rho_nn=so2nn_reshape(rho_so,Nspin,Norb)
       !
       !#    <S>(iorb,jorb)   #
       do iorb=1,Norb
          do jorb=1,Norb
             Simp(1,iorb,jorb) = 0.5d0*( rho_nn(1,2,iorb,jorb) + rho_nn(2,1,iorb,jorb) )
             Simp(2,iorb,jorb) = 0.5d0*( rho_nn(2,1,iorb,jorb) - rho_nn(1,2,iorb,jorb) )*xi
             Simp(3,iorb,jorb) = 0.5d0*( rho_nn(1,1,iorb,jorb) - rho_nn(2,2,iorb,jorb) )
          enddo
       enddo
       !
       !#   <L>(ispin,jspin)   #
       do ispin=1,Nspin
          do jspin=1,Nspin
             Limp(1,ispin,jspin) = ( rho_nn(ispin,jspin,3,2) - rho_nn(ispin,jspin,2,3) )*xi
             Limp(2,ispin,jspin) = ( rho_nn(ispin,jspin,1,3) - rho_nn(ispin,jspin,3,1) )*xi
             Limp(3,ispin,jspin) = ( rho_nn(ispin,jspin,2,1) - rho_nn(ispin,jspin,1,2) )*xi
          enddo
       enddo
       !
       !#          <J>         #
       jimp(1) = trace(matmul(rho_so,atomic_j("x")))
       jimp(2) = trace(matmul(rho_so,atomic_j("y")))
       jimp(3) = trace(matmul(rho_so,atomic_j("z")))
       !
       !#          <Jz>        #
       jimp_sq(1) = trace(matmul(rho_so,matmul(atomic_j("x"),atomic_j("x"))))
       jimp_sq(2) = trace(matmul(rho_so,matmul(atomic_j("y"),atomic_j("y"))))
       jimp_sq(3) = trace(matmul(rho_so,matmul(atomic_j("z"),atomic_j("z"))))
       !
       Jsq=sum(jimp_sq)
       !#          <LS>        #
       LSimp = trace(matmul(rho_so,atomic_SOC()))
    endif
    !
    !   IMPURITY SPIN OPERATOR - (S)
    !
    unit_ = free_unit()
    open(unit=unit_,file='S_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Sx_(iorb,jorb)}, Im{Sx_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(Simp(1,iorb,jorb)),jorb=1,Norb),(aimag(Simp(1,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sy_(iorb,jorb)}, Im{Sy_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(Simp(2,iorb,jorb)),jorb=1,Norb),(aimag(Simp(2,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sz_(iorb,jorb)}, Im{Sz_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(Simp(3,iorb,jorb)),jorb=1,Norb),(aimag(Simp(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(Simp(1,:,:))),aimag(trace(Simp(1,:,:))),&
                                   real(trace(Simp(2,:,:))),aimag(trace(Simp(2,:,:))),&
                                   real(trace(Simp(3,:,:))),aimag(trace(Simp(3,:,:)))
    close(unit_)
    !
    !   IMPURITY ORBITAL ANGULAR MOMENTUM OPERATOR - (L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='L_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Lx_(ipin,jspin)}, Im{Lx_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(Limp(1,ispin,jspin)),jspin=1,Nspin),(aimag(Limp(1,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Ly_(ipin,jspin)}, Im{Ly_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(Limp(2,ispin,jspin)),jspin=1,Nspin),(aimag(Limp(2,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Lz_(ipin,jspin)}, Im{Lz_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(Limp(3,ispin,jspin)),jspin=1,Nspin),(aimag(Limp(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(Limp(1,:,:))),aimag(trace(Limp(1,:,:))),&
                                   real(trace(Limp(2,:,:))),aimag(trace(Limp(2,:,:))),&
                                   real(trace(Limp(3,:,:))),aimag(trace(Limp(3,:,:)))
    close(unit_)
    !
    !   IMPURITY TOTAL ANGULAR MOMENTUM OPERATOR - (J = S + L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='J_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "# 1-Re{jx}   "," 2-Im{jx}   "," 3-Re{jy}   "," 4-Im{jy}   "," 5-Re{jz}   "," 6-Im{jz}   ", &
                                 " 7-Re{jx_sq}"," 8-Im{jx_sq}"," 9-Re{jy_sq}","10-Im{jy_sq}","11-Re{jz_sq}","12-Im{jz_sq}", &
                                 "13-Re{L.S}","14-Im{L.S}"
    write(unit_,'(30(F20.12,1X))') real(jimp(1))   ,aimag(jimp(1))   ,real(jimp(2))   ,aimag(jimp(2))   ,real(jimp(3))   ,aimag(jimp(3))   , &
                                   real(jimp_sq(1)),aimag(jimp_sq(1)),real(jimp_sq(2)),aimag(jimp_sq(2)),real(jimp_sq(3)),aimag(jimp_sq(3)), &
                                   real(LSimp),aimag(LSimp),(real(LSbth(io)),io=1,Nbath)
    close(unit_)
    !
  end subroutine ed_get_quantum_SOC_operators



END MODULE ED_IO
