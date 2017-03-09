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





  !DEBUG>>
  subroutine ed_get_density_matrix(dm_,dm_eig_,dm_rot_,dm_custom_rot_)
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    real(8),allocatable,intent(out),optional     :: dm_eig_(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
    complex(8),allocatable,intent(in),optional   :: dm_custom_rot_(:,:)
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
    if(present(dm_custom_rot_))then
       dm_diag=matmul(transpose(conjg(dm_custom_rot_)),matmul(dm_,dm_custom_rot_))
    endif
    !
    dm_eig_=0.0d0;dm_rot_=zero;dm_rot_=dm_
    call matrix_diagonalize(dm_rot_,dm_eig_,'V','U')
    !
    unit = free_unit()
    open(unit,file="imp_density_matrix.dat",action="write",position="rewind",status='unknown')
    !
    write(unit,"(A10)")"# Re{rho}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    !
    write(unit,"(A10)")"# Im{rho}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    !
    write(unit,"(A10)")"# rho_tilda"
    write(unit,'(10F22.12)') dm_eig_
    write(unit,"(A10)")
    !
    write(unit,"(A100)")"# Re{theta}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    !
    write(unit,"(A100)")"# Im{theta}:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,"(A10)")
    !
    if(present(dm_custom_rot_))then
       !
       write(unit,"(A30)")"# Re{rho tilda CUSTOM}"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_diag(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       !
       write(unit,"(A30)")"# Im{rho tilda CUSTOM}"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_diag(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,"(A10)")
       write(unit,"(A30)")"# J basis densities"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_diag(io,io)),jo=1,Nspin*Norb)
       enddo
       !
    endif
    close(unit)
    !
  !  if(ed_verbose<1) then
  !     Tr=zero;Tr=trace(dm_)
  !     write(LOGfile,'(A25,6F15.7)') 'test #1: Tr[rho]:    ',real(Tr),aimag(Tr)
  !     Tr=zero;Tr=trace(matmul(dm_,dm_))
  !     write(LOGfile,'(A25,6F15.7)') 'test #2: Tr[rho^2]:  ',real(Tr),aimag(Tr)
  !     Tr=zero;Tr=sum(dm_eig_)
  !     write(LOGfile,'(A25,6F15.7)') 'test #3: Tr[rot(rho)]:',Tr
  !  endif
  end subroutine ed_get_density_matrix


  subroutine ed_get_quantum_SOC_operators(Simp,Limp,jimp)
    complex(8),optional,intent(inout) ::  Simp(:,:,:)
    complex(8),optional,intent(inout) ::  Limp(:,:,:)
    complex(8),optional,intent(inout) ::  Jimp(:)
    integer                           ::  unit_
    integer                           ::  iorb,ispin,jorb,jspin
    !real(8)                           ::  Lxsq,Lysq,Lzsq,Lsq
    !real(8)                           ::  Sxsq,Sysq,Szsq,Ssq
    real(8)                           ::  jxsq,jysq,jzsq,Jsq
    if(Norb/=3)stop "SOC_operators implemented for 3 orbitals"
    if(present(Simp).and.((size(Simp,dim=1)/=3).or.(size(Simp,dim=2)/=3).or.(size(Simp,dim=3)/=3)))stop "wrong S size (3,3,3)"
    if(present(Limp).and.((size(Limp,dim=1)/=3).or.(size(Limp,dim=2)/=2).or.(size(Limp,dim=3)/=2)))stop "wrong L size (3,2,2)"
    if(present(jimp).and.(size(jimp)/=3))stop "wrong j size (3)"
    !
    if(present(Simp))        Simp=impStot
    if(present(Limp))        Limp=impLtot
    if(present(jimp))        jimp=impj_aplha
    !
    !Sxsq = (trace(impStot(1,:,:)))*conjg(trace(impStot(1,:,:)))
    !Sysq = (trace(impStot(2,:,:)))*conjg(trace(impStot(2,:,:)))
    !Szsq = (trace(impStot(3,:,:)))*conjg(trace(impStot(3,:,:)))
    !Ssq  = Sxsq + Sysq + Szsq
    !
    !Lxsq = (trace(impLtot(1,:,:)))*conjg(trace(impStot(1,:,:)))
    !Lysq = (trace(impLtot(2,:,:)))*conjg(trace(impLtot(2,:,:)))
    !Lzsq = (trace(impLtot(3,:,:)))*conjg(trace(impLtot(3,:,:)))
    !Lsq  = Lxsq + Lysq + Lzsq
    !
    jxsq = impj_aplha_sq(1)
    jysq = impj_aplha_sq(2)
    jzsq = impj_aplha_sq(3)
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
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(impStot(1,:,:))),aimag(trace(impStot(1,:,:))),&
                                   real(trace(impStot(2,:,:))),aimag(trace(impStot(2,:,:))),&
                                   real(trace(impStot(3,:,:))),aimag(trace(impStot(3,:,:)))
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
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(impLtot(1,:,:))),aimag(trace(impLtot(1,:,:))),&
                                   real(trace(impLtot(2,:,:))),aimag(trace(impLtot(2,:,:))),&
                                   real(trace(impLtot(3,:,:))),aimag(trace(impLtot(3,:,:)))
    close(unit_)
    !
    !   IMPURITY TOTAL ANGULAR MOMENTUM OPERATOR - (J = S + L)
    !
    unit_ = free_unit()
    open(unit=unit_,file='J_imp.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "#1-Re{jx}","2-Im{jx}","3-Re{jy}","4-Im{jy}","5-Re{jz}","6-Im{jz}","7-|jx|^2","8-|jy|^2","9-|jz|^2","10-|j|^2","11-Re{L.S}","12-Im{L.S}"
    write(unit_,'(30(F20.12,1X))') real(impj_aplha(1)),aimag(impj_aplha(1)), &
                                   real(impj_aplha(2)),aimag(impj_aplha(2)), &
                                   real(impj_aplha(3)),aimag(impj_aplha(3)), &
                                   jxsq,jysq,jzsq,Jsq                      , &
                                   real(impLdotS),aimag(impLdotS)
    close(unit_)
    !
  end subroutine ed_get_quantum_SOC_operators



END MODULE ED_IO
