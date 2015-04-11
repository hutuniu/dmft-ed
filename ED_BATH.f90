!########################################################################
! The dimensions of the bath components are:
! e = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites per orbital]
! v = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites per orbital]
! N = Nspin*(2*Norb)*Nbath
!
! e = Nspin[# of spins] * 1 * Nbath[# of bath sites]
! v = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites]
! N = Nspin*(Norb+1)*Nbath
!########################################################################
MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none

  private

  !\Delta NORMAL hybridization function Matsubara
  interface delta_bath_mats
     module procedure delta_bath_mats_1
     module procedure delta_bath_mats_2
  end interface delta_bath_mats

  !\Delta ANOMALOUS hybridization function Matsubara
  interface fdelta_bath_mats
     module procedure fdelta_bath_mats_1
     module procedure fdelta_bath_mats_2
  end interface fdelta_bath_mats

  !\Nabla\Delta NORMAL hybridization function gradient Matsubara
  interface grad_delta_bath_mats
     module procedure grad_delta_bath_mats_1
     module procedure grad_delta_bath_mats_2
  end interface grad_delta_bath_mats

  !\Nabla\Delta ANOMALOUS hybridization function gradient Matsubara
  interface grad_fdelta_bath_mats
     module procedure grad_fdelta_bath_mats_1
     module procedure grad_fdelta_bath_mats_2
  end interface grad_fdelta_bath_mats

  !Weiss NORMAL non-interacting Green's function Matsubara
  interface weiss_bath_mats
     module procedure weiss_bath_mats_1
     module procedure weiss_bath_mats_2
  end interface weiss_bath_mats

  !Weiss ANOMALOUS non-interacting Green's function Matsubara
  interface fweiss_bath_mats
     module procedure fweiss_bath_mats_1
     module procedure fweiss_bath_mats_2
  end interface fweiss_bath_mats


  !\Delta NORMAL hybridization function Real
  interface delta_bath_real
     module procedure delta_bath_real_1
     module procedure delta_bath_real_2
  end interface delta_bath_real

  !\Delta ANOMALOUS hybridization function Real
  interface fdelta_bath_real
     module procedure fdelta_bath_real_1
     module procedure fdelta_bath_real_2
  end interface fdelta_bath_real

  public :: allocate_bath              !INTERNAL (for effective_bath type)
  public :: deallocate_bath            !INTERNAL (for effective_bath type)
  public :: get_bath_size              !PUBLIC
  public :: check_bath_dimension       !PUBLIC
  public :: init_bath_ed               !INTERNAL (for effective_bath type)
  public :: write_bath                 !INTERNAL (for effective_bath type)
  public :: save_bath                  !INTERNAL (for effective_bath type)
  public :: set_bath                   !INTERNAL (for effective_bath type)
  public :: copy_bath                  !INTERNAL (for effective_bath type)
  public :: break_symmetry_bath        !PUBLIC
  public :: spin_symmetrize_bath       !PUBLIC
  public :: ph_symmetrize_bath         !PUBLIC
  public :: ph_trans_bath              !PUBLIC
  public :: enforce_normal_bath        !PUBLIC
  public :: set_nonsu2_components      !INTERNAL (for effective_bath type)
  public :: get_array_bath_dimension   !INTERNAL (for array bath type)
  public :: check_array_bath_dimension !INTERNAL (for array bath type)
  public :: dmft_bath2array            !INTERNAL (for effective_bath type)
  public :: array2dmft_bath            !INTERNAL (for effective_bath type)
  !
  public :: delta_bath_mats            !DELTA FUNCTION
  public :: fdelta_bath_mats           !DELTA FUNCTION
  !
  public :: weiss_bath_mats            !WEISS FUNCTION
  public :: fweiss_bath_mats           !WEISS FUNCTION
  !
  public :: grad_delta_bath_mats       !GRAD_DELTA FUNCTION
  public :: grad_fdelta_bath_mats      !GRAD_DELTA FUNCTION
  !
  public :: delta_bath_real            !DELTA FUNCTION
  public :: fdelta_bath_real           !DELTA FUNCTION


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(dmft_bath_%status)call deallocate_bath(dmft_bath_)
    allocate(dmft_bath_%v(Nhel,Nspin,Norb,Nbath))   !hybridization from impurities to bath (up->up,dw->dw)
    select case(bath_type)
    case default
       allocate(dmft_bath_%e(Nspin,Norb,Nbath))        !local energies of the bath
       select case(ed_mode)
       case default                                    !normal [N,Sz]
       case ("superc")                                 !superc [Sz] 
          allocate(dmft_bath_%d(Nspin,Norb,Nbath))     !local SC order parameters the bath  
          dmft_bath_%superc_status=.true.
       case ("nonsu2")                                 !nonsu2 [N]
          allocate(dmft_bath_%w(Nhel,Nhel,Norb,Nbath)) !hybridization W^{ss`} (up->up,up->dw,dw->up,dw->dw)
          allocate(dmft_bath_%h(Nhel,Norb,Nbath))     !local energies E^{s} (Nhel=2)
          dmft_bath_%nonsu2_status=.true.
       end select
    case('hybrid')
       allocate(dmft_bath_%v(Nhel,Nspin,Norb,Nbath))   !hybridization from impurities to bath (up->up,dw->dw)
       allocate(dmft_bath_%e(Nspin,1   ,Nbath))        !local energies of the bath
       select case(ed_mode)
       case default                                    !normal  [N,Sz]
       case ("superc")                                 !superc  [Sz]
          allocate(dmft_bath_%d(Nspin,1,Nbath))        !local SC order parameters the bath
          dmft_bath_%superc_status=.true.
       case ("nonsu2")                                 !nonsu2 case [N] qn
          allocate(dmft_bath_%w(Nhel,Nhel,Norb,Nbath)) !hybridization from impurities to bath (up->dw,dw->up)
          allocate(dmft_bath_%h(Nhel,1   ,Nbath))     !local energies E^{s} (Nhel=2)
          dmft_bath_%nonsu2_status=.true.
       end select
    end select
    dmft_bath_%status=.true.
  end subroutine allocate_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%v))deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%d))deallocate(dmft_bath_%d)
    if(allocated(dmft_bath_%w))deallocate(dmft_bath_%w)
    if(allocated(dmft_bath_%h))deallocate(dmft_bath_%h)
    dmft_bath_%superc_status=.false.
    dmft_bath_%nonsu2_status=.false.
    dmft_bath_%status=.false.
  end subroutine deallocate_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  ! in the nonSU2 channel the auxiliary arrays W and EH are not 
  ! taken into account for the calculation of the size
  !+-------------------------------------------------------------------+
  function get_bath_size() result(bath_size) !result(dims)
    integer :: dims(2)
    integer :: bath_size
    dims(1)=Nspin
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          !(e [][Norb][Nbath] + v [Nhel][][Norb][Nbath])
          dims(2) = (1+Nhel)*Norb*Nbath
       case ("superc")
          !(e [][Norb][Nbath] + v [Nhel][][Norb][Nbath] + d [][Norb][Nbath])
          dims(2) = (2+Nhel)*Norb*Nbath
       case ("nonsu2")
          !(e [][Norb][Nbath] + v [Nhel][][Norb][Nbath])
          dims(2) = (1+Nhel)*Norb*Nbath
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          !(e [][1][Nbath] + v [Nhel][][Norb][Nbath])
          dims(2) = (Nhel*Norb+1)*Nbath
       case ("superc")
          !(e [][1][Nbath] + v [Nhel][][Norb][Nbath] + d [][1][Nbath])
          dims(2) = (Nhel*Norb+2)*Nbath
       case ("nonsu2")
          !(e [][1][Nbath] + v [Nhel][][Norb][Nbath])
          dims(2) = (Nhel*Norb+1)*Nbath
       end select
    end select
    bath_size=dims(1)*dims(2)
  end function get_bath_size










  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    real(8),dimension(:) :: bath_
    integer              :: N1_,N2_,Ntrue
    logical              :: bool
    ! N1_=size(bath_,1)
    ! N2_=size(bath_,2)
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          Ntrue = (1+Nhel)*Norb*Nbath !(2*Norb)*Nbath
       case ("superc")
          Ntrue = (2+Nhel)*Norb*Nbath !(3*Norb)*Nbath
       case ("nonsu2")
          Ntrue = (1+Nhel)*Norb*Nbath !(3*Norb)*Nbath
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          Ntrue = (Nhel*Norb+1)*Nbath !(Norb+1)*Nbath
       case ("superc")
          Ntrue = (Nhel*Norb+2)*Nbath !(Norb+2)*Nbath
       case ("nonsu2")
          Ntrue = (Nhel*Norb+1)*Nbath !(2*Norb+1)*Nbath
       endif
    end select
    !bool = (N1_ == Nspin).AND.(N2_ == Ntrue)
    bool = (size(bath) == Nspin*Ntrue)
  end function check_bath_dimension







  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed(dmft_bath_,hwband_)
    type(effective_bath) :: dmft_bath_
    real(8)              :: hwband_,wband_
    integer              :: i,iorb,ispin,ihel,unit,flen,Nh
    logical              :: IOfile
    real(8)              :: de
    if(.not.dmft_bath_%status)stop "init_bath: bath not allocated"
    !Generating the bath anyway, then you may want to read it to 
    !update some entries. This way you can restart even 
    !from different Ns calculation, this is better than 
    !start from a complete guess
    dmft_bath_%e(:,:,1)    =-hwband_
    dmft_bath_%e(:,:,Nbath)= hwband_
    Nh=Nbath/2
    if(mod(Nbath,2)==0)then
       de=hwband_/max(Nh-1,1)
       dmft_bath_%e(:,:,Nh)  = -1.d-4
       dmft_bath_%e(:,:,Nh+1)=  1.d-4
       do i=2,Nh-1
          dmft_bath_%e(:,:,i)   =-hwband_ + (i-1)*de 
          dmft_bath_%e(:,:,Nbath-i+1)= hwband_ - (i-1)*de
       enddo
    else
       de=hwband_/Nh
       dmft_bath_%e(:,:,Nh+1)= 0.d0
       do i=2,Nh
          dmft_bath_%e(:,:,i)        =-hwband_ + (i-1)*de
          dmft_bath_%e(:,:,Nbath-i+1)= hwband_ - (i-1)*de
       enddo
    endif
    do i=1,Nbath
       dmft_bath_%v(:,:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))
    enddo
    if(ed_mode=="superc")then
       dmft_bath_%d(:,:,:)=deltasc
    endif
    !
    !Read from file if exist:
    !
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")'Reading bath from file '//trim(Hfile)//trim(ed_file_suffix)//".restart"
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
       read(unit,*)
       select case(bath_type)
       case default
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%d(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
       case ('hybrid')
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     dmft_bath_%d(ispin,1,i),&
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
       end select
       close(unit)
    endif
  end subroutine init_bath_ed





  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i,ispin,iorb,ihel,
    if(ED_MPI_ID==0)then
       unit_=LOGfile;if(present(unit))unit_=unit
       if(.not.dmft_bath_%status)stop "WRITE_BATH: bath not allocated"
       select case(bath_type)
       case default
          select case(ed_mode)
          case default
             write(unit,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Norb),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             write(unit,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%d(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("nonsu2")
             write(unit,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb), ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     (dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
       case('hybrid')
          select case(ed_mode)
          case default
             write(unit,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  (("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),&
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             write(unit,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  "Dk_s"//reg(txtfy(ispin)) ,&
                  (("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),&
                     dmft_bath_%d(ispin,1,i),&
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("nonsu2")
             write(unit,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  (("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_h"//reg(txtfy(ihel)),&
                  ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),    &
                     ((dmft_bath_%v(ihel,ispin,iorb,i),&
                     ihel=1,Nhel),iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
       end select
    endif
  end subroutine write_bath





  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given file with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine save_bath(dmft_bath_,file)
    type(effective_bath) :: dmft_bath_
    integer              :: unit_
    character(len=*)     :: file
    integer              :: i,ispin,iorb
    if(ED_MPI_ID==0)then
       if(.not.dmft_bath_%status)stop "save_bath error: bath not allocated"
       unit_=free_unit()
       open(unit_,file=trim(file)//".used")
       call write_bath(dmft_bath_,unit_)
       close(unit_)
    endif
  end subroutine save_bath





  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath_,dmft_bath_)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    integer              :: i,iorb,ispin,ihel,stride
    logical              :: check
    if(.not.dmft_bath_%status)stop "SET_BATH: bath not allocated"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_bath: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath_(io)
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%d(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath_(io)
                   enddo
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath_(io)
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    case ('hybrid')
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath(io)
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%d(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath_(io)
                   enddo
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      dmft_bath_%v(ihel,ispin,iorb,i) = bath(io)
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    end select
  end subroutine set_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine copy_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:,:) :: bath_
    integer                :: iorb,ispin
    integer                :: stride_spin,stride_orb
    integer                :: ei(2),vi(2),di(2)
    logical                :: check
    if(.not.dmft_bath_%status)stop "COPY_BATH: bath not allocated"
    check=check_bath_dimension(bath_)
    if(.not.check)stop "copy_bath: wrong bath dimensions"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%e(ispin,iorb,i) 
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      bath_(io) = dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath_%d(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      bath_(io) =  dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb + (ihel-1)*Nbath*Norb*Nspin
                      bath_(io) =  dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    case ('hybrid')
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      bath_(io) =  dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%d(ispin,1,i)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      bath_(io) =  dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) = dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ihel=1,Nhel
             do ispin=1,Nspin
                do iorb=1,Norb
                   do i=1,Nbath
                      io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath + (ihel-1)*Nspin*Norb*Nbath
                      bath_(io) =  dmft_bath_%v(ihel,ispin,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    end select
  end subroutine copy_bath
















  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array apply a specific transformation or 
  ! impose a given symmetry:
  ! - break spin symmetry by applying a symmetry breaking field
  ! - given a bath array set both spin components to have 
  !    the same bath, i.e. impose non-magnetic solution
  ! - given a bath array enforces the particle-hole symmetry 
  !    by setting the positive energies in modulo identical to the negative
  !    ones.
  ! - given a bath enforce normal (i.e. non superconducting) solution
  !+-------------------------------------------------------------------+
  subroutine break_symmetry_bath(bath_,field,sign)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    real(8)              :: field
    real(8)              :: sign
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine break_symmetry_bath

  subroutine spin_symmetrize_bath(bath_)
    real(8),dimension(:,:) :: bath_
    type(effective_bath)   :: dmft_bath_
    if(Nspin==1)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
    dmft_bath_%v(:,Nspin,:,:)=dmft_bath_%v(:,1,:,:)
    if(ed_mode=="superc")dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath

  subroutine ph_symmetrize_bath(bath_)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    integer              :: i
    if(ed_mode=="nonsu2")stop "PH symmetry not yet implemented for non-SU(2)..."
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    if(mod(Nbath,2)==0)then
       do i=1,Nbath/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_supercond)dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
    else
       do i=1,(Nbath-1)/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_supercond)dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
       dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
    endif
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath

  subroutine ph_trans_bath(bath_)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    type(effective_bath) :: tmp_dmft_bath
    integer              :: i
    if(ed_mode=="nonsu2")stop "PH symmetry not yet implemented for non-SU(2)..."
    call allocate_bath(dmft_bath_)
    call allocate_bath(tmp_dmft_bath)
    call set_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    do i=1,Nbath
       select case(Norb)
       case default
          ! do nothing
          dmft_bath_%e(:,:,i)= dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)= dmft_bath_%v(:,:,i)
          if(ed_supercond)dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)
       case(1)
          dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
          if(ed_supercond)dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)
       case(2)
          tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
          tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
          dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
          tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
          tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
          dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
          if(ed_supercond)dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)          
       end select
    end do
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine ph_trans_bath

  subroutine enforce_normal_bath(bath_)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    integer              :: i
    if(ed_mode=="nonsu2")stop "PH symmetry not yet implemented for non-SU(2)..."
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    do i=1,Nbath
       if(ed_supercond)dmft_bath_%d(:,:,i)=0d0
    enddo
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine enforce_normal_bath





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: given a dmft bath for the ed_mode=nonsu2 case set the hidden 
  ! components W^{ss'}_\a(l) and H^{s}_\a(l) used in the calculation of 
  !+-----------------------------------------------------------------------------+!
  subroutine pull_nonsu2_components(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(ed_mode/="nonsu2")then
       write(LOGunit,"(A)")"set_nonsu2_components: called with ed_mode != nonsu2. return."
       return
    endif
    if(.not.dmft_bath_%status)stop "set_nonsu2_components: bath not allocated"
    select case(Nspin)          !Nhel=2 for ed_mode=nonsu2
    case (1)
       dmft_bath_%w(1,1,:,:) = dmft_bath_%v(1,1,:,:) ![Nhel][Nspin][Norb][Nbath]
       dmft_bath_%w(1,2,:,:) = dmft_bath_%v(2,1,:,:)
       dmft_bath_%w(2,1,:,:) = dmft_bath_%v(2,1,:,:)
       dmft_bath_%w(2,2,:,:) = dmft_bath_%v(1,1,:,:)
       !
       dmft_bath_%h(1,:,:) = dmft_bath_%e(1,:,:)    ![Nspin][1/Norb][Nbath]
       dmft_bath_%h(2,:,:) = dmft_bath_%e(1,:,:)    ![Nspin][1/Norb][Nbath]
    case (2)
       dmft_bath_%w(1,1,:,:) = dmft_bath_%v(1,1,:,:)
       dmft_bath_%w(1,2,:,:) = dmft_bath_%v(2,1,:,:)
       dmft_bath_%w(2,1,:,:) = dmft_bath_%v(1,2,:,:)
       dmft_bath_%w(2,2,:,:) = dmft_bath_%v(2,2,:,:)
       !
       dmft_bath_%h(1,:,:) = dmft_bath_%e(1,:,:)    ![Nspin][1/Norb][Nbath]
       dmft_bath_%h(2,:,:) = dmft_bath_%e(2,:,:)    ![Nspin][1/Norb][Nbath]
    end select
  end subroutine pull_nonsu2_components

  !+-----------------------------------------------------------------------------+!
  !PURPOSE: given a dmft bath for the ed_mode=nonsu2 case set the hidden 
  ! components W^{ss'}_\a(l) and Eh^{s}_\a(l) used in the calculation of 
  !+-----------------------------------------------------------------------------+!
  subroutine push_nonsu2_components(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(ed_mode/="nonsu2")then
       write(LOGunit,"(A)")"set_nonsu2_components: called with ed_mode != nonsu2. return."
       return
    endif
    if(.not.dmft_bath_%status)stop "set_nonsu2_components: bath not allocated"
    select case(Nspin)          !Nhel=2 for ed_mode=nonsu2
    case (1)
       dmft_bath_%v(1,1,:,:) = dmft_bath_%w(1,1,:,:) !+ dmft_bath_%w(2,2,:,:))/2
       dmft_bath_%v(2,1,:,:) = dmft_bath_%w(1,2,:,:) !+ dmft_bath_%w(2,1,:,:))/2
       !
       dmft_bath_%e(1,:,:) = dmft_bath_%h(1,:,:)     !+ dmft_bath_%h(2,:,:))/2
    case (2)
       dmft_bath_%v(1,1,:,:) = dmft_bath_%w(1,1,:,:)
       dmft_bath_%v(2,1,:,:) = dmft_bath_%w(1,2,:,:)
       dmft_bath_%v(1,2,:,:) = dmft_bath_%w(2,1,:,:)
       dmft_bath_%v(2,2,:,:) = dmft_bath_%w(2,2,:,:)
       !
       dmft_bath_%e(1,:,:) = dmft_bath_%h(1,:,:)
       dmft_bath_%e(2,:,:) = dmft_bath_%h(2,:,:)
    end select
  end subroutine push_nonsu2_components




  !+-------------------------------------------------------------------+
  !PURPOSE  : Return the correct dimension of the 
  ! array used to store the bath in the chi2 fit procedure 
  ! fixed spin ed_mode=normal,superc
  ! both spins ed_mode=nonsu2
  !+-------------------------------------------------------------------+
  function get_array_bath_dimension() result(Ntrue)
    integer              :: Ntrue
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          !E_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][  \a][Nbath]
          !Nhel = 1
          Ntrue = 2*Nbath 
       case ("superc")
          !E_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !D_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][  \a][Nbath]
          !Nhel = 1
          Ntrue = 3*Nbath 
       case ("nonsu2")
          !H_{h,\a}(k)          [Nhel][  \a][Nbath]
          !W_{h,h`,\a}(k) [Nhel][Nhel][  \a][Nbath]
          !Nhel = 2
          Ntrue = (1+Nhel)*Nhel*Nbath
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          !E_{ 1, 1}(k)         [   1 ][  1 ][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][Norb][Nbath]
          !Nhel = 1
          Ntrue = (Norb+1)*Nbath
       case ("superc")
          !E_{\s, 1}(k)         [   1 ][  1 ][Nbath]
          !D_{\s, 1}(k)         [   1 ][  1 ][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][Norb][Nbath]
          !Nhel = 1
          Ntrue = (Norb+2)*Nbath
       case ("nonsu2")
          !H_{h,1}(k)           [Nhel][  1 ][Nbath]
          !W_{h,h`,\a}(k) [Nhel][Nhel][Norb][Nbath]
          !Nhel = 2
          Ntrue = (Nhel*Norb+1)*Nhel*Nbath
       endif
    end select
  end function get_array_bath_dimension

  !+-------------------------------------------------------------------+
  !PURPOSE  : check the correct dimension of the array used to store 
  ! the bath in the chi2 fit procedure 
  !+-------------------------------------------------------------------+
  function check_array_bath_dimension(a) result(bool)
    real(8),dimension(:),intent(in) :: a
    integer                         :: Ntrue
    logical                         :: bool
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          !E_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][  \a][Nbath]
          !Nhel = 1
          Ntrue = 2*Nbath 
       case ("superc")
          !E_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !D_{ 1,\a}(k)         [   1 ][  \a][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][  \a][Nbath]
          !Nhel = 1
          Ntrue = 3*Nbath 
       case ("nonsu2")
          !H_{h,\a}(k)          [Nhel][  \a][Nbath]
          !W_{h,h`,\a}(k) [Nhel][Nhel][  \a][Nbath]
          !Nhel = 2
          Ntrue = (1+Nhel)*Nhel*Nbath
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          !E_{ 1, 1}(k)         [   1 ][  1 ][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][Norb][Nbath]
          !Nhel = 1
          Ntrue = (Norb+1)*Nbath
       case ("superc")
          !E_{\s, 1}(k)         [   1 ][  1 ][Nbath]
          !D_{\s, 1}(k)         [   1 ][  1 ][Nbath]
          !V_{1, 1,\a}(k) [  1 ][   1 ][Norb][Nbath]
          !Nhel = 1
          Ntrue = (Norb+2)*Nbath
       case ("nonsu2")
          !H_{h,1}(k)           [Nhel][  1 ][Nbath]
          !W_{h,h`,\a}(k) [Nhel][Nhel][Norb][Nbath]
          !Nhel = 2
          Ntrue = (Nhel*Norb+1)*Nhel*Nbath
       endif
    end select
    bool = (size(a)==Ntrue)
  end function check_array_bath_dimension



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: A procedure that transforms an effective_bath type 
  ! into a simple array as used in the Chi**2 fit procedure.
  !+-----------------------------------------------------------------------------+!
  subroutine dmft_bath2array(dmft_bath,array,ispin,iorb)
    type(effective_bath) :: dmft_bath
    real(8),dimension(:) :: array
    integer,optional     :: ispin,iorb
    integer              :: ispin_,iorb_,i_,io_,ihel_,jhel_,stride_
    logical              :: check
    ispin_=1;if(present(ispin))ispin_=ispin
    iorb_ =1;if(present(iorb)) iorb_ =iorb
    check = check_array_bath_dimension(array)
    if(.not.check) stop "dmft_bath2array error: size[array] is incorrect"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          !2*Nbath == Nbath + Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%e(ispin_,iorb_,i_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%v(1, ispin_,iorb_,i_)
          enddo
          !
       case ("superc")
          !3*Nbath == Nbath + Nbath + Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%e(ispin_,iorb_,i_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%d(ispin_,iorb_,i_)
          enddo
          !Nbath:
          stride_=2*Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%v(1, ispin_,iorb_,i_)
          enddo
          !
       case ("nonsu2")
          !retrieve the nonSU2 components of the bath which go into the array
          call pull_nonsu2_components(dmft_bath) 
          !(1+Nhel)*Nhel*Nbath = Nhel*Nbath [H_h(k)] + Nhel*Nhel*Nbath [W_{ss`}(k)]
          !Nspin*Nbath
          stride_=0
          do ihel_=1,Nhel
             do i_=1,Nbath
                io_ = stride_ + i_ + (ihel_-1)*Nbath
                array(io_) =  dmft_bath%h(ihel_,iorb_,i_)
             enddo
          enddo
          !Nhel*Nbath
          stride_=Nhel*Nbath
          do ihel_=1,Nhel
             do jhel_=1,Nhel
                do i_=1,Nbath
                   io_ = stride_ +  i_ + (jhel_-1)*Nbath + (ihel_-1)*Nhel*Nbath
                   array(io_) = dmft_bath%w(ihel_,jhel_,iorb_,i_)
                enddo
             enddo
          enddo
          !
       end select
    case ("hybrid")
       select case(ed_mode)
       case default
          !(Norb+1)*Nbath = Nbath + Norb*Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%e(ispin_,1,i_)
          enddo
          !Norb*Nbath:
          stride_=Nbath
          do iorb_=1,Norb
             do i_=1,Nbath
                io_ = stride_ + i_ + (iorb_-1)*Nbath
                array(io_) = dmft_bath%v(1,ispin_,iorb_,i_)
             enddo
          enddo
          !
       case ("superc")
          !you should not reach this point.
          !stop "effective_bath2array error: ed_mode=superc, bath_type=hybrid not yet implemented"
          !(Norb+2)*Nbath = Nbath + Nbath + Norb*Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%e(ispin_,1,i_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             array(io_) = dmft_bath%d(ispin_,1,i_)
          enddo
          !Norb*Nbath
          stride_=2*Nbath
          do iorb_=1,Norb
             do i_=1,Nbath
                io_ = stride_ + i_ + (iorb_-1)*Nbath
                array(io_) = dmft_bath%v(1,ispin_,iorb_,i_)
             enddo
          enddo
          !
       case ("nonsu2")
          !retrieve the nonSU2 components of the bath which go into the array
          call pull_nonsu2_components(dmft_bath) 
          ! (Nhel*Norb+1)*Nhel*Nbath = Nhel*Nbath + Nhel*Nhel*Norb*Nbath
          !Nhel*Nbath
          stride_=0
          do ihel_=1,Nhel
             do i_=1,Nbath
                io_ = stride_ + i_ + (ihel_-1)*Nbath
                array(io_) = dmft_bath%h(ihel_,1,i_)
             enddo
          enddo
          !Nhel*Nhel*Norb*Nbath
          stride_ = Nhel*Nbath
          do ihel_=1,Nhel
             do jhel_=1,Nhel
                do iorb_=1,Norb
                   do i_=1,Nbath
                      io_ = stride_ + i_ + (iorb_-1)*Nbath + (jhel_-1)*Norb*Nbath + (ihel_-1)*Nhel*Norb*Nbath
                      array(io_) =  dmft_bath_%w(ihel_,jhel_,iorb_,i_)
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end subroutine dmft_bath2array








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: A procedure that transforms a simple array as used in the Chi**2 
  ! fit procedure into an effective_bath type.
  !+-----------------------------------------------------------------------------+!
  subroutine array2dmft_bath(array,dmft_bath,ispin,iorb)
    real(8),dimension(:) :: array
    type(effective_bath) :: dmft_bath
    integer,optional     :: ispin,iorb
    integer              :: ispin_,iorb_,i_,io_,ihel_,jhel_,stride_
    logical              :: check
    ispin_=1;if(present(ispin))ispin_=ispin
    iorb_ =1;if(present(iorb)) iorb_ =iorb
    check = check_array_bath_dimension(array)
    if(.not.check) stop "dmft_bath2array error: size[array] is incorrect"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          !2*Nbath == Nbath + Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%e(ispin_,iorb_,i_) = array(io_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%v(1, ispin_,iorb_,i_) = array(io_)
          enddo
          !
       case ("superc")
          !3*Nbath == Nbath + Nbath + Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%e(ispin_,iorb_,i_) = array(io_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%d(ispin_,iorb_,i_) = array(io_)
          enddo
          !Nbath:
          stride_=2*Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%v(1, ispin_,iorb_,i_) = array(io_)
          enddo
          !
       case ("nonsu2")
          !(1+Nhel)*Nhel*Nbath = Nhel*Nbath [H_h(k)] + Nhel*Nhel*Nbath [W_{ss`}(k)]
          !Nspin*Nbath
          stride_=0
          do ihel_=1,Nspin
             do i_=1,Nbath
                io_ = stride_ + i_ + (ihel_-1)*Nbath
                dmft_bath%h(ihel_,iorb_,i_) = array(io_)
             enddo
          enddo
          !Nhel*Nhel*Nbath
          stride_=Nhel*Nbath
          do ihel_=1,Nhel
             do jhel_=1,Nhel
                do i_=1,Nbath
                   io_ = stride_ +  i_ + (jhel_-1)*Nbath + (ihel_-1)*Nhel*Nbath
                   dmft_bath%w(ihel_,jhel_,iorb_,i_) = array(io_)
                enddo
             enddo
          enddo
          call push_nonsu2_components(dmft_bath)
          !
          !
       end select
    case ("hybrid")
       select case(ed_mode)
       case default
          !(Norb+1)*Nbath = Nbath + Norb*Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%e(ispin_,1,i_) = array(io_) 
          enddo
          !Norb*Nbath:
          stride_=Nbath
          do iorb_=1,Norb
             do i_=1,Nbath
                io_ = stride_ + i_ + (iorb_-1)*Nbath
                dmft_bath%v(1,ispin_,iorb_,i_) = array(io_) 
             enddo
          enddo
          !
       case ("superc")
          !you should not reach this point.
          !stop "effective_bath2array error: ed_mode=superc, bath_type=hybrid not yet implemented"
          !(Norb+2)*Nbath = Nbath + Nbath + Norb*Nbath
          !Nbath:
          stride_=0
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%e(ispin_,1,i_) = array(io_)
          enddo
          !Nbath:
          stride_=Nbath
          do i_=1,Nbath
             io_ = stride_ + i_
             dmft_bath%d(ispin_,1,i_) = array(io_) 
          enddo
          !Norb*Nbath
          stride_=2*Nbath
          do iorb_=1,Norb
             do i_=1,Nbath
                io_ = stride_ + i_ + (iorb_-1)*Nbath
                dmft_bath%v(1,ispin_,iorb_,i_) = array(io_)
             enddo
          enddo
          !
       case ("nonsu2")
          ! (Nhel*Norb+1)*Nhel*Nbath = Nhel*Nbath + Nhel*Nhel*Norb*Nbath
          !Nhel*Nbath
          stride_=0
          do ihel_=1,Nhel
             do i_=1,Nbath
                io_ = stride_ + i_ + (ihel_-1)*Nbath
                dmft_bath%h(ihel_,1,i_) = array(io_) 
             enddo
          enddo
          !Nhel*Nhel*Norb*Nbath
          stride_ = Nhel*Nbath
          do ihel_=1,Nhel
             do jhel_=1,Nhel
                do iorb_=1,Norb
                   do i_=1,Nbath
                      io_ = stride_ + i_ + (iorb_-1)*Nbath + (jhel_-1)*Norb*Nbath + (ihel_-1)*Nhel*Norb*Nbath
                      dmft_bath_%w(ihel_,jhel_,iorb_,i_) = array(io_)
                   enddo
                enddo
             enddo
          enddo
          call push_nonsu2_components(dmft_bath)
          !
       end select
       !
    end select
  end subroutine array2dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  ! DELTA_BATH_MATS:
  include "ed_delta_bath_mats.f90"


  !+-------------------------------------------------------------------+
  !PURPOSE  :  compute the gradient of the hybridization function
  ! at a point x
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  ! GRAD_DELTA_BATH_MATS:
  include "ed_grad_delta_bath_mats.f90"


  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  ! WEISS_BATH_MATS:
  include "ed_weiss_bath_mats.f90"


  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  ! DELTA_BATH_REAL:
  include "ed_delta_bath_real.f90"



















  ! function delta_bath_irred_mats(ispin,iorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   complex(8),intent(in) :: x
  !   integer,intent(in)    :: iorb,ispin
  !   complex(8)            :: fg
  !   if(.not.ed_supercond)then
  !      fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
  !   else
  !      fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
  !           (dimag(x)**2+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  !   endif
  ! end function delta_bath_irred_mats
  ! !
  ! function fdelta_bath_irred_mats(ispin,iorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   complex(8),intent(in) :: x
  !   integer,intent(in)    :: iorb,ispin
  !   complex(8)            :: fg
  !   fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
  !        (dimag(x)**2+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  ! end function fdelta_bath_irred_mats
  !  
  ! !NORMAL/IRREDUCIBLE BATH:
  ! !Real axis:
  ! function delta_bath_irred_real(ispin,iorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   integer,intent(in)    :: iorb,ispin
  !   complex(8),intent(in) :: x
  !   complex(8)            :: fg
  !   if(.not.ed_supercond)then
  !      fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
  !   else
  !      fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
  !           ( x*(-x)+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  !   endif
  ! end function delta_bath_irred_real
  ! !
  ! function fdelta_bath_irred_real(ispin,iorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   integer,intent(in)    :: iorb,ispin
  !   complex(8),intent(in) :: x
  !   complex(8)            :: fg
  !   fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
  !        ( x*(-x) + dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  ! end function fdelta_bath_irred_real
  !
  !
  ! !<ACTHUNG/TODO: extend the irred expressions for the Delta to the hybrid case.
  ! !HYBRIDIZED BATH:
  ! !Matsubara:
  ! function delta_bath_hybrd_mats(ispin,iorb,jorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   complex(8),intent(in) :: x
  !   integer,intent(in)    :: iorb,jorb,ispin
  !   complex(8)            :: fg
  !   if(.not.ed_supercond)then
  !      fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/(x - dmft_bath_%e(ispin,1,1:Nbath)))
  !   else
  !      fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x + dmft_bath_%e(ispin,1,1:Nbath))/&
  !           (dimag(x)**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  !   endif
  ! end function delta_bath_hybrd_mats
  ! !
  ! function fdelta_bath_hybrd_mats(ispin,iorb,jorb,x,dmft_bath_) result(fg)
  !   type(effective_bath)  :: dmft_bath_
  !   complex(8),intent(in) :: x
  !   integer,intent(in)    :: iorb,jorb,ispin
  !   complex(8)            :: fg
  !   fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
  !        ( dimag(x)**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  ! end function fdelta_bath_hybrd_mats
  !
  !   !HYBRIDIZED BATH:
  !   !Real-axis:
  !   function delta_bath_hybrd_real(ispin,iorb,jorb,x,dmft_bath_) result(fg)
  !     type(effective_bath)  :: dmft_bath_
  !     integer,intent(in)    :: iorb,jorb,ispin
  !     complex(8)            :: fg,x
  !     if(.not.ed_supercond)then
  !        fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)&
  !             /(x-dmft_bath_%e(ispin,1,1:Nbath)))
  !     else
  !        fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x+dmft_bath_%e(ispin,1,1:Nbath))/&
  !             ( -x**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  !     endif
  !   end function delta_bath_hybrd_real
  !   !
  !   function fdelta_bath_hybrd_real(ispin,iorb,jorb,x,dmft_bath_) result(fg)
  !     type(effective_bath)  :: dmft_bath_
  !     integer,intent(in)    :: iorb,jorb,ispin
  !     complex(8)            :: fg,x
  !     fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
  !          ( -x**2  + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  !   end function fdelta_bath_hybrd_real









END MODULE ED_BATH
