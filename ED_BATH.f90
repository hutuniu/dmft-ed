MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none

  private

  interface save_bath
     module procedure save_dmft_bath
     module procedure save_user_bath
  end interface save_bath

  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !procedures:
  public :: allocate_bath              !INTERNAL (for effective_bath)
  public :: deallocate_bath            !INTERNAL (for effective_bath)
  public :: get_bath_size              !PUBLIC   (for user_bath)
  public :: check_bath_dimension       !PUBLIC   (for user_bath)
  public :: init_bath_ed               !INTERNAL (for effective_bath)
  public :: write_bath                 !INTERNAL (for effective_bath)
  public :: save_bath                  !INTERNAL (for effective_bath)
  public :: set_bath                   !INTERNAL (for effective_bath)
  public :: copy_bath                  !INTERNAL (for effective_bath)
  public :: break_symmetry_bath        !PUBLIC   (for user_bath)
  public :: spin_symmetrize_bath       !PUBLIC   (for user_bath)
  public :: ph_symmetrize_bath         !PUBLIC   (for user_bath)
  public :: ph_trans_bath              !PUBLIC   (for user_bath)
  public :: enforce_normal_bath        !PUBLIC   (for user_bath)



  interface get_Whyb_matrix
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix
  public :: get_Whyb_matrix


  !functions:
  !
  !\Delta hybridization function Matsubara
  !NORMAL
  interface delta_bath_mats
     module procedure delta_bath_mats_main
     module procedure delta_bath_mats_ispin_jspin
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb
     module procedure delta_bath_mats_main_
     module procedure delta_bath_mats_ispin_jspin_
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb_
  end interface delta_bath_mats
  !
  !ANOMALOUS
  interface fdelta_bath_mats
     module procedure fdelta_bath_mats_main
     module procedure fdelta_bath_mats_ispin_jspin
     module procedure fdelta_bath_mats_ispin_jspin_iorb_jorb
     module procedure fdelta_bath_mats_main_
     module procedure fdelta_bath_mats_ispin_jspin_
     module procedure fdelta_bath_mats_ispin_jspin_iorb_jorb_
  end interface fdelta_bath_mats
  !
  public :: delta_bath_mats
  public :: fdelta_bath_mats



  !\Delta hybridization function Real
  !NORMAL
  interface delta_bath_real
     module procedure delta_bath_real_main
     module procedure delta_bath_real_ispin_jspin
     module procedure delta_bath_real_ispin_jspin_iorb_jorb
     module procedure delta_bath_real_main_
     module procedure delta_bath_real_ispin_jspin_
     module procedure delta_bath_real_ispin_jspin_iorb_jorb_
  end interface delta_bath_real
  !
  !ANOMALOUS
  interface fdelta_bath_real
     module procedure fdelta_bath_real_main
     module procedure fdelta_bath_real_ispin_jspin
     module procedure fdelta_bath_real_ispin_jspin_iorb_jorb
     module procedure fdelta_bath_real_main_
     module procedure fdelta_bath_real_ispin_jspin_
     module procedure fdelta_bath_real_ispin_jspin_iorb_jorb_
  end interface fdelta_bath_real
  !
  public :: delta_bath_real
  public :: fdelta_bath_real




  !Non-interacting Green's function Matsubara
  !NORMAL
  interface g0and_bath_mats
     module procedure g0and_bath_mats_main
     module procedure g0and_bath_mats_ispin_jspin
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb
     module procedure g0and_bath_mats_main_
     module procedure g0and_bath_mats_ispin_jspin_
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface g0and_bath_mats
  !
  !ANOMALOUS
  interface f0and_bath_mats
     module procedure f0and_bath_mats_main
     module procedure f0and_bath_mats_ispin_jspin
     module procedure f0and_bath_mats_ispin_jspin_iorb_jorb
     module procedure f0and_bath_mats_main_
     module procedure f0and_bath_mats_ispin_jspin_
     module procedure f0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface f0and_bath_mats
  !
  public :: g0and_bath_mats
  public :: f0and_bath_mats




  !Inverse Non-interacting Green's function Matsubara
  !NORMAL
  interface invg0_bath_mats
     module procedure invg0_bath_mats_main
     module procedure invg0_bath_mats_ispin_jspin
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb
     module procedure invg0_bath_mats_main_
     module procedure invg0_bath_mats_ispin_jspin_
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb_
  end interface invg0_bath_mats
  !
  !ANOMALOUS
  interface invf0_bath_mats
     module procedure invf0_bath_mats_main
     module procedure invf0_bath_mats_ispin_jspin
     module procedure invf0_bath_mats_ispin_jspin_iorb_jorb
     module procedure invf0_bath_mats_main_
     module procedure invf0_bath_mats_ispin_jspin_
     module procedure invf0_bath_mats_ispin_jspin_iorb_jorb_
  end interface invf0_bath_mats
  !
  public :: invg0_bath_mats
  public :: invf0_bath_mats





  !Non-interacting Green's function Real-axis
  !NORMAL
  interface g0and_bath_real
     module procedure g0and_bath_real_main
     module procedure g0and_bath_real_ispin_jspin
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb
     module procedure g0and_bath_real_main_
     module procedure g0and_bath_real_ispin_jspin_
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb_
  end interface g0and_bath_real
  !
  !ANOMALOUS
  interface f0and_bath_real
     module procedure f0and_bath_real_main
     module procedure f0and_bath_real_ispin_jspin
     module procedure f0and_bath_real_ispin_jspin_iorb_jorb
     module procedure f0and_bath_real_main_
     module procedure f0and_bath_real_ispin_jspin_
     module procedure f0and_bath_real_ispin_jspin_iorb_jorb_
  end interface f0and_bath_real
  !
  public :: g0and_bath_real
  public :: f0and_bath_real




  !Inverse Non-interacting Green's function Real-axis
  !NORMAL
  interface invg0_bath_real
     module procedure invg0_bath_real_main
     module procedure invg0_bath_real_ispin_jspin
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb
     module procedure invg0_bath_real_main_
     module procedure invg0_bath_real_ispin_jspin_
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb_
  end interface invg0_bath_real
  !
  !ANOMALOUS
  interface invf0_bath_real
     module procedure invf0_bath_real_main
     module procedure invf0_bath_real_ispin_jspin
     module procedure invf0_bath_real_ispin_jspin_iorb_jorb
     module procedure invf0_bath_real_main_
     module procedure invf0_bath_real_ispin_jspin_
     module procedure invf0_bath_real_ispin_jspin_iorb_jorb_
  end interface invf0_bath_real
  !
  public :: invg0_bath_real
  public :: invf0_bath_real



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(dmft_bath_%status)call deallocate_bath(dmft_bath_)
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default                                 !normal [N,Sz]
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       case ("superc")                              !superc [Sz] 
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%d(Nspin,Norb,Nbath))  !local SC order parameters the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       case ("nonsu2")                              !nonsu2 [N]
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
          allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default                                 !normal  [N,Sz]
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       case ("superc")                              !superc  [Sz]
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%d(Nspin,1,Nbath))     !local SC order parameters the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       case ("nonsu2")                              !nonsu2 case [N] qn
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
          allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    end select
    dmft_bath_%status=.true.
  end subroutine allocate_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%d))deallocate(dmft_bath_%d)
    if(allocated(dmft_bath_%v))deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%u))deallocate(dmft_bath_%u)
    dmft_bath_%status=.false.
  end subroutine deallocate_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program:
  ! - user bath
  ! - chi2 bath
  !+-------------------------------------------------------------------+
  function get_bath_size() result(bath_size)
    integer :: bath_size
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nspin*Norb*Nbath + Nspin*Norb*Nbath
       case ("superc")
          !( e [Nspin][Norb][Nbath] + d [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nspin*Norb*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
       case ("nonsu2")
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
          bath_size = Nspin*Norb*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nspin*Nbath + Nspin*Norb*Nbath
       case ("superc")
          !(e [Nspin][1][Nbath] + d [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nspin*Nbath + Nspin*Nbath + Nspin*Norb*Nbath
       case ("nonsu2")
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
          bath_size = Nspin*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
       end select
       !
    end select
    !
  end function get_bath_size








  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  ! - user bath
  ! - chi2 bath
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    real(8),dimension(:) :: bath_
    integer              :: Ntrue
    logical              :: bool
    !
    Ntrue = get_bath_size()
    !
    bool  = (size(bath_) == Ntrue)
    !
  end function check_bath_dimension









  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed(dmft_bath_,hwband_)
    type(effective_bath) :: dmft_bath_
    real(8)              :: hwband_,wband_
    integer              :: i,iorb,ispin,unit,flen,Nh
    logical              :: IOfile
    real(8)              :: de
    if(.not.dmft_bath_%status)stop "init_bath: bath not allocated"
    !Generating the bath anyway, then you may want to read it to 
    !update some entries. This way you can restart even 
    !from different Ns calculation, this is better than 
    !start from a complete guess
    !
    !This is valid for all ed_modes: normal,superc,nonsu2
    !
    !Get energies:
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
    !Get spin-keep yhbridizations
    do i=1,Nbath
       dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))
       if(bath_type=="hybrid") then
          dmft_bath_%v(:,:,i)=0.1d0
       endif
    enddo
    !Get SC amplitudes
    if(ed_mode=="superc")dmft_bath_%d(:,:,:) = deltasc
    !Get spin-flip hybridizations
    if(ed_mode=="nonsu2")dmft_bath_%u(:,:,:) = dmft_bath_%v(:,:,:)
    !
    !Read from file if exist:
    !
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
       !
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
       read(unit,*)
       !
       select case(bath_type)
       case default
          !
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%d(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     dmft_bath_%u(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
          !
       case ('hybrid')
          !
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     dmft_bath_%d(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     dmft_bath_%u(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          end select
          !
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
    integer              :: i,ispin,iorb
    if(ED_MPI_ID==0)then
       unit_=LOGfile;if(present(unit))unit_=unit
       if(.not.dmft_bath_%status)stop "WRITE_BATH: bath not allocated"
       select case(bath_type)
       case default
          !
          select case(ed_mode)
          case default
             write(unit_,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit_,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             write(unit_,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)) ,&
                  "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit_,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%d(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("nonsu2")
             write(unit_,"(90(A21,1X))")&
                  ((&
                  "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  iorb=1,Norb), ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     dmft_bath_%u(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
          !
       case('hybrid')
          !
          select case(ed_mode)
          case default
             write(unit_,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
                  ispin=1,Nspin)
             do i=1,Nbath
                write(unit_,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),&
                     (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("superc")
             write(unit_,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  "Dk_s"//reg(txtfy(ispin)) ,&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
                  ispin=1,Nspin)
             do i=1,Nbath
                write(unit_,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),&
                     dmft_bath_%d(ispin,1,i),&
                     (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("nonsu2")
             write(unit_,"(90(A21,1X))")(&
                  "#Ek_s"//reg(txtfy(ispin)),&
                  (&
                  "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  iorb=1,Norb),&
                  ispin=1,Nspin)
             do i=1,Nbath
                write(unit_,"(90(F21.12,1X))")(&
                     dmft_bath_%e(ispin,1,i),    &
                     (dmft_bath_%v(ispin,iorb,i),dmft_bath_%u(ispin,iorb,i),iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          end select
          !
       end select
    endif
  end subroutine write_bath







  !+-------------------------------------------------------------------+
  !PURPOSE  : save the bath to a given file using the write bath
  ! procedure and formatting: 
  !+-------------------------------------------------------------------+
  subroutine save_dmft_bath(dmft_bath_,file,used)
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
    if(ED_MPI_ID==0)then
       if(.not.dmft_bath_%status)stop "save_used_bath error: bath is not allocated"
       used_=.false.;if(present(used))used_=used
       extension=".restart";if(used_)extension=".used"
       file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
       if(present(file))file_=reg(file)
       unit_=free_unit()
       open(unit_,file=reg(file_))
       call write_bath(dmft_bath_,unit_)
       close(unit_)
    endif
  end subroutine save_dmft_bath

  subroutine save_user_bath(bath_,file,used)
    real(8),dimension(:)      :: bath_
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_,check
    character(len=16)         :: extension
    integer                   :: unit_
    if(ED_MPI_ID==0)then
       check= check_bath_dimension(bath_)
       if(.not.check)stop "save_used_bath error: wrong bath dimensions"
       call allocate_bath(dmft_bath_)
       call set_bath(bath_,dmft_bath_)
       used_=.false.;if(present(used))used_=used
       extension=".restart";if(used_)extension=".used"
       file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
       if(present(file))file_=reg(file)
       unit_=free_unit()
       open(unit_,file=reg(file_))
       call write_bath(dmft_bath_,unit_)
       close(unit_)
       call deallocate_bath(dmft_bath_)
    endif
  end subroutine save_user_bath







  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath_,dmft_bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i,iorb,jorb,ispin,jspin,stride,io
    logical                :: check
    if(.not.dmft_bath_%status)stop "SET_BATH: bath not allocated"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_bath: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
          !
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
    case ('hybrid')
       !
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
    end select
  end subroutine set_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine copy_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:)   :: bath_
    integer                :: iorb,jorb,ispin,jspin,stride,io,i
    logical                :: check
    if(.not.dmft_bath_%status)stop "COPY_BATH: bath not allocated"
    check=check_bath_dimension(bath_)
    if(.not.check)stop "copy_bath: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
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
                   bath_(io) = dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
    case ('hybrid')
       !
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
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
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath_%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
    end select
  end subroutine copy_bath







  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build up the all-spin hybridization matrix W_{ss`}
  !+-----------------------------------------------------------------------------+!
  function get_Whyb_matrix_1orb(v,u) result(w)
    real(8),dimension(Nspin,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Nbath) :: w
    integer                              :: ispin
    !
    if(ed_spin_sym)then
       do ispin=1,Nspin
          w(ispin,ispin,:) = v(1,:)
       enddo
       w(1,Nspin,:) = u(1,:)
       w(Nspin,1,:) = u(1,:)
    else
       do ispin=1,Nspin
          w(ispin,ispin,:) = v(ispin,:)
       enddo
       w(1,Nspin,:) = u(1,:)
       w(Nspin,1,:) = u(2,:)
    endif
    !
  end function get_Whyb_matrix_1orb

  function get_Whyb_matrix_Aorb(v,u) result(w)
    real(8),dimension(Nspin,Norb,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    if(ed_spin_sym)then
       do ispin=1,Nspin
          w(ispin,ispin,:,:) = v(1,:,:)
       enddo
       w(1,Nspin,:,:) = u(1,:,:)
       w(Nspin,1,:,:) = u(1,:,:)
    else
       do ispin=1,Nspin
          w(ispin,ispin,:,:) = v(ispin,:,:)
       enddo
       w(1,Nspin,:,:) = u(1,:,:)
       w(Nspin,1,:,:) = u(2,:,:)
    endif
    !
  end function get_Whyb_matrix_Aorb

  function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
    type(effective_bath)                      :: dmft_bath_
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    if(ed_spin_sym)then
       do ispin=1,Nspin
          w(ispin,ispin,:,:) = dmft_bath_%v(1,:,:)
       enddo
       w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
       w(Nspin,1,:,:) = dmft_bath_%u(1,:,:)
    else
       do ispin=1,Nspin
          w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
       enddo
       w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
       w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
    endif
    !
  end function get_Whyb_matrix_dmft_bath







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
  ! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization 
  !    matrix
  ! - given a dmft bath pull/push the nonsu2 components
  !+-------------------------------------------------------------------+
  subroutine break_symmetry_bath(bath_,field,sign,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    save_=.false.;if(present(save))save_=save
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    if(save_)call save_bath(dmft_bath_)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine break_symmetry_bath

  subroutine spin_symmetrize_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    save_=.false.;if(present(save))save_=save
    if(Nspin==1)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
    dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
    if(ed_mode=="superc")dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    if(save_)call save_bath(dmft_bath_)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath

  subroutine ph_symmetrize_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.false.;if(present(save))save_=save
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    if(mod(Nbath,2)==0)then
       do i=1,Nbath/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
    else
       do i=1,(Nbath-1)/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
       dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
    endif
    if(save_)call save_bath(dmft_bath_)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath

  subroutine ph_trans_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    type(effective_bath)   :: tmp_dmft_bath
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.false.;if(present(save))save_=save
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
          if(ed_mode=="superc")dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)
       case(1)
          dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)
       case(2)
          tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
          tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
          dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
          tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
          tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
          dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,i)=dmft_bath_%d(:,:,i)          
       end select
    end do
    if(save_)call save_bath(dmft_bath_)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine ph_trans_bath

  subroutine enforce_normal_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.false.;if(present(save))save_=save
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    if(ed_mode=="superc")dmft_bath_%d(:,:,:)=0.d0
    if(save_)call save_bath(dmft_bath_)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine enforce_normal_bath








  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  include "ed_bath_delta_bath_mats.f90"  ! DELTA_BATH_MATS:
  include "ed_bath_delta_bath_real.f90"  ! DELTA_BATH_REAL:


  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  include "ed_bath_g0and_bath_mats.f90"  ! G0and_BATH_MATS:
  include "ed_bath_g0and_bath_real.f90"  ! G0and_BATH_REAL:


  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the inverse G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  include "ed_bath_invg0_bath_mats.f90"  ! invG0_BATH_MATS:
  include "ed_bath_invg0_bath_real.f90"  ! invG0_BATH_REAL:












  ! function check_chi2_bath_dimension(a) result(bool)
  !   real(8),dimension(:),intent(in) :: a
  !   integer                         :: Ntrue
  !   logical                         :: bool
  !   !
  !   Ntrue = get_chi2_bath_size()
  !   !
  !   bool  = (size(a) == Ntrue)
  !   !
  ! end function check_chi2_bath_dimension



  ! function get_chi2_bath_size() result(bath_size)
  !   integer :: bath_size
  !   !
  !   select case(bath_type)
  !   case default
  !      !
  !      select case(ed_mode)
  !      case default             !per spin/per orb
  !         !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !         !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !         bath_size = Nbath + Nbath
  !      case ("superc")          !per spin/per orb
  !         !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !         !D_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !         !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !         bath_size = Nbath + Nbath + Nbath
  !      case ("nonsu2")          !per orb
  !         !H_{:,\a}(:)   [Nspin][ 1 ][Nbath]
  !         !W_{:,:,\a}(:) [Nspin][Nspin][ 1 ][Nbath]
  !         bath_size = Nspin*Nbath + Nspin*Nspin*Nbath
  !      end select
  !      !
  !   case('hybrid')
  !      !
  !      select case(ed_mode)
  !      case default
  !         !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !         !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
  !         bath_size = Nbath + Norb*Nbath
  !      case ("superc")
  !         !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !         !D_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !         !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
  !         bath_size = Nbath + Nbath + Norb*Nbath
  !      case ("nonsu2")
  !         !H_{:,1}(:)   [Nspin][ 1 ][Nbath]
  !         !W_{:,:,:}(:) [Nspin][Nspin][Norb][Nbath]
  !         bath_size = Nspin*Nbath + Nspin*Nspin*Norb*Nbath
  !      end select
  !      !
  !   end select
  !   !
  ! end function get_chi2_bath_size



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  :  compute the gradient of the hybridization function
  ! ! at a point x
  ! ! type(effective_bath) :: dmft_bath
  ! ! OR
  ! ! real(8),dimension(:) :: bath_array
  ! !+-------------------------------------------------------------------+
  ! ! GRAD_DELTA_BATH_MATS:
  ! include "ed_bath_grad_delta_bath_mats.f90"
  ! !\Nabla\Delta NORMAL hybridization function gradient Matsubara
  ! !NORMAL
  ! interface grad_delta_bath_mats
  !    module procedure grad_delta_bath_mats_main
  !    module procedure grad_delta_bath_mats_ispin_jspin
  !    module procedure grad_delta_bath_mats_ispin_jspin_iorb_jorb
  !    module procedure grad_delta_bath_mats_main_
  !    module procedure grad_delta_bath_mats_ispin_jspin_
  !    module procedure grad_delta_bath_mats_ispin_jspin_iorb_jorb_
  ! end interface grad_delta_bath_mats
  ! !
  ! !ANOMALOUS
  ! interface grad_fdelta_bath_mats
  !    module procedure grad_fdelta_bath_mats_main
  !    module procedure grad_fdelta_bath_mats_ispin_jspin
  !    module procedure grad_fdelta_bath_mats_ispin_jspin_iorb_jorb
  !    module procedure grad_fdelta_bath_mats_main_
  !    module procedure grad_fdelta_bath_mats_ispin_jspin_
  !    module procedure grad_fdelta_bath_mats_ispin_jspin_iorb_jorb_
  ! end interface grad_fdelta_bath_mats
  ! !
  ! public :: grad_delta_bath_mats
  ! public :: grad_fdelta_bath_mats


  ! !+-----------------------------------------------------------------------------+!
  ! !PURPOSE: A procedure that transforms an effective_bath type 
  ! ! into a simple array as used in the Chi**2 fit procedure.
  ! !+-----------------------------------------------------------------------------+!
  ! subroutine dmft_bath2chi2_bath(dmft_bath,array,ispin,iorb)
  !   type(effective_bath) :: dmft_bath
  !   real(8),dimension(:) :: array
  !   integer,optional     :: ispin,iorb
  !   integer              :: ispin_,jspin_,iorb_,i_,io_,stride_
  !   logical              :: check
  !   !
  !   ispin_=1;if(present(ispin))ispin_=ispin
  !   iorb_ =1;if(present(iorb)) iorb_ =iorb
  !   !
  !   check = check_chi2_bath_dimension(array)
  !   if(.not.check) stop "dmft_bath2chi2_bath error: size[array] is incorrect"
  !   !
  !   select case(bath_type)
  !   case default
  !      select case(ed_mode)
  !      case default             !per spin/ per orb
  !         !2*Nbath == Nbath + Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%e(ispin_,iorb_,i_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%v(ispin_,iorb_,i_)
  !         enddo
  !         !
  !      case ("superc")          !per spin/ per orb
  !         !3*Nbath == Nbath + Nbath + Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%e(ispin_,iorb_,i_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%d(ispin_,iorb_,i_)
  !         enddo
  !         stride_=2*Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%v(ispin_,iorb_,i_)
  !         enddo
  !         !
  !      end select
  !      !
  !   case ("hybrid")
  !      !
  !      select case(ed_mode)     !per spin
  !      case default
  !         !Nbath + Norb*Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%e(ispin_,1,i_)
  !         enddo
  !         stride_=Nbath
  !         do iorb_=1,Norb
  !            do i_=1,Nbath
  !               io_ = stride_ + i_ + (iorb_-1)*Nbath
  !               array(io_) = dmft_bath%v(ispin_,iorb_,i_)
  !            enddo
  !         enddo
  !         !
  !      case ("superc")          !per spin
  !         !Nbath + Nbath + Norb*Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%e(ispin_,1,i_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            array(io_) = dmft_bath%d(ispin_,1,i_)
  !         enddo
  !         stride_=2*Nbath
  !         do iorb_=1,Norb
  !            do i_=1,Nbath
  !               io_ = stride_ + i_ + (iorb_-1)*Nbath
  !               array(io_) = dmft_bath%v(ispin_,iorb_,i_)
  !            enddo
  !         enddo
  ! case ("nonsu2")          !none
  !    !Nspin*Nbath + Nspin*Nspin*Norb*Nbath
  !   stride_=0
  !   do ispin_=1,Nspin
  !      do i_=1,Nbath
  !         io_ = stride_ + i_ + (ispin-1)*Nbath
  !         array(io_) = dmft_bath%e(ispin_,1,i_)
  !      enddo
  !   enddo
  !   stride_ = Nspin*Nbath
  !   do ispin_=1,Nspin
  !      do jspin_=1,Nspin
  !         do iorb_=1,Norb
  !            do i_=1,Nbath
  !               io_ = stride_ + i_ + (iorb_-1)*Nbath + (jspin_-1)*Norb*Nbath + (ispin_-1)*Nspin*Norb*Nbath
  !               array(io_) =  dmft_bath%w(ispin_,jspin_,iorb_,i_)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !         !
  !      end select
  !      !
  !   end select
  ! end subroutine dmft_bath2chi2_bath


  ! !+-----------------------------------------------------------------------------+!
  ! !PURPOSE: A procedure that transforms a simple array as used in the Chi**2 
  ! ! fit procedure into an effective_bath type.
  ! !+-----------------------------------------------------------------------------+!
  ! subroutine chi2_bath2dmft_bath(array,dmft_bath,ispin,iorb)
  !   real(8),dimension(:) :: array
  !   type(effective_bath) :: dmft_bath
  !   integer,optional     :: ispin,iorb
  !   integer              :: ispin_,jspin_,iorb_,i_,io_,stride_
  !   logical              :: check
  !   !
  !   ispin_=1;if(present(ispin))ispin_=ispin
  !   iorb_ =1;if(present(iorb)) iorb_ =iorb
  !   !
  !   check = check_chi2_bath_dimension(array)
  !   if(.not.check) stop "chi2_bath2dmft_bath error: size[array] is incorrect"
  !   !
  !   select case(bath_type)
  !   case default
  !      !
  !      select case(ed_mode)     !per spin / per orb
  !      case default
  !         !Nbath + Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%e(ispin_,iorb_,i_) = array(io_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%v(ispin_,iorb_,i_) = array(io_)
  !         enddo
  !         !
  !      case ("superc")          !per spin / per orb
  !         !Nbath + Nbath + Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%e(ispin_,iorb_,i_) = array(io_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%d(ispin_,iorb_,i_) = array(io_)
  !         enddo
  !         stride_=2*Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%v(ispin_,iorb_,i_) = array(io_)
  !         enddo
  !         !
  !      end select
  !      !
  !   case ("hybrid")
  !      !
  !      select case(ed_mode)     !per spin
  !      case default
  !         !Nbath + Norb*Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%e(ispin_,1,i_) = array(io_) 
  !         enddo
  !         stride_=Nbath
  !         do iorb_=1,Norb
  !            do i_=1,Nbath
  !               io_ = stride_ + i_ + (iorb_-1)*Nbath
  !               dmft_bath%v(ispin_,iorb_,i_) = array(io_) 
  !            enddo
  !         enddo
  !         !
  !      case ("superc")    !per spin
  !         !Nbath + Nbath + Norb*Nbath
  !         stride_=0
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%e(ispin_,1,i_) = array(io_)
  !         enddo
  !         stride_=Nbath
  !         do i_=1,Nbath
  !            io_ = stride_ + i_
  !            dmft_bath%d(ispin_,1,i_) = array(io_) 
  !         enddo
  !         stride_=2*Nbath
  !         do iorb_=1,Norb
  !            do i_=1,Nbath
  !               io_ = stride_ + i_ + (iorb_-1)*Nbath
  !               dmft_bath%v(ispin_,iorb_,i_) = array(io_)
  !            enddo
  !         enddo
  !         
  !    case ("nonsu2")          !none
  ! !Nspin*Nbath + Nspin*Nspin*Norb*Nbath
  ! stride_=0
  ! do ispin_=1,Nspin
  !    do i_=1,Nbath
  !        io_ = stride_ + i_ + (ispin_-1)*Nbath
  !        dmft_bath%e(ispin_,1,i_) = array(io_) 
  !     enddo
  !  enddo
  !  stride_ = Nspin*Nbath
  !  do ispin_=1,Nspin
  !     do jspin_=1,Nspin
  !        do iorb_=1,Norb
  !           do i_=1,Nbath
  !              io_ = stride_ + i_ + (iorb_-1)*Nbath + (jspin_-1)*Norb*Nbath + (ispin_-1)*Nspin*Norb*Nbath
  !              dmft_bath%w(ispin_,jspin_,iorb_,i_) = array(io_)
  !           enddo
  !        enddo
  !     enddo
  !  enddo
  !
  !      end select
  !      !
  !   end select
  ! end subroutine chi2_bath2dmft_bath


END MODULE ED_BATH
