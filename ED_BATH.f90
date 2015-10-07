MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none

  private

  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !
  !DMFT BATH procedures:
  public :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)
  !
  !USER BATH procedures
  public :: get_size_bath                    !PUBLIC   (for user_bath)
  public :: get_component_size_bath          !PUBLIC   (for user_bath)
  public :: get_spin_component_size_bath     !PUBLIC   (for user_bath)
  public :: get_orb_component_size_bath      !PUBLIC   (for user_bath)
  public :: get_spin_orb_component_size_bath !PUBLIC   (for user_bath)
  public :: check_size_bath                  !PUBLIC   (for user_bath)
  public :: get_component_bath               !PUBLIC   (for user_bath)
  public :: get_spin_component_bath          !PUBLIC   (for user_bath)
  public :: get_orb_component_bath           !PUBLIC   (for user_bath)
  public :: get_spin_orb_component_bath      !PUBLIC   (for user_bath)
  public :: set_component_bath               !PUBLIC   (for user_bath)
  public :: set_spin_component_bath          !PUBLIC   (for user_bath)
  public :: set_orb_component_bath           !PUBLIC   (for user_bath)
  public :: set_spin_orb_component_bath      !PUBLIC   (for user_bath)
  public :: copy_component_bath              !PUBLIC   (for user_bath)
  public :: copy_spin_component_bath         !PUBLIC   (for user_bath)
  public :: copy_orb_component_bath          !PUBLIC   (for user_bath)
  public :: copy_spin_orb_component_bath     !PUBLIC   (for user_bath)
  public :: save_bath                        !PUBLIC   (for user_bath)
  !explicit symmetries:
  public :: break_symmetry_bath              !PUBLIC   (for user_bath)
  public :: spin_symmetrize_bath             !PUBLIC   (for user_bath)
  public :: ph_symmetrize_bath               !PUBLIC   (for user_bath)
  public :: ph_trans_bath                    !PUBLIC   (for user_bath)
  public :: enforce_normal_bath              !PUBLIC   (for user_bath)


  interface get_bath_size
     module procedure get_size_bath
  end interface get_bath_size
  public :: get_bath_size

  interface check_bath_dimension
     module procedure check_size_bath
  end interface check_bath_dimension
  public :: check_bath_dimension

  interface get_Whyb_matrix
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix
  public :: get_Whyb_matrix



contains





  !##################################################################
  !
  !     DMFT BATH PROCEDURES:
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
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
  end subroutine allocate_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%d))deallocate(dmft_bath_%d)
    if(allocated(dmft_bath_%v))deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%u))deallocate(dmft_bath_%u)
    dmft_bath_%status=.false.
  end subroutine deallocate_dmft_bath




  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_dmft_bath(dmft_bath_,hwband_)
    type(effective_bath) :: dmft_bath_
    real(8)              :: hwband_
    integer              :: i,iorb,ispin,unit,flen,Nh
    logical              :: IOfile
    real(8)              :: de
    if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
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
    enddo
    !Get SC amplitudes
    if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
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
  end subroutine init_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_dmft_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i,ispin,iorb
    if(ED_MPI_ID==0)then
       unit_=LOGfile;if(present(unit))unit_=unit
       if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
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
  end subroutine write_dmft_bath




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
       if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
       used_=.false.;if(present(used))used_=used
       extension=".restart";if(used_)extension=".used"
       file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
       if(present(file))file_=reg(file)
       unit_=free_unit()
       open(unit_,file=reg(file_))
       call write_dmft_bath(dmft_bath_,unit_)
       close(unit_)
    endif
  end subroutine save_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_dmft_bath(bath_,dmft_bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i,iorb,ispin,stride,io
    logical                :: check
    if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
    check = check_size_bath(bath_)
    if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
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
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
    end select
  end subroutine set_dmft_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine get_dmft_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:)   :: bath_
    integer                :: iorb,ispin,stride,io,i
    logical                :: check
    if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
    check=check_size_bath(bath_)
    if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
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
                   bath_(io) = dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
    end select
  end subroutine get_dmft_bath






  !##################################################################
  !
  !     USER BATH PROCEDURES:
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !
  ! Get size of each dimension of the component array. 
  ! The Result is an rank 1 integer array Ndim with dimension:
  ! 3 for get_component_size_bath
  ! 2 for get_spin_component_size_bath & get_orb_component_size_bath
  ! 1 for get_spin_orb_component_size_bath
  !+-------------------------------------------------------------------+
  function get_size_bath(ispin) result(bath_size)
    integer :: bath_size
    integer,optional :: ispin
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Norb*Nbath + Norb*Nbath
       case ("superc")
          !( e [Nspin][Norb][Nbath] + d [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
       case ("nonsu2")
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nbath + Norb*Nbath
       case ("superc")
          !(e [Nspin][1][Nbath] + d [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
          bath_size = Nbath + Nbath + Norb*Nbath
       case ("nonsu2")
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
          bath_size = Nbath + Norb*Nbath + Norb*Nbath
       end select
       !
    end select
    !
    if(.not.present(ispin))bath_size=Nspin*bath_size
    !
  end function get_size_bath

  function get_component_size_bath(itype) result(Ndim)
    integer                  :: itype
    integer                  :: Ndim(3)
    call check_type_bath(itype)
    select case(bath_type)
    case default
       Ndim=[Nspin,Norb,Nbath]
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       case ("superc")
          if(itype==1.OR.itype==2)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       case ("nonsu2")
          if(itype==1)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       end select
    end select
  end function get_component_size_bath
  !
  function get_spin_component_size_bath(itype) result(Ndim) 
    integer                  :: itype
    integer                  :: Ndim(2)
    call check_type_bath(itype)
    select case(bath_type)
    case default
       Ndim=[Norb,Nbath]
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       case ("superc")
          if(itype==1.OR.itype==2)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       case ("nonsu2")
          if(itype==1)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       end select
    end select
  end function get_spin_component_size_bath
  !
  function get_orb_component_size_bath(itype) result(Ndim)
    integer                  :: itype
    integer                  :: Ndim(2)
    Ndim=[Nspin,Nbath]
  end function get_orb_component_size_bath
  !
  function get_spin_orb_component_size_bath(itype) result(Ndim)
    integer                  :: itype
    integer                  :: Ndim
    Ndim=Nbath
  end function get_spin_orb_component_size_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_size_bath(bath_) result(bool)
    real(8),dimension(:) :: bath_
    integer              :: Ntrue
    logical              :: bool
    Ntrue = get_size_bath()
    bool  = ( size(bath_) == Ntrue )
  end function check_size_bath


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  check if the specified itype is consistent with the input parameters.
  !+-----------------------------------------------------------------------------+!
  subroutine check_type_bath(itype)
    integer :: itype
    select case(ed_mode)
    case default
       if(itype<1.OR.itype>2)stop "check_type_bath error: ed_mode=normal, itype!=1,2"
    case ("superc")
       if(itype<1.OR.itype>3)stop "check_type_bath error: ed_mode=superc, itype!=1,2,3"
    case ("nonsu2")
       if(itype<1.OR.itype>3)stop "check_type_bath error: ed_mode=nonsu2, itype!=1,2,3"
    end select
    return
  end subroutine check_type_bath



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: check that the input array hsa the correct dimensions specified 
  ! for the choice of itype and possiblty ispin and/or iorb.
  !+-----------------------------------------------------------------------------+!
  subroutine assert_component_size_bath(array,itype,string1,string2)
    real(8),dimension(:,:,:) :: array
    integer                  :: itype
    character(len=*)         :: string1,string2
    integer                  :: Ndim(3)
    call check_type_bath(itype)
    select case(bath_type)
    case default
       Ndim=[Nspin,Norb,Nbath]
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       case ("superc")
          if(itype==1.OR.itype==2)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       case ("nonsu2")
          if(itype==1)then
             Ndim=[Nspin,1,Nbath]
          else
             Ndim=[Nspin,Norb,Nbath]
          endif
       end select
    end select
    call assert_shape(Array,Ndim,reg(string1),reg(string2))
  end subroutine assert_component_size_bath
  !
  subroutine assert_spin_component_size_bath(array,itype,string1,string2)
    real(8),dimension(:,:) :: array
    integer                  :: itype
    character(len=*)         :: string1,string2
    integer                  :: Ndim(2)
    call check_type_bath(itype)
    select case(bath_type)
    case default
       Ndim=[Norb,Nbath]
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       case ("superc")
          if(itype==1.OR.itype==2)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       case ("nonsu2")
          if(itype==1)then
             Ndim=[1,Nbath]
          else
             Ndim=[Norb,Nbath]
          endif
       end select
    end select
    call assert_shape(Array,Ndim,reg(string1),reg(string2))
  end subroutine assert_spin_component_size_bath
  !
  subroutine assert_orb_component_size_bath(array,itype,string1,string2)
    real(8),dimension(:,:) :: array
    integer                  :: itype
    character(len=*)         :: string1,string2
    integer                  :: Ndim(2)
    Ndim=[Nspin,Nbath]
    call assert_shape(Array,Ndim,reg(string1),reg(string2))
  end subroutine assert_orb_component_size_bath
  !
  subroutine assert_spin_orb_component_size_bath(array,itype,string1,string2)
    real(8),dimension(:) :: array
    integer                  :: itype
    character(len=*)         :: string1,string2
    integer                  :: Ndim
    Ndim=Nbath
    if(size(array)/=Nbath)stop "assert_spin_orb_component_size_bath error: size(array)!=Ndim"
  end subroutine assert_spin_orb_component_size_bath




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Get a specified itype,ispin,iorb component of the user bath.
  ! The component is returned into an Array of rank D
  ! get_component_bath         : return the entire itype component (D=3)
  ! get_spin_component_bath    : return the itype component for the select ispin (D=2)
  ! get_orb_component_bath     : return the itype component for the select iorb (D=2)
  ! get_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
  !+-----------------------------------------------------------------------------+!
  subroutine get_component_bath(array,bath_,itype)
    real(8),dimension(:,:,:) :: array
    real(8),dimension(:)     :: bath_
    integer                  :: itype
    logical                  :: check
    type(effective_bath)     :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "get_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_component_size_bath(array,itype,"get_component_bath","Array")
    call check_type_bath(itype)
    select case(ed_mode)
    case default
       if(itype==1)then
          Array = dmft_bath_%e(:,:,:)
       else
          Array = dmft_bath_%v(:,:,:)
       endif
    case ("superc")
       if(itype==1)then
          Array = dmft_bath_%e(:,:,:)           
       elseif(itype==2)then
          Array = dmft_bath_%d(:,:,:)
       else
          Array = dmft_bath_%v(:,:,:)
       endif
    case ("nonsu2")
       if(itype==1)then
          Array = dmft_bath_%e(:,:,:)           
       elseif(itype==2)then
          Array = dmft_bath_%v(:,:,:)
       else
          Array = dmft_bath_%u(:,:,:)
       endif
    end select
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine get_component_bath

  subroutine get_spin_component_bath(array,bath_,itype,ispin)
    real(8),dimension(:,:) :: array
    real(8),dimension(:)   :: bath_
    integer                :: itype,ispin
    logical                :: check
    type(effective_bath)   :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "get_spin_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_spin_component_size_bath(array,itype,"get_spin_component_bath","Array")
    call check_type_bath(itype)
    if(ispin>Nspin)stop "get_spin_component_bath error: ispin > Nspin"
    select case(ed_mode)
    case default
       if(itype==1)then
          Array = dmft_bath_%e(ispin,:,:)
       else
          Array = dmft_bath_%v(ispin,:,:)
       endif
    case ("superc")
       if(itype==1)then
          Array = dmft_bath_%e(ispin,:,:)           
       elseif(itype==2)then
          Array = dmft_bath_%d(ispin,:,:)
       else
          Array = dmft_bath_%v(ispin,:,:)
       endif
    case ("nonsu2")
       if(itype==1)then
          Array = dmft_bath_%e(ispin,:,:)           
       elseif(itype==2)then
          Array = dmft_bath_%v(ispin,:,:)
       else
          Array = dmft_bath_%u(ispin,:,:)
       endif
    end select
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine get_spin_component_bath

  subroutine get_orb_component_bath(array,bath_,itype,iorb)
    real(8),dimension(:,:) :: array
    real(8),dimension(:)   :: bath_
    integer                :: itype,iorb
    logical                :: check
    type(effective_bath)   :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "get_orb_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_orb_component_size_bath(array,itype,"get_orb_component_bath","Array")
    call check_type_bath(itype)
    if(iorb>Norb)stop "get_orb_component_bath error: iorb > Norb"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(itype==1)then
             Array = dmft_bath_%e(:,iorb,:)
          else
             Array = dmft_bath_%v(:,iorb,:)
          endif
       case ("superc")
          if(itype==1)then
             Array = dmft_bath_%e(:,iorb,:)           
          elseif(itype==2)then
             Array = dmft_bath_%d(:,iorb,:)
          else
             Array = dmft_bath_%v(:,iorb,:)
          endif
       case ("nonsu2")
          if(itype==1)then
             Array = dmft_bath_%e(:,iorb,:)           
          elseif(itype==2)then
             Array = dmft_bath_%v(:,iorb,:)
          else
             Array = dmft_bath_%u(:,iorb,:)
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Array = dmft_bath_%e(:,1,:)
          else
             Array = dmft_bath_%v(:,iorb,:)
          endif
       case ("superc")
          if(itype==1)then
             Array = dmft_bath_%e(:,1,:)           
          elseif(itype==2)then
             Array = dmft_bath_%d(:,1,:)
          else
             Array = dmft_bath_%v(:,iorb,:)
          endif
       case ("nonsu2")
          if(itype==1)then
             Array = dmft_bath_%e(:,1,:)           
          elseif(itype==2)then
             Array = dmft_bath_%v(:,iorb,:)
          else
             Array = dmft_bath_%u(:,iorb,:)
          endif
       end select
    end select
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine get_orb_component_bath

  subroutine get_spin_orb_component_bath(array,bath_,itype,ispin,iorb)
    real(8),dimension(:) :: array
    real(8),dimension(:) :: bath_
    integer              :: itype,ispin,iorb
    logical              :: check
    type(effective_bath) :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "get_spin_orb_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_spin_orb_component_size_bath(array,itype,"get_spin_orb_component_bath","Array")
    call check_type_bath(itype)
    if(ispin>Nspin)stop "get_spin_orb_component_bath error: ispin > Nspin"
    if(iorb>Norb)stop "get_spin_orb_component_bath error: iorb > Norb"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(itype==1)then
             Array = dmft_bath_%e(ispin,iorb,:)
          else
             Array = dmft_bath_%v(ispin,iorb,:)
          endif
       case ("superc")
          if(itype==1)then
             Array = dmft_bath_%e(ispin,iorb,:)           
          elseif(itype==2)then
             Array = dmft_bath_%d(ispin,iorb,:)
          else
             Array = dmft_bath_%v(ispin,iorb,:)
          endif
       case ("nonsu2")
          if(itype==1)then
             Array = dmft_bath_%e(ispin,iorb,:)           
          elseif(itype==2)then
             Array = dmft_bath_%v(ispin,iorb,:)
          else
             Array = dmft_bath_%u(ispin,iorb,:)
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             Array = dmft_bath_%e(ispin,1,:)
          else
             Array = dmft_bath_%v(ispin,iorb,:)
          endif
       case ("superc")
          if(itype==1)then
             Array = dmft_bath_%e(ispin,1,:)           
          elseif(itype==2)then
             Array = dmft_bath_%d(ispin,1,:)
          else
             Array = dmft_bath_%v(ispin,iorb,:)
          endif
       case ("nonsu2")
          if(itype==1)then
             Array = dmft_bath_%e(ispin,1,:)           
          elseif(itype==2)then
             Array = dmft_bath_%v(ispin,iorb,:)
          else
             Array = dmft_bath_%u(ispin,iorb,:)
          endif
       end select
    end select
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine get_spin_orb_component_bath





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Set a specified itype,ispin,iorb component of the user bath.
  ! The component is set from an Array of rank D
  ! set_component_bath         : return the entire itype component (D=3)
  ! set_spin_component_bath    : return the itype component for the select ispin (D=2)
  ! set_orb_component_bath     : return the itype component for the select iorb (D=2)
  ! set_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
  !+-----------------------------------------------------------------------------+!
  subroutine set_component_bath(array,bath_,itype)
    real(8),dimension(:,:,:) :: array
    real(8),dimension(:)     :: bath_
    integer                  :: itype
    logical                  :: check
    type(effective_bath)     :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "set_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_component_size_bath(array,itype,"set_component_bath","Array")
    call check_type_bath(itype)
    select case(ed_mode)
    case default
       if(itype==1)then
          dmft_bath_%e(:,:,:)  = Array
       else
          dmft_bath_%v(:,:,:)  = Array
       endif
    case ("superc")
       if(itype==1)then
          dmft_bath_%e(:,:,:)  = Array           
       elseif(itype==2)then
          dmft_bath_%d(:,:,:)  = Array
       else
          dmft_bath_%v(:,:,:)  = Array
       endif
    case ("nonsu2")
       if(itype==1)then
          dmft_bath_%e(:,:,:)  = Array           
       elseif(itype==2)then
          dmft_bath_%v(:,:,:)  = Array
       else
          dmft_bath_%u(:,:,:)  = Array
       endif
    end select
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine set_component_bath

  subroutine set_spin_component_bath(array,bath_,itype,ispin)
    real(8),dimension(:,:) :: array
    real(8),dimension(:)   :: bath_
    integer                :: itype,ispin
    logical                :: check
    type(effective_bath)   :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "set_spin_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_spin_component_size_bath(array,itype,"set_spin_component_bath","Array")
    call check_type_bath(itype)
    if(ispin>Nspin)stop "set_spin_component_bath error: ispin > Nspin"
    select case(ed_mode)
    case default
       if(itype==1)then
          dmft_bath_%e(ispin,:,:)  = Array
       else
          dmft_bath_%v(ispin,:,:)  = Array
       endif
    case ("superc")
       if(itype==1)then
          dmft_bath_%e(ispin,:,:)  = Array           
       elseif(itype==2)then
          dmft_bath_%d(ispin,:,:)  = Array
       else
          dmft_bath_%v(ispin,:,:)  = Array
       endif
    case ("nonsu2")
       if(itype==1)then
          dmft_bath_%e(ispin,:,:)  = Array           
       elseif(itype==2)then
          dmft_bath_%v(ispin,:,:)  = Array
       else
          dmft_bath_%u(ispin,:,:)  = Array
       endif
    end select
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine set_spin_component_bath

  subroutine set_orb_component_bath(array,bath_,itype,iorb)
    real(8),dimension(:,:) :: array
    real(8),dimension(:)   :: bath_
    integer                :: itype,iorb
    logical                :: check
    type(effective_bath)   :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "set_orb_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_orb_component_size_bath(array,itype,"set_orb_component_bath","Array")
    call check_type_bath(itype)
    if(iorb>Norb)stop "set_orb_component_bath error: iorb > Norb"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(itype==1)then
             dmft_bath_%e(:,iorb,:)  = Array
          else
             dmft_bath_%v(:,iorb,:)  = Array
          endif
       case ("superc")
          if(itype==1)then
             dmft_bath_%e(:,iorb,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%d(:,iorb,:)  = Array
          else
             dmft_bath_%v(:,iorb,:)  = Array
          endif
       case ("nonsu2")
          if(itype==1)then
             dmft_bath_%e(:,iorb,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%v(:,iorb,:)  = Array
          else
             dmft_bath_%u(:,iorb,:)  = Array
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             dmft_bath_%e(:,1,:)  = Array
          else
             dmft_bath_%v(:,iorb,:)  = Array
          endif
       case ("superc")
          if(itype==1)then
             dmft_bath_%e(:,1,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%d(:,1,:)  = Array
          else
             dmft_bath_%v(:,iorb,:)  = Array
          endif
       case ("nonsu2")
          if(itype==1)then
             dmft_bath_%e(:,1,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%v(:,iorb,:)  = Array
          else
             dmft_bath_%u(:,iorb,:)  = Array
          endif
       end select
    end select
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine set_orb_component_bath

  subroutine set_spin_orb_component_bath(array,bath_,itype,ispin,iorb)
    real(8),dimension(:) :: array
    real(8),dimension(:) :: bath_
    integer              :: itype,ispin,iorb
    logical              :: check
    type(effective_bath) :: dmft_bath_
    !
    check= check_size_bath(bath_)
    if(.not.check)stop "set_spin_orb_component_bath error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call assert_spin_orb_component_size_bath(array,itype,"set_spin_orb_component_bath","Array")
    call check_type_bath(itype)
    if(ispin>Nspin)stop "set_spin_orb_component_bath error: ispin > Nspin"
    if(iorb>Norb)stop "set_spin_orb_component_bath error: iorb > Norb"
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(itype==1)then
             dmft_bath_%e(ispin,iorb,:)  = Array
          else
             dmft_bath_%v(ispin,iorb,:)  = Array
          endif
       case ("superc")
          if(itype==1)then
             dmft_bath_%e(ispin,iorb,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%d(ispin,iorb,:)  = Array
          else
             dmft_bath_%v(ispin,iorb,:)  = Array
          endif
       case ("nonsu2")
          if(itype==1)then
             dmft_bath_%e(ispin,iorb,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%v(ispin,iorb,:)  = Array
          else
             dmft_bath_%u(ispin,iorb,:)  = Array
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(itype==1)then
             dmft_bath_%e(ispin,1,:)  = Array
          else
             dmft_bath_%v(ispin,iorb,:)  = Array
          endif
       case ("superc")
          if(itype==1)then
             dmft_bath_%e(ispin,1,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%d(ispin,1,:)  = Array
          else
             dmft_bath_%v(ispin,iorb,:)  = Array
          endif
       case ("nonsu2")
          if(itype==1)then
             dmft_bath_%e(ispin,1,:)  = Array           
          elseif(itype==2)then
             dmft_bath_%v(ispin,iorb,:)  = Array
          else
             dmft_bath_%u(ispin,iorb,:)  = Array
          endif
       end select
    end select
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine set_spin_orb_component_bath



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Copy a specified component of IN bath to the OUT bath.
  ! copy_component_bath         : copy the entire itype component
  ! copy_spin_component_bath    : copy ispin to jspin component
  ! copy_orb_component_bath     : copy iorb to jorb component
  ! copy_spin_orb_component_bath: copy ispin.iorb to jspin.jorb components
  !+-----------------------------------------------------------------------------+!
  subroutine copy_component_bath(bathIN,bathOUT,itype)
    real(8),dimension(:)     :: bathIN,bathOUT
    integer                  :: itype
    logical                  :: check
    type(effective_bath)     :: dIN,dOUT
    !
    check= check_size_bath(bathIN)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
    check= check_size_bath(bathOUT)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
    call allocate_dmft_bath(dIN)
    call allocate_dmft_bath(dOUT)
    call set_dmft_bath(bathIN,dIN)
    call set_dmft_bath(bathOUT,dOUT)
    call check_type_bath(itype)
    select case(ed_mode)
    case default
       if(itype==1)then
          dOUT%e(:,:,:)  = dIN%e(:,:,:)
       else
          dOUT%v(:,:,:)  = dIN%v(:,:,:)
       endif
    case ("superc")
       if(itype==1)then
          dOUT%e(:,:,:)  = dIN%e(:,:,:)
       elseif(itype==2)then
          dOUT%d(:,:,:)  = dIN%d(:,:,:)
       else
          dOUT%v(:,:,:)  = dIN%v(:,:,:)
       endif
    case ("nonsu2")
       if(itype==1)then
          dOUT%e(:,:,:)  = dIN%e(:,:,:)
       elseif(itype==2)then
          dOUT%v(:,:,:)  = dIN%v(:,:,:)
       else
          dOUT%u(:,:,:)  = dIN%u(:,:,:)
       endif
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_component_bath

  subroutine copy_spin_component_bath(bathIN,ispin,bathOUT,jspin,itype)
    real(8),dimension(:)     :: bathIN,bathOUT
    integer                  :: ispin,jspin
    integer,optional         :: itype
    logical                  :: check
    type(effective_bath)     :: dIN,dOUT
    !
    check= check_size_bath(bathIN)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
    check= check_size_bath(bathOUT)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
    call allocate_dmft_bath(dIN)
    call allocate_dmft_bath(dOUT)
    call set_dmft_bath(bathIN,dIN)
    call set_dmft_bath(bathOUT,dOUT)
    if(present(itype))call check_type_bath(itype)
    if(ispin>Norb.OR.jspin>Nspin)stop "copy_spin_component_bath error: ispin/jspin > Nspin"
    select case(ed_mode)
    case default
       if(present(itype))then
          if(itype==1)then
             dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          else
             dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
          endif
       else
          dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
       endif
    case ("superc")
       if(present(itype))then
          if(itype==1)then
             dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          elseif(itype==2)then
             dOUT%d(ispin,:,:)  = dIN%d(jspin,:,:)
          else
             dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
          endif
       else
          dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          dOUT%d(ispin,:,:)  = dIN%d(jspin,:,:)
          dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
       endif
    case ("nonsu2")
       if(present(itype))then
          if(itype==1)then
             dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          elseif(itype==2)then
             dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
          else
             dOUT%u(ispin,:,:)  = dIN%u(jspin,:,:)
          endif
       else
          dOUT%e(ispin,:,:)  = dIN%e(jspin,:,:)
          dOUT%v(ispin,:,:)  = dIN%v(jspin,:,:)
          dOUT%u(ispin,:,:)  = dIN%u(jspin,:,:)
       endif
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_spin_component_bath

  subroutine copy_orb_component_bath(bathIN,iorb,bathOUT,jorb,itype)
    real(8),dimension(:)     :: bathIN,bathOUT
    integer                  :: iorb,jorb
    integer,optional         :: itype
    logical                  :: check
    type(effective_bath)     :: dIN,dOUT
    !
    check= check_size_bath(bathIN)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
    check= check_size_bath(bathOUT)
    if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
    call allocate_dmft_bath(dIN)
    call allocate_dmft_bath(dOUT)
    call set_dmft_bath(bathIN,dIN)
    call set_dmft_bath(bathOUT,dOUT)
    if(present(itype))call check_type_bath(itype)
    if(iorb>Norb.OR.jorb>Norb)stop "copy_orb_component_bath error: iorb/jorb > Norb"
    !
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             else
                dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             endif
          else
             dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             elseif(itype==2)then
                dOUT%d(:,iorb,:)  = dIN%d(:,jorb,:)
             else
                dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             endif
          else
             dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             dOUT%d(:,iorb,:)  = dIN%d(:,jorb,:)
             dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             elseif(itype==2)then
                dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             else
                dOUT%u(:,iorb,:)  = dIN%u(:,jorb,:)
             endif
          else
             dOUT%e(:,iorb,:)  = dIN%e(:,jorb,:)
             dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             dOUT%u(:,iorb,:)  = dIN%u(:,jorb,:)
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,1,:)    = dIN%e(:,1,:)
             else
                dOUT%v(:,iorb,:) = dIN%v(:,jorb,:)
             endif
          else
             dOUT%e(:,1,:)     = dIN%e(:,1,:)
             dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,1,:)    = dIN%e(:,1,:)
             elseif(itype==2)then
                dOUT%d(:,1,:)    = dIN%d(:,1,:)
             else
                dOUT%v(:,iorb,:) = dIN%v(:,jorb,:)
             endif
          else
             dOUT%e(:,1,:)    = dIN%e(:,1,:)
             dOUT%d(:,1,:)    = dIN%d(:,1,:)
             dOUT%v(:,iorb,:) = dIN%v(:,jorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(:,1,:)     = dIN%e(:,1,:)
             elseif(itype==2)then
                dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             else
                dOUT%u(:,iorb,:)  = dIN%u(:,jorb,:)
             endif
          else
             dOUT%e(:,1,:)     = dIN%e(:,1,:)
             dOUT%v(:,iorb,:)  = dIN%v(:,jorb,:)
             dOUT%u(:,iorb,:)  = dIN%u(:,jorb,:)
          endif
       end select
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_orb_component_bath

  subroutine copy_spin_orb_component_bath(bathIN,ispin,iorb,bathOUT,jspin,jorb,itype)
    real(8),dimension(:)     :: bathIN,bathOUT
    integer                  :: ispin,jspin
    integer                  :: iorb,jorb
    integer,optional         :: itype
    logical                  :: check
    type(effective_bath)     :: dIN,dOUT
    !
    check= check_size_bath(bathIN)
    if(.not.check)stop "copy_spin_orb_component_bath error: wrong bath dimensions IN"
    check= check_size_bath(bathOUT)
    if(.not.check)stop "copy_spin_orb_component_bath error: wrong bath dimensions OUT"
    call allocate_dmft_bath(dIN)
    call allocate_dmft_bath(dOUT)
    call set_dmft_bath(bathIN,dIN)
    call set_dmft_bath(bathOUT,dOUT)
    if(present(itype))call check_type_bath(itype)
    if(ispin>Norb.OR.jspin>Nspin)stop "copy_spin_orb_component_bath error: ispin/jspin > Nspin"
    if(iorb>Norb.OR.jorb>Norb)stop "copy_spin_orb_component_bath error: iorb/jorb > Norb"
    !
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             else
                dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             elseif(itype==2)then
                dOUT%d(ispin,iorb,:)  = dIN%d(jspin,jorb,:)
             else
                dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             dOUT%d(ispin,iorb,:)  = dIN%d(jspin,jorb,:)
             dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             elseif(itype==2)then
                dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             else
                dOUT%u(ispin,iorb,:)  = dIN%u(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,iorb,:)  = dIN%e(jspin,jorb,:)
             dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             dOUT%u(ispin,iorb,:)  = dIN%u(jspin,jorb,:)
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,1,:)    = dIN%e(jspin,1,:)
             else
                dOUT%v(ispin,iorb,:) = dIN%v(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,1,:)     = dIN%e(jspin,1,:)
             dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,1,:)    = dIN%e(jspin,1,:)
             elseif(itype==2)then
                dOUT%d(ispin,1,:)    = dIN%d(jspin,1,:)
             else
                dOUT%v(ispin,iorb,:) = dIN%v(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,1,:)    = dIN%e(jspin,1,:)
             dOUT%d(ispin,1,:)    = dIN%d(jspin,1,:)
             dOUT%v(ispin,iorb,:) = dIN%v(jspin,jorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(ispin,1,:)     = dIN%e(jspin,1,:)
             elseif(itype==2)then
                dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             else
                dOUT%u(ispin,iorb,:)  = dIN%u(jspin,jorb,:)
             endif
          else
             dOUT%e(ispin,1,:)     = dIN%e(jspin,1,:)
             dOUT%v(ispin,iorb,:)  = dIN%v(jspin,jorb,:)
             dOUT%u(ispin,iorb,:)  = dIN%u(jspin,jorb,:)
          endif
       end select
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_spin_orb_component_bath




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: save the user defined bath to a file
  !+-----------------------------------------------------------------------------+!
  subroutine save_bath(bath_,file,used)
    real(8),dimension(:)      :: bath_
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_,check
    character(len=16)         :: extension
    integer                   :: unit_
    if(ED_MPI_ID==0)then
       check= check_size_bath(bath_)
       if(.not.check)stop "save_bath error: wrong bath dimensions"
       call allocate_dmft_bath(dmft_bath_)
       call set_dmft_bath(bath_,dmft_bath_)
       used_=.false.;if(present(used))used_=used
       extension=".restart";if(used_)extension=".used"
       file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
       if(present(file))file_=reg(file)
       unit_=free_unit()
       open(unit_,file=reg(file_))
       call write_dmft_bath(dmft_bath_,unit_)
       close(unit_)
       call deallocate_dmft_bath(dmft_bath_)
    endif
  end subroutine save_bath















  !##################################################################
  !
  !     USER BATH PREDEFINED SYMMETRIES:
  !
  !##################################################################

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
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine break_symmetry_bath

  subroutine spin_symmetrize_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    if(Nspin==1)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
    dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
    if(ed_mode=="superc")dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath

  subroutine ph_symmetrize_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
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
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath

  subroutine ph_trans_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    type(effective_bath)   :: tmp_dmft_bath
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call allocate_dmft_bath(tmp_dmft_bath)
    call set_dmft_bath(bath_,dmft_bath_)
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
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_trans_bath

  subroutine enforce_normal_bath(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(ed_mode=="superc")dmft_bath_%d(:,:,:)=0.d0
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine enforce_normal_bath







  !##################################################################
  !
  !     W_hyb PROCEDURES
  !
  !##################################################################
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



END MODULE ED_BATH
