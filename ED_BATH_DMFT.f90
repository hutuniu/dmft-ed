MODULE ED_BATH_DMFT
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private

  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !
  !DMFT BATH procedures:
  public :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public :: init_dmft_bath_mask              !INTERNAL (for effective_bath)
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)
  !

  interface get_Whyb_matrix
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix
  public :: get_Whyb_matrix



contains


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
    case('replica')
       !
       select case(ed_mode)
       case default                                               !normal  [N,Sz]
          allocate(dmft_bath_%h(Nspin,Nspin,Norb,Norb,Nbath))     !replica hamilt of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))                !hybridization 
          allocate(dmft_bath_%mask(Nspin,Nspin,Norb,Norb,2))      !mask on components 
       case ("superc")                                            !superc  [Sz]

       case ("nonsu2")                                            !nonsu2 case [N] qn
          allocate(dmft_bath_%h(Nspin,Nspin,Norb,Norb,Nbath))     !replica hamilt of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))                !hybridization 
          allocate(dmft_bath_%mask(Nspin,Nspin,Norb,Norb,2))      !mask on components 1=Re,2=Im
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
    if(allocated(dmft_bath_%h))deallocate(dmft_bath_%h)
    if(allocated(dmft_bath_%mask))deallocate(dmft_bath_%mask)
    dmft_bath_%status=.false.
  end subroutine deallocate_dmft_bath




  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_dmft_bath(dmft_bath_,hwband_)
    type(effective_bath) :: dmft_bath_
    real(8)              :: hwband_
    complex(8)           :: himp_aux(Nspin*Norb,Nspin*Norb)
    real(8)              :: hybr_aux(Nspin*Norb)
    real(8)              :: himp_aux_R(Nspin*Norb,Nspin*Norb)
    real(8)              :: himp_aux_I(Nspin*Norb,Nspin*Norb)
    integer              :: i,unit,flen,Nh
    integer              :: io,jo,iorb,ispin,jorb,jspin
    logical              :: IOfile
    real(8)              :: de
    real(8),allocatable  :: noise(:),noise2(:,:)
    if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
    select case(bath_type)
    !
    case default
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
       allocate(noise(Nbath));noise=0.d0 
       call random_number(noise)
       noise=noise*ed_bath_noise_thr
       !Get spin-keep yhbridizations
       do i=1,Nbath
          dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))+noise(i)
       enddo
       !Get SC amplitudes
       if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
       !Get spin-flip hybridizations
       if(ed_mode=="nonsu2")then
          noise=0.d0
          call random_number(noise)
          noise=noise*ed_bath_noise_thr
          do i=1,Nbath
             dmft_bath_%u(:,:,i) = dmft_bath_%v(:,:,i)*ed_vsf_ratio+noise(i)
          enddo
       endif
       deallocate(noise)
    !
    case('replica')
       !
       do i=1,Nbath
          dmft_bath_%h(:,:,:,:,i)=impHloc
          dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath*Norb)))
       enddo
   !    if (ed_mode=="normal")then
          Nh=Nbath/2
          if(mod(Nbath,2)==0)de=hwband_/max(Nh-1,1)
          if(mod(Nbath,2)/=0)de=hwband_/Nh
          do ispin=1,Nspin
             do iorb=1,Norb
                dmft_bath_%h(ispin,ispin,iorb,iorb,1)    =-hwband_
                dmft_bath_%h(ispin,ispin,iorb,iorb,Nbath)= hwband_
                if(mod(Nbath,2)==0)then
                   dmft_bath_%h(ispin,ispin,iorb,iorb,Nh)  = -1.d-4
                   dmft_bath_%h(ispin,ispin,iorb,iorb,Nh+1)=  1.d-4
                   do i=2,Nh-1
                      dmft_bath_%h(ispin,ispin,iorb,iorb,i)   =-hwband_ + (i-1)*de 
                      dmft_bath_%h(ispin,ispin,iorb,iorb,Nbath-i+1)= hwband_ - (i-1)*de
                   enddo
                else
                   dmft_bath_%h(ispin,ispin,iorb,iorb,Nh+1)= 0.d0
                   do i=2,Nh
                      dmft_bath_%h(ispin,ispin,iorb,iorb,i)        =-hwband_ + (i-1)*de
                      dmft_bath_%h(ispin,ispin,iorb,iorb,Nbath-i+1)= hwband_ - (i-1)*de
                   enddo
                endif
             enddo
          enddo
    !   else
    !      allocate(noise2(Nspin*Norb,Nspin*Norb));noise2=0.d0
    !      call random_number(noise2)
    !      noise2=noise2*ed_bath_noise_thr
    !      do i=1,Nbath
    !         dmft_bath_%h(:,:,:,:,i)=impHloc+so2nn_reshape(noise2,Nspin,Norb)
    !      enddo
    !      deallocate(noise2)
    !   endif

       call init_dmft_bath_mask(dmft_bath_)
       !
    end select
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
       case ('replica')
          !
          select case(ed_mode)
          case default
             do i=1,Nbath
                himp_aux_R=0.0d0
                himp_aux_I=0.0d0
                hybr_aux=0.0d0
                himp_aux=zero
                do io=1,Nspin*Norb
                   read(unit,*) hybr_aux(io),(himp_aux_R(io,jo),jo=1,Nspin*Norb),(himp_aux_I(io,jo),jo=1,Nspin*Norb)
                enddo
                read(unit,*)
                himp_aux=cmplx(himp_aux_R,himp_aux_I)
                dmft_bath_%h(:,:,:,:,i)=so2nn_reshape(himp_aux,Nspin,Norb)
                do ispin=1,Nspin
                   do iorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      dmft_bath_%v(ispin,iorb,i)=hybr_aux(io)
                   enddo
                enddo
             enddo
          case ("superc")

          case ("nonsu2")
             do i=1,Nbath
                himp_aux_R=0.0d0
                himp_aux_I=0.0d0
                hybr_aux=0.0d0
                himp_aux=zero
                do io=1,Nspin*Norb
                   read(unit,*) hybr_aux(io),(himp_aux_R(io,jo),jo=1,Nspin*Norb),(himp_aux_I(io,jo),jo=1,Nspin*Norb)
                enddo
                read(unit,*)
                himp_aux=cmplx(himp_aux_R,himp_aux_I)
                dmft_bath_%h(:,:,:,:,i)=so2nn_reshape(himp_aux,Nspin,Norb)
                do ispin=1,Nspin
                   do iorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      dmft_bath_%v(ispin,iorb,i)=hybr_aux(io)
                   enddo
                enddo
             enddo
          end select
          !
       end select
       close(unit)
    endif
  end subroutine init_dmft_bath


  subroutine init_dmft_bath_mask(dmft_bath_)
    type(effective_bath),intent(inout) :: dmft_bath_
    integer              :: io,jo,iorb,ispin,jorb,jspin
    !
    if(.not.(allocated(impHloc))) then
       stop "impHloc not allocated on mask initialization"
    endif
    dmft_bath_%mask=.false.
    select case (ed_mode)
    case default
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             !diagonal elements always present
             dmft_bath_%mask(ispin,ispin,iorb,iorb,1)=.true.
             !off-diagonal elements
             do jspin=ispin+1,Nspin
                do jorb=iorb+1,Norb
                   !Re
                   if( abs(real(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                      dmft_bath_%mask(ispin,jspin,iorb,jorb,1)=.true.
                      dmft_bath_%mask(jspin,ispin,jorb,iorb,1)=.true.
                   endif
                   !Im
                   if(abs(aimag(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                      dmft_bath_%mask(ispin,jspin,iorb,jorb,2)=.true.
                      dmft_bath_%mask(jspin,ispin,jorb,iorb,2)=.true.
                      if(ed_mode=="d") stop "complex impHloc and ed_mode='d' are not compatible"
                   endif
                enddo
             enddo
           enddo
        enddo
        !
     case("normal")
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             !diagonal elements always present
             dmft_bath_%mask(ispin,ispin,iorb,iorb,1)=.true.
             !off-diagonal elements
             do jorb=1,Norb
                !Re
                if( abs(real(impHloc(ispin,ispin,iorb,jorb))).gt.1e-6)then
                   dmft_bath_%mask(ispin,ispin,iorb,jorb,1)=.true.
                   dmft_bath_%mask(ispin,ispin,jorb,iorb,1)=.true.
                endif
                !Im
                if(abs(aimag(impHloc(ispin,ispin,iorb,jorb))).gt.1e-6)then
                   dmft_bath_%mask(ispin,ispin,iorb,jorb,2)=.true.
                   dmft_bath_%mask(ispin,ispin,jorb,iorb,2)=.true.
                   if(ed_mode=="d") stop "complex impHloc and ed_mode='d' are not compatible"
                endif
             enddo
           enddo
        enddo
        !
     end select
     !
  end subroutine init_dmft_bath_mask




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_dmft_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i,flen,Nh
    integer              :: io,jo,iorb,ispin,jorb,jspin
    real(8)              :: hybr_aux(Nspin*Norb)
    complex(8)           :: himp_aux(Nspin*Norb,Nspin*Norb)
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
       case ('replica')
          !
          select case(ed_mode)
          case default
             do i=1,Nbath
                himp_aux=zero;himp_aux=nn2so_reshape(dmft_bath_%h(:,:,:,:,i),Nspin,Norb)
                do ispin=1,Nspin
                   do iorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      hybr_aux(io)=dmft_bath_%v(ispin,iorb,i)
                   enddo
                enddo
                do io=1,Nspin*Norb
                   write(unit_,"(90(F21.12,1X))") hybr_aux(io),(real(himp_aux(io,jo)),jo=1,Nspin*Norb),(aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                enddo
                write(unit_,*)
             enddo
          case ("superc")

          case ("nonsu2")
             do i=1,Nbath
                himp_aux=zero;himp_aux=nn2so_reshape(dmft_bath_%h(:,:,:,:,i),Nspin,Norb)
                do ispin=1,Nspin
                   do iorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      hybr_aux(io)=dmft_bath_%v(ispin,iorb,i)
                   enddo
                enddo
                do io=1,Nspin*Norb
                   write(unit_,"(90(F21.12,1X))") hybr_aux(io),(real(himp_aux(io,jo)),jo=1,Nspin*Norb),(aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                enddo
                write(unit_,*)
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
       !<DEBUG
       write(*,*) reg(file_)
       !DEBUG>
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
    integer                :: stride,io,i
    integer                :: iorb,ispin,jorb,jspin
    logical                :: check
    real(8)                :: element_R,element_I
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
    case ('replica')
       !
       select case(ed_mode)
       case default
          !
          dmft_bath_%v=zero
          dmft_bath_%h=zero
          !
          io = 0
          do ispin=1,Nspin
             !all hybrd
             do iorb=1,Norb
                do i=1,Nbath
                   io=io+1
                   dmft_bath_%v(ispin,iorb,i)=bath_(io)
                enddo
             enddo
             !all non-vanishing terms in imploc - diagonal spin
             do i=1,Nbath
                do iorb=1,Norb
                   do jorb=1,Norb
                      element_R=0.0d0;element_I=0.0d0
                      if(dmft_bath_%mask(ispin,ispin,iorb,jorb,1)) then
                         io=io+1
                         element_R=bath_(io)
                      endif
                      if(dmft_bath_%mask(ispin,ispin,iorb,jorb,2)) then
                         io=io+1
                         element_I=bath_(io)
                      endif
                      dmft_bath_%h(ispin,ispin,iorb,jorb,i)=cmplx(element_R,element_I)
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")

          !
       case("nonsu2")
          !
          dmft_bath_%v=zero
          dmft_bath_%h=zero
          !
          io = 0
          do ispin=1,Nspin
             !all hybrd
             do iorb=1,Norb
                do i=1,Nbath
                   io=io+1
                   dmft_bath_%v(ispin,iorb,i)=bath_(io)
                enddo
             enddo
             !all non-vanishing terms in imploc - all spin
             do jspin=1,Nspin
                do i=1,Nbath
                   do iorb=1,Norb
                      do jorb=1,Norb
                         element_R=0.0d0;element_I=0.0d0
                         if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                            io=io+1
                            element_R=bath_(io)
                         endif
                         if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                            io=io+1
                            element_I=bath_(io)
                         endif
                         dmft_bath_%h(ispin,jspin,iorb,jorb,i)=cmplx(element_R,element_I)
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
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
    integer                :: stride,io,i
    integer                :: iorb,ispin,jorb,jspin
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
    case ('replica')
       !
       select case(ed_mode)
       case default
          !
          io = 0
          do ispin=1,Nspin
             !all hybrd
             do iorb=1,Norb
                do i=1,Nbath
                   io=io+1
                   bath_(io)=dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
             !all non-vanishing terms in imploc - diagonal spin
             do i=1,Nbath
                do iorb=1,Norb
                   do jorb=1,Norb
                      if(dmft_bath_%mask(ispin,ispin,iorb,jorb,1)) then
                         io=io+1
                         bath_(io)=real(dmft_bath_%h(ispin,ispin,iorb,jorb,i))
                       endif
                       if(dmft_bath_%mask(ispin,ispin,iorb,jorb,2)) then
                         io=io+1
                         bath_(io)=aimag(dmft_bath_%h(ispin,ispin,iorb,jorb,i))
                      endif
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0

          !
       case("nonsu2")
          io = 0
          do ispin=1,Nspin
             !all hybrd
             do iorb=1,Norb
                do i=1,Nbath
                   io=io+1
                   bath_(io)=dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
             !all non-vanishing terms in imploc - all spin
             do i=1,Nbath
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                            io=io+1
                            bath_(io)=real(dmft_bath_%h(ispin,jspin,iorb,jorb,i))
                         endif
                         if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                            io=io+1
                            bath_(io)=aimag(dmft_bath_%h(ispin,jspin,iorb,jorb,i))
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end subroutine get_dmft_bath



  !##################################################################
  !
  !     CHECK USER BATH SIZE (this is here for self-consistence)
  !     for bath_type='replica' this subroutine is used when impHloc 
  !     is already allocated. Hence no intercafe is needed as in the 
  !     case of the function in ED_BATH_USER
  !
  !##################################################################
  function check_size_bath(bath_) result(bool)
    real(8),dimension(:) :: bath_
    integer              :: bath_size
    integer              :: ndx,ispin,iorb,jspin,jorb
    logical              :: bool
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          bath_size = Norb*Nbath + Norb*Nbath
       case ("superc")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
       case ("nonsu2")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
       end select
       bath_size = Nspin*bath_size
    case('hybrid')
       select case(ed_mode)
       case default
          bath_size = Nbath + Norb*Nbath
       case ("superc")
          bath_size = Nbath + Nbath + Norb*Nbath
       case ("nonsu2")
          bath_size = Nbath + Norb*Nbath + Norb*Nbath
       end select
       bath_size = Nspin*bath_size
    case('replica')
       !
       if(.not.allocated(impHloc))stop "ERROR: bath_type='replica' but impHloc_nn not provided to get_size_bath"
       ndx=0
       !off-diagonal non-vanishing elements
       do ispin=1,Nspin
          do jspin=ispin+1,Nspin
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   if( abs(real(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)ndx=ndx+1
                   if(abs(aimag(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       ndx=ndx*2
       !Real diagonal elements (always assumed)
       ndx= ndx + Nspin * Norb
       !complex diagonal elements checked
       do ispin=1,Nspin
          do iorb=1,Norb
             if(abs(aimag(impHloc(ispin,ispin,iorb,iorb))).gt.1e-6)ndx=ndx+1
          enddo
       enddo
       select case(ed_mode)
       case default
          bath_size = ndx * Nbath + Nspin * Norb * Nbath
       case ("superc")
          !
       case ("nonsu2")
          bath_size = ndx * Nbath + Nspin * Norb * Nbath
       end select
    end select
    bool  = ( size(bath_) == bath_size )
  end function check_size_bath



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
    if(ed_para)then
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
    if(ed_para)then
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
    if(ed_para)then
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



END MODULE ED_BATH_DMFT
