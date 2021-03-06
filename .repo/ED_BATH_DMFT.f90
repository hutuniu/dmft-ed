MODULE ED_BATH_DMFT
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,matrix_diagonalize
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
       allocate(dmft_bath_%h(Nspin,Nspin,Norb,Norb,Nbath))     !replica hamilt of the bath
       allocate(dmft_bath_%vr(Nbath))                          !hybridization 
       allocate(dmft_bath_%mask(Nspin,Nspin,Norb,Norb,2))      !mask on components 
       allocate(dmft_bath_%LS(Nspin,Nspin,Norb,Norb,2))        !rotations
       !
    end select
    dmft_bath_%status=.true.
  end subroutine allocate_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%d))   deallocate(dmft_bath_%d)
    if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%u))   deallocate(dmft_bath_%u)
    if(allocated(dmft_bath_%vr))  deallocate(dmft_bath_%vr)
    if(allocated(dmft_bath_%h))   deallocate(dmft_bath_%h)
    if(allocated(dmft_bath_%mask))deallocate(dmft_bath_%mask)
    if(allocated(dmft_bath_%LS))  deallocate(dmft_bath_%LS)
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
    real(8)              :: hybr_aux_R,hybr_aux_I
    real(8)              :: himp_aux_R(Nspin*Norb,Nspin*Norb)
    real(8)              :: himp_aux_I(Nspin*Norb,Nspin*Norb)
    real(8)              :: re,im
    integer              :: i,unit,flen,Nh
    integer              :: io,jo,iorb,ispin,jorb,jspin
    logical              :: IOfile
    real(8)              :: de,noise_tot
    real(8),allocatable  :: noise_b(:),noise_s(:),noise_o(:)
    real(8)              :: phi_i,phi_j,A_ij,B_ij
    real(8)              :: Re1,Im1,Re2,Im2
    character(len=21)    :: space
    if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
    allocate(noise_b(Nbath));noise_b=0.d0 
    allocate(noise_s(Nspin));noise_s=0.d0 
    allocate(noise_o(Norb)); noise_o=0.d0 
    call random_number(noise_b)
    call random_number(noise_s)
    call random_number(noise_o)
    noise_b=noise_b*ed_bath_noise_thr
    noise_s=noise_s*ed_bath_noise_thr
    noise_o=noise_o*ed_bath_noise_thr

    select case(bath_type)
    case default
       !Get energies:
       dmft_bath_%e(:,:,1)    =-hwband_ + noise_b(1)
       dmft_bath_%e(:,:,Nbath)= hwband_ + noise_b(Nbath)
       Nh=Nbath/2
       if(mod(Nbath,2)==0.and.Nbath>=4)then
          de=hwband_/max(Nh-1,1)
          dmft_bath_%e(:,:,Nh)  = -1.d-3 + noise_b(Nh)
          dmft_bath_%e(:,:,Nh+1)=  1.d-3 + noise_b(Nh+1)
          do i=2,Nh-1
             dmft_bath_%e(:,:,i)   =-hwband_ + (i-1)*de  + noise_b(i)
             dmft_bath_%e(:,:,Nbath-i+1)= hwband_ - (i-1)*de + noise_b(Nbath-i+1)
          enddo
       elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
          de=hwband_/Nh
          dmft_bath_%e(:,:,Nh+1)= 0.0d0 + noise_b(Nh+1)
          do i=2,Nh
             dmft_bath_%e(:,:,i)        =-hwband_ + (i-1)*de + noise_b(i)
             dmft_bath_%e(:,:,Nbath-i+1)= hwband_ - (i-1)*de + noise_b(Nbath-i+1)
          enddo
       endif
       !Get spin-keep yhbridizations
       do i=1,Nbath
          dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))+noise_b(i)
       enddo
       !Get SC amplitudes
       if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
       !Get spin-flip hybridizations
       if(ed_mode=="nonsu2")then
          do i=1,Nbath
             dmft_bath_%u(:,:,i) = dmft_bath_%v(:,:,i)*ed_vsf_ratio+noise_b(i)
          enddo
       endif
       !
    case('replica')
       !BATH INITIALIZATION
       dmft_bath_%h=zero
       do i=1,Nbath
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      re=0.0d0;im=0.0d0
                      noise_tot=noise_b(i)+noise_s(ispin)+noise_s(jspin)+noise_o(iorb)+noise_o(jorb)
                      if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1))re=1.0d0+noise_tot
                      if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2))im=0.1d0+noise_tot
                      io = iorb + (ispin-1)*Norb                           
                      jo = jorb + (jspin-1)*Norb
                      if(io==jo)then
                         dmft_bath_%h(ispin,jspin,iorb,jorb,i)=cmplx(re,im)
                      else
                         dmft_bath_%h(ispin,jspin,iorb,jorb,i)=cmplx(re,im)
                      endif
                   enddo
                enddo
             enddo
          enddo
          if(ed_para)then
             dmft_bath_%h(:,:,:,:,i)=zero
             dmft_bath_%h(:,:,:,:,i)=-0.6d0+noise_b(i)
          endif
       enddo
       !HYBR. INITIALIZATION
       dmft_bath_%vr=zero
       do i=1,Nbath
          noise_tot=noise_b(i)
          dmft_bath_%vr(i)=cmplx(0.2d0+noise_tot,0.0d0)!*(-1)**(i-1)
       enddo
       !
       deallocate(noise_b,noise_s,noise_o)
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
       !
       select case(bath_type)
       case default
          !
          read(unit,*)
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
          read(unit,*)
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
          case ("normal","nonsu2")
             do i=1,Nbath
                himp_aux_R=0.0d0;himp_aux_I=0.0d0
                hybr_aux_R=0.0d0;hybr_aux_I=0.0d0
                himp_aux=zero
                do io=1,Nspin*Norb
                   if(io==1)read(unit,"(90(F21.12,1X))")     hybr_aux_R,hybr_aux_I,(himp_aux_R(io,jo),jo=1,Nspin*Norb),(himp_aux_I(io,jo),jo=1,Nspin*Norb)
                   if(io/=1)read(unit,"(2a21,90(F21.12,1X))")   space  ,   space  ,(himp_aux_R(io,jo),jo=1,Nspin*Norb),(himp_aux_I(io,jo),jo=1,Nspin*Norb)
                enddo
                read(unit,*)
                himp_aux=cmplx(himp_aux_R,himp_aux_I)
                dmft_bath_%h(:,:,:,:,i)=so2nn_reshape(himp_aux,Nspin,Norb)
                dmft_bath_%vr(i)=cmplx(hybr_aux_R,hybr_aux_I)
             enddo
          !
          case ("superc")
          !
          end select
          !
       end select
       close(unit)
    endif
  end subroutine init_dmft_bath

  !+-------------------------------------------------------------------+
  !PURPOSE  : set the mask based on impHloc in the replica bath topology
  !+-------------------------------------------------------------------+
  subroutine init_dmft_bath_mask(dmft_bath_)
    type(effective_bath),intent(inout) :: dmft_bath_
    integer                            :: iorb,ispin,jorb,jspin
    integer                            :: io,jo
    complex(8)                         :: LS(Nspin*Norb,Nspin*Norb)
    complex(8)                         :: LS_rot(Nspin*Norb,Nspin*Norb)
    !
    if(.not.(allocated(impHloc))) then
       stop "impHloc not allocated on mask initialization"
    endif
    dmft_bath_%mask=.false.
    dmft_bath_%LS=zero
    !
    ! MASK INITIALIZATION
    !
    do ispin=1,Nspin
       do iorb=1,Norb
         !Re-diagonal elements always present
         dmft_bath_%mask(ispin,ispin,iorb,iorb,1)=.true.
         !Im-diagonal elements checked
         if(abs(aimag(impHloc(ispin,ispin,iorb,iorb))).gt.1e-6)stop"impHloc is not Hermitian"
         !off-diagonal elements
         do jspin=1,Nspin
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                if(io/=jo)then
                   !Re
                   if( abs(real(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                      dmft_bath_%mask(ispin,jspin,iorb,jorb,1)=.true.
                   endif
                   !Im
                   if(abs(aimag(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                      dmft_bath_%mask(ispin,jspin,iorb,jorb,2)=.true.
                      if(ed_type=="d") stop "complex impHloc and ed_mode='d' are not compatible"
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
    !
    ! LS
    !
    LS=zero
    LS(1:2,3:4)= -Xi * pauli_z ! / 2.
    LS(1:2,5:6)= +Xi * pauli_y ! / 2.
    LS(3:4,5:6)= -Xi * pauli_x ! / 2.
    do io=1,Nspin*Norb
       do jo=io+1,Nspin*Norb
          LS(jo,io)=conjg(LS(io,jo))
       enddo
    enddo
    LS=so2os_reshape(LS,Nspin,Norb)
    !
    ! LS eigenvectors
    !
    LS_rot=zero
    !J=1/2 jz=-1/2
    LS_rot(1,1)=-Xi
    LS_rot(3,1)=-1.0d0
    LS_rot(6,1)=+Xi
    LS_rot(:,1)=LS_rot(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    LS_rot(2,2)=-Xi
    LS_rot(4,2)=+1.0d0
    LS_rot(5,2)=-Xi
    LS_rot(:,2)=LS_rot(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    LS_rot(2,3)=-Xi
    LS_rot(4,3)=+1.0d0
    LS_rot(5,3)=+2.0d0*Xi
    LS_rot(:,3)=LS_rot(:,3)/sqrt(6.)
    !J=3/2 jz=-1/2
    LS_rot(1,4)=+Xi
    LS_rot(3,4)=-1.0d0
    LS_rot(:,4)=LS_rot(:,4)/sqrt(2.)
    !J=3/2 jz=+1/2
    LS_rot(2,5)=-Xi 
    LS_rot(4,5)=-1.0d0
    LS_rot(:,5)=LS_rot(:,5)/sqrt(2.)
    !J=3/2 jz=+3/2
    LS_rot(1,6)=+Xi
    LS_rot(3,6)=+1.0d0
    LS_rot(6,6)=+2.0d0*Xi
    LS_rot(:,6)=LS_rot(:,6)/sqrt(6.)
    LS_rot=so2os_reshape(LS_rot,Nspin,Norb)
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                dmft_bath_%LS(ispin,jspin,iorb,jorb,1)=LS(io,jo)
                dmft_bath_%LS(ispin,jspin,iorb,jorb,2)=LS_rot(io,jo)
             enddo
          enddo
       enddo
    enddo
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
    complex(8)           :: hybr_aux
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
          case ("normal","nonsu2")
             do i=1,Nbath
                himp_aux=zero;himp_aux=nn2so_reshape(dmft_bath_%h(:,:,:,:,i),Nspin,Norb)
                hybr_aux=dmft_bath_%vr(i)
                do io=1,Nspin*Norb
                   if(unit_==LOGfile)then
                      if(ed_type=="d") then
                         if(io==1)write(unit_,"(F8.3,a5,90(F8.3,1X))")  real(hybr_aux),"|",(real(himp_aux(io,jo)),jo=1,Nspin*Norb)
                         if(io/=1)write(unit_,"(a8,a5,90(F8.3,1X))")        "  "      ,"|",(real(himp_aux(io,jo)),jo=1,Nspin*Norb)
                      endif
                      if(ed_type=="c") then
                         if(io==1) write(unit_,"(2F8.3,a5,90(F8.3,1X))") real(hybr_aux),aimag(hybr_aux),"|",( real(himp_aux(io,jo)),jo=1,Nspin*Norb),&
                                                                                                              (aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                         if(io/=1) write(unit_,"(2a8,a5,90(F8.3,1X))")        "  "     ,      "  "     ,"|",( real(himp_aux(io,jo)),jo=1,Nspin*Norb),&
                                                                                                              (aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                      endif
                   else
                      if(io==1)write(unit_,"(90(F21.12,1X))")      real(hybr_aux),aimag(hybr_aux),(real(himp_aux(io,jo)),jo=1,Nspin*Norb),(aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                      if(io/=1)write(unit_,"(2a21,90(F21.12,1X))")      "  "     ,     "  "      ,(real(himp_aux(io,jo)),jo=1,Nspin*Norb),(aimag(himp_aux(io,jo)),jo=1,Nspin*Norb)
                   endif
                enddo
                write(unit_,*)
             enddo
             !
          case ("superc")
             !
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
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check
    real(8)                :: element_R,element_I
    real(8)                :: phi_i,phi_j,A_ij,B_ij
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: LSmatrix
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
       case("normal","nonsu2")
          !
          dmft_bath_%h=zero
          dmft_bath_%vr=zero
          i = 0
          if(ed_para)then
             do ibath=1,Nbath
                i=i+1
                element_R=0.0d0;element_I=0.0d0
                !diagonal eps_k
                element_R=bath_(i)
                do ispin=1,Nspin
                   do iorb=1,Norb
                      dmft_bath_%h(ispin,ispin,iorb,iorb,ibath)=cmplx(element_R,element_I)
                   enddo
                enddo
                i=i+1
                element_R=0.0d0;element_I=0.0d0
                !off-diagonal lambda_k
                element_R=bath_(i)
                !I can perform this sum since LS is purely non-diagonal while eps_k only diagonal
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1)*Norb
                            jo = jorb + (jspin-1)*Norb
                            dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)+element_R*dmft_bath_%LS(ispin,jspin,iorb,jorb,1)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          else
             !all non-vanishing terms in imploc - all spin
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do ibath=1,Nbath
                            io = iorb + (ispin-1)*Norb
                            jo = jorb + (jspin-1)*Norb
                            if(io.gt.jo)cycle!only diagonal and upper triangular are saved for hermiticity
                            element_R=0.0d0;element_I=0.0d0
                            if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                               i=i+1
                               element_R=bath_(i)
                            endif
                            if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                               i=i+1
                               element_I=bath_(i)
                            endif
                            dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)=cmplx(element_R,element_I)
                            !hermiticity
                            if((ispin==jspin).and.(iorb/=jorb))dmft_bath_%h(ispin,ispin,jorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,ispin,iorb,jorb,ibath))
                            if((ispin/=jspin).and.(iorb==jorb))dmft_bath_%h(jspin,ispin,iorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,jspin,iorb,iorb,ibath))
                            if((ispin/=jspin).and.(iorb/=jorb))dmft_bath_%h(jspin,ispin,jorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          endif
          !
          !all Re[Hybr]
          do ibath=1,Nbath
             element_R=0.0d0;element_I=0.0d0
             i=i+1
             element_R=bath_(i)
             dmft_bath_%vr(ibath)=cmplx(element_R,element_I)
          enddo
          !
       case ("superc")

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
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
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
       case("normal","nonsu2")
          i = 0
          if(ed_para)then
             do ibath=1,Nbath
                !all diagonal per bath *all equal*
                i=i+1
                bath_(i)=real(dmft_bath_%h(1,1,1,1,ibath))
                !search for first non-vanishing off-diagonal per bath *all equal*
                loop1: do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1)*Norb
                            jo = jorb + (jspin-1)*Norb
                            if(io>=jo)cycle!only upper triangular is checked
                            if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) exit loop1
                         enddo
                      enddo
                   enddo
                enddo loop1 
                i=i+1
                bath_(i)=real(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
             enddo
          else
             !all non-vanishing terms in imploc - all spin
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do ibath=1,Nbath
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         if(io.gt.jo)cycle!only diagonal and upper triangular are saved for hermiticity
                            if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                               i=i+1
                               bath_(i)=real(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                            endif
                            if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                               i=i+1
                               bath_(i)=aimag(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          endif
          !
          !all Re[Hybr]
          do ibath=1,Nbath
             i=i+1
             bath_(i)=real(dmft_bath_%vr(ibath))
          enddo
          !
       case ("superc")
          stride = 0

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
    integer              :: ndx,ispin,iorb,jspin,jorb,io,jo,off_im_ndx
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
       ndx=0
       !Re/Im off-diagonal non-vanishing elements
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io.lt.jo)then
                      if( abs(real(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)ndx=ndx+1
                      if(abs(aimag(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                         ndx=ndx+1
                         if(ed_mode=="d")stop "complex impHloc and ed_mode='d' are not compatible"
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       !Real diagonal elements (always assumed)
       ndx= ndx + Nspin * Norb
       !complex diagonal elements checked
       do ispin=1,Nspin
          do iorb=1,Norb
             if(abs(aimag(impHloc(ispin,ispin,iorb,iorb))).gt.1e-6)stop"impHloc is not Hermitian"
          enddo
       enddo
       !number of non vanishing elements for each replica
       ndx = ndx * Nbath
       !real diagonal hybridizations
       ndx = ndx + Nbath
       !
       if(ed_para)then
          bath_size = ( 1+1+1 ) * Nbath
       else
          bath_size = ndx
       endif
       !
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
