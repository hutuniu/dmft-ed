MODULE ED_BATH_USER
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH_DMFT
  implicit none

  private

  !Ensure backward compatibility: [yeah I am lazy...]
  interface get_bath_size
     module procedure get_size_bath
  end interface get_bath_size

  interface get_size_bath
     module procedure get_size_bath!_normal_hybrid
     !module procedure get_size_bath_replica
  end interface get_size_bath
  !
  interface check_bath_dimension
     module procedure check_size_bath
  end interface check_bath_dimension

  !Interface for user bath I/O operations: get,set,copy
  interface get_component_bath
     module procedure get_full_component_bath
     module procedure get_spin_component_bath
     module procedure get_spin_orb_component_bath
  end interface get_component_bath

  interface set_component_bath
     module procedure set_full_component_bath
     module procedure set_spin_component_bath
     module procedure set_spin_orb_component_bath
  end interface set_component_bath

  interface copy_component_bath
     module procedure copy_full_component_bath
     module procedure copy_spin_component_bath
     module procedure copy_spin_orb_component_bath
  end interface copy_component_bath

  public :: get_size_bath

  public :: check_size_bath                  
  public :: get_bath_size
  public :: check_bath_dimension
  public :: get_component_size_bath
  public :: get_spin_component_size_bath
  public :: get_orb_component_size_bath
  public :: get_spin_orb_component_size_bath
  public :: get_component_bath
  public :: set_component_bath
  public :: copy_component_bath
  public :: save_bath
  !explicit symmetries:
  public :: break_symmetry_bath              
  public :: spin_symmetrize_bath             
  public :: SOC_symmetrize_bath             
  public :: ph_symmetrize_bath               
  public :: ph_trans_bath                    
  public :: enforce_normal_bath              


contains


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
  function get_size_bath(Hloc_nn,ispin_) result(bath_size)
    integer :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo
    integer,optional :: ispin_
    complex(8),allocatable,optional,intent(in) :: Hloc_nn(:,:,:,:)
    select case(bath_type)
    case default
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
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    case('hybrid')
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
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    case('replica')
       !
       if(.not.present(Hloc_nn))stop "ERROR: bath_type='replica' but impHloc_nn not provided to get_size_bath"
       ndx=0
       !off-diagonal non-vanishing elements
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io.lt.jo)then
                      if(real_Hrepl)then
                         if( abs((Hloc_nn(ispin,jspin,iorb,jorb))).gt.1e-12)ndx=ndx+1
                         if(abs(aimag(Hloc_nn(ispin,jspin,iorb,jorb))).gt.1e-6)then
                            if(ed_mode=="d")stop "complex impHloc and ed_mode='d' are not compatible"
                         endif
                      else
                         if( abs(real(Hloc_nn(ispin,jspin,iorb,jorb))).gt.1e-6)ndx=ndx+1
                         if(abs(aimag(Hloc_nn(ispin,jspin,iorb,jorb))).gt.1e-6)then
                            ndx=ndx+1
                            if(ed_mode=="d")stop "complex impHloc and ed_mode='d' are not compatible"
                         endif
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
             if(abs(aimag(Hloc_nn(ispin,ispin,iorb,iorb))).gt.1e-6)stop"impHloc is not Hermitian"
          enddo
       enddo
       select case(ed_mode)
       case default
          !
          !if(ed_type=="d")bath_size = ndx * Nbath + Nspin * Norb * Nbath
          if(ed_type=="d")bath_size = ndx * Nbath + Nbath
          if(ed_type=="c")then
             if(real_hybr)then
                bath_size = ndx * Nbath + Nbath
             else
                bath_size = ndx * Nbath + Nbath * 2
             endif
          endif
       case ("superc")
          !
       case ("nonsu2")
          if(ed_para)then
           !  if((ed_type=="d").or.(real_Hrepl.and.(ed_type=="c")))then!tutto reale
                bath_size = 2 * Nbath
                if(real_hybr)      bath_size = bath_size + Nbath 
                if(.not.real_hybr) bath_size = bath_size + Nbath * 2
           !  elseif((.not.real_Hrepl).and.(ed_type=="c"))then
           !     bath_size = 4 * Nbath
           !     if(real_hybr)      bath_size = bath_size + Nbath 
           !     if(.not.real_hybr) bath_size = bath_size + Nbath * 2
           !  endif
          else
             if(ed_type=="d")bath_size = ndx * Nbath + Nbath
             if(ed_type=="c")then
                if(real_hybr)then
                   bath_size = ndx * Nbath + Nbath
                else
                   bath_size = ndx * Nbath + Nbath * 2
                endif
             endif
          endif
       end select
       !
    end select
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









  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Get a specified itype,ispin,iorb component of the user bath.
  ! The component is returned into an Array of rank D
  ! get_full_component_bath    : return the entire itype component (D=3)
  ! get_spin_component_bath    : return the itype component for the select ispin (D=2)
  ! get_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
  !+-----------------------------------------------------------------------------+!
  subroutine get_full_component_bath(array,bath_,itype)
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
  end subroutine get_full_component_bath

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
  !This is super-seeded by the spin_orb procedure, called for each ispin value.
  !In this way it is possible to interface all the procedures to a single one.
  !   subroutine get_orb_component_bath(array,bath_,itype,iorb)
  !     real(8),dimension(:,:) :: array
  !     real(8),dimension(:)   :: bath_
  !     integer                :: itype,iorb
  !     logical                :: check
  !     type(effective_bath)   :: dmft_bath_
  !     !
  !     check= check_size_bath(bath_)
  !     if(.not.check)stop "get_orb_component_bath error: wrong bath dimensions"
  !     call allocate_dmft_bath(dmft_bath_)
  !     call set_dmft_bath(bath_,dmft_bath_)
  !     call assert_orb_component_size_bath(array,itype,"get_orb_component_bath","Array")
  !     call check_type_bath(itype)
  !     if(iorb>Norb)stop "get_orb_component_bath error: iorb > Norb"
  !     select case(bath_type)
  !     case default
  !        select case(ed_mode)
  !        case default
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,iorb,:)
  !           else
  !              Array = dmft_bath_%v(:,iorb,:)
  !           endif
  !        case ("superc")
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,iorb,:)           
  !           elseif(itype==2)then
  !              Array = dmft_bath_%d(:,iorb,:)
  !           else
  !              Array = dmft_bath_%v(:,iorb,:)
  !           endif
  !        case ("nonsu2")
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,iorb,:)           
  !           elseif(itype==2)then
  !              Array = dmft_bath_%v(:,iorb,:)
  !           else
  !              Array = dmft_bath_%u(:,iorb,:)
  !           endif
  !        end select
  !     case('hybrid')
  !        select case(ed_mode)
  !        case default
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,1,:)
  !           else
  !              Array = dmft_bath_%v(:,iorb,:)
  !           endif
  !        case ("superc")
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,1,:)           
  !           elseif(itype==2)then
  !              Array = dmft_bath_%d(:,1,:)
  !           else
  !              Array = dmft_bath_%v(:,iorb,:)
  !           endif
  !        case ("nonsu2")
  !           if(itype==1)then
  !              Array = dmft_bath_%e(:,1,:)           
  !           elseif(itype==2)then
  !              Array = dmft_bath_%v(:,iorb,:)
  !           else
  !              Array = dmft_bath_%u(:,iorb,:)
  !           endif
  !        end select
  !     end select
  !     call deallocate_dmft_bath(dmft_bath_)
  !   end subroutine get_orb_component_bath
  ! !
  !   subroutine get_so_component_bath(array,bath_,itype,spin,orb)
  !     real(8),dimension(:,:) :: array
  !     real(8),dimension(:)   :: bath_
  !     integer                :: itype,ispin,iorb
  !     integer,optional       :: spin,orb
  !     if(present(spin).AND..not.present(orb))then
  !        call get_spin_component_bath(array,bath_,itype,spin)
  !     elseif(.not.present(spin).AND.present(orb))then
  !        call get_orb_component_bath(array,bath_,itype,orb)
  !     else
  !        stop "get_so_component_bath error: either spin or orb input"
  !     endif
  !   end subroutine get_so_component_bath





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Set a specified itype,ispin,iorb component of the user bath.
  ! The component is set from an Array of rank D
  ! set_full_component_bath    : return the entire itype component (D=3)
  ! set_spin_component_bath    : return the itype component for the select ispin (D=2)
  ! set_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
  !+-----------------------------------------------------------------------------+!
  subroutine set_full_component_bath(array,bath_,itype)
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
  end subroutine set_full_component_bath

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
  !
  !This is super-seeded by the spin_orb procedure, called for each ispin value.
  !In this way it is possible to interface all the procedures to a single one.
  ! subroutine set_orb_component_bath(array,bath_,itype,iorb)
  !   real(8),dimension(:,:) :: array
  !   real(8),dimension(:)   :: bath_
  !   integer                :: itype,iorb
  !   logical                :: check
  !   type(effective_bath)   :: dmft_bath_
  !   !
  !   check= check_size_bath(bath_)
  !   if(.not.check)stop "set_orb_component_bath error: wrong bath dimensions"
  !   call allocate_dmft_bath(dmft_bath_)
  !   call set_dmft_bath(bath_,dmft_bath_)
  !   call assert_orb_component_size_bath(array,itype,"set_orb_component_bath","Array")
  !   call check_type_bath(itype)
  !   if(iorb>Norb)stop "set_orb_component_bath error: iorb > Norb"
  !   select case(bath_type)
  !   case default
  !      select case(ed_mode)
  !      case default
  !         if(itype==1)then
  !            dmft_bath_%e(:,iorb,:)  = Array
  !         else
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         endif
  !      case ("superc")
  !         if(itype==1)then
  !            dmft_bath_%e(:,iorb,:)  = Array           
  !         elseif(itype==2)then
  !            dmft_bath_%d(:,iorb,:)  = Array
  !         else
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         endif
  !      case ("nonsu2")
  !         if(itype==1)then
  !            dmft_bath_%e(:,iorb,:)  = Array           
  !         elseif(itype==2)then
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         else
  !            dmft_bath_%u(:,iorb,:)  = Array
  !         endif
  !      end select
  !   case('hybrid')
  !      select case(ed_mode)
  !      case default
  !         if(itype==1)then
  !            dmft_bath_%e(:,1,:)  = Array
  !         else
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         endif
  !      case ("superc")
  !         if(itype==1)then
  !            dmft_bath_%e(:,1,:)  = Array           
  !         elseif(itype==2)then
  !            dmft_bath_%d(:,1,:)  = Array
  !         else
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         endif
  !      case ("nonsu2")
  !         if(itype==1)then
  !            dmft_bath_%e(:,1,:)  = Array           
  !         elseif(itype==2)then
  !            dmft_bath_%v(:,iorb,:)  = Array
  !         else
  !            dmft_bath_%u(:,iorb,:)  = Array
  !         endif
  !      end select
  !   end select
  !   call get_dmft_bath(dmft_bath_,bath_)
  !   call deallocate_dmft_bath(dmft_bath_)
  ! end subroutine set_orb_component_bath
  ! subroutine set_so_component_bath(array,bath_,itype,spin,orb)
  !   real(8),dimension(:,:) :: array
  !   real(8),dimension(:)   :: bath_
  !   integer                :: itype,ispin,iorb
  !   integer,optional       :: spin,orb
  !   if(present(spin).AND..not.present(orb))then
  !      call set_spin_component_bath(array,bath_,itype,spin)
  !   elseif(.not.present(spin).AND.present(orb))then
  !      call set_orb_component_bath(array,bath_,itype,orb)
  !   else
  !      stop "set_so_component_bath error: either spin or orb input"
  !   endif
  ! end subroutine set_so_component_bath



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Copy a specified component of IN bath to the OUT bath.
  ! copy_full_component_bath    : copy the entire itype component
  ! copy_spin_component_bath    : copy ispin to jspin component
  ! copy_spin_orb_component_bath: copy ispin.iorb to jspin.jorb components
  !+-----------------------------------------------------------------------------+!
  subroutine copy_full_component_bath(bathIN,bathOUT,itype)
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
  end subroutine copy_full_component_bath

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
             dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          else
             dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
          endif
       else
          dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
       endif
    case ("superc")
       if(present(itype))then
          if(itype==1)then
             dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          elseif(itype==2)then
             dOUT%d(jspin,:,:)  = dIN%d(ispin,:,:)
          else
             dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
          endif
       else
          dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          dOUT%d(jspin,:,:)  = dIN%d(ispin,:,:)
          dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
       endif
    case ("nonsu2")
       if(present(itype))then
          if(itype==1)then
             dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          elseif(itype==2)then
             dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
          else
             dOUT%u(jspin,:,:)  = dIN%u(ispin,:,:)
          endif
       else
          dOUT%e(jspin,:,:)  = dIN%e(ispin,:,:)
          dOUT%v(jspin,:,:)  = dIN%v(ispin,:,:)
          dOUT%u(jspin,:,:)  = dIN%u(ispin,:,:)
       endif
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_spin_component_bath

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
                dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             else
                dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             elseif(itype==2)then
                dOUT%d(jspin,jorb,:)  = dIN%d(ispin,iorb,:)
             else
                dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             dOUT%d(jspin,jorb,:)  = dIN%d(ispin,iorb,:)
             dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             elseif(itype==2)then
                dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             else
                dOUT%u(jspin,jorb,:)  = dIN%u(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,jorb,:)  = dIN%e(ispin,iorb,:)
             dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             dOUT%u(jspin,jorb,:)  = dIN%u(ispin,iorb,:)
          endif
       end select
    case('hybrid')
       select case(ed_mode)
       case default
          if(present(itype))then
             if(itype==1)then
                dOUT%e(jspin,1,:)    = dIN%e(ispin,1,:)
             else
                dOUT%v(jspin,jorb,:) = dIN%v(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,1,:)     = dIN%e(ispin,1,:)
             dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
          endif
       case ("superc")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(jspin,1,:)    = dIN%e(ispin,1,:)
             elseif(itype==2)then
                dOUT%d(jspin,1,:)    = dIN%d(ispin,1,:)
             else
                dOUT%v(jspin,jorb,:) = dIN%v(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,1,:)    = dIN%e(ispin,1,:)
             dOUT%d(jspin,1,:)    = dIN%d(ispin,1,:)
             dOUT%v(jspin,jorb,:) = dIN%v(ispin,iorb,:)
          endif
       case ("nonsu2")
          if(present(itype))then
             if(itype==1)then
                dOUT%e(jspin,1,:)     = dIN%e(ispin,1,:)
             elseif(itype==2)then
                dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             else
                dOUT%u(jspin,jorb,:)  = dIN%u(ispin,iorb,:)
             endif
          else
             dOUT%e(jspin,1,:)     = dIN%e(ispin,1,:)
             dOUT%v(jspin,jorb,:)  = dIN%v(ispin,iorb,:)
             dOUT%u(jspin,jorb,:)  = dIN%u(ispin,iorb,:)
          endif
       end select
    end select
    call get_dmft_bath(dOUT,bathOUT)
    call deallocate_dmft_bath(dIN)
    call deallocate_dmft_bath(dOUT)
  end subroutine copy_spin_orb_component_bath
  !
  !This is super-seeded by the spin_orb procedure, called for each ispin value.
  !In this way it is possible to interface all the procedures to a single one.
  ! subroutine copy_orb_component_bath(bathIN,iorb,bathOUT,jorb,itype)
  !   real(8),dimension(:)     :: bathIN,bathOUT
  !   integer                  :: iorb,jorb
  !   integer,optional         :: itype
  !   logical                  :: check
  !   type(effective_bath)     :: dIN,dOUT
  !   !
  !   check= check_size_bath(bathIN)
  !   if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
  !   check= check_size_bath(bathOUT)
  !   if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
  !   call allocate_dmft_bath(dIN)
  !   call allocate_dmft_bath(dOUT)
  !   call set_dmft_bath(bathIN,dIN)
  !   call set_dmft_bath(bathOUT,dOUT)
  !   if(present(itype))call check_type_bath(itype)
  !   if(iorb>Norb.OR.jorb>Norb)stop "copy_orb_component_bath error: iorb/jorb > Norb"
  !   !
  !   select case(bath_type)
  !   case default
  !      select case(ed_mode)
  !      case default
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            else
  !               dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !         endif
  !      case ("superc")
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            elseif(itype==2)then
  !               dOUT%d(:,jorb,:)  = dIN%d(:,iorb,:)
  !            else
  !               dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            dOUT%d(:,jorb,:)  = dIN%d(:,iorb,:)
  !            dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !         endif
  !      case ("nonsu2")
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            elseif(itype==2)then
  !               dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            else
  !               dOUT%u(:,jorb,:)  = dIN%u(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,jorb,:)  = dIN%e(:,iorb,:)
  !            dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            dOUT%u(:,jorb,:)  = dIN%u(:,iorb,:)
  !         endif
  !      end select
  !   case('hybrid')
  !      select case(ed_mode)
  !      case default
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,1,:)    = dIN%e(:,1,:)
  !            else
  !               dOUT%v(:,jorb,:) = dIN%v(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,1,:)     = dIN%e(:,1,:)
  !            dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !         endif
  !      case ("superc")
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,1,:)    = dIN%e(:,1,:)
  !            elseif(itype==2)then
  !               dOUT%d(:,1,:)    = dIN%d(:,1,:)
  !            else
  !               dOUT%v(:,jorb,:) = dIN%v(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,1,:)    = dIN%e(:,1,:)
  !            dOUT%d(:,1,:)    = dIN%d(:,1,:)
  !            dOUT%v(:,jorb,:) = dIN%v(:,iorb,:)
  !         endif
  !      case ("nonsu2")
  !         if(present(itype))then
  !            if(itype==1)then
  !               dOUT%e(:,1,:)     = dIN%e(:,1,:)
  !            elseif(itype==2)then
  !               dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            else
  !               dOUT%u(:,jorb,:)  = dIN%u(:,iorb,:)
  !            endif
  !         else
  !            dOUT%e(:,1,:)     = dIN%e(:,1,:)
  !            dOUT%v(:,jorb,:)  = dIN%v(:,iorb,:)
  !            dOUT%u(:,jorb,:)  = dIN%u(:,iorb,:)
  !         endif
  !      end select
  !   end select
  !   call get_dmft_bath(dOUT,bathOUT)
  !   call deallocate_dmft_bath(dIN)
  !   call deallocate_dmft_bath(dOUT)
  ! end subroutine copy_orb_component_bath




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
    integer                :: bath_size,shift,ibath
    complex(8),dimension(Norb):: dum
    save_=.true.;if(present(save))save_=save
    if(Nspin==1)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if (bath_type/="replica") then
       dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
       dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
       if(ed_mode=="superc")dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    else
       stop"spin symmetrize not implemented for replica"
    endif
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath

  subroutine SOC_symmetrize_bath(bath_,save,rot_)
    real(8),dimension(:),intent(inout)          :: bath_
    logical,optional                            :: save
    complex(8),allocatable,intent(in)           :: rot_(:,:)
    type(effective_bath)                        :: dmft_bath_
    logical                                     :: save_
    integer                                     :: ispin,jspin,iorb,jorb,ibath,io,jo
    complex(8)                                  :: Hrepl(Nspin*Norb,Nspin*Norb,Nbath)
    complex(8)                                  :: Hrepl_tilda(Nspin*Norb,Nspin*Norb,Nbath)
    !
    write(LOGfile,"(A)")"SOC_symmetrize_bath"
    save_=.true.;if(present(save))save_=save
    if(bath_type/="replica")then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"SOC_symmetrize_bath: not in replica bath mode, no SOC"
       return
    endif
    if(Nspin==1)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"SOC_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    do ibath=1,Nbath
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Hrepl(io,jo,ibath)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    do ibath=1,Nbath
       Hrepl_tilda(:,:,ibath)=matmul(transpose(conjg(rot_)),matmul(Hrepl(:,:,ibath),rot_))
       Hrepl(:,:,ibath)=zero
    enddo
    !
    do ibath=1,Nbath
       Hrepl_tilda(1,1,ibath)=(Hrepl_tilda(1,1,ibath)+Hrepl_tilda(2,2,ibath))/2.d0
       Hrepl_tilda(2,2,ibath)=(Hrepl_tilda(1,1,ibath)+Hrepl_tilda(2,2,ibath))/2.d0
       Hrepl_tilda(3,3,ibath)=(Hrepl_tilda(3,3,ibath)+Hrepl_tilda(4,4,ibath)+Hrepl_tilda(5,5,ibath)+Hrepl_tilda(6,6,ibath))/2.d0
       Hrepl_tilda(4,4,ibath)=(Hrepl_tilda(3,3,ibath)+Hrepl_tilda(4,4,ibath)+Hrepl_tilda(5,5,ibath)+Hrepl_tilda(6,6,ibath))/2.d0
       Hrepl_tilda(5,5,ibath)=(Hrepl_tilda(3,3,ibath)+Hrepl_tilda(4,4,ibath)+Hrepl_tilda(5,5,ibath)+Hrepl_tilda(6,6,ibath))/2.d0
       Hrepl_tilda(6,6,ibath)=(Hrepl_tilda(3,3,ibath)+Hrepl_tilda(4,4,ibath)+Hrepl_tilda(5,5,ibath)+Hrepl_tilda(6,6,ibath))/2.d0
       do io=1,Nspin*Norb
          do jo=1,Nspin*Norb
             if(io/=jo)Hrepl_tilda(:,:,ibath)=zero
          enddo
       enddo
    enddo
    !
    do ibath=1,Nbath
       Hrepl(:,:,ibath)=matmul(rot_,matmul(Hrepl_tilda(:,:,ibath),transpose(conjg(rot_))))
    enddo
  end subroutine SOC_symmetrize_bath

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
  !     USER BATH CHECKS:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_size_bath(bath_,Hloc_nn) result(bool)
    real(8),dimension(:) :: bath_
    integer              :: Ntrue
    logical              :: bool
    complex(8),optional,allocatable,intent(in) :: Hloc_nn(:,:,:,:)
    if (present(Hloc_nn))then
       Ntrue = get_size_bath(Hloc_nn)
    else
       Ntrue = get_size_bath()
    endif
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




END MODULE ED_BATH_USER
