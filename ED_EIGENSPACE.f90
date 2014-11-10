module ED_EIGENSPACE
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none
  private

  type sparse_estate
     integer                         :: sector        !index of the sector
     real(8)                         :: e             !energy of the eigen-state
     real(8),dimension(:),pointer    :: vec           !double precision eigen-vector
     complex(8),dimension(:),pointer :: cvec          !double complex eigen-vector
     logical                         :: isreal=.true. !real*8 flag
     logical                         :: itwin=.false. !twin sector label
     type(sparse_estate),pointer     :: twin=>null()   !link to twin box 
     type(sparse_estate),pointer     :: next          !link to next box (chain)
  end type sparse_estate

  type sparse_espace
     integer                     :: size
     real(8)                     :: emax,emin
     logical                     :: status=.false.
     type(sparse_estate),pointer :: root !head/root of the list\== list itself
  end type sparse_espace


  interface es_insert_state
     module procedure es_insert_state_d,es_insert_state_c
  end interface es_insert_state

  interface es_add_state
     module procedure es_add_state_d,es_add_state_c
  end interface es_add_state

  public :: sparse_estate
  public :: sparse_espace
  !
  public :: es_init_espace      !init the espace                 !checked
  public :: es_delete_espace    !del the espace                  !checked
  public :: es_free_espace      !free the espace                 !checked
  public :: es_print_espace     !print the espace                !checked
  !
  public :: es_insert_state     !insert a state                  !checked
  public :: es_add_state        !add a state w/ costraint        !checked
  public :: es_pop_state        !pop a state                     !checked
  !
  public :: es_return_type         !get the type of a state         !checked
  public :: es_return_sector       !get the sector of a state       !checked
  public :: es_return_energy       !get the energy of a state       !checked
  public :: es_return_vector       !get the vector of a state       !checked
  public :: es_return_cvector      !get the vector of a state       !checked
  public :: es_return_gs_degeneracy!get the number of degenerate GS !checked
  !
  type(sparse_espace),public    :: state_list
  !
contains        !some routine to perform simple operation on the lists


  !+------------------------------------------------------------------+
  !PURPOSE  : initialize the list of states
  !+------------------------------------------------------------------+
  function es_init_espace() result(space)
    type(sparse_espace) :: space
    allocate(space%root)
    space%status=.true.
    space%root%next => null()
    space%size=0
    space%emax=-huge(1.d0)
    space%emin=huge(1.d0)
  end function es_init_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : destroy the list of states
  !+------------------------------------------------------------------+
  subroutine es_delete_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       if(associated(c%vec))deallocate(c%vec)
       if(associated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    deallocate(space%root)
    space%status=.false.
    p=>null()
    c=>null()
  end subroutine es_delete_espace




  !+------------------------------------------------------------------+
  !PURPOSE  : empty the list of states
  !+------------------------------------------------------------------+
  subroutine es_free_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       if(associated(c%vec))deallocate(c%vec)
       if(associated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    space%size=0
    space%emax=-huge(1.d0)
    space%emin=huge(1.d0)
    p=>null()
    c=>null()
  end subroutine es_free_espace







  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_add_state_d(espace,e,vec,sector,twin,size,verbose)
    type(sparse_espace),intent(inout) :: espace
    real(8),intent(in)                :: e
    real(8),dimension(:),intent(in)   :: vec
    integer,intent(in)                :: sector
    integer,intent(in),optional       :: size
    logical,intent(in),optional       :: verbose
    logical,intent(in),optional       :: twin
    logical                           :: twin_
    twin_=.false.;if(present(twin))twin_=twin
    if(present(size))then     !if present size add respecting the size constraint.
       if(espace%size<size)then
          call es_insert_state_d(espace,e,vec,sector,twin_)
       else
          if(e < es_return_energy(espace))then
             if(present(verbose).AND.(verbose.eqv..true.))print*,"found a new state:"
             call es_pop_state(espace)
             call es_insert_state_d(espace,e,vec,sector,twin_)
          endif
       endif
    else                      !else add normally
       call es_insert_state_d(espace,e,vec,sector,twin_)
    endif
  end subroutine es_add_state_d

  subroutine es_add_state_c(espace,e,cvec,sector,twin,size,verbose)
    type(sparse_espace),intent(inout) :: espace
    real(8),intent(in)                :: e
    complex(8),dimension(:),intent(in):: cvec
    integer,intent(in)                :: sector
    integer,intent(in),optional       :: size
    logical,intent(in),optional       :: verbose
    logical,intent(in),optional       :: twin
    logical                           :: twin_
    twin_=.false.;if(present(twin))twin_=twin
    if(present(size))then !if present size add respecting the size constraint.
       if(espace%size<size)then
          call es_insert_state_c(espace,e,cvec,sector,twin_)
       else
          if(e < es_return_energy(espace))then
             if(present(verbose).AND.(verbose.eqv..true.))print*,"found a new state:"
             call es_pop_state(espace)
             call es_insert_state_c(espace,e,cvec,sector,twin_)
          endif
       endif
    else                      !else add normally
       call es_insert_state_c(espace,e,cvec,sector,twin_)
    endif
  end subroutine es_add_state_c




  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_insert_state_d(space,e,vec,sector,twin)
    type(sparse_espace),intent(inout) :: space
    real(8),intent(in)                :: e
    real(8),dimension(:),intent(in)   :: vec
    integer,intent(in)                :: sector
    logical                           :: twin
    type(sparse_estate),pointer       :: p,c,pn,pt
    integer                           :: step
    p => space%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    if(e > space%emax)space%emax=e !update the max energy (corresponds to the top entry)
    if(e < space%emin)space%emin=e !update the min energy (corresponds to the first entry)
    p%next%e = e    
    allocate(p%next%vec(size(vec)))
    p%next%vec = vec
    p%next%isreal=.true.
    p%next%itwin=.false.
    p%next%sector=sector
    space%size = space%size+1
    if(twin)then               !Create a twin element in the list with same energy, no vector and twin flag T
       allocate(p%next%next)
       p%next%next%e = e    
       p%next%next%isreal=.true.
       p%next%next%itwin=.true.
       p%next%next%sector=get_twin_sector(sector)
       p%next%next%twin => p%next
       space%size = space%size+1
    endif
    if(.not.associated(c))then !end of the list special case
       if(.not.twin)then
          p%next%next  => null()
       else
          p%next%next%next  => null()
       endif
    else
       if(.not.twin)then
          p%next%next  => c      !the %next of the new node come to current
       else
          p%next%next%next  => c 
       endif
    end if
    p=>null()
    c=>null()
  end subroutine es_insert_state_d

  subroutine es_insert_state_c(space,e,vec,sector,twin)
    type(sparse_espace),intent(inout)  :: space
    real(8),intent(in)                 :: e
    complex(8),dimension(:),intent(in) :: vec
    integer,intent(in)                 :: sector
    logical                            :: twin
    type(sparse_estate),pointer        :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until e < value (ordered list)
       if(.not.associated(c))exit
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    p%next%e = e
    if(e > space%emax)space%emax=e !update the max energy (corresponds to the top entry)
    if(e < space%emin)space%emin=e !update the min energy (corresponds to the first entry)
    allocate(p%next%cvec(size(vec)))
    p%next%cvec = vec
    p%next%isreal=.false.
    p%next%itwin=.false.
    p%next%sector=sector
    space%size = space%size+1
    if(twin)then !Create a twin element in the list with same energy, no vector and twin flag T
       allocate(p%next%next)
       p%next%next%e = e    
       p%next%next%isreal=.false.
       p%next%next%itwin=.true.
       p%next%next%sector=get_twin_sector(sector)
       p%next%next%twin => p%next
       space%size = space%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       if(.not.twin)then
          p%next%next  => null()
       else
          p%next%next%next  => null()
       endif
    else
       if(.not.twin)then
          p%next%next  => c      !the %next of the new node come to current
       else
          p%next%next%next  => c
       endif
    end if
    p=>null()
    c=>null()
  end subroutine es_insert_state_c






  !+------------------------------------------------------------------+
  !PURPOSE  : remove last element from the list, if +n is given remove 
  ! the n-th element, if +e is given remove the state with state%e=e
  !+------------------------------------------------------------------+
  subroutine es_pop_state(space,n)
    type(sparse_espace),intent(inout) :: space
    integer,optional,intent(in)       :: n
    integer                           :: i,pos
    type(sparse_estate),pointer       :: p,c
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)stop "es_pop_state: pos > espace.size"
    if(space%size==0)stop "es_pop_state: empty list"
    c => space%root
    do i=1,pos
       p => c
       c => c%next
       if(.not.associated(c))return !empty or end of the list
    end do
    p%next => c%next !reallocate and skip the deleted link
    if(associated(c%vec))deallocate(c%vec)
    if(associated(c%cvec))deallocate(c%cvec)
    if(associated(c%twin))c%twin=>null()
    deallocate(c)           !free link
    if(pos==space%size)then     !pop last term carrying e=emax, update emax
       space%emax = p%e
    elseif(pos==1)then          !pop first term carrying e=emin, update emin
       space%emin = p%e
    endif
    space%size=space%size-1
    p=>null()
    c=>null()
  end subroutine es_pop_state






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_gs_degeneracy(space,gsthreshold) result(numzero)
    type(sparse_espace),intent(in) :: space
    real(8),optional               :: gsthreshold
    real(8)                        :: gsthreshold_
    integer                        :: numzero,pos
    type(sparse_estate),pointer    :: p,c
    real(8)                        :: oldzero,enemin
    gsthreshold_=1.d-9;if(present(gsthreshold))gsthreshold_=gsthreshold
    if(.not.space%status) stop "es_return_gs_degeneracy: espace not allocated"
    oldzero=1000.d0
    numzero=0
    c => space%root
    pos=0
    do 
       c => c%next
       pos=pos+1
       if(.not.associated(c))exit !end of the list
       enemin = c%e
       if (enemin < oldzero-10.d0*gsthreshold_) then
          numzero=1
          oldzero=enemin
       elseif(abs(enemin-oldzero) <= gsthreshold_)then
          numzero=numzero+1
          oldzero=min(oldzero,enemin)
       endif
    end do
  end function es_return_gs_degeneracy




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_sector(space,n) result(sector)
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n
    integer                        :: sector
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_sector: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_sector: n > espace.size"
    sector=0
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    sector = c%sector
    c=>null()
  end function es_return_sector



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_type(space,n) result(ctype)
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n
    logical                        :: ctype
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_type: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)stop "es_return_type: n > espace.size"
    if(space%size==0)stop "es_return_type: espace empty"
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(.not.c%itwin)then
       ctype = c%isreal
    else
       ctype = c%twin%isreal
    endif
  end function es_return_type



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_energy(space,n) result(egs)
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n
    real(8)                        :: egs
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_energy: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)    stop "es_return_energy: n > espace.size"
    c => space%root
    egs=space%emax
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    if(.not.c%itwin)then
       egs = c%e
    else
       egs = c%twin%e
    endif
  end function es_return_energy



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_vector(space,n) result(vector)
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n
    real(8),dimension(:),pointer     :: vector
    type(sparse_estate),pointer      :: c
    integer                          :: i,pos
    integer                          :: dim,sector
    integer,dimension(:),allocatable :: order
    if(.not.space%status) stop "es_return_vector: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_vector: n > espace.size"
    if(space%size==0)stop "es_return_vector: espace emtpy"
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(.not.c%itwin)then
       if(.not.c%isreal)stop "pop vector isreal=F: can not associate to complex vector"
       vector => c%vec
    else
       if(.not.c%twin%isreal)stop "pop vector isreal=F: can not associate to complex vector"
       dim = getdim(c%sector)
       allocate(Order(dim))
       call twin_sector_order(c%twin%sector,Order)
       allocate(vector(dim))
       do i=1,dim
          vector(i) = c%twin%vec(Order(i))
       enddo
       deallocate(order)
    endif
  end function es_return_vector

  function es_return_cvector(space,n) result(vector)
    type(sparse_espace),intent(in)  :: space
    integer,optional,intent(in)    :: n
    complex(8),dimension(:),pointer :: vector
    type(sparse_estate),pointer     :: c
    integer                         :: i,pos
    integer                          :: dim,sector
    integer,dimension(:),allocatable :: order
    if(.not.space%status) stop "es_return_cvector: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_cvector: n > espace.size"
    if(space%size==0)stop "es_return_cvector: espace emtpy"
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(.not.c%itwin)then
       if(c%isreal)stop "pop cvector isreal=T: can not associate to real vector"
       vector => c%cvec
    else
       if(c%twin%isreal)stop "pop cvector isreal=T: can not associate to real vector"
       dim = getdim(c%sector)
       allocate(Order(dim))
       call twin_sector_order(c%twin%sector,Order)
       allocate(vector(dim))
       do i=1,dim
          vector(i) = c%twin%vec(Order(i))
       enddo
       deallocate(order)
    endif
  end function es_return_cvector







  !+------------------------------------------------------------------+
  !PURPOSE  : pretty print the list of states
  !+------------------------------------------------------------------+
  subroutine es_print_espace(space,unit,wvec)
    type(sparse_espace),intent(in) :: space
    type(sparse_estate),pointer    :: c
    integer                        :: counter,i
    integer,optional               :: unit
    integer                        :: unit_
    logical,optional               :: wvec
    logical                        :: wvec_
    wvec_=.false.;if(present(wvec))wvec_=wvec
    unit_=6;if(present(unit))unit_=unit
    write(*,"(A,I3)")"Print sparse espace unit ->",unit_
    c => space%root%next   !assume is associated,ie list exists
    counter = 0
    if(space%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          write(unit_,"(A10,I5)")   "Index   : ",counter
          write(unit_,"(A10,I5)")   "Sector  : ",c%sector
          write(unit_,"(A10,3L)")   "Twin    : ",c%itwin,associated(c%vec),associated(c%cvec)
          write(unit_,"(A10,I5)")   "Size    : ",getdim(c%sector)!size(c%vec)
          write(unit_,"(A10,f18.9)")"Energy  : ",c%e
          if(wvec_)then
             if(c%isreal)then
                write(unit_,"(A10)")"Vec     : "
                write(unit_,*)c%vec
             else
                write(unit_,"(A10)")"Vec     : "
                write(unit_,*)c%cvec
             endif
          endif
          c => c%next  !traverse list
          write(unit_,*)""
       end do
    else
       write(unit_,*)"Empty space"
       return
    endif
    c=>null()
  end subroutine es_print_espace



end module ED_EIGENSPACE
