!+-----------------------------------------------------------------------------+!
!PURPOSE:  G0 and F0 non-interacting Green's functions on the Matsubara axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
!NORMAL:
function weiss_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_) result(fg)
  integer,intent(in)              :: iorb,jorb,ispin,jspin,k
  type(effective_bath)            :: dmft_bath_
  complex(8),intent(in)           :: x
  integer                         :: l,m,s,r,k
  complex(8)                      :: det
  complex(8)                      :: fg,delta,ff,fdelta
  complex(8),dimension(Nhel,Nhel) :: fghel,dhel
  complex(8),dimension(Norb,Norb) :: fgorb,dorb
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
        !
     case default
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        fg    = one/fg
        !
        !
     case ("superc")
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        ff    = -fdelta
        det   = abs(fg)**2 + ff**2
        fg    = conjg(fg)/det
        !
        !
     case ("nonsu2")
        !
        stop "Weiss_bath_mats error: called with ed_mode=nonsu2, bath_type=normal. THIS IS NOT YET CHECKED"
        !
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        do l=1,Norb
           fgorb(l,l) = x + xmu - impHloc(ispin,ispin,l,l) - delta_bath_mats(ispin,ispin,l,l,x,dmft_bath_)
           do m=l+1,Norb
              fgorb(l,m) = -impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              fgorb(m,l) = -impHloc(ispin,ispin,m,l) - delta_bath_mats(ispin,ispin,m,l,x,dmft_bath_)
           enddo
        enddo
        call inv(fgorb)
        fg = fgorb(iorb,jorb)
        !
        !
     case ("superc")
        !
        stop "Weiss_bath_mats error: called with ed_mode=superc, bath_type=hybrid. THIS IS NOT YET CHECKED"
        !
     case ("nonsu2")
        !
        stop "Weiss_bath_mats error: called with ed_mode=nonsu2, bath_type=hybrid. THIS IS NOT YET CHECKED"
        !
     end select
  end select
end function weiss_bath_mats_1

function weiss_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_) result(fg)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: fg
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "weiss_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fg = weiss_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function weiss_bath_mats_2


!ANOMALous:
function fweiss_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_) result(fg)
  type(effective_bath)            :: dmft_bath_
  complex(8),intent(in)           :: x
  integer,intent(in)              :: iorb,ispin,jorb,jspin,k
  integer                         :: l,m,s,r
  complex(8)                      :: det
  complex(8)                      :: fg,delta,ff,fdelta
  complex(8),dimension(Norb,Norb) :: fgorb,dorb
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default
        !
        stop "Fweiss_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        ff    = -fdelta
        det   = abs(fg)**2 + ff**2
        fg    = ff/det
     case ("nonsu2")
        !
        stop "Fweiss_bath_mats error: called with ed_mode=nonsu2, bath_type=normal"
        !
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        !
        stop "Fweiss_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
        !
     case ("superc")
        !
        stop "Fweiss_bath_mats error: called with ed_mode=superc, bath_type=hybrid. THIS IS NOT YET CHECKED"
        !
     case ("nonsu2")
        !
        stop "Fweiss_bath_mats error: called with ed_mode=nonsu2, bath_type=hybrid"
        !
     end select
  end select
end function fweiss_bath_mats_1

function fweiss_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_) result(fg)
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  integer,intent(in)    :: iorb,ispin,jorb,jspin
  complex(8)            :: fg
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "weiss_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fg = fweiss_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fweiss_bath_mats_2
