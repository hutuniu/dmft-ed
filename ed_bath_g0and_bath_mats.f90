!+-----------------------------------------------------------------------------+!
!PURPOSE:  G0 and F0 non-interacting Green's functions on the Matsubara axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
!NORMAL:
function g0and_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_) result(fg)
  integer,intent(in)                  :: iorb,jorb,ispin,jspin
  type(effective_bath)                :: dmft_bath_
  complex(8),intent(in)               :: x
  integer                             :: l,m,s,r,k
  complex(8)                          :: det
  complex(8)                          :: fg,delta,ff,fdelta
  complex(8),dimension(Norb,Norb)     :: fgorb,zeta,dorb
  complex(8),dimension(2*Norb,2*Norb) :: szeta,sfgorb
  complex(8),dimension(Nhel,Nhel) :: fghel,dhel
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        fg    = one/fg
        !
     case ("superc")
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        ff    =                                          - fdelta
        det   = abs(fg)**2 + ff**2
        fg    = conjg(fg)/det
        !
     case ("nonsu2")
        !
        stop "G0And_bath_mats error: called with ed_mode=nonsu2, bath_type=normal. THIS IS NOT YET CHECKED"
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        zeta = zero
        do l=1,Norb
           zeta(l,l) = x + xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              fgorb(l,m) = zeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(fgorb)
        fg = fgorb(iorb,jorb)
        !
     case ("superc")
        szeta = zero
        do l=1,Norb
           szeta(l,l)           = x + xmu
           szeta(l+Norb,l+Norb) = x - xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              sfgorb(l,m)           = szeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l,m+Norb)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m+Norb) = szeta(l,m) + impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(sfgorb)
        fg = sfgorb(iorb,jorb)
        !
     case ("nonsu2")
        !
        stop "G0AND_bath_mats error: called with ed_mode=nonsu2, bath_type=hybrid. THIS IS NOT YET CHECKED"
        !
     end select
  end select
end function g0and_bath_mats_1

function g0and_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_) result(fg)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: fg
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fg = g0and_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_2

function g0and_bath_mats_3(ispin,jspin,x,dmft_bath_) result(fgorb)
  integer,intent(in)                  :: ispin,jspin
  type(effective_bath)                :: dmft_bath_
  complex(8),intent(in)               :: x
  integer                             :: l,m,s,r,k,iorb,jorb
  complex(8)                          :: det
  complex(8)                          :: fg,delta,ff,fdelta
  complex(8),dimension(Norb,Norb)     :: fgorb,zeta
  complex(8),dimension(2*Norb,2*Norb) :: szeta,sfgorb
  !
  fgorb=zero
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        do iorb=1,Norb
           delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
           fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
           fgorb(iorb,iorb)    = one/fg
        enddo
        !
     case ("superc")
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        ff    =                                          - fdelta
        det   = abs(fg)**2 + ff**2
        fgorb(iorb,iorb)    = conjg(fg)/det
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        zeta = zero
        do l=1,Norb
           zeta(l,l) = x + xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              fgorb(l,m) = zeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(fgorb)
        !fg = fgorb(iorb,jorb)
        !
     case ("superc")
        szeta = zero
        do l=1,Norb
           szeta(l,l)           = x + xmu
           szeta(l+Norb,l+Norb) = x - xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              sfgorb(l,m)           = szeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l,m+Norb)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m+Norb) = szeta(l,m) + impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(sfgorb)
        fgorb(:,:) = sfgorb(:Norb,:Norb)
        !
     end select
  end select
end function g0and_bath_mats_3

function g0and_bath_mats_4(ispin,jspin,x,bath_) result(fgorb)
  integer,intent(in)    :: ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: fgorb(Norb,Norb)
  real(8),dimension(:)  :: bath_
  logical               :: check
  integer               :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fgorb = g0and_bath_mats_3(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_4




!ANOMALous:
function f0and_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_) result(fg)
  type(effective_bath)                :: dmft_bath_
  complex(8),intent(in)               :: x
  integer,intent(in)                  :: iorb,ispin,jorb,jspin
  integer                             :: l,m,s,r,k
  complex(8)                          :: det
  complex(8)                          :: fg,delta,ff,fdelta
  complex(8),dimension(Norb,Norb)     :: fgorb
  complex(8),dimension(2*Norb,2*Norb) :: szeta,sfgorb
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        stop "F0And_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
        fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
        ff    = -fdelta
        det   = abs(fg)**2 + ff**2
        fg    = ff/det
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        stop "F0And_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
        !
     case ("superc")
        szeta = zero
        do l=1,Norb
           szeta(l,l)           = x + xmu
           szeta(l+Norb,l+Norb) = x - xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              sfgorb(l,m)           = szeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l,m+Norb)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m+Norb) = szeta(l,m) + impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(sfgorb)
        fg = sfgorb(iorb+Norb,jorb+Norb)
        !
     end select
  end select
end function f0and_bath_mats_1

function f0and_bath_mats_2(ispin,jspin,iorb,jorb,x,bath_) result(fg)
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  integer,intent(in)    :: iorb,ispin,jorb,jspin
  complex(8)            :: fg
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fg = f0and_bath_mats_1(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_2


function f0and_bath_mats_3(ispin,jspin,x,dmft_bath_) result(fgorb)
  type(effective_bath)                :: dmft_bath_
  complex(8),intent(in)               :: x
  integer,intent(in)                  :: ispin,jspin
  integer                             :: l,m,s,r,k,iorb,jorb
  complex(8)                          :: det
  complex(8)                          :: fg,delta,ff,fdelta
  complex(8),dimension(Norb,Norb)     :: fgorb
  complex(8),dimension(2*Norb,2*Norb) :: szeta,sfgorb
  !
  fgorb=zero
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        stop "F0And_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        do iorb=1,Norb
           delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
           fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
           fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
           ff    = -fdelta
           det   = abs(fg)**2 + ff**2
           fgorb(iorb,iorb)    = ff/det
        enddo
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        stop "F0And_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
        !
     case ("superc")
        szeta = zero
        do l=1,Norb
           szeta(l,l)           = x + xmu
           szeta(l+Norb,l+Norb) = x - xmu
        enddo
        do l=1,Norb
           do m=1,Norb
              sfgorb(l,m)           = szeta(l,m) - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l,m+Norb)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m)      =                                       - fdelta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
              sfgorb(l+Norb,m+Norb) = szeta(l,m) + impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,ispin,l,m,x,dmft_bath_)
           enddo
        enddo
        call inv(sfgorb)
        fgorb(:,:) = sfgorb(:Norb,1+Norb:2*Norb)
        !
     end select
  end select
end function f0and_bath_mats_3

function f0and_bath_mats_4(ispin,jspin,x,bath_) result(fgorb)
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  integer,intent(in)    :: ispin,jspin
  complex(8)            :: fgorb(Norb,Norb)
  real(8),dimension(:)  :: bath_
  logical               :: check
  integer               :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  fgorb = f0and_bath_mats_3(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_4
