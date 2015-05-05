!+-----------------------------------------------------------------------------+!
!PURPOSE:  G0 and F0 non-interacting Green's functions on the Matsubara axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
!NORMAL:
function g0and_bath_mats_main(x,dmft_bath_) result(G0and)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: G0and
  integer                                     :: iorb,jorb,ispin,jspin
  complex(8)                                  :: det
  complex(8)                                  :: fg,delta,ff,fdelta
  complex(8),dimension(:,:),allocatable       :: fgorb,zeta
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        do ispin=1,Nspin
           do iorb=1,Norb
              delta = delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
              fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) - delta
              G0and(ispin,ispin,iorb,iorb) = one/fg
           enddo
        enddo
        !
     case ("superc")
        do ispin=1,Nspin
           do iorb=1,Norb
              delta =  delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
              fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
              fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) -  delta
              ff    =                                          - fdelta
              det   = abs(fg)**2 + ff**2
              G0and(ispin,ispin,iorb,iorb) = conjg(fg)/det
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
        G0and=zero
        do ispin=1,Nspin         !Spin diagonal
           zeta = zero
           fgorb= zero
           do iorb=1,Norb
              zeta(iorb,iorb) = x + xmu
           enddo
           do iorb=1,Norb
              do jorb=1,Norb
                 fgorb(iorb,jorb) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-delta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
              enddo
           enddo
           call inv(fgorb)
           G0and(ispin,ispin,:,:)=fgorb
        enddo
        deallocate(fgorb,zeta)
        !
     case ("superc")
        allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
        G0and = zero
        do ispin=1,Nspin
           zeta = zero
           fgorb= zero
           do iorb=1,Norb
              zeta(iorb,iorb)           = x + xmu
              zeta(iorb+Norb,iorb+Norb) = x - xmu
           enddo
           do iorb=1,Norb
              do jorb=1,Norb
                 fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - delta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - fdelta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - fdelta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg(delta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_))
              enddo
           enddo
           call inv(fgorb)
           G0and(ispin,ispin,:,:) = fgorb(1:Norb,1:Norb)
        enddo
        deallocate(fgorb,zeta)
        !
     end select
  end select
end function g0and_bath_mats_main


function g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: G0and
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  G0out = G0and(ispin,jspin,:,:)
end function g0and_bath_mats_ispin_jspin


function g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
  integer,intent(in)                          :: iorb,jorb,ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8)                                  :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: G0and
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  G0out = G0and(ispin,jspin,iorb,jorb)
end function g0and_bath_mats_ispin_jspin_iorb_jorb


function g0and_bath_mats_main_(x,bath_) result(G0and)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: G0and
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_main_


function g0and_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: G0out
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  integer                                     :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_ispin_jspin_

function g0and_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: G0out
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_ispin_jspin_iorb_jorb_










!ANOMALous:
function f0and_bath_mats_main(x,dmft_bath_) result(F0and)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: F0and
  integer                                     :: iorb,jorb,ispin,jspin
  complex(8)                                  :: det
  complex(8)                                  :: fg,delta,ff,fdelta
  complex(8),dimension(:,:),allocatable       :: fgorb,zeta
  !
  F0and=zero
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        stop "F0and_bath_mats error: called with ed_mode=normal, bath_type=normal"
        !
     case ("superc")
        do ispin=1,Nspin
           do iorb=1,Norb
              delta =  delta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
              fdelta= fdelta_bath_mats(ispin,ispin,iorb,iorb,x,dmft_bath_)
              fg    = x + xmu - impHloc(ispin,ispin,iorb,iorb) -  delta
              ff    =                                          - fdelta
              det   = abs(fg)**2 + ff**2
              F0and(ispin,ispin,iorb,iorb) = ff/det
           enddo
        enddo
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        stop "F0and_bath_mats error: called with ed_mode=normal, bath_type=hybrid"
        !
     case ("superc")
        allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
        do ispin=1,Nspin
           zeta = zero
           fgorb= zero
           do iorb=1,Norb
              zeta(iorb,iorb)           = x + xmu
              zeta(iorb+Norb,iorb+Norb) = x - xmu
           enddo
           do iorb=1,Norb
              do jorb=1,Norb
                 fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - delta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - fdelta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - fdelta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_)
                 fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg(delta_bath_mats(ispin,ispin,iorb,jorb,x,dmft_bath_))
              enddo
           enddo
           call inv(fgorb)
           F0and(ispin,ispin,:,:) = fgorb(iorb+Norb,jorb+Norb)
        enddo
        deallocate(fgorb,zeta)
        !
     end select
     !
  end select
end function f0and_bath_mats_main

function f0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: F0and
  F0and = F0and_bath_mats_main(x,dmft_bath_)
  F0out = F0and(ispin,jspin,:,:)
end function f0and_bath_mats_ispin_jspin

function f0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
  integer,intent(in)                          :: iorb,jorb,ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8)                                  :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: F0and
  F0and = f0and_bath_mats_main(x,dmft_bath_)
  F0out = F0and(ispin,jspin,iorb,jorb)
end function f0and_bath_mats_ispin_jspin_iorb_jorb

function f0and_bath_mats_main_(x,bath_) result(F0and)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: F0and
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0and = f0and_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_main_

function f0and_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: F0out
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  integer                                     :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_ispin_jspin_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = f0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_ispin_jspin_

function f0and_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: F0out
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = f0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_ispin_jspin_iorb_jorb_
