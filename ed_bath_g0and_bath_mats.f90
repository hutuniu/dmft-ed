!+-----------------------------------------------------------------------------+!
!PURPOSE:  G0 and F0 non-interacting Green's functions on the Matsubara axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
!NORMAL:
function g0and_bath_mats_main(x,dmft_bath_) result(G0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
  real(8),dimension(size(x))                          :: det
  complex(8),dimension(size(x))                       :: fg,ff
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  !
  G0and = zero
  !
  L=size(x)
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        !
        Delta = delta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin
           do iorb=1,Norb
              fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
              G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
           enddo
        enddo
        !
     case ("superc")
        !
        Delta =  delta_bath_mats(x,dmft_bath_)
        Fdelta= fdelta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin
           do iorb=1,Norb
              fg(:)  = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
              ff(:)  =                                             - Fdelta(ispin,ispin,iorb,iorb,:)
              det(:) = abs(fg(:))**2 + ff(:)*ff(:)
              G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(:))/det(:)
           enddo
        enddo
        !
     case ("nonsu2")
        !
        allocate(fgorb(Nspin,Nspin),zeta(Nspin,Nspin))
        Delta = delta_bath_mats(x,dmft_bath_)
        do i=1,L
           zeta  = (x(i) + xmu)*eye(Nspin)
           fgorb = zero
           do iorb=1,Norb
              do ispin=1,Nspin
                 do jspin=1,Nspin
                    fgorb(ispin,jspin) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
                 enddo
              enddo
              call inv(fgorb)
              do ispin=1,Nspin
                 do jspin=1,Nspin
                    G0and(ispin,jspin,iorb,iorb,i) = fgorb(ispin,jspin)
                 enddo
              enddo
           enddo
        enddo
        deallocate(fgorb,zeta)
        !
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     !
     select case(ed_mode)
     case default
        !
        allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
        Delta = delta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin         !Spin diagonal
           do i=1,L
              fgorb= zero
              zeta = (x(i)+xmu)*eye(Norb)
              do iorb=1,Norb
                 do jorb=1,Norb
                    fgorb(iorb,jorb) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-Delta(ispin,ispin,iorb,jorb,i)
                 enddo
              enddo
              call inv(fgorb)
              G0and(ispin,ispin,:,:,i)=fgorb
           enddo
        enddo
        deallocate(fgorb,zeta)
        !
     case ("superc")
        !
        allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
        Delta  = delta_bath_mats(x,dmft_bath_)
        Fdelta = fdelta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin
           do i=1,L
              zeta = zero
              fgorb= zero
              do iorb=1,Norb
                 zeta(iorb,iorb)           = x(i) + xmu
                 zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
              enddo
              do iorb=1,Norb
                 do jorb=1,Norb
                    fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,i) )
                 enddo
              enddo
              call inv(fgorb)
              G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
           enddo
        enddo
        deallocate(fgorb,zeta)
        !
     case ("nonsu2")
        !
        Nso=Nspin*Norb
        allocate(fgorb(Nso,Nso),zeta(Nso,Nso))
        Delta = delta_bath_mats(x,dmft_bath_)
        do i=1,L
           zeta  = (x(i) + xmu)*eye(Nso)
           fgorb = zero
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       fgorb(io,jo) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                    enddo
                 enddo
              enddo
           enddo
           call inv(fgorb)
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        deallocate(fgorb,zeta)
        !
     end select
  end select
end function g0and_bath_mats_main


function g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
  integer,intent(in)                                  :: ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x))             :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  G0out = G0and(ispin,jspin,:,:,:)
end function g0and_bath_mats_ispin_jspin


function g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
  integer,intent(in)                                  :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8)                                          :: G0out(size(x))
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  G0out = G0and(ispin,jspin,iorb,jorb,:)
end function g0and_bath_mats_ispin_jspin_iorb_jorb


function g0and_bath_mats_main_(x,bath_) result(G0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0and = g0and_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_main_


function g0and_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
  integer,intent(in)                      :: ispin,jspin
  complex(8),dimension(:),intent(in)      :: x
  type(effective_bath)                    :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x)) :: G0out
  real(8),dimension(:)                    :: bath_
  logical                                 :: check
  integer                                 :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_ispin_jspin_

function g0and_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
  integer,intent(in)                 :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in) :: x
  type(effective_bath)               :: dmft_bath_
  complex(8)                         :: G0out(size(x))
  real(8),dimension(:)               :: bath_
  logical                            :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function g0and_bath_mats_ispin_jspin_iorb_jorb_










!ANOMALous:
function f0and_bath_mats_main(x,dmft_bath_) result(F0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,jspin,i,L
  real(8),dimension(size(x))                          :: det
  complex(8),dimension(size(x))                       :: fg,ff
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  !
  F0and=zero
  !
  L = size(x)
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        stop "F0and_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
        !
     case ("superc")
        Delta =  delta_bath_mats(x,dmft_bath_)
        Fdelta= fdelta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin
           do iorb=1,Norb
              fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
              ff(:) =                                              - Fdelta(ispin,ispin,iorb,iorb,:)
              det(:)= abs(fg(:))**2 + ff(:)*ff(:)
              F0and(ispin,ispin,iorb,iorb,:) = ff(:)/det(:)
           enddo
        enddo
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     select case(ed_mode)
     case default
        stop "F0and_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
        !
     case ("superc")
        allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
        Delta =  delta_bath_mats(x,dmft_bath_)
        Fdelta= fdelta_bath_mats(x,dmft_bath_)
        do ispin=1,Nspin
           do i=1,L
              zeta = zero
              fgorb= zero
              do iorb=1,Norb
                 zeta(iorb,iorb)           = x(i) + xmu
                 zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
              enddo
              do iorb=1,Norb
                 do jorb=1,Norb
                    fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,i) )
                 enddo
              enddo
              call inv(fgorb)
              F0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1+Norb:Norb+Norb)
           enddo
        enddo
        deallocate(fgorb,zeta)
        !
     end select
     !
  end select
end function f0and_bath_mats_main

function f0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
  integer,intent(in)                                  :: ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x))             :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  F0and = F0and_bath_mats_main(x,dmft_bath_)
  F0out = F0and(ispin,jspin,:,:,:)
end function f0and_bath_mats_ispin_jspin

function f0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
  integer,intent(in)                                  :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8)                                          :: F0out(size(x))
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  F0and = f0and_bath_mats_main(x,dmft_bath_)
  F0out = F0and(ispin,jspin,iorb,jorb,:)
end function f0and_bath_mats_ispin_jspin_iorb_jorb

function f0and_bath_mats_main_(x,bath_) result(F0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0and = f0and_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_main_

function f0and_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
  integer,intent(in)                      :: ispin,jspin
  complex(8),dimension(:),intent(in)      :: x
  type(effective_bath)                    :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x)) :: F0out
  real(8),dimension(:)                    :: bath_
  logical                                 :: check
  integer                                 :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_ispin_jspin_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = f0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_ispin_jspin_

function f0and_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
  integer,intent(in)                 :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in) :: x
  type(effective_bath)               :: dmft_bath_
  complex(8)                         :: F0out(size(x))
  real(8),dimension(:)               :: bath_
  logical                            :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "f0and_bath_mats_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = f0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function f0and_bath_mats_ispin_jspin_iorb_jorb_
