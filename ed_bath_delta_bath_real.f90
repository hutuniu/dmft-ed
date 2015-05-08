!+-----------------------------------------------------------------------------+!
!PURPOSE:  Delta and Fdelta functions on the Real axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
function delta_bath_real_main(x,dmft_bath_) result(Delta)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Delta
  integer                                     :: iorb,jorb,ispin,jspin,k,ih
  complex(8),dimension(Nbath)                 :: den
  real(8),dimension(Nbath)                    :: eps,dps,vps
  real(8),dimension(Norb,Nbath)               :: vops
  real(8),dimension(Nspin,Nbath)              :: hps
  real(8),dimension(Nspin,Nspin,Nbath)        :: wps
  real(8),dimension(Nspin,Nspin,Norb,Nbath)   :: wops
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        !
        !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(w+i\h - E_{a}(k)) ]
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              Delta(ispin,ispin,iorb,iorb) = sum( vps(:)*vps(:)/(x - eps(:)) )
           enddo
        enddo
        !
     case ("superc")
        !
        !\Delta_{aa}^{ss} = - \sum_k [ V_{a}(k) * V_{a}(k) * (w+i\h + E_{a}(k)) / Den(k) ]
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              forall(k=1:Nbath)den(k) =  x*(-x) + eps(k)**2 + dps(k)**2 !den(k) = (w+i\h)*(-w-i\h) + E_{a}(k)**2 + \Delta_{a}(k)**2
              Delta(ispin,ispin,iorb,iorb) = -sum( vps(:)*vps(:)*(x + eps(:))/den(:) )
           enddo
        enddo
        !
     case ("nonsu2")
        !
        !\Delta_{aa}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{a}^{hs`}(k)/(w+i\h - H_{a}^{h}(k))]
        !we assume that the %e and %w are correctly popolated in this channel
        hps = dmft_bath_%e(        1:Nspin,iorb,1:Nbath)
        wps = dmft_bath_%w(1:Nspin,1:Nspin,iorb,1:Nbath)
        Delta(ispin,ispin,iorb,iorb) = zero
        do ih=1,Nspin
           Delta(ispin,ispin,iorb,iorb) = Delta(ispin,ispin,iorb,iorb) + sum( wps(ispin,ih,:)*wps(ih,jspin,:)/(x - hps(ih,:)) )
        enddo
        !
     end select
     !
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     !
     select case(ed_mode)
     case default
        !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(w+i\h  - E(k)) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           do iorb=1,Norb
              do jorb=1,Norb
                 Delta(ispin,ispin,iorb,jorb) = sum( vops(iorb,:)*vops(jorb,:)/(x - eps(:)) )
              enddo
           enddo
        enddo
        !
     case ("superc")
        !
        !\Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (w+i\h + E(k)) / Den(k) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1      ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1      ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           forall(k=1:Nbath)den(k) =  x*(-x) + eps(k)**2 + dps(k)**2 !den(k) = ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2
           do iorb=1,Norb
              do jorb=1,Norb
                 Delta(ispin,ispin,iorb,jorb) = -sum( vops(iorb,:)*vops(jorb,:)*(x + eps(:))/Den(:) )
              enddo
           enddo
        enddo
     case ("nonsu2")
        !
        !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{hs`}(k)/(w+i\h - H^{h}(k))]
        !we assume that the %e and %w are correctly popolated in this channel
        hps  = dmft_bath_%e(        1:Nspin,1     ,1:Nbath)
        wops = dmft_bath_%w(1:Nspin,1:Nspin,1:Norb,1:Nbath)
        Delta(ispin,ispin,iorb,iorb) = zero
        do ih=1,Nspin
           Delta(ispin,ispin,iorb,iorb) = Delta(ispin,ispin,iorb,iorb) + sum( wops(ispin,ih,iorb,:)*wops(ih,jspin,jorb,:)/(x - hps(ih,:)) )
        enddo
        !
     end select
     !
  end select
end function delta_bath_real_main


function delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Delta
  Delta = delta_bath_real_main(x,dmft_bath_)
  G0out = Delta(ispin,jspin,:,:)
end function delta_bath_real_ispin_jspin


function delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
  integer,intent(in)                          :: iorb,jorb,ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8)                                  :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Delta
  Delta = delta_bath_real_main(x,dmft_bath_)
  G0out = Delta(ispin,jspin,iorb,jorb)
end function delta_bath_real_ispin_jspin_iorb_jorb



function delta_bath_real_main_(x,bath_) result(Delta)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Delta
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_real_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  Delta = delta_bath_real_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_real_main_


function delta_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: G0out
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  integer                                     :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_real_ispin_jspin_

function delta_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: G0out
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_real_ispin_jspin_iorb_jorb_



















!ANOMALOUS:
function fdelta_bath_real_main(x,dmft_bath_) result(Fdelta)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Fdelta
  integer                                     :: iorb,ispin,jorb,jspin,ih
  real(8),dimension(Norb,Norb)                :: delta_orb
  complex(8),dimension(Nbath)                 :: den
  real(8),dimension(Nbath)                    :: eps,dps,vps
  real(8),dimension(Norb,Nbath)               :: vops
  integer                                     :: k
  !
  Fdelta=zero
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     select case(ed_mode)
     case default
        stop "Fdelta_bath_real error: called with ed_mode=normal/nonsu2, bath_type=normal"
        !
     case ("superc")
        !
        !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              forall(k=1:Nbath)den(k) =  x*(-x) + eps(k)**2 + dps(k)**2 !den(k) = (w+i\h)*(-w-i\h) + E_{a}(k)**2 + \D_{a}(k)**2
              Fdelta(ispin,ispin,iorb,iorb) = sum( dps(:)*vps(:)*vps(:)/Den(:) )
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        stop "Fdelta_bath_real error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
        !
     case ("superc")
        !
        !\FDelta_{ab} = - \sum_k [ \Delta(k) * V_{a}(k) * V_{b}(k) / Den(k) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           forall(k=1:Nbath)den(k) =  x*(-x) + eps(k)**2 + dps(k)**2 !den(k) = (w+i\h)*(-w-i\h) + E_{a}(k)**2 + \D_{a}(k)**2
           do iorb=1,Norb
              do jorb=1,Norb
                 Fdelta(ispin,ispin,iorb,jorb) = -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/Den(:) )
              enddo
           enddo
        enddo
        !
     end select
     !
  end select
end function fdelta_bath_real_main

function fdelta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Fdelta
  Fdelta = fdelta_bath_real_main(x,dmft_bath_)
  F0out  = Fdelta(ispin,jspin,:,:)
end function fdelta_bath_real_ispin_jspin

function fdelta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
  integer,intent(in)                          :: iorb,jorb,ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8)                                  :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Fdelta
  Fdelta = fdelta_bath_real_main(x,dmft_bath_)
  F0out = Fdelta(ispin,jspin,iorb,jorb)
end function fdelta_bath_real_ispin_jspin_iorb_jorb

function fdelta_bath_real_main_(x,bath_) result(Fdelta)
  complex(8),intent(in)                       :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Fdelta
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_real_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  Fdelta = fdelta_bath_real_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_real_main_

function fdelta_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
  integer,intent(in)                          :: ispin,jspin
  type(effective_bath)                        :: dmft_bath_
  complex(8),intent(in)                       :: x
  complex(8),dimension(Norb,Norb)             :: F0out
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  integer                                     :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_real_ispin_jspin_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = fdelta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_real_ispin_jspin_

function fdelta_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
  integer,intent(in)    :: iorb,jorb,ispin,jspin
  type(effective_bath)  :: dmft_bath_
  complex(8),intent(in) :: x
  complex(8)            :: F0out
  real(8),dimension(:)  :: bath_
  logical               :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_real_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = fdelta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_real_ispin_jspin_iorb_jorb_
