!+-----------------------------------------------------------------------------+!
!PURPOSE:  Delta and Fdelta functions on the Matsubara axis:
! _1 : input type(effective_bath) dmft_bath
! _2 : input array bath
! Delta_ : normal
! Fdelta_: anomalous
!+-----------------------------------------------------------------------------+!
!NORMAL:
function delta_bath_mats_main(x,dmft_bath_) result(Delta)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  integer                                             :: i,iorb,jorb,ispin,jspin,ih,k,L
  real(8),dimension(Nbath)                            :: eps,dps,vps
  real(8),dimension(Norb,Nbath)                       :: vops
  real(8),dimension(Nspin,Nbath)                      :: hps
  real(8),dimension(Nspin,Nspin,Nbath)                :: wps
  real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wops
  !
  Delta=zero
  !
  L = size(x)
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        !
        !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              do i=1,L
                 Delta(ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
              enddo
           enddo
        enddo
        !
     case ("superc")
        !
        !\Delta_{aa} = - \sum_k [ V_{a}(k) * V_{a}(k) * (iw_n + E_{a}(k)) / Den(k) ]
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              do i=1,L
                 Delta(ispin,ispin,iorb,iorb,i) = -sum( vps(:)*vps(:)*(x(i) + eps(:))/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
              enddo
           enddo
        enddo
        !
     case ("nonsu2")
        !
        !\Delta_{aa}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{a}^{hs`}(k)/(iw_n - H_{a}^{h}(k))]
        do iorb=1,Norb
           hps = dmft_bath_%e(1:Nspin,iorb,1:Nbath)
           wps = get_Whyb_matrix(dmft_bath_%v(1:Nspin,iorb,1:Nbath),dmft_bath_%u(1:Nspin,iorb,1:Nbath))
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do ih=1,Nspin
                    do i=1,L
                       Delta(ispin,jspin,iorb,iorb,i) = Delta(ispin,jspin,iorb,iorb,i) + sum( wps(ispin,ih,:)*wps(ih,jspin,:)/(x(i) - hps(ih,:)) )
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        !
        !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           do iorb=1,Norb
              do jorb=1,Norb
                 do i=1,L
                    Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                 enddo
              enddo
           enddo
        enddo
        !
     case ("superc")
        !
        !\Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (iw_n + E(k)) / Den(k) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           do iorb=1,Norb
              do jorb=1,Norb
                 do i=1,L
                    Delta(ispin,ispin,iorb,jorb,i) = -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                 enddo
              enddo
           enddo
        enddo
        !
     case ("nonsu2")
        !
        !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{hs`}(k)/(iw_n - H^{h}(k))]
        hps  = dmft_bath_%e(1:Nspin,1     ,1:Nbath)
        wops = get_Whyb_matrix(dmft_bath_%v(1:Nspin,1:Norb,1:Nbath),dmft_bath_%u(1:Nspin,1:Norb,1:Nbath))
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,Nspin
                 do jspin=1,Nspin
                    do ih=1,Nspin
                       do i=1,L
                          Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                               sum( wops(ispin,ih,iorb,:)*wops(ih,jspin,jorb,:)/(x(i) - hps(ih,:)) )
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !
     end select
     !
  end select
end function delta_bath_mats_main


function delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
  integer,intent(in)                                  :: ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x))             :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  Delta = delta_bath_mats_main(x,dmft_bath_)
  G0out = Delta(ispin,jspin,:,:,:)
end function delta_bath_mats_ispin_jspin


function delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
  integer,intent(in)                                  :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(size(x))                       :: G0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  Delta = delta_bath_mats_main(x,dmft_bath_)
  G0out = Delta(ispin,jspin,iorb,jorb,:)
end function delta_bath_mats_ispin_jspin_iorb_jorb


function delta_bath_mats_main_(x,bath_) result(Delta)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  Delta = delta_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_mats_main_


function delta_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
  integer,intent(in)                      :: ispin,jspin
  type(effective_bath)                    :: dmft_bath_
  complex(8),dimension(:),intent(in)      :: x
  complex(8),dimension(Norb,Norb,size(x)) :: G0out
  real(8),dimension(:)                    :: bath_
  logical                                 :: check
  integer                                 :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_mats_ispin_jspin_

function delta_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
  integer,intent(in)                 :: iorb,jorb,ispin,jspin
  type(effective_bath)               :: dmft_bath_
  complex(8),dimension(:),intent(in) :: x
  complex(8),dimension(size(x))      :: G0out
  real(8),dimension(:)               :: bath_
  logical                            :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  G0out = delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function delta_bath_mats_ispin_jspin_iorb_jorb_









!ANOMALous:
function fdelta_bath_mats_main(x,dmft_bath_) result(Fdelta)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  integer                                             :: iorb,ispin,jorb,jspin
  real(8),dimension(Norb,Norb)                        :: delta_orb
  real(8),dimension(Nbath)                            :: eps,dps,vps
  real(8),dimension(Norb,Nbath)                       :: vops
  integer                                             :: i,k,L
  !
  Fdelta=zero
  !
  L = size(x)
  !
  select case(bath_type)
  case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
     !
     select case(ed_mode)
     case default
        stop "Fdelta_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
        !
     case ("superc")
        !
        !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
        !
        do ispin=1,Nspin
           do iorb=1,Norb
              eps = dmft_bath_%e(ispin,iorb,1:Nbath)
              dps = dmft_bath_%d(ispin,iorb,1:Nbath)
              vps = dmft_bath_%v(ispin,iorb,1:Nbath)
              do i=1,L
                 Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
              enddo
           enddo
        enddo
        !
     end select
     !
  case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
     !
     select case(ed_mode)
     case default
        stop "Fdelta_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
        !
     case ("superc")
        !
        !\FDelta_{ab} = - \sum_k [ \Delta(k) * V_{a}(k) * V_{b}(k) / Den(k) ]
        do ispin=1,Nspin
           eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
           dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
           vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
           do iorb=1,Norb
              do jorb=1,Norb
                 do i=1,L
                    Fdelta(ispin,ispin,iorb,jorb,i) = -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                 enddo
              enddo
           enddo
        enddo
        !
     end select
  end select
end function fdelta_bath_mats_main

function fdelta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
  integer,intent(in)                                  :: ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x))             :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  Fdelta = fdelta_bath_mats_main(x,dmft_bath_)
  F0out  = Fdelta(ispin,jspin,:,:,:)
end function fdelta_bath_mats_ispin_jspin

function fdelta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
  integer,intent(in)                                  :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(size(x))                       :: F0out
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  Fdelta = fdelta_bath_mats_main(x,dmft_bath_)
  F0out = Fdelta(ispin,jspin,iorb,jorb,:)
end function fdelta_bath_mats_ispin_jspin_iorb_jorb

function fdelta_bath_mats_main_(x,bath_) result(Fdelta)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  Fdelta = fdelta_bath_mats_main(x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_mats_main_

function fdelta_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
  integer,intent(in)                          :: ispin,jspin
  complex(8),dimension(:),intent(in)          :: x
  type(effective_bath)                        :: dmft_bath_
  complex(8),dimension(Norb,Norb,size(x))     :: F0out
  real(8),dimension(:)                        :: bath_
  logical                                     :: check
  integer                                     :: iorb,jorb
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_ispin_jspin_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = fdelta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_mats_ispin_jspin_

function fdelta_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
  integer,intent(in)                 :: iorb,jorb,ispin,jspin
  complex(8),dimension(:),intent(in) :: x
  type(effective_bath)               :: dmft_bath_
  complex(8),dimension(size(x))      :: F0out
  real(8),dimension(:)               :: bath_
  logical                            :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "fdelta_bath_mats_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
  call allocate_bath(dmft_bath_)
  call set_bath(bath_,dmft_bath_)
  F0out = fdelta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
  call deallocate_bath(dmft_bath_)
end function fdelta_bath_mats_ispin_jspin_iorb_jorb_
