MODULE ED_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,zeye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH_DMFT
  USE ED_BATH_USER
  USE ED_AUX_FUNX
  implicit none

  private


  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface delta_bath_mats
     module procedure delta_bath_mats_main
     module procedure delta_bath_mats_ispin_jspin
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb
     module procedure delta_bath_mats_main_
     module procedure delta_bath_mats_ispin_jspin_
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb_
  end interface delta_bath_mats
  !
  !ANOMALOUS
  interface fdelta_bath_mats
     module procedure fdelta_bath_mats_main
     module procedure fdelta_bath_mats_ispin_jspin
     module procedure fdelta_bath_mats_ispin_jspin_iorb_jorb
     module procedure fdelta_bath_mats_main_
     module procedure fdelta_bath_mats_ispin_jspin_
     module procedure fdelta_bath_mats_ispin_jspin_iorb_jorb_
  end interface fdelta_bath_mats




  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION REAL
  !
  !##################################################################
  !NORMAL
  interface delta_bath_real
     module procedure delta_bath_real_main
     module procedure delta_bath_real_ispin_jspin
     module procedure delta_bath_real_ispin_jspin_iorb_jorb
     module procedure delta_bath_real_main_
     module procedure delta_bath_real_ispin_jspin_
     module procedure delta_bath_real_ispin_jspin_iorb_jorb_
  end interface delta_bath_real
  !
  !ANOMALOUS
  interface fdelta_bath_real
     module procedure fdelta_bath_real_main
     module procedure fdelta_bath_real_ispin_jspin
     module procedure fdelta_bath_real_ispin_jspin_iorb_jorb
     module procedure fdelta_bath_real_main_
     module procedure fdelta_bath_real_ispin_jspin_
     module procedure fdelta_bath_real_ispin_jspin_iorb_jorb_
  end interface fdelta_bath_real



  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_mats
     module procedure g0and_bath_mats_main
     module procedure g0and_bath_mats_ispin_jspin
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb
     module procedure g0and_bath_mats_main_
     module procedure g0and_bath_mats_ispin_jspin_
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface g0and_bath_mats
  !
  !ANOMALOUS
  interface f0and_bath_mats
     module procedure f0and_bath_mats_main
     module procedure f0and_bath_mats_ispin_jspin
     module procedure f0and_bath_mats_ispin_jspin_iorb_jorb
     module procedure f0and_bath_mats_main_
     module procedure f0and_bath_mats_ispin_jspin_
     module procedure f0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface f0and_bath_mats



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_mats
     module procedure invg0_bath_mats_main
     module procedure invg0_bath_mats_ispin_jspin
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb
     module procedure invg0_bath_mats_main_
     module procedure invg0_bath_mats_ispin_jspin_
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb_
  end interface invg0_bath_mats
  !
  !ANOMALOUS
  interface invf0_bath_mats
     module procedure invf0_bath_mats_main
     module procedure invf0_bath_mats_ispin_jspin
     module procedure invf0_bath_mats_ispin_jspin_iorb_jorb
     module procedure invf0_bath_mats_main_
     module procedure invf0_bath_mats_ispin_jspin_
     module procedure invf0_bath_mats_ispin_jspin_iorb_jorb_
  end interface invf0_bath_mats



  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_real
     module procedure g0and_bath_real_main
     module procedure g0and_bath_real_ispin_jspin
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb
     module procedure g0and_bath_real_main_
     module procedure g0and_bath_real_ispin_jspin_
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb_
  end interface g0and_bath_real
  !
  !ANOMALOUS
  interface f0and_bath_real
     module procedure f0and_bath_real_main
     module procedure f0and_bath_real_ispin_jspin
     module procedure f0and_bath_real_ispin_jspin_iorb_jorb
     module procedure f0and_bath_real_main_
     module procedure f0and_bath_real_ispin_jspin_
     module procedure f0and_bath_real_ispin_jspin_iorb_jorb_
  end interface f0and_bath_real



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_real
     module procedure invg0_bath_real_main
     module procedure invg0_bath_real_ispin_jspin
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb
     module procedure invg0_bath_real_main_
     module procedure invg0_bath_real_ispin_jspin_
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb_
  end interface invg0_bath_real
  !
  !ANOMALOUS
  interface invf0_bath_real
     module procedure invf0_bath_real_main
     module procedure invf0_bath_real_ispin_jspin
     module procedure invf0_bath_real_ispin_jspin_iorb_jorb
     module procedure invf0_bath_real_main_
     module procedure invf0_bath_real_ispin_jspin_
     module procedure invf0_bath_real_ispin_jspin_iorb_jorb_
  end interface invf0_bath_real




  public :: delta_bath_mats
  public :: fdelta_bath_mats
  public :: delta_bath_real
  public :: fdelta_bath_real
  public :: g0and_bath_mats
  public :: f0and_bath_mats
  public :: invg0_bath_mats
  public :: invf0_bath_mats
  public :: g0and_bath_real
  public :: f0and_bath_real
  public :: invg0_bath_real
  public :: invf0_bath_real




contains




  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
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
    integer                                             :: i,ih,k,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    real(8),dimension(Nspin,Nbath)                      :: hps
    real(8),dimension(Nspin,Nspin,Nbath)                :: wps
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    real(8),dimension(Nspin*Norb,Nspin*Norb)            :: V_k
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: invH_k
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: Delta_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)   :: invH_knn
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
          !\Delta_{aa}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{a}^{s`h}(k)/(iw_n - H_{a}^{h}(k))]
          do iorb=1,Norb
             ehel = dmft_bath_%e(1:Nspin,iorb,1:Nbath)
             whel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,iorb,1:Nbath),dmft_bath_%u(1:Nspin,iorb,1:Nbath))
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do i=1,L
                      do ih=1,Nspin
                         Delta(ispin,jspin,iorb,iorb,i) = Delta(ispin,jspin,iorb,iorb,i) + &
                              sum( whel(ispin,ih,:)*whel(jspin,ih,:)/(x(i) - ehel(ih,:)) )
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
          !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{s`h}(k)/(iw_n - H^{h}(k))]
          ehel  = dmft_bath_%e(1:Nspin,1,1:Nbath)
          wohel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,1:Norb,1:Nbath),dmft_bath_%u(1:Nspin,1:Norb,1:Nbath))
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do i=1,L
                         do ih=1,Nspin
                            Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                                 sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(x(i) - ehel(ih,:)) )
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    case ("replica")
       !
       select case(ed_mode)
       case default
          !
          Delta_so=zero
          invH_k=zero
          do i=1,L
             !
             do ibath=1,Nbath
                !
                V_k=0.0d0
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1) * Norb
                            jo = jorb + (jspin-1) * Norb
                            invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                            V_k(io,io)=dmft_bath_%v(ispin,iorb,ibath)
                         enddo
                      enddo
                   enddo
                enddo
                !
                invH_k(:,:,i) = zeye(Nspin*Norb) * x(i) - invH_k(:,:,i)
                call inv(invH_k(:,:,i))
                !
                Delta_so(:,:,i)=Delta_so(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
                !
             enddo
             !
             Delta(:,:,:,:,i)=so2nn_reshape(Delta_so(:,:,i),Nspin,Norb)
             !
          enddo
          !
       case ("superc")
          !

          !
       case ("nonsu2")
          !
          Delta_so=zero
          invH_k=zero
          do i=1,L
             !VERSIONE 1 ===>
             !do ibath=1,Nbath
             !   !
             !   V_k=0.0d0
             !   do ispin=1,Nspin
             !      do jspin=1,Nspin
             !         do iorb=1,Norb
             !            do jorb=1,Norb
             !               io = iorb + (ispin-1) * Norb
             !               jo = jorb + (jspin-1) * Norb
             !               invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
             !               V_k(io,io)=dmft_bath_%v(ispin,iorb,ibath)
             !            enddo
             !         enddo
             !      enddo
             !   enddo
             !   !
             !   invH_k(:,:,i) = zeye(Nspin*Norb) * x(i) - invH_k(:,:,i)
             !   call inv(invH_k(:,:,i))
             !   !
             !   Delta_so(:,:,i)=Delta_so(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
             !   !
             !enddo
             !!
             !Delta(:,:,:,:,i)=so2nn_reshape(Delta_so(:,:,i),Nspin,Norb)
             !===> VERSIONE 1

             !VERSIONE 2 ===>
             invH_knn=zero
             do ibath=1,Nbath
                !
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1) * Norb
                            jo = jorb + (jspin-1) * Norb
                            invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                         enddo
                      enddo
                   enddo
                enddo
                !
                invH_k(:,:,i) = zeye(Nspin*Norb) * x(i) - invH_k(:,:,i)
                call inv(invH_k(:,:,i))
                invH_knn(:,:,:,:,ibath)=so2nn_reshape(invH_k(:,:,i),Nspin,Norb)
                !
             enddo
             !
             do ibath=1,Nbath
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            Delta(ispin,jspin,iorb,jorb,i)=Delta(ispin,jspin,iorb,jorb,i)+ &
                            dmft_bath_%v(ispin,iorb,ibath)*invH_knn(ispin,jspin,iorb,jorb,ibath)*dmft_bath_%v(jspin,jorb,ibath)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             !===> VERSIONE 2

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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Delta = delta_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Fdelta = fdelta_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = fdelta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = fdelta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function fdelta_bath_mats_ispin_jspin_iorb_jorb_



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta and Fdelta functions on the Real axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  function delta_bath_real_main(x,dmft_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                             :: i,ih,k,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    real(8),dimension(Nspin,Nbath)                      :: hps
    real(8),dimension(Nspin,Nspin,Nbath)                :: wps
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    real(8),dimension(Nspin*Norb,Nspin*Norb)            :: V_k
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: invH_k
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: Delta_so
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          !
          !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(w+i\h - E_{a}(k)) ]
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
          !\Delta_{aa}^{ss} = - \sum_k [ V_{a}(k) * V_{a}(k) * (w+i\h + E_{a}(k)) / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
          do ispin=1,Nspin
             do iorb=1,Norb
                eps = dmft_bath_%e(ispin,iorb,1:Nbath)
                dps = dmft_bath_%d(ispin,iorb,1:Nbath)
                vps = dmft_bath_%v(ispin,iorb,1:Nbath)
                do i=1,L
                   Delta(ispin,ispin,iorb,iorb,i) = -sum( vps(:)*vps(:)*(x(i) + eps(:))/(x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
                enddo
             enddo
          enddo
          !
       case ("nonsu2")
          !
          !\Delta_{aa}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{a}^{s`h}(k)/(w+i\h - H_{a}^{h}(k))]
          do iorb=1,Norb
             ehel = dmft_bath_%e(1:Nspin,iorb,1:Nbath)
             whel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,iorb,1:Nbath),dmft_bath_%u(1:Nspin,iorb,1:Nbath))
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do i=1,L
                      do ih=1,Nspin
                         Delta(ispin,jspin,iorb,iorb,i) = Delta(ispin,jspin,iorb,iorb,i) + &
                              sum( whel(ispin,ih,:)*whel(jspin,ih,:)/(x(i) - ehel(ih,:)) )
                      enddo
                   enddo
                enddo
             enddo
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
                   do i=1,L
                      Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          !
          !\Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (w+i\h + E(k)) / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
          do ispin=1,Nspin
             eps  = dmft_bath_%e(ispin,1      ,1:Nbath)
             dps  = dmft_bath_%d(ispin,1      ,1:Nbath)
             vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
             do iorb=1,Norb
                do jorb=1,Norb
                   do i=1,L
                      Delta(ispin,ispin,iorb,jorb,i) = -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/((x(i)*(-x(i)) + eps(:)**2 + dps(:)**2)) )
                   enddo
                enddo
             enddo
          enddo
       case ("nonsu2")
          !
          !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{s`h}(k)/(w+i\h - H^{h}(k))]
          ehel  = dmft_bath_%e(1:Nspin,1,1:Nbath)
          wohel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,1:Norb,1:Nbath),dmft_bath_%u(1:Nspin,1:Norb,1:Nbath))
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do i=1,L
                         do ih=1,Nspin
                            Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                                 sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(x(i) - ehel(ih,:)) )
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    case ("replica")
       !
       select case(ed_mode)
       case default
          !
          Delta_so=zero
          invH_k=zero
          do i=1,L
             !
             do ibath=1,Nbath
                !
                V_k=0.0d0
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1) * Norb
                            jo = jorb + (jspin-1) * Norb
                            invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                            V_k(io,io)=dmft_bath_%v(ispin,iorb,ibath)
                         enddo
                      enddo
                   enddo
                enddo
                !
                invH_k(:,:,i) = eye(Nspin*Norb) * x(i) - invH_k(:,:,i)
                call inv(invH_k(:,:,i))
                !
                Delta_so(:,:,i)=Delta_so(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
                !
             enddo
             !
             Delta(:,:,:,:,i)=so2nn_reshape(Delta_so(:,:,i),Nspin,Norb)
             !
          enddo
          !
       case ("superc")
          !

          !
       case ("nonsu2")
          !
          Delta_so=zero
          invH_k=zero
          do i=1,L
             !
             do ibath=1,Nbath
                !
                V_k=0.0d0
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            io = iorb + (ispin-1) * Norb
                            jo = jorb + (jspin-1) * Norb
                            invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                            V_k(io,io)=dmft_bath_%v(ispin,iorb,ibath)
                         enddo
                      enddo
                   enddo
                enddo
                !
                invH_k(:,:,i) = eye(Nspin*Norb) * x(i) - invH_k(:,:,i)
                call inv(invH_k(:,:,i))
                !
                Delta_so(:,:,i)=Delta_so(:,:,i)+matmul(V_k,matmul(invH_k(:,:,i),V_k))
                !
             enddo
             !
             Delta(:,:,:,:,i)=so2nn_reshape(Delta_so(:,:,i),Nspin,Norb)
             !
          enddo
          !
       end select
       !
    end select
  end function delta_bath_real_main


  function delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_real_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,:,:,:)
  end function delta_bath_real_ispin_jspin


  function delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(size(x))                       :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_real_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,iorb,jorb,:)
  end function delta_bath_real_ispin_jspin_iorb_jorb


  function delta_bath_real_main_(x,bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Delta = delta_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_main_


  function delta_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_ispin_jspin_

  function delta_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8),dimension(size(x))      :: G0out
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_ispin_jspin_iorb_jorb_



















  !ANOMALOUS:
  function fdelta_bath_real_main(x,dmft_bath_) result(Fdelta)
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
       select case(ed_mode)
       case default
          stop "Fdelta_bath_real error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          !
          !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
          do ispin=1,Nspin
             do iorb=1,Norb
                eps = dmft_bath_%e(ispin,iorb,1:Nbath)
                dps = dmft_bath_%d(ispin,iorb,1:Nbath)
                vps = dmft_bath_%v(ispin,iorb,1:Nbath)
                do i=1,L
                   Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/( x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
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
          stop "Fdelta_bath_real error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          !
          !\FDelta_{ab} = - \sum_k [ \Delta(k) * V_{a}(k) * V_{b}(k) / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
          do ispin=1,Nspin
             eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
             dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
             vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
             do iorb=1,Norb
                do jorb=1,Norb
                   do i=1,L
                      Fdelta(ispin,ispin,iorb,jorb,i) = -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end function fdelta_bath_real_main

  function fdelta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: F0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
    Fdelta = fdelta_bath_real_main(x,dmft_bath_)
    F0out  = Fdelta(ispin,jspin,:,:,:)
  end function fdelta_bath_real_ispin_jspin

  function fdelta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(size(x))                       :: F0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
    Fdelta = fdelta_bath_real_main(x,dmft_bath_)
    F0out = Fdelta(ispin,jspin,iorb,jorb,:)
  end function fdelta_bath_real_ispin_jspin_iorb_jorb

  function fdelta_bath_real_main_(x,bath_) result(Fdelta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "fdelta_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Fdelta = fdelta_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function fdelta_bath_real_main_

  function fdelta_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: F0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "fdelta_bath_real_ispin_jspin_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = fdelta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function fdelta_bath_real_ispin_jspin_

  function fdelta_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8),dimension(size(x))      :: F0out
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "fdelta_bath_real_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = fdelta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function fdelta_bath_real_ispin_jspin_iorb_jorb_






  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
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
             zeta  = (x(i) + xmu)*zeye(Nspin)
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
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
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
             zeta  = (x(i) + xmu)*zeye(Nso)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = g0and_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0and = f0and_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = f0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
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
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = f0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function f0and_bath_mats_ispin_jspin_iorb_jorb_



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 and F0 non-interacting Green's functions on the real-axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function g0and_bath_real_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    complex(8),dimension(size(x))                       :: det,fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    !
    G0and = zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       select case(ed_mode)
       case default
          !
          Delta = delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                fg(:)    = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
                G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
             enddo
          enddo
          !
       case ("superc")
          !
          Delta  =  delta_bath_real(x,dmft_bath_)
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                fg(:)  =  dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                ff(:)  =                                              - Fdelta(ispin,ispin,iorb,iorb,:)
                det(:) = -fg(:)*conjg(fg(L:1:-1)) - ff(:)*ff(:)
                G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(L:1:-1))/det(:)
             enddo
          enddo
          !
       case ("nonsu2")
          !
          Delta = delta_bath_real(x,dmft_bath_)
          allocate(fgorb(Nspin,Nspin),zeta(Nspin,Nspin))
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nspin)
             fgorb = zero
             !
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
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       !
       select case(ed_mode)
       case default
          !
          allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
          Delta = delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
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
          Delta  =  delta_bath_real(x,dmft_bath_)
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = zero
                fgorb= zero
                do iorb=1,Norb
                   zeta(iorb,iorb)           =        x(i)     + xmu
                   zeta(iorb+Norb,iorb+Norb) = -conjg(x(L-i+1) + xmu)
                enddo
                do iorb=1,Norb
                   do jorb=1,Norb
                      fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
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
          Delta = delta_bath_real(x,dmft_bath_)
          allocate(fgorb(Nso,Nso),zeta(Nso,Nso))
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nso)
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
  end function g0and_bath_real_main


  function g0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_real_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,:,:,:)
  end function g0and_bath_real_ispin_jspin


  function g0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_real_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function g0and_bath_real_ispin_jspin_iorb_jorb


  function g0and_bath_real_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = g0and_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_main_


  function g0and_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_ispin_jspin_


  function g0and_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_ispin_jspin_iorb_jorb_










  !ANOMALous:
  function f0and_bath_real_main(x,dmft_bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,i,L
    complex(8),dimension(size(x))                       :: det,fg,ff
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
          stop "F0and_bath_real error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          Delta  = delta_bath_real(x,dmft_bath_)
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                fg(:)  = dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                ff(:)  =                                                    - Fdelta(ispin,ispin,iorb,iorb,:)
                det(:) = fg(:)*conjg(fg(L:1:-1)) + ff(:)*ff(:)
                F0and(ispin,ispin,iorb,iorb,:) = ff(:)/det(:)
             enddo
          enddo
       end select
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       select case(ed_mode)
       case default
          stop "F0and_bath_real error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
          Delta  = delta_bath_real( x,dmft_bath_)
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = zero
                fgorb= zero
                do iorb=1,Norb
                   zeta(iorb,iorb)           =        x(i)     + xmu
                   zeta(iorb+Norb,iorb+Norb) = -conjg(x(L-i+1) + xmu)
                enddo
                do iorb=1,Norb
                   do jorb=1,Norb
                      fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                      fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
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
  end function f0and_bath_real_main

  function f0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: F0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = F0and_bath_real_main(x,dmft_bath_)
    F0out = F0and(ispin,jspin,:,:,:)
  end function f0and_bath_real_ispin_jspin

  function f0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: F0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = f0and_bath_real_main(x,dmft_bath_)
    F0out = F0and(ispin,jspin,iorb,jorb,:)
  end function f0and_bath_real_ispin_jspin_iorb_jorb

  function f0and_bath_real_main_(x,bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "f0and_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0and = f0and_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function f0and_bath_real_main_

  function f0and_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: F0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "f0and_bath_real_ispin_jspin_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = f0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function f0and_bath_real_ispin_jspin_

  function f0and_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: F0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "f0and_bath_real_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = f0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function f0and_bath_real_ispin_jspin_iorb_jorb_



  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the inverse G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0^{-1} and F0{-1} non-interacting Green's functions on the Matsubara axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function invg0_bath_mats_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
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
                G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
          !
       case ("superc")
          !
          Delta =  delta_bath_mats(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                G0and(ispin,ispin,iorb,iorb,:)  =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
          !
       case ("nonsu2")
          !
          Delta = delta_bath_mats(x,dmft_bath_)
          allocate(zeta(Nspin,Nspin))
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nspin)
             do iorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      G0and(ispin,jspin,iorb,iorb,i) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       end select
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       select case(ed_mode)
       case default
          !
          allocate(zeta(Norb,Norb))
          Delta = delta_bath_mats(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = (x(i)+xmu)*eye(Norb)
                do iorb=1,Norb
                   do jorb=1,Norb
                      G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-Delta(ispin,ispin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       case ("superc")
          !
          allocate(zeta(2*Norb,2*Norb))
          Delta  = delta_bath_mats(x,dmft_bath_)
          Fdelta = fdelta_bath_mats(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = zero
                do iorb=1,Norb
                   zeta(iorb,iorb)           = x(i) + xmu
                   zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
                enddo
                do iorb=1,Norb
                   do jorb=1,Norb
                      G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb) - Delta(ispin,ispin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       case ("nonsu2")
          !
          Nso=Nspin*Norb
          allocate(zeta(Nso,Nso))
          Delta = delta_bath_mats(x,dmft_bath_)
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nso)
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         G0and(ispin,jspin,iorb,jorb,i) = zeta(io,jo) -impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                      enddo
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       end select
       !
    end select
    !
  end function invg0_bath_mats_main


  function invg0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,:,:,:)
  end function invg0_bath_mats_ispin_jspin


  function invg0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function invg0_bath_mats_ispin_jspin_iorb_jorb


  function invg0_bath_mats_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_main_

  function invg0_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = invg0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_ispin_jspin_

  function invg0_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = invg0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_ispin_jspin_iorb_jorb_










  !ANOMALous:
  function invf0_bath_mats_main(x,dmft_bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,i,L
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
          !
          stop "Invf0_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          !
          Fdelta= fdelta_bath_mats(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                F0and(ispin,ispin,iorb,iorb,:) = -Fdelta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
       end select
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       select case(ed_mode)
       case default
          !
          stop "Invf0_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          !
          Fdelta= fdelta_bath_mats(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   F0and(ispin,ispin,iorb,jorb,:) = -Fdelta(ispin,ispin,iorb,jorb,:)
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end function invf0_bath_mats_main

  function invf0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: F0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = Invf0_bath_mats_main(x,dmft_bath_)
    F0out = F0and(ispin,jspin,:,:,:)
  end function invf0_bath_mats_ispin_jspin

  function invf0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: F0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = invf0_bath_mats_main(x,dmft_bath_)
    F0out = F0and(ispin,jspin,iorb,jorb,:)
  end function invf0_bath_mats_ispin_jspin_iorb_jorb

  function invf0_bath_mats_main_(x,bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0and = invf0_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_mats_main_

  function invf0_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: F0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_mats_ispin_jspin_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = invf0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_mats_ispin_jspin_

  function invf0_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: F0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_mats_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = invf0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_mats_ispin_jspin_iorb_jorb_




  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 and F0 non-interacting Green's functions on the real-axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function invg0_bath_real_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    complex(8),dimension(size(x))                       :: det,fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    !
    G0and = zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       select case(ed_mode)
       case default
          !
          Delta = delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                G0and(ispin,ispin,iorb,iorb,:) =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
          !
       case ("superc")
          !
          Delta  =  delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                G0and(ispin,ispin,iorb,iorb,:) = dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
          !
       case ("nonsu2")
          !
          Delta = delta_bath_real(x,dmft_bath_)
          allocate(zeta(Nspin,Nspin))
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nspin)
             do iorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      G0and(ispin,jspin,iorb,iorb,i) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       end select
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       select case(ed_mode)
       case default
          !
          allocate(zeta(Norb,Norb))
          Delta = delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = (x(i)+xmu)*eye(Norb)
                do iorb=1,Norb
                   do jorb=1,Norb
                      G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb) - Delta(ispin,ispin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       case ("superc")
          !
          allocate(zeta(Norb,Norb))
          Delta = delta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                zeta = (dreal(x(i))  + xmu)*eye(Norb)
                do iorb=1,Norb
                   do jorb=1,Norb
                      G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       case ("nonsu2")
          !
          Nso=Nspin*Norb
          allocate(zeta(Nso,Nso))
          Delta = delta_bath_real(x,dmft_bath_)
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nso)
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         G0and(ispin,jspin,iorb,jorb,i) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                      enddo
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
          !
       end select
       !
    end select
    !
  end function invg0_bath_real_main


  function invg0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_real_main(x,dmft_bath_)
    G0out = zero
    G0out = G0and(ispin,jspin,:,:,:)
  end function invg0_bath_real_ispin_jspin

  function invg0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_real_main(x,dmft_bath_)
    G0out = zero
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function invg0_bath_real_ispin_jspin_iorb_jorb


  function invg0_bath_real_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = invg0_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_main_


  function invg0_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = zero
    G0out = invg0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_ispin_jspin_

  function invg0_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = zero
    G0out = invg0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_ispin_jspin_iorb_jorb_










  !ANOMALous:
  function invf0_bath_real_main(x,dmft_bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin
    !
    F0and=zero
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          stop "Invf0_bath_real error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                F0and(ispin,ispin,iorb,iorb,:) = -Fdelta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
       end select
       !
       !
    case ("hybrid")
       select case(ed_mode)
       case default
          stop "Invf0_bath_real error: called with ed_mode=normal, bath_type=hybrid"
          !
       case ("superc")
          Fdelta = fdelta_bath_real(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   F0and(ispin,ispin,iorb,jorb,:) = -Fdelta(ispin,ispin,iorb,jorb,:)
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end function invf0_bath_real_main

  function invf0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: F0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = Invf0_bath_real_main(x,dmft_bath_)
    F0out = zero
    F0out = F0and(ispin,jspin,:,:,:)
  end function invf0_bath_real_ispin_jspin

  function invf0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(F0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: F0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    F0and = invf0_bath_real_main(x,dmft_bath_)
    F0out = zero
    F0out = F0and(ispin,jspin,iorb,jorb,:)
  end function invf0_bath_real_ispin_jspin_iorb_jorb

  function invf0_bath_real_main_(x,bath_) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0and = invf0_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_real_main_

  function invf0_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(F0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: F0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    integer                                 :: iorb,jorb
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_real_ispin_jspin_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = zero
    F0out = invf0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_real_ispin_jspin_

  function invf0_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(F0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: F0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invf0_bath_real_ispin_jspin_iorb_jorb_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    F0out = zero
    F0out = invf0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invf0_bath_real_ispin_jspin_iorb_jorb_



END MODULE ED_BATH_FUNCTIONS
