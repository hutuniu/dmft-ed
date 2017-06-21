!##################################################################
! THE CALCULATION OF THE \chi^2 FUNCTIONS USE PROCEDURES FURTHER 
! BELOW TO EVALUATE INDEPENDENTLY THE ANDERSON MODEL:
!  - DELTA, 
!  -\GRAD DELTA
!  - G0
! THE LATTER ARE ADAPTED FROM THE PROCEDURES:
! DELTA_BATH_MATS
! GRAD_DELTA_BATH_MATS
! G0 BATH_MATS
! FOR, YOU NEED TO DECOMPOSE THE a INPUT ARRAY INTO ELEMENTS.
!##################################################################


!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath nonSU2 phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_nonsu2(fg,bath_)
  complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: iter,i,j,io,stride,iorb,ispin,jspin,cspin,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=20)                  :: suffix
  integer                            :: unit
  !
  if(size(fg,1)/=Nspin)stop "chi2_fitgf_normal_nonsu2 error: size[fg,1]!=Nspin"
  if(size(fg,2)/=Nspin)stop "chi2_fitgf_normal_nonsu2 error: size[fg,2]!=Nspin"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_nonsu2 error: size[fg,3]!=Norb"
  if(size(fg,4)/=Norb)stop "chi2_fitgf_normal_nonsu2 error: size[fg,4]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_nonsu2 error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
  !
  totNspin = Nspin*(Nspin+1)/2
  allocate(getIspin(totNspin),getJspin(totNspin))
  cspin=0
  do ispin=1,Nspin
     do jspin=ispin,Nspin
        cspin=cspin+1
        getIspin(cspin)=ispin
        getJspin(cspin)=jspin
     enddo
  enddo
  if(cspin/=totNspin)stop "chi2_fitgf_normal_nonsu2: error counting the spins"
  !
  allocate(Gdelta(totNspin,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  select case(cg_weight)
  case default
     Wdelta=1d0*Ldelta
  case(1)
     Wdelta=1d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  select case(ed_para)
  case (.true.)
     ! E_{:,\a}(:) [1][ 1 ][Nbath]
     ! V_{:,\a}(:) [1][ 1 ][Nbath]
     ! U_{:,\a}(:) [1][ 1 ][Nbath]
     Asize = Nbath + Nbath + Nbath
     allocate(array_bath(Asize))
  case (.false.)
     ! E_{:,\a}(:) [Nspin][ 1 ][Nbath]
     ! V_{:,\a}(:) [Nspin][ 1 ][Nbath]
     ! U_{:,\a}(:) [Nspin][ 1 ][Nbath]
     Asize = Nspin*Nbath + Nspin*Nbath + Nspin*Nbath
     allocate(array_bath(Asize))
  end select
  !
  do iorb=1,Norb
     Orb_indx=iorb
     !
     do i=1,totNspin
        Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),iorb,iorb,1:Ldelta)
     enddo
     !
     select case(ed_para)
     case (.true.)
        ! size = Nbath + Nbath + Nbath
        stride = 0
        ispin  = 1 
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           array_bath(io) = dmft_bath%e(ispin,iorb,i)
        enddo
        stride = Nbath
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           array_bath(io) = dmft_bath%v(ispin,iorb,i)
        enddo
        stride = Nbath + Nbath
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           array_bath(io) = dmft_bath%u(ispin,iorb,i)
        enddo
     case (.false.)
        ! Nspin*Nbath + Nspin*Nbath + Nspin*Nbath
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              array_bath(io) = dmft_bath%e(ispin,iorb,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              array_bath(io) = dmft_bath%v(ispin,iorb,i)
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              array_bath(io) = dmft_bath%u(ispin,iorb,i)
           enddo
        enddo
     end select
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_nonsu2 error: cg_scheme != [weiss,delta]"
        end select
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,chi2_weiss_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,chi2_delta_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_nonsu2 error: cg_scheme != [weiss,delta]"
        end select
        !
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(array_bath,chi2_weiss_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,chi2_delta_normal_nonsu2,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_nonsu2 error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     !
     write(LOGfile,"(A,ES18.9,A,I5,A)")&
          "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(iorb))//" All spins"
     !
     suffix="_orb"//reg(txtfy(iorb))//"_ALLspins_"//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
     !
     select case(ed_para)
     case (.true.)
        ! size  = Nbath + Nbath + Nbath
        stride = 0
        ispin  = 1
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           dmft_bath%e(ispin,iorb,i) = array_bath(io)
        enddo
        stride = Nbath
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           dmft_bath%v(ispin,iorb,i)  = array_bath(io)
        enddo
        stride = Nbath + Nbath
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           dmft_bath%u(ispin,iorb,i)  = array_bath(io)
        enddo
        dmft_bath%e(Nspin,iorb,:) = dmft_bath%e(1,iorb,:)
        dmft_bath%v(Nspin,iorb,:) = dmft_bath%v(1,iorb,:)
        dmft_bath%u(Nspin,iorb,:) = dmft_bath%u(1,iorb,:)
     case (.false.)
        ! size = Nspin*Nbath + Nspin*Nbath + Nspin*Nbath
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath%e(ispin,iorb,i) = array_bath(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath%v(ispin,iorb,i)  = array_bath(io)
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath%u(ispin,iorb,i)  = array_bath(io)
           enddo
        enddo
     end select
     !
  enddo
  !
  call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  call write_fit_result(ispin)
  !
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  deallocate(getIspin,getJspin)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(Nspin,Nspin,Norb,Norb,Ldelta)
    integer           :: i,j,s,iorb,ispin,jspin
    real(8)           :: w
    if(cg_scheme=='weiss')then
       fgand(:,:,:,:,:) = g0and_bath_mats(xi*Xdelta(:),dmft_bath)
    else
       fgand(:,:,:,:,:) = delta_bath_mats(xi*Xdelta(:),dmft_bath)
    endif
    !
    do iorb=1,Norb
       do s=1,totNspin
          ispin = getIspin(s)
          jspin = getJspin(s)
          suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
          unit=free_unit()
          if(cg_scheme=='weiss')then
             open(unit,file="fit_weiss"//reg(suffix)//".ed")
          else
             open(unit,file="fit_delta"//reg(suffix)//".ed")
          endif
          !
          do i=1,Ldelta
             w = Xdelta(i)
             write(unit,"(5F24.15)")Xdelta(i),&
                  dimag(fg(ispin,jspin,iorb,iorb,i)),dimag(fgand(ispin,jspin,iorb,iorb,i)),&
                  dreal(fg(ispin,jspin,iorb,iorb,i)),dreal(fgand(ispin,jspin,iorb,iorb,i))
          enddo
          close(unit)
       enddo
    enddo
  end subroutine write_fit_result
  !
end subroutine chi2_fitgf_normal_nonsu2






!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_normal_nonsu2(a) result(chi2)
  real(8),dimension(:)                     ::  a
  real(8)                                  ::  chi2
  real(8),dimension(totNspin)              ::  chi2_spin
  complex(8),dimension(Nspin,Nspin,Ldelta) ::  Delta
  integer                                  ::  s,ispin,jspin
  !
  Delta = delta_normal_nonsu2(a)
  !
  do s=1,totNspin
     ispin=getIspin(s)
     jspin=getJspin(s)
     chi2_spin(s) = sum(abs(Gdelta(s,:)-Delta(ispin,jspin,:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi2_spin)
  !
end function chi2_delta_normal_nonsu2


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_normal_nonsu2(a) result(chi2)
  real(8),dimension(:)                     :: a
  real(8),dimension(totNspin)              :: chi2_spin
  complex(8),dimension(Nspin,Nspin,Ldelta) :: g0and
  real(8)                                  :: chi2
  real(8)                                  :: w
  integer                                  :: i,s,ispin,jspin
  !
  g0and(:,:,:) = g0and_normal_nonsu2(a)
  !
  do s=1,totNspin
     ispin=getIspin(s)
     jspin=getJspin(s)
     chi2_spin(s) = sum(abs(Gdelta(s,:)-g0and(ispin,jspin,:))**2/Wdelta(:))
  enddo
  !
  chi2=sum(chi2_spin)
  !
end function chi2_weiss_normal_nonsu2






!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_normal_nonsu2(a) result(Delta)
  real(8),dimension(:)                     :: a
  complex(8),dimension(Nspin,Nspin,Ldelta) :: Delta
  integer                                  :: ispin,jspin
  integer                                  :: i,io,stride,ih
  real(8),dimension(Nspin,Nbath)           :: hps
  real(8),dimension(Nspin,Nbath)           :: vps
  real(8),dimension(Nspin,Nbath)           :: ups
  real(8),dimension(Nhel,Nbath)            :: ehel
  real(8),dimension(Nhel,Nhel,Nbath)       :: whel

  select case(ed_para)
  case (.true.)
     stride = 0
     ispin  = 1
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        hps(ispin,i) = a(io)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        vps(ispin,i)  = a(io)
     enddo
     stride = Nbath + Nbath
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        ups(ispin,i)  = a(io)
     enddo
     hps(Nspin,:)=hps(ispin,:)
     vps(Nspin,:)=vps(ispin,:)
     ups(Nspin,:)=ups(ispin,:)
  case (.false.)
     stride = 0
     do ispin=1,Nspin
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           hps(ispin,i) = a(io)
        enddo
     enddo
     stride = Nspin*Nbath
     do ispin=1,Nspin
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           vps(ispin,i)  = a(io)
        enddo
     enddo
     stride = Nspin*Nbath + Nspin*Nbath
     do ispin=1,Nspin
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           ups(ispin,i)  = a(io)
        enddo
     enddo
  end select
  !
  whel = get_Whyb_matrix(vps(1:Nspin,1:Nbath),ups(1:Nspin,1:Nbath))
  !
  Delta=zero
  do i=1,Ldelta
     do ispin=1,Nspin
        do jspin=1,Nspin
           do ih=1,Nspin
              Delta(ispin,jspin,i) = Delta(ispin,jspin,i) + sum( whel(ispin,ih,:)*whel(jspin,ih,:)/(xi*Xdelta(i) - hps(ih,:)) )
           enddo
        enddo
     enddo
  enddo
  !
end function delta_normal_nonsu2

function g0and_normal_nonsu2(a) result(G0and)
  integer                                  :: iorb,ispin,jspin
  real(8),dimension(:)                     :: a
  complex(8),dimension(Nspin,Nspin,Ldelta) :: G0and,Delta
  complex(8),dimension(Nspin,Nspin)        :: zeta,fgorb
  integer                                  :: i
  !
  iorb  = Orb_indx
  !
  Delta(:,:,:) = delta_normal_nonsu2(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta = (xi*Xdelta(i)+xmu)*zeye(Nspin)
     fgorb(:,:)   = zeta(:,:) - impHloc(:,:,iorb,iorb) - Delta(:,:,i)
     call inv(fgorb)
     G0and(:,:,i) = fgorb
  enddo
  !
end function g0and_normal_nonsu2
