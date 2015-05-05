!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)        :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  integer                            :: ispin
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  complex(8)                         :: fgand
  real(8)                            :: w
  character(len=20)                  :: suffix
  integer                            :: unit
  !
  if(size(fg,1)/=Norb)stop"chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(Gdelta(1,Ldelta))
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
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
  !
  call allocate_bath(chi2_bath)
  call set_bath(bath_,chi2_bath)
  !
  Asize = get_chi2_bath_size()
  allocate(array_bath(Asize))
  !
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     !
     Gdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
     !
     call dmft_bath2chi2_bath(dmft_bath,array_bath,ispin,iorb)
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,&
                chi2_weiss_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,&
                chi2_delta_normal_normal,&
                grad_chi2_delta_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,&
                chi2_weiss_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,&
                chi2_delta_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(array_bath,&
                chi2_weiss_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,&
                chi2_delta_normal_normal,&
                grad_chi2_delta_normal_normal,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     !
     if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)")&
          "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
     !
     if(ed_verbose<2)then
        suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
        unit=free_unit()
        open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
        write(unit,"(ES18.9,1x,I5)") chi,iter
        close(unit)
     endif
     !
     call chi2_bath2dmft_bath(array_bath,dmft_bath,ispin,iorb)
     !
  enddo
  !
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  !
  ! unit=free_unit()
  ! open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  ! call write_bath(dmft_bath,unit)
  ! close(unit)
  call save_bath(dmft_bath)
  !
  if(ed_verbose<3)call write_fit_result(ispin)
  call copy_bath(dmft_bath,bath_)
  call deallocate_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Gdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
       fgand=zero
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          w = Xdelta(i)
          if(cg_scheme=='weiss')then
             fgand = g0and_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          else
             fgand = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          endif
          write(unit,"(5F24.15)")Xdelta(i),dimag(Gdelta(1,i)),dimag(fgand),dreal(Gdelta(1,i)),dreal(fgand)
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_normal







!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0
  real(8)                      ::  chi2,w
  integer                      ::  i,iorb,ispin
  type(effective_bath)         ::  dmft_bath
  iorb=Orb_indx
  ispin=Spin_indx
  call chi2_bath2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w = xdelta(i)
     g0(i)   = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
  enddo
  !
  chi2=sum(abs(Gdelta(1,:)-g0(:))**2/Wdelta(:))
  !
end function chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_normal(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(size(a))           :: df
  complex(8),dimension(Ldelta)         :: g0
  complex(8),dimension(Ldelta,size(a)) :: dg0
  integer                              :: i,j,iorb,ispin
  real(8)                              :: w
  df=0d0
  iorb=Orb_indx
  ispin=Spin_indx
  !push the array into a dmft_bath
  call chi2_bath2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w        = Xdelta(i)
     g0(i)    = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     dg0(i,:) = grad_delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath,size(a))
  enddo
  !
  do j=1,size(a)
     df(j)=sum( dreal(Gdelta(1,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
          sum(  dimag(Gdelta(1,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
  enddo
  !
  dchi2 = -2.d0*df
  !
end function grad_chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0
  real(8)                      ::  chi2,w
  integer                      ::  i,iorb,ispin
  chi2 = 0d0 
  iorb=Orb_indx
  ispin=Spin_indx
  !push the array into a dmft_bath
  call chi2_bath2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w      = Xdelta(i)
     g0(i)  = g0and_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
  enddo
  !
  chi2=sum(abs(Gdelta(1,:)-g0(:))**2/Wdelta(:))
  !
end function chi2_weiss_normal_normal

