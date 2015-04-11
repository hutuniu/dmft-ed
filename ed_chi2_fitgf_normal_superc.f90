!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_superc(fg,bath_,ispin)
  complex(8),dimension(:,:,:,:)        :: fg ![2][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)   :: bath_
  integer                              :: ispin
  real(8),dimension(:),allocatable     :: array_bath
  integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
  real(8)                              :: chi
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  real(8)                              :: w
  character(len=20)                    :: suffix
  integer                              :: unit
  if(size(fg,1)/=2)stop"chi2_fitgf_normal_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
  if(size(fg,3)/=Norb)stop"chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
  !
  allocate(Fdelta(2,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  select case(Cg_weight)
  case default
     Wdelta=dble(Ldelta)
  case(1)
     Wdelta=1.d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  !
  call allocate_bath(dmft_bath)
  call set_bath(bath_,dmft_bath)
  !
  call allocate_bath(chi2_bath)
  call set_bath(bath_,chi2_bath)
  !
  Asize = get_array_bath_dimension()
  allocate(array_bath(Asize))
  !
  do iorb=1,Norb
     Orb_indx=iorb
     Spin_indx=ispin
     !
     Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
     Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
     !
     call dmft_bath2array(dmft_bath,array_bath,ispin,iorb)
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,&
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case ("delta")
           call fmin_cg(array_bath,&
                chi2_delta_normal_superc,&
                grad_chi2_delta_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol  ,&
                istop=cg_stop ,&
                eps=cg_eps)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,&
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgminimize(array_bath,&
                chi2_delta_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     case (2)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgplus(array_bath,&
                chi2_weiss_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case ("delta")
           call fmin_cgplus(array_bath,&
                chi2_delta_normal_superc,&
                grad_chi2_delta_normal_superc,&
                iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol)
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)")&
          'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
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
     call array2dmft_bath(array_bath,dmft_bath,ispin,iorb)
     !
  enddo
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  unit=free_unit()
  open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  call write_bath(dmft_bath,unit)
  close(unit)
  if(ed_verbose<3)call write_fit_result(ispin)
  call copy_bath(dmft_bath,bath_)
  call deallocate_bath(dmft_bath)
  deallocate(Fdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    complex(8)        :: fgand(2),det
    integer           :: i,j,iorb,ispin
    real(8)           :: w
    do iorb=1,Norb
       suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
       Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
       fgand=zero
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          w = Xdelta(i)
          if(cg_scheme=='weiss')then
             fgand(1) = weiss_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
             fgand(2) = fweiss_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          else
             fgand(1) = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
             fgand(2) = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,dmft_bath)
          endif
          write(unit,"(10F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand(1)),dreal(Fdelta(1,i)),dreal(fgand(1)), &
               dimag(Fdelta(2,i)),dimag(fgand(2)),dreal(Fdelta(2,i)),dreal(fgand(2))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_superc



!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
! The gradient is not evaluated.
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_delta_normal_superc(a) result(chi2_sc)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0
  real(8)                        ::  chi2_sc,w
  integer                        ::  i,iorb
  chi2_sc = 0.d0 
  iorb=Orb_indx
  ispin=Spin_indx
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w = xdelta(i)
     g0(1,i) = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     g0(2,i) = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
  enddo
  !
  chi2_sc=sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:))
  !
end function chi2_delta_normal_superc
! the analytic GRADIENT of \chi^2
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_superc(a) result(dchi2_sc)
  real(8),dimension(:)                   ::  a
  complex(8),dimension(2,Ldelta)         ::  g0
  real(8),dimension(size(a))             ::  dchi2_sc,df
  complex(8),dimension(2,Ldelta,size(a)) ::  dg0
  integer                                ::  i,j,iorb
  real(8)                                ::  w
  dchi2_sc = 0.d0 
  df = 0.d0
  iorb=Orb_indx
  ispin=Spin_indx
  !push the array into a dmft_bath
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w      = Xdelta(i)
     g0(1,i)  = delta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     g0(2,i)  = fdelta_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     dg0(1,i) = grad_delta_bath_math(ispin,ispin,iorb,iorb,xi*2,chi2_bat,size(a))
     dg0(2,i) = grad_fdelta_bath_math(ispin,ispin,iorb,iorb,xi*2,chi2_bat,size(a))
  enddo
  !
  do j=1,size(a)
     df(j) = &
          sum( dreal(Fdelta(1,:)-g0(1,:))*dreal(dg0(1,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(1,:)-g0(1,:))*dimag(dg0(1,:,j))/Wdelta(:) ) + &
          sum( dreal(Fdelta(2,:)-g0(2,:))*dreal(dg0(2,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(2,:)-g0(2,:))*dimag(dg0(2,:,j))/Wdelta(:) )
  enddo
  !
  dchi2_sc=-2.d0*df
  !
end function grad_chi2_delta_normal_superc


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0 function 
!         in the SUPERCONDUCTING case.
! The gradient is not evaluated.
! NORMAL bath.
! SPIN & ORBITAL DIAGONAL
!+-------------------------------------------------------------+
function chi2_weiss_normal_superc(a) result(chi2_sc)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0
  real(8)                        ::  chi2_sc
  integer                        ::  i,iorb
  chi2_sc = 0.d0 
  iorb=Orb_indx
  ispin=Spin_indx
  !push the array into a dmft_bath
  call array2dmft_bath(a,chi2_bath,ispin,iorb)
  do i=1,Ldelta
     w   = Xdelta(i)
     g0(1,i)  = weiss_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
     g0(2,i)  = fweiss_bath_mats(ispin,ispin,iorb,iorb,xi*w,chi2_bath)
  enddo
  !
  chi2_sc = sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) 
  chi2_sc = chi2_sc + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:))
  !
end function chi2_weiss_normal_superc
