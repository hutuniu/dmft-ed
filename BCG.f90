MODULE BCG
  implicit none
  private

  interface bcg_solve
     module procedure bcg_solve_dd,bcg_solve_cc
  end interface bcg_solve
  public :: bcg_solve

contains

  !--------------------------------------------------------------------!
  !PURPOSE: solve a sparse linear system using biconjugate gradient 
  ! method. this routine is adapted from Num. Rec. F90 adding direct 
  ! interface to A*v and M*b=r with pre-conditioned M matrix operations. 
  !NOTE: routine is restricted to symmetric (possibly positive definite)
  ! matrices A, so that matrix-vector multiplication for A and A^T are 
  ! identical which for LL sparse matrices seems to involved.
  !--------------------------------------------------------------------!
  subroutine bcg_solve_dd(atimes,asolve,b,x,itol,tol,itmax,iter,err,iverbose)
    real(8), dimension(:), intent(in)          :: b
    real(8), dimension(size(b)), intent(inout) :: x
    integer, intent(in)                        :: itol,itmax
    real(8), intent(in)                        :: tol
    integer, intent(out)                       :: iter
    real(8), intent(out)                       :: err
    real(8), parameter                         :: eps=1.d-14
    integer                                    :: n
    real(8)                                    :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
    real(8), dimension(size(b))                :: p,pp,r,rr,z,zz
    logical :: iverbose
    interface 
       subroutine atimes(vin,vout)
         real(8),dimension(:)                  :: vin
         real(8),dimension(size(vin))          :: vout
       end subroutine atimes
       subroutine asolve(vin,vout)
         real(8),dimension(:)                  :: vin
         real(8),dimension(size(vin))          :: vout
       end subroutine asolve
    end interface
    n=size(b)
    iter=0
    call atimes(x,r)
    r=b-r
    rr=r
    select case(itol)
    case(1)
       bnrm=dnrm(b,itol)
       call asolve(r,z)
    case(2)
       call asolve(b,z)
       bnrm=dnrm(z,itol)
       call asolve(r,z)
    case(3:4)
       call asolve(b,z)
       bnrm=dnrm(z,itol)
       call asolve(r,z)
       znrm=dnrm(z,itol)
    case default
       stop 'illegal itol in linbcg'
    end select
    do
       if (iter > itmax) exit
       iter=iter+1
       call asolve(rr,zz)
       bknum=dot_product(z,rr)
       if (iter == 1) then
          p=z
          pp=zz
       else
          bk=bknum/bkden
          p=bk*p+z
          pp=bk*pp+zz
       end if
       bkden=bknum
       call atimes(p,z)
       akden=dot_product(z,pp)
       ak=bknum/akden
       call atimes(pp,zz)
       x=x+ak*p
       r=r-ak*z
       rr=rr-ak*zz
       call asolve(r,z)
       select case(itol)
       case(1)
          err=dnrm(r,itol)/bnrm
       case(2)
          err=dnrm(z,itol)/bnrm
       case(3:4)
          zm1nrm=znrm
          znrm=dnrm(z,itol)
          if (abs(zm1nrm-znrm) > EPS*znrm) then
             dxnrm=abs(ak)*dnrm(p,itol)
             err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
             err=znrm/bnrm
             cycle
          end if
          xnrm=dnrm(x,itol)
          if (err <= 0.5d0*xnrm) then
             err=err/xnrm
          else
             err=znrm/bnrm
             cycle
          end if
       end select
       if(iverbose)write (*,*) ' iter=',iter,' err=',err
       if (err <= tol) exit
    end do
  end subroutine bcg_solve_dd


  subroutine bcg_solve_cc(atimes,asolve,b,x,itol,tol,itmax,iter,err,iverbose)
    complex(8), dimension(:), intent(in)          :: b
    complex(8), dimension(size(b)),intent(inout) :: x
    integer, intent(in)                        :: itol,itmax
    real(8), intent(in)                        :: tol
    integer, intent(out)                       :: iter
    real(8), intent(out)                       :: err
    real(8), parameter                         :: eps=1.d-14
    integer                                    :: n
    complex(8)                                    :: ak,akden,bk,bkden,bknum
    real(8) :: bnrm,dxnrm,xnrm,zm1nrm,znrm
    complex(8), dimension(size(b))                :: p,pp,r,rr,z,zz
    logical :: iverbose
    interface 
       subroutine atimes(vin,vout)
         complex(8),dimension(:)                  :: vin
         complex(8),dimension(size(vin))          :: vout
       end subroutine atimes
       subroutine asolve(vin,vout)
         complex(8),dimension(:)                  :: vin
         complex(8),dimension(size(vin))          :: vout
       end subroutine asolve
    end interface
    n=size(b)
    iter=0
    call atimes(x,r)
    r=b-r
    rr=r
    select case(itol)
    case(1)
       bnrm=cnrm(b,itol)
       call asolve(r,z)
    case(2)
       call asolve(b,z)
       bnrm=cnrm(z,itol)
       call asolve(r,z)
    case(3:4)
       call asolve(b,z)
       bnrm=cnrm(z,itol)
       call asolve(r,z)
       znrm=cnrm(z,itol)
    case default
       stop 'illegal itol in linbcg'
    end select
    do
       if (iter > itmax) exit
       iter=iter+1
       call asolve(rr,zz)
       bknum=dot_product(z,rr)
       if (iter == 1) then
          p=z
          pp=zz
       else
          bk=bknum/bkden
          p=bk*p+z
          pp=bk*pp+zz
       end if
       bkden=bknum
       call atimes(p,z)
       akden=dot_product(z,pp)
       ak=bknum/akden
       call atimes(pp,zz)
       x=x+ak*p
       r=r-ak*z
       rr=rr-ak*zz
       call asolve(r,z)
       select case(itol)
       case(1)
          err=cnrm(r,itol)/bnrm
       case(2)
          err=cnrm(z,itol)/bnrm
       case(3:4)
          zm1nrm=znrm
          znrm=cnrm(z,itol)
          if (abs(zm1nrm-znrm) > EPS*znrm) then
             dxnrm=abs(ak)*cnrm(p,itol)
             err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
             err=znrm/bnrm
             cycle
          end if
          xnrm=cnrm(x,itol)
          if (err <= 0.5d0*xnrm) then
             err=err/xnrm
          else
             err=znrm/bnrm
             cycle
          end if
       end select
       if(iverbose)write (*,*) ' iter=',iter,' err=',err
       if (err <= tol) exit
    end do
  end subroutine bcg_solve_cc


  function dnrm(sx,itol)  result(snrm)
    real(8),dimension(:) :: sx
    integer              :: n,itol,i,isamax 
    real(8)              :: snrm
    n=size(sx)
    if (itol<=3)then
       snrm=0.d0
       do i=1,n
          snrm=snrm+sx(i)**2 
       enddo
       snrm=sqrt(snrm) 
    else
       isamax=1
       do i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i 
       enddo
       snrm=abs(sx(isamax)) 
    endif
    return 
  end function dnrm


  function cnrm(sx,itol) result(snrm)
    complex(8),dimension(:) :: sx
    integer              :: n,itol,i,isamax 
    real(8)              :: snrm
    n=size(sx)
    if (itol<=3)then
       snrm = sqrt(dot_product(sx,sx))
    else
       isamax=1
       do i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i 
       enddo
       snrm=abs(sx(isamax)) 
    endif
    return 
  end FUNCTION cnrm


END MODULE BCG
