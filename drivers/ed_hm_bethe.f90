!include "MIXING.f90"
!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  character(len=16)      :: finput
  integer                :: M
  real(8)                :: alpha,K,wmixing,Eout(2)
  real(8),allocatable    :: He(:),Wte(:)
  real(8),allocatable    :: Gtau(:)

#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  ! call parse_input_variable(M,"M",finput,default=0)
  ! call parse_input_variable(alpha,"ALPHA",finput,default=1.d0)
  !
  call ed_read_input(trim(finput))
  !Hloc=zero

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe()

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(delta,bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     !Check convergence (if required change chemical potential)

     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu,converged)
     call end_loop
  enddo

  ! allocate(Wte(500),He(500))
  ! call bethe_lattice(Wte,He,500,1d0)
  ! Eout = ed_kinetic_energy(one*He,Wte,impSmats(1,1,1,1,:))

contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                     :: i,j,iorb
    complex(8)                  :: iw,zita,g0loc
    complex(8),dimension(Lmats) :: gloc,sigma,Tiw
    complex(8),dimension(Lreal) :: grloc
    real(8)                     :: wm(Lmats),wr(Lreal),tau(0:Lmats),C0,C1,n0
    real(8),dimension(0:Lmats)  :: sigt,gtau,Ttau
    real(8),dimension(3)  :: Scoeff

    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)

    do iorb=1,Norb
       do i=1,Lmats
          iw = xi*wm(i)
          zita    = iw + xmu - impSmats(1,1,iorb,iorb,i)
          gloc(i) = gfbethe(wm(i),zita,Wband)
          if(cg_scheme=='weiss')then
             delta(iorb,iorb,i)= one/(one/gloc(i) + impSmats(1,1,iorb,iorb,i))
          else
             delta(iorb,iorb,i)= iw + xmu - impSmats(1,1,iorb,iorb,i) - one/gloc(i)
          endif
       enddo

       do i=1,Lreal
          iw=cmplx(wr(i),eps)
          zita     = iw + xmu - impSreal(1,1,iorb,iorb,i)
          grloc(i) = gfbether(wr(i),zita,Wband)
       enddo
       if(ED_MPI_ID==0)then
          call splot("Gloc_"//reg(txtfy(iorb))//"_iw.ed",wm,gloc)
          call splot("Gloc_"//reg(txtfy(iorb))//"_realw.ed",wr,grloc)
          call splot("DOS"//reg(txtfy(iorb))//".ed",wr,-dimag(grloc)/pi)
          call splot("Delta_"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,iorb,:))
       endif

       ! tau(0:) = linspace(0.d0,beta,Lmats+1)

       ! C0=Uloc(1)*(ed_dens_up(1)-0.5d0)
       ! C1=Uloc(1)**2*ed_dens_up(1)*(1.d0-ed_dens_dw(1))
       ! Tiw=dcmplx(C0,-C1/wm)
       ! call splot("Tail_"//reg(txtfy(iorb))//"_iw.ed",wm,Tiw)

       ! Ttau = -C1/2.d0
       ! Sigma = impSmats(1,1,iorb,iorb,:)  - Tiw
       ! call fftgf_iw2tau(Sigma,Sigt(0:),beta,notail=.true.)
       ! Sigt=Sigt + Ttau
       ! call splot("Sigma_"//reg(txtfy(iorb))//"_tau.ed",tau,sigt)

       ! Sigt=Sigt !+ Ttau
       ! call fftgf_tau2iw(sigt(0:),sigma,beta)
       ! Sigma=Sigma !+ Tiw
       ! call splot("Sigma_"//reg(txtfy(iorb))//"_iw.ed",wm,sigma)


       ! call fftgf_iw2tau(gloc,gtau(0:),beta)
       ! call splot("Gloc_tau.ed",tau(0:),gtau(0:))
    enddo

  end subroutine get_delta_bethe
  !+----------------------------------------+







  ! !<DEBUG:
  ! ! subroutine get_ed_energy(Lk) 
  ! !   integer               :: Lk
  ! !   real(8),dimension(Lk) :: ek
  ! !   real(8)               :: de
  ! !   real(8),dimension(Lk) :: Wtk
  ! !   ek  = linspace(-Wband,Wband,Lk,mesh=de)
  ! !   Wtk = dens_bethe(ek,wband)*de
  ! !   call ed_kinetic_energy(impSmats(1,1,1,1,:),ek,wtk)
  ! ! end subroutine get_ed_energy


  ! function get_energy(Lk) result(H0)
  !   integer                     :: Lk
  !   complex(8),dimension(Lk)    :: Hk
  !   complex(8),dimension(Lmats) :: Sigma
  !   real(8)                     :: H0
  !   real(8),dimension(Lk)       :: Wtk
  !   real(8)                     :: Tail0,Tail1
  !   real(8)                     :: Sigma_HF,Ak,Bk
  !   complex(8)                  :: Ck,Dk,Zk
  !   complex(8)                  :: Zeta,Gk,Tk
  !   integer                     :: i,ik,iorb
  !   real(8),dimension(Lmats)    :: wm
  !   real(8)                     :: de
  !   !
  !   wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !   !
  !   Hk  = one*linspace(-Wband,Wband,Lk,mesh=de)
  !   Wtk = dens_bethe(dreal(Hk(:)),wband)*de
  !   Sigma = impSmats(1,1,1,1,:) 

  !   Sigma_HF = dreal(Sigma(Lmats))
  !   !
  !   H0=0.d0
  !   do ik=1,Lk
  !      Ak = Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      do i=1,Lmats
  !         Gk = one/(xi*wm(i) + xmu - Hk(ik) - Sigma(i) )
  !         Tk = one/(xi*wm(i)) - Bk/(xi*wm(i))**2
  !         Ck = Ak*(Gk - Tk)
  !         H0 = H0 + Ck*Wtk(ik)
  !      enddo
  !   enddo
  !   H0=H0/beta*4d0
  !   !
  !   Tail0=zero
  !   Tail1=zero
  !   do ik=1,Lk
  !      Ak= Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      Ck= Ak*Bk
  !      Tail0 = Tail0 + 0.5d0*Ak*Wtk(ik)
  !      Tail1 = Tail1 + 0.25d0*Ck*Wtk(ik)*beta
  !   enddo
  !   Tail0=2d0*Tail0
  !   Tail1=2d0*Tail1

  !   H0 = H0 + Tail0 + Tail1
  ! end function get_energy
  ! !>DEBUG


end program lancED



