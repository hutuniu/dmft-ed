!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE CONSTANTS
  USE PARSE_INPUT
  USE IOTOOLS
  USE TOOLS
  USE FUNCTIONS
  USE ERROR
  USE ARRAYS
  implicit none
  integer                :: iloop,Nb(2),Le
  logical                :: converged

  !Bath:
  real(8),allocatable    :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:)
  character(len=16)      :: finput
  real(8)                :: wband
  real(8)                :: wmix
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(wmix,"wmix",finput,default=0.5d0)
  call parse_input_variable(Le,"Le",finput,default=1000)
  !
  call ed_read_input(trim(finput))

  !Allocate Weiss Field:
  allocate(delta(Nspin,Norb,Norb,Lmats))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  allocate(bath_(Nb(1),Nb(2)))
  call init_ed_solver(bath)
  call break_symmetry_bath(bath,sb_field,1.d0)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe()

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta(1,:,:,:),bath,ispin=1)
     call chi2_fitgf(delta(2,:,:,:),bath,ispin=2)

     if(iloop>1)bath = wmix*bath + (1.d0-wmix)*Bath_
     Bath_=bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu,converged)
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                       :: i,j,ie
    complex(8)                    :: iw,zita(2),zeta
    complex(8),dimension(2,Lmats) :: gloc
    complex(8),dimension(2,Lreal) :: grloc
    real(8)                       :: wm(Lmats),wr(Lreal)
    real(8)                       :: epsi(Le),dos(Le)
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    call bethe_lattice(dos,epsi,Le,wband)
    do i=1,Lmats
       iw      = xi*wm(i)
       zita(1) = iw + xmu - impSmats(1,1,1,1,i)
       zita(2) = iw + xmu - impSmats(2,2,1,1,i)
       zeta    = zita(1)*zita(2)
       gloc(:,i)    = zero
       do ie=1,Le
          gloc(1,i) = gloc(1,i) + dos(ie)/(zeta - epsi(ie)**2)
       enddo
       gloc(2,i) = zita(1)*gloc(1,i)
       gloc(1,i) = zita(2)*gloc(1,i)
       if(cg_scheme=='weiss')then
          delta(1,1,1,i)= one/(one/gloc(1,i) + impSmats(1,1,1,1,i))
          delta(2,1,1,i)= one/(one/gloc(2,i) + impSmats(2,2,1,1,i))
       else
          delta(1,1,1,i)= iw + xmu - impSmats(1,1,1,1,i) - one/gloc(1,i)
          delta(2,1,1,i)= iw + xmu - impSmats(2,2,1,1,i) - one/gloc(2,i)
       endif
    enddo

    do i=1,Lreal
       iw=cmplx(wr(i),eps)
       zita(1) = iw + xmu - impSreal(1,1,1,1,i)
       zita(2) = iw + xmu - impSreal(2,2,1,1,i)
       zeta    = zita(1)*zita(2)
       grloc(:,i)    = zero
       do ie=1,Le
          grloc(1,i) = grloc(1,i) + dos(ie)/(zeta - epsi(ie)**2)
       enddo
       grloc(2,i) = zita(1)*grloc(1,i)
       grloc(1,i) = zita(2)*grloc(1,i)
    enddo

    call splot("Gloc_iw.ed",wm,gloc(1,:),gloc(2,:))
    call splot("Gloc_realw.ed",wr,-dimag(grloc(1,:))/pi,dreal(grloc(1,:)),-dimag(grloc(2,:))/pi,dreal(grloc(2,:)))
    call splot("Delta_iw.ed",wm,delta(1,1,1,:),delta(2,1,1,:))

  end subroutine get_delta_bethe
  !+----------------------------------------+

end program lancED



