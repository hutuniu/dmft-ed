!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                :: iloop,Nb,Ne,ie
  logical                :: converged
  real(8),allocatable    :: wm(:),wr(:)
  real(8)                :: wband,ts,de
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  real(8),allocatable    :: epsik(:),wt(:)

  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)

  call parse_input_variable(ts,"TS","inputED.in",default=1.d0)
  call parse_input_variable(Ne,"NE","inputED.in",default=2000)
  call ed_read_input("inputED.in")

  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr = linspace(wini,wfin,Lreal)

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))

  allocate(wt(Ne),epsik(Ne))
  wband=4.d0*ts
  epsik = linspace(-wband,wband,Ne,mesh=de)
  do ie=1,Ne
     wt(ie)=dens_2dsquare(epsik(ie),ts)
  enddo
  if(ED_MPI_ID==0)call splot("DOS2d.ed",epsik,wt)
  wt = wt*de

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call ed_init_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     !if(nread/=0.d0)call search_mu(nimp(1),converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     if(ED_MPI_ID==0)call end_loop
  enddo


contains

  !+----------------------------------------+
  subroutine get_delta
    integer                   :: i,j,ie
    complex(8)                :: iw,zita,g0and,g0loc,gg
    complex(8),dimension(Lmats)  :: self,gloc
    complex(8),dimension(Lreal)  :: selfr,grloc

    do i=1,Lmats
       iw = xi*wm(i)
       zita    = iw + xmu - impSmats(1,1,1,1,i)
       gloc(i) = zero
       do ie=1,Ne
          gloc(i)=gloc(i)+wt(ie)/(zita-epsik(ie))
       enddo
       delta(1,1,i)= iw + xmu - impSmats(1,1,1,1,i) - one/gloc(i)
    enddo

    do i=1,Lreal
       iw=cmplx(wr(i),eps)
       zita     = iw + xmu - impSreal(1,1,1,1,i)
       grloc(i) = zero
       do ie=1,Ne
          grloc(i)=grloc(i)+wt(ie)/(zita-epsik(ie))
       enddo
    enddo
    if(ED_MPI_ID==0)then
       call splot("Gloc_iw.ed",wm,gloc)
       call splot("Gloc_realw.ed",wr,grloc)
       call splot("DOS.ed",wr,-dimag(grloc)/pi)
       call splot("Delta_iw.ed",wm,delta(1,1,:))
    endif
  end subroutine get_delta
  !+----------------------------------------+

end program lancED



