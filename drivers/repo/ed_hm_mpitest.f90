program hm_mpitest
  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: Nb(2)
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  character(len=16)      :: finput
  real(8)                :: wmixing
  real(8),allocatable    :: Hk(:),wt(:)

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
  call ed_read_input(trim(finput))

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call ed_init_solver(bath)

  !DMFT loop
  if(ED_MPI_ID==0)call start_loop(1,1,"DMFT-loop")
  !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
  call ed_solve(bath) 
  !Get the Weiss field/Delta function to be fitted (user defined)
  call get_delta_bethe()
  !Perform the SELF-CONSISTENCY by fitting the new bath
  call ed_chi2_fitgf(delta,bath,ispin=1)
  if(ED_MPI_ID==0)call end_loop

contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                     :: i,j,iorb
    complex(8)                  :: iw,zita,g0loc
    complex(8),dimension(Lmats) :: gloc,sigma,Tiw
    complex(8),dimension(Lreal) :: grloc
    real(8)                     :: wm(Lmats),wr(Lreal),tau(0:Lmats),C0,C1,n0
    real(8),dimension(0:Lmats)  :: sigt,gtau,Ttau
    real(8),dimension(3)        :: Scoeff
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
    enddo
  end subroutine get_delta_bethe
  !+----------------------------------------+

end program hm_mpitest



