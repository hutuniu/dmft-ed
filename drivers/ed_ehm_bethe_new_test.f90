program ed_ehm
  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR

  implicit none
  integer                                       :: iloop
  logical                                       :: converged
  real(8)                                       :: wband,W0
  character(len=16)                             :: finput
  real(8)                                       :: wmixing
  integer                                       :: Nk
  real(8)                                       :: dens,Eout(2,2)
  !real(8)					:: xmu0
  !Bath:
  integer                                       :: Nb(2)
  real(8),allocatable                           :: Bath(:,:,:),Bath_prev(:,:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta ![Nlat][Nspin][Nspin][Norb][Norb][Nfreq]
  !the green functions
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_loc
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_loc
  !Hamiltonian
  real(8),allocatable                           :: BetheDOS(:)
  real(8),allocatable                           :: BetheEk(:,:,:)
  complex(8),allocatable                        :: Hloc(:,:,:,:,:)
  character(len=5)                              :: tmp_suffix

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Nk,"Nk",finput,default=500)
  call parse_input_variable(dens,"dens",finput,default=1.0d0)
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(W0,"W0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))

  !xmu0=xmu+W0

  Nlat=1                        !only one lattice
  if(Norb/=1)stop "This drivers requires Norb==1"
  if(Nspin/=1)stop "This drivers requires Nspin==1"

  allocate(BetheDOS(Nk),BetheEk(Nlat,Nlat,Nk))
  call build_Hbethe2x2
  !<DEBUG 
  call splot("testDOSbethe.dat",BetheEk(1,1,:),BetheDOS(:))
  !>DEBUG


  !Allocate Weiss (Delta) Field:
  allocate(Delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_loc(Nlat,Nspin,Nspin,Norb,Norb,Lmats))

  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_loc(Nlat,Nspin,Nspin,Norb,Norb,Lreal))

  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))  

  !setup solver
  Nb=get_bath_size()
  print*,"Nb = ",Nb(1)," ",Nb(2)
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_prev(Nlat,Nb(1),Nb(2)))
  call ed_init_solver(Bath(1,:,:))
  
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath(1,:,:)) 

     dens = ed_dens(1)
     print'(A6,F10.8,A6,F10.8)',"dens = ",dens," xmu = ",xmu

        call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(1,:,:,:,:,:))
        call ed_get_gimp_matsubara(Gmats(1,:,:,:,:,:))
        call ed_get_gimp_real(Greal(1,:,:,:,:,:))


        !xmu=xmu0-W0*dens
        call ed_get_gloc(one*BetheEk(1,1,:),BetheDOS,&
             Gmats_loc(1,1,1,1,1,:),&
             Greal_loc(1,1,1,1,1,:),&
             Smats(1,1,1,1,1,:),&
             Sreal(1,1,1,1,1,:),iprint=1)




     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe()

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(delta(1,1,1,:,:,:),bath(1,:,:),ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     !Check convergence (if required change chemical potential)

     converged = check_convergence(delta(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu,converged)
     call end_loop
  enddo







contains

  subroutine build_Hbethe2x2
    real(8),dimension(Nk) :: ek,wtk
    real(8)               :: de
    ek  = linspace(-Wband,Wband,Nk,mesh=de)
    Wtk = dens_bethe(ek,wband)*de
    BetheDOS        = zero
    BetheEk         = 0d0
    BetheDOS(:)     = Wtk(:)
    BetheEk(1,1,:)  = ek
 !   BetheEk(2,2,:)  = ek
  end subroutine build_Hbethe2x2

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
             delta(1,1,1,iorb,iorb,i)= one/(one/gloc(i) + impSmats(1,1,iorb,iorb,i))
          else
             delta(1,1,1,iorb,iorb,i)= iw + xmu - impSmats(1,1,iorb,iorb,i) - one/gloc(i)
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
          call splot("Delta_"//reg(txtfy(iorb))//"_iw.ed",wm,delta(1,1,1,iorb,iorb,:))
       endif
    enddo

  end subroutine get_delta_bethe
  !+----------------------------------------+
 

end program ed_ehm



