program ed_cdw
  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR


  implicit none
  integer                                       :: iloop,ineq
  logical                                       :: converged
  real(8)                                       :: wband,W0
  character(len=16)                             :: finput
  real(8)                                       :: wmixing
  integer                                       :: Nk
  real(8)                                       :: n0(2),xmu0,Eout(2,2)
  !Bath:
  integer                                       :: Nb(2)
  real(8),allocatable                           :: Bath(:,:,:),Bath_prev(:,:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta ![Nlat][Nspin][Nspin][Norb][Norb][Nfreq]
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
  call parse_input_variable(n0,"n0",finput,default=[0.5d0,0.5d0])
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(W0,"W0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))

  xmu0=xmu+W0

  Nlat=2                        !the two ineq. sublattices A,B
  if(Norb/=1)stop "This drivers requires Norb==1"
  if(Nspin/=1)stop "This drivers requires Nspin==1"

  allocate(BetheDOS(Nk),BetheEk(Nlat,Nlat,Nk))
  call build_Hbethe2x2
  !<DEBUG 
  call splot("testDOSbethe.dat",BetheEk(1,1,:),BetheDOS(:))
  !>DEBUG


  !Allocate Weiss Field:
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
  do ineq=1,2
     write(tmp_suffix,'(I4.4)') ineq
     ed_file_suffix="_site"//trim(tmp_suffix)
     call ed_init_solver(Bath(ineq,:,:))
  enddo


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     bath_prev=bath

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     ! solve the impurities on each inequivalent site:
     do ineq=1,2
        xmu=xmu0-W0*n0(3-ineq)
        write(tmp_suffix,'(I4.4)') ineq
        ed_file_suffix="_site"//trim(tmp_suffix)
        call ed_solve(bath(ineq,:,:))
        n0(ineq) = ed_dens(1)
        call ed_get_sigma_matsubara(Smats(ineq,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(ineq,:,:,:,:,:))
        call ed_get_gimp_matsubara(Gmats(ineq,:,:,:,:,:))
        call ed_get_gimp_real(Greal(ineq,:,:,:,:,:))
     enddo
     print*,""
     print*,"n1 = ",n0(1)," n2= ",n0(2)
     print*,""

     do ineq=1,2
        xmu=xmu0-W0*n0(3-ineq)
        write(tmp_suffix,'(I4.4)') ineq
        ed_file_suffix="_site"//trim(tmp_suffix)
        call ed_get_gloc(one*BetheEk(1,1,:),BetheDOS,&
             Gmats_loc(ineq,1,1,1,1,:),&
             Greal_loc(ineq,1,1,1,1,:),&
             Smats(ineq,1,1,1,1,:),&
             Sreal(ineq,1,1,1,1,:),iprint=1)
        !Get Delta and fit it
        ! call ed_get_weiss(Gmats_loc(ineq,:,:,:,:,:),Smats(ineq,:,:,:,:,:),Delta(ineq,:,:,:,:,:),iprint=1)
        Delta(ineq,:,:,:,:,:) = (wband**2d0)/4d0*Gmats(3-ineq,:,:,:,:,:)!Gmats_loc(3-ineq,:,:,:,:,:)
        call ed_chi2_fitgf(Delta(ineq,1,1,:,:,:),bath(ineq,:,:),ispin=1)
     enddo


     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,1,1,:),dmft_error,nsuccess,nloop)

     !if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu0,converged)
     call end_loop
  enddo

  do ineq=1,2
     xmu = xmu0 - W0*n0(3-ineq)
     write(tmp_suffix,'(I4.4)') ineq
     ed_file_suffix="_site"//trim(tmp_suffix)
     call ed_kinetic_energy(one*BetheEk(1,1,:),BetheDOS,Smats(ineq,1,1,1,1,:))
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
    BetheEk(2,2,:)  = ek
  end subroutine build_Hbethe2x2



end program ed_cdw



