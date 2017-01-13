program hm_2bands_dos
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso
  logical                                     :: converged
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !
  real(8),dimension(2)                        :: Wband
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0,de
  real(8)                                     :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,Eout(2)
  character(len=10)                           :: dos_model


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wband,"WBAND",finput,default=[1d0,0.5d0])
  call parse_input_variable(delta,"DELTA",finput,default=0d0)
  call parse_input_variable(dos_model,"DOS_MODEL",finput,"bethe")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=2)stop "Wrong setup from input file: Nspin=1; Norb=2"
  Nso=Nspin*Norb

  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(de(Nso))
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  select case(dos_model)
  case ("bethe")
     Dbands(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
     Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)
  case ("flat")
     Dbands(1,:) = dens_flat(Ebands(1,:),Wband(1))*de(1)
     Dbands(2,:) = dens_flat(Ebands(2,:),Wband(2))*de(2)
  case default
     stop "error: dos_model not in {bethe,flat}. Add your own if needed"
  end select

  allocate(H0(Nso))
  H0=[-Delta/2,Delta/2]
  call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc(1,1,:,:)=diag(H0)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Ebands,Dbands,H0,Gmats,Smats,iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc,iprint=1)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc,iprint=1)
     endif
     !
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss,bath,ispin=1)
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     !
     call end_loop
  enddo


  call dmft_gloc_realaxis(Ebands,Dbands,H0,Greal,Sreal,iprint=1)
  Eout = dmft_kinetic_energy(Ebands,Dbands,H0,Smats(1,1,:,:,:))


contains



  function dens_flat(ebands,wband) result(dens)
    real(8),dimension(:)            :: ebands
    real(8)                         :: wband
    real(8),dimension(size(ebands)) :: dens
    integer                         :: i
    real(8)                         :: e
    do i=1,size(ebands)
       e=ebands(i)
       dens(i)= step(wband-abs(e))/(2*wband)
    enddo
  end function dens_flat



end program hm_2bands_dos



