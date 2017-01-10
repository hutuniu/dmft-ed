program ed_SOC
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  !
  !#########   VARIABLEs DECLARATION   #########
  !
  integer                                       :: iloop,i,j
  integer                                       :: Nlat,ilat
  integer                                       :: Nso,io,jo
  integer                                       :: iorb,jorb,ispin,jspin
  logical                                       :: converged
  real(8)                                       :: wmixing
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  !Mpi:
  integer                                       :: comm,rank
  logical                                       :: master
  !Bath:
  integer                                       :: Nb,unit
  real(8),allocatable                           :: Bath(:)
  real(8),allocatable                           :: Bath_old(:)
  !Local functions:
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  !Weiss&Hybridization functions
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss_old
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Delta_old
  !Hmiltonian input:
  integer                                       :: Nk
  integer                                       :: Nkpath
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  real(8),allocatable,dimension(:)              :: Wtk
  complex(8),allocatable,dimension(:,:)         :: SOC_Hloc
  complex(8),allocatable,dimension(:,:,:,:)     :: SOC_Hloc_nn
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sigma_correction
  !Variables for the model:
  real(8)                                       :: soc,ivb
  !custom variables for rotations:
  logical                                       :: surface
  logical                                       :: Hk_test
  logical                                       :: rotateG0loc
  !custom variables for convergence test:
  complex(8),allocatable,dimension(:)           :: delta_conv
  !custom variables for chempot search:
  logical                                       :: converged_n,upprshft
  integer                                       :: conv_n_loop=0
  integer                                       :: shift_n_loop=0
  real(8)                                       :: Alvl=0.d0
  real(8)                                       :: bottom,top,shift
  real(8)                                       :: dw,sumdens,xmu_old
  real(8),allocatable,dimension(:)              :: w,orb_dens
  !custom variables for density matrix:
  real(8),allocatable,dimension(:)              :: dm_eig
  complex(8),allocatable,dimension(:,:)         :: dm,dm_rot
  !custom variables for SOC expectations:
  complex(8),allocatable,dimension(:,:,:)       :: Stot
  complex(8),allocatable,dimension(:,:,:)       :: Ltot
  complex(8),allocatable,dimension(:)           :: jz
  !
  !#########   MPI INITIALIZATION   #########
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !#########    VARIABLE PARSING    #########
  !
  call parse_cmd_variable(finput,       "FINPUT",             default='inputED_SOC.in')
  call parse_input_variable(hkfile,     "HKFILE",finput,      default="hkfile.in")
  call parse_input_variable(nk,         "NK",finput,          default=10)
  call parse_input_variable(nkpath,     "NKPATH",finput,      default=500)
  call parse_input_variable(wmixing,    "WMIXING",finput,     default=0.5d0)
  call parse_input_variable(soc,        "SOC",finput,         default=0.0d0)
  call parse_input_variable(ivb,        "IVB",finput,         default=0.0d0)
  call parse_input_variable(surface,    "SURFACE",finput,     default=.false.)
  call parse_input_variable(Hk_test,    "HK_TEST",finput,     default=.true.)
  call parse_input_variable(upprshft,   "upprshft",finput,    default=.false.)
  call parse_input_variable(rotateG0loc,"ROTATEG0loc",finput, default=.false.)
  !
  call ed_read_input(trim(finput),comm)
  !
  Nso=Nspin*Norb
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !#########       ALLOCATION       #########
  !
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));            Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));            Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));            Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));            Greal=zero
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));            Weiss=zero
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats));            delta=zero
  !
  allocate(weiss_old(Nspin,Nspin,Norb,Norb,Lmats));        weiss_old=zero
  allocate(delta_old(Nspin,Nspin,Norb,Norb,Lmats));        delta_old=zero
  allocate(Sigma_correction(Nspin,Nspin,Norb,Norb,Lmats)); Sigma_correction=zero
  !
  allocate(delta_conv(Lmats));                             delta_conv=zero
  !
  allocate(dm(Nspin*Norb,Nspin*Norb));                     dm=zero
  allocate(dm_eig(Nspin*Norb));                            dm_eig=zero
  allocate(dm_rot(Nspin*Norb,Nspin*Norb));                 dm_rot=zero
  !
  allocate(Stot(3,Norb,Norb));                             Stot=zero
  allocate(Ltot(3,Nspin,Nspin));                           Ltot=zero
  allocate(jz(3));                                         jz=zero
  !
  !#########        BUILD Hk        #########
  !
 ! call build_hk(trim(hkfile))
  !
  !#########          BATH          #########
  !
  if (bath_type/="replica") then
     Nb=get_bath_dimension()
  else
     Nb=get_bath_dimension(SOC_Hloc_nn)
  endif
  if(master)write(LOGfile,*)"Bath_size:",Nb
  allocate(Bath(Nb));     Bath=0.0d0
  allocate(Bath_old(Nb)); Bath_old=0.0d0
  !
  !#########      INIT SOLVER       #########
  !
  call ed_init_solver(Comm,Bath,SOC_Hloc_nn)
  !
  !#########          DMFT          #########
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve impurity
     call ed_solve(comm,Bath)
     !
     !get sigmas
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     !get local Gf's
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats,iprint=4)
     call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal,iprint=4)
     !
     !get new Weiss/Delta
     if(cg_scheme=='weiss')then
        weiss_old=weiss
        call dmft_weiss(Gmats,Smats,Weiss,SOC_Hloc_nn,iprint=4)
     else
        delta_old=delta
        call dmft_delta(Gmats,Smats,Delta,SOC_Hloc_nn,iprint=4)
     endif
     !
     !mixing Weiss/Delta
     if(master)write(LOGfile,'(a10,F7.2,a10,i3,2a10)') "wmixing",wmixing,"cg_weight",cg_weight,"scheme",cg_scheme
     if(iloop>1)then
        if(cg_scheme=='weiss')then
           weiss = wmixing*weiss + (1.d0-wmixing)*weiss_old
        else
           delta = wmixing*delta + (1.d0-wmixing)*delta_old
        endif
     endif
     !
     !fit for new Anderson parameters
     bath_old  = bath
     if (ed_mode=="normal") then
        call ed_chi2_fitgf(Comm,delta,bath,ispin=1)
        call spin_symmetrize_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(Comm,delta,bath)
     endif
     !
     !each loop operations
     if(bath_type=="replica")then
        Alvl=0.8d0
        if(master)then
           call ed_get_density_matrix(dm,dm_eig,dm_rot)
           call ed_get_quantum_SOC_operators(Stot,Ltot,jz)
        endif
  !      call Jz_rotate(Greal,"Gw","A",bottom,top,pi*Alvl)
  !      call Jz_rotate(Smats,"Sw","A")
     endif
     !
     !chemical potential find
     converged_n=.true.
     xmu_old=xmu
     allocate(orb_dens(Norb));orb_dens=0.d0
     call ed_get_dens(orb_dens);sumdens=sum(orb_dens)
     deallocate(orb_dens)
     if(master)write(*,'(3(a10,F10.5))') "sumdens",sumdens,"diffdens",abs(nread-sumdens),"nread",nread
     if(nread/=0.d0)then
        converged_n=.false.
        if(iloop>=3)call search_chempot(xmu,sumdens,converged_n)
        if(master)write(*,'(2(a10,F10.5))') "xmu_old",xmu_old,"xmu_new",xmu
     endif
     if(converged_n)then
        conv_n_loop=conv_n_loop+1
     else
        conv_n_loop=0
     endif
     !
     !convergence
     do i=1,Lmats
        delta_conv(i)=sum(nn2so_reshape(delta(:,:,:,:,i),Nspin,Norb))
     enddo
     if(master) then
        converged = check_convergence(delta_conv,dmft_error,nsuccess,nloop)
        write(LOGfile,'(a35,L3)') "sigma converged",converged
        write(LOGfile,'(a35,L3)') "dens converged",converged_n
        converged = converged .and. converged_n
        write(LOGfile,'(a35,L3)') "total converged",converged
        write(LOGfile,'(a35,I3)') "global iloop",iloop
        write(LOGfile,'(a35,I3)') "times dens is ok",conv_n_loop
        write(LOGfile,'(a35,I3)') "times rigid shift",shift_n_loop
     endif
     call Bcast_MPI(Comm,converged)
     !call MPI_Barrier(Comm)
     !
     !final mu shift
     if(converged_n.and.upprshft.and.((nread==5.d0).or.(nread==2.d0)))then
        shift_n_loop=shift_n_loop+1
        if(bath_type/="replica")then
           if(allocated(w))deallocate(w);allocate(w(Lreal));w=0.0d0
           w = linspace(wini,wfin,Lreal,mesh=dw)
           loop1: do i=1,Lreal
              if(abs(aimag(Greal(1,1,1,1,i))).gt.0.8d0)then
                 bottom=w(i)
                 exit loop1
              endif
           enddo loop1
           loop2: do i=1,Lreal
              if(abs(aimag(Greal(1,1,1,1,Lreal-i+1))).gt.0.8d0)then
                 top=w(Lreal-i+1)
                 exit loop2
              endif
           enddo loop2
        endif
        if(master)write(LOGfile,*)"top",top,"bottom",bottom
        shift      = bottom + ( top - bottom ) / 2.d0
        xmu_old    = xmu
        if(abs(shift)>=0.005)then
           xmu        = xmu_old + shift
           converged  = .false.
           nread  = 0.0d0!con questo una volta che comincio a shiftare rigidamente la densità non la ricontrollo piu
        endif
        if(master)then
           write(LOGfile,'(5(a10,F10.5))') "shift",shift,"xmu_old",xmu_old,"xmu_new",xmu
           unit=free_unit()
           open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
           write(unit,*)xmu,sumdens,sumdens-nerr,"shift"
           close(unit)
        endif
     endif
     !
     if(master)call end_loop
  enddo
  !
  !#########    BUILD Hk ON PATH    #########
  !
 ! call build_eigenbands()
  !
  !
  call finalize_MPI()
  !
  !
contains






end program ed_SOC
