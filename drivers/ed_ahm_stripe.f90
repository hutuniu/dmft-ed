!########################################################
!PURPOSE  :solve the attractive (A) Hubbard model
!          on a stripe geometry with periodic Hubbard moduation
!AUTHOR   :G Mazza & A Amaricci
!########################################################
program ed_stripe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  ! USE RDMFT
  ! USE RANDOM,    only:nrand,init_random_number
  ! USE ERROR
  ! USE TOOLS
  ! USE IOTOOLS
  ! USE ARRAYS
  ! USE STATISTICS
  ! USE SQUARE_LATTICE
  ! USE DMFT_ED
  implicit none
  ! complex(8),allocatable,dimension(:,:,:) :: Smats,Sreal !self_energies  !(2,Nsites,Nspin,Nspin,Norb,Norb,Lmats) 
  ! complex(8),allocatable,dimension(:,:,:) :: Smats_,Sreal_ !self_energies
  ! !WARNING:  in get_sigma separate normal and anomalous components !!!
  ! complex(8),allocatable,dimension(:,:,:) :: Gmats,Greal !local green's functions  ! Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats) 
  ! complex(8),allocatable,dimension(:,:,:) :: Gmats_,Greal_ !local green's functions
  ! complex(8),allocatable,dimension(:,:,:) :: Delta          !(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
  ! complex(8),allocatable,dimension(:,:,:) :: Delta_      
  complex(8),allocatable                :: Hloc(:,:,:,:,:)       ![Nlat][Nspin][Nspin][Norb][Norb]
  complex(8),allocatable                :: Hloc_(:,:,:,:,:)      ![Nindep][Nspin][Nspin][Norb][Norb]

  complex(8),allocatable                :: Smats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable                :: Sreal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable                :: Gmats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable                :: Greal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable                :: Delta(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]

  complex(8),allocatable                :: Smats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable                :: Sreal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable                :: Gmats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable                :: Greal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable                :: Delta_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]


  !<DEBUG
  complex(8),allocatable,dimension(:,:,:) :: Gmats_1b          !(2,Nlat,Lmats)
  complex(8),allocatable,dimension(:,:,:) :: Gmats_1b_         !(2,Nindep,Lmats)         
  complex(8),allocatable,dimension(:,:,:) :: Smats_1b          !(2,Nlat,Lmats)
  complex(8),allocatable,dimension(:,:,:) :: Smats_1b_         !(2,Nindep,Lmats)         
  complex(8),allocatable,dimension(:,:,:) :: Delta_1b          !(2,Nlat,Lmats)
  complex(8),allocatable,dimension(:,:,:) :: Delta_1b_         !(2,Nindep,Lmats)       
  !DEBUG>



  real(8),allocatable,dimension(:,:,:)    :: bath,bath_old
  real(8),allocatable,dimension(:,:,:)    :: bath_,bath_old_

  real(8),allocatable,dimension(:)        :: nii,dii,pii,eii
  real(8),allocatable,dimension(:)        :: nii_,dii_,pii_,eii_

  real(8),dimension(:),allocatable          :: wm,wr  
  complex(8),allocatable :: Hk(:,:,:)              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]

  real(8) :: ts
  real(8),allocatable,dimension(:,:)        :: Usite
  real(8),allocatable,dimension(:,:)      :: Uij,Unodes

  logical                                 :: converged
  real(8)                                 :: Uamplitude
  real(8)                                 :: r,wmixing
  real(8),allocatable,dimension(:)        :: epsik,wt,k_grid
  logical,allocatable,dimension(:)        :: hk_symm  
  integer                                 :: Uperiod,Nperiod
  integer                                 :: i,is,iloop,ik
  integer                                 :: Nb(2),Nx,Lk
  integer                                 :: Nrow,Ncol
  integer                                 :: row,col,ilat,i_ind

  integer                                 :: unit
  logical                                 :: pbc_row,pbc_col
  logical                                 :: symmetry_flag
  integer                                 :: symmetry_type

  !+---------+!
  ! START MPI !
  !+---------+!
#ifdef _MPI_INEQ
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif
  !+------------------------+!
  !+- READ INPUT VARIABLES +-!
  !+------------------------+!
  call parse_input_variable(ts,"TS","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(Nrow,"Nrow","inputRDMFT.in",default=4)
  call parse_input_variable(Uperiod,"U_PERIOD","inputRDMFT.in",default=10,comment='Period of the Interaction modulation (lattice sites units)')
  call parse_input_variable(Nperiod,"NU_PERIOD","inputRDMFT.in",default=1,comment='Number of periods contained in the stripe')
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=1.d0)
  call parse_input_variable(Uamplitude,"Uamplitude","inputRDMFT.in",default=0.d0)
  call parse_input_variable(pbc_row,"PBC_ROW","inputRDMFT.in",default=.true.)
  call parse_input_variable(pbc_col,"PBC_COL","inputRDMFT.in",default=.false.)
  call parse_input_variable(symmetry_flag,"REFLECTION_SYMM","inputRDMFT.in",default=.true.,comment='Set lattice reflection symmetry')  
  !
  !call rdmft_read_input("inputRDMFT.in")
  call ed_read_input("inputRDMFT.in")
  call set_store_size(1024)

  !+-----------------------------+!
  !+- BUILD LATTICE HAMILTONIAN -+!
  !+-----------------------------+!
  Ncol = Uperiod
  Nlat = Nrow*Ncol
  !
  !call get_lattice_hamiltonian(Nrow,Ncol,pbc_row=pbc_row,pbc_col=pbc_col)
  !
  unit=free_unit()
  open(unit,file='k_grid.lattice')
  Lk=Nperiod
  allocate(wt(Lk),epsik(Lk),k_grid(Lk),hk_symm(Lk))
  do ik=1,Lk
     k_grid(ik) = 2*pi/dble(Lk)*dble(ik-1)
     epsik(ik) = 2.d0*cos(k_grid(ik))
     wt(ik) = 1.d0/dble(Lk)
     write(unit,'(6(F18.10))') dble(ik),k_grid(ik),epsik(ik),wt(ik)
  end do
  close(unit)
  call get_k_hamiltonian_stripe(Nrow,Ncol,pbc_row,pbc_col,k_grid)

  !+-------------------------+!
  !+- BUILD LATTICE DETAILS -+!
  !+-------------------------+!
  allocate(Usite(Nlat,Norb),Uij(Ncol,Nrow))
  Usite=Uloc(1)
  if(mpiID==0) open(unit,file='Unodes.lattice')
  do row=0,Nrow-1
     do col=0,Ncol-1
        ilat = col + row*Ncol + 1
        if(Uperiod.gt.1) Usite(ilat,:) = Usite(ilat,:) + Uamplitude*dsin(2.d0*pi*dble(col)/dble(Uperiod))
        Uij(col+1,row+1) = Usite(ilat,1)
        if(abs(dsin(2.d0*pi*dble(col)/dble(Uperiod)))<1.d-10) then
           if(mpiID==0) write(unit,*) dble(col+1),dble(row+1)
        end if
     end do
  end do
  if(mpiID==0) call splot3d("Ustripe.ed",(/(dble(i),i=1,Ncol)/),(/(dble(i),i=1,Nrow)/),Uij)  
  if(mpiID==0) close(unit)  


  !+----------------------------------+!  
  !+- ALLOCATE GF & INITIALIZE BATHS -+!
  !+----------------------------------+!  
  ! Matsubara and Real freq
  allocate(wm(Lmats),wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
  ! Independent sites baths
  Nb=get_bath_size()
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_old(Nlat,Nb(1),Nb(2)))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  ! Observables
  allocate(nii(Nlat))
  allocate(dii(Nlat))
  allocate(pii(Nlat))
  allocate(eii(Nlat))
  ! Self energies
  allocate(Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! Green function
  allocate(Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! Impurity-bath hybritizations
  allocate(Delta(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))

  !<DEBUG
  allocate(Gmats_1b(2,Nlat,Lmats))
  allocate(Smats_1b(2,Nlat,Lmats))
  allocate(Delta_1b(2,Nlat,Lmats))
  !DEBUG>


  !
  call  ed_init_solver_lattice(bath)
  Hloc=0.d0
  !call init_lattice_baths(bath)

  if(symmetry_flag) then     
     Nsymm=1
     call get_independent_sites(reflect)
     !
     allocate(bath_(Nindep,Nb(1),Nb(2)))
     allocate(bath_old_(Nindep,Nb(1),Nb(2)))
     allocate(Hloc_(Nindep,Nspin,Nspin,Norb,Norb))
     ! Observables
     allocate(nii_(Nindep))
     allocate(dii_(Nindep))
     allocate(pii_(Nindep))
     allocate(eii_(Nindep))
     ! Self energies
     allocate(Smats_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))
     allocate(Sreal_(2,Nindep,Nspin,Nspin,Norb,Norb,Lreal))
     ! Green function
     allocate(Gmats_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))
     allocate(Greal_(2,Nindep,Nspin,Nspin,Norb,Norb,Lreal))
     ! Impurity-bath hybritizations
     allocate(Delta_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))

     !<DEBUG
     allocate(Gmats_1b_(2,Nindep,Lmats))
     allocate(Smats_1b_(2,Nindep,Lmats))
     allocate(Delta_1b_(2,Nindep,Lmats))
     !DEBUG>

     do i_ind=1,Nindep
        bath_(i_ind,:,:) = bath(indep_list(i_ind),:,:)
        Hloc_(i_ind,:,:,:,:) = Hloc(indep_list(i_ind),:,:,:,:)
     end do
  end if




  !+-------------+!
  !+- DMFT LOOP -+!
  !+-------------+!
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop")

     if(symmetry_flag) then
        !+- LOOP WITH SYMMETRIES -+!
        if(rdmft_phsym)then
           do i_ind=1,Nindep
              call ph_symmetrize_bath(bath_(i_ind,:,:))
           enddo
        endif
        bath_old_=bath_
        !call ed_solve_impurity(bath_,Smats_,Sreal_,nii_,dii_,pii_,eii_,usite=usite)
        call ed_solve_lattice(bath_,Hloc_,Uloc_ii=Usite)
        nii_ = ed_get_dens_lattice(Nindep,1)
        dii_ = ed_get_docc_lattice(Nindep,1)
        eii_ = ed_get_epot_lattice(Nindep)
        pii_ = ed_get_phisc_lattice(Nindep,1)
        call ed_get_sigma_matsubara_lattice(Smats_(1,:,:,:,:,:,:),Nindep)
        call ed_get_sigma_real_lattice(Sreal_(1,:,:,:,:,:,:),Nindep)
        call ed_get_self_matsubara_lattice(Smats_(2,:,:,:,:,:,:),Nindep)
        call ed_get_self_real_lattice(Sreal_(2,:,:,:,:,:,:),Nindep)
        !
        do ilat=1,Nlat
           i_ind=map_lat2ind(ilat)
           Smats(:,ilat,:,:,:,:,:) = Smats_(:,i_ind,:,:,:,:,:)
           Sreal(:,ilat,:,:,:,:,:) = Sreal_(:,i_ind,:,:,:,:,:)
           nii(ilat) = nii_(i_ind)
           dii(ilat) = dii_(i_ind)
           pii(ilat) = pii_(i_ind)
           eii(ilat) = eii_(i_ind)
        end do
        !
        !call rdmft_get_gloc(Hk,wt,Gmats,Greal,Smats,Sreal,hk_symm=hk_symm)
        call ed_get_gloc_lattice(Hk,Wt,Gmats,Greal,Smats,Sreal,1,hk_symm=hk_symm)
        !
        do i_ind=1,Nindep
           ilat=indep_list(i_ind)
           Gmats_(:,i_ind,:,:,:,:,:)=Gmats(:,ilat,:,:,:,:,:)
           Greal_(:,i_ind,:,:,:,:,:)=Greal(:,ilat,:,:,:,:,:)
        end do
        
        !<DEBUG
        do i_ind=1,Nindep
           do i=1,Lmats
              Gmats_1b_(:,i_ind,i)=Gmats_(:,i_ind,1,1,1,1,i)
              Smats_1b_(:,i_ind,i)=Smats_(:,i_ind,1,1,1,1,i)
           end do
        end do
        ilat=1
        do i=1,Lmats
           write(500,'(10(F18.10))') wm(i),Smats_1b_(1,ilat,i),Smats_(1,ilat,1,1,1,1,i)
           write(501,'(10(F18.10))') wm(i),Gmats_1b_(1,ilat,i),Gmats_(1,ilat,1,1,1,1,i)
        end do

        Delta_1b_=zero

        !<DEBUG
        ! write(*,*)
        ! write(*,*)
        ! write(*,*) "eloc-hloc Nmats=10"
        ! write(*,*)
        ! write(*,*)        
        !DEBUG>

        !<DEBUG
        ! write(*,*) 
        ! write(*,*) "input-Smats  hloc "
        ! write(*,'(4(F28.20))') dreal(Smats_(:,1,1,1,1,1,10))
        ! write(*,'(4(F28.20))') dimag(Smats_(:,1,1,1,1,1,10))
        ! write(*,*)
        !DEBUG>
        call ed_get_weiss_lattice(Nindep,Gmats_,Smats_,Delta_,Hloc_)                        
        !<DEBUG
        ! write(*,*) 
        ! write(*,*) "output-Smats  hloc"
        ! write(*,'(4(F28.20))') dreal(Smats_(:,1,1,1,1,1,10))
        ! write(*,'(4(F28.20))') dimag(Smats_(:,1,1,1,1,1,10))
        ! write(*,*)
        !DEBUG>

        !<DEBUG
        ! write(*,*) 
        ! write(*,*) "input-Smats  hloc "
        ! write(*,'(4(F28.20))') dreal(Smats_1b_(:,1,10))
        ! write(*,'(4(F28.20))') dimag(Smats_1b_(:,1,10))
        ! write(*,*)
        !DEBUG>
        call ed_get_weiss_lattice(Nindep,Gmats_1b_,Smats_1b_,Delta_1b_)        
        !<DEBUG
        ! write(*,*) 
        ! write(*,*) "output-Smats  eloc "
        ! write(*,'(4(F28.20))') dreal(Smats_1b_(:,1,10))
        ! write(*,'(4(F28.20))') dimag(Smats_1b_(:,1,10))
        ! write(*,*)
        !DEBUG>

        !<DEBUG
        ! do i=1,Lmats
        !    write(600,'(10(F18.10))') wm(i),Delta_1b_(1,ilat,i),Delta_(1,ilat,1,1,1,1,i)
        !    write(601,'(10(F18.10))') wm(i),Delta_1b_(2,ilat,i),Delta_(2,ilat,1,1,1,1,i)
        ! end do
        !DEBUG>

 !       write(*,'(10(F28.20))') Delta_(1,1,1,1,1,1,10)

        !<DEBUG
        ! write(*,*)
        ! write(*,*)
        ! write(*,*) "hloc-eloc Nmats=10"
        ! write(*,*)
        ! write(*,*)        
        !DEBUG>


        ! call ed_get_weiss_lattice_hloc(Nindep,Gmats_,Smats_,Delta_,Hloc_)       
        ! call ed_get_weiss_lattice_eloc(Nindep,Gmats_1b_,Smats_1b_,Delta_1b_)        


!        write(*,'(10(F28.20))') Delta_(1,1,1,1,1,1,10)

        !Delta_(1:2,1:Nindep,1,1,1,1,1:Lmats) = Delta_1b_(1:2,1:Nindep,1:Lmats)
        !
!        ilat=1
        !DEBUG>

        !<DEBUG
        ! do i=1,Lmats
        !    write(700,'(10(F18.10))') wm(i),Delta_1b_(1,ilat,i),Delta_(1,ilat,1,1,1,1,i)
        !    write(701,'(10(F18.10))') wm(i),Delta_1b_(2,ilat,i),Delta_(2,ilat,1,1,1,1,i)
        ! end do
        !DEBUG>
        call ed_chi2_fitgf_lattice(bath_,Delta_,Hloc_)

        !
        bath_=wmixing*bath_ + (1.d0-wmixing)*bath_old_
        !
     else
        bath_old=bath
        !call ed_solve_impurity(bath,Smats,Sreal,nii,dii,pii,eii,usite=usite)
        call ed_solve_lattice(bath,Hloc)
        nii = ed_get_dens_lattice(Nlat,1)
        dii = ed_get_docc_lattice(Nlat,1)
        eii = ed_get_epot_lattice(Nlat)
        pii = ed_get_phisc_lattice(Nlat,1)        
        call ed_get_sigma_matsubara_lattice(Smats(1,:,:,:,:,:,:),Nlat)
        call ed_get_sigma_real_lattice(Sreal(1,:,:,:,:,:,:),Nlat)
        call ed_get_self_matsubara_lattice(Smats(2,:,:,:,:,:,:),Nlat)
        call ed_get_self_real_lattice(Sreal(2,:,:,:,:,:,:),Nlat)
        !
        !call rdmft_get_gloc(Hk,wt,Gmats,Greal,Smats,Sreal,hk_symm=hk_symm)
        call ed_get_gloc_lattice(Hk,Wt,Gmats,Greal,Smats,Sreal,1,hk_symm=hk_symm)
        !
        !call rdmft_get_weiss_field(Nlat,Gmats,Smats,Delta)
        !call ed_fit_bath(bath,Delta)
        call ed_get_weiss_lattice(Nlat,Gmats,Smats,Delta,Hloc)        
        call ed_chi2_fitgf_lattice(bath,Delta,Hloc)
        !
        bath=wmixing*bath + (1.d0-wmixing)*bath_old
        if(rdmft_phsym)then
           do i=1,Nlat
              call ph_symmetrize_bath(bath(i,:,:))
           enddo
        endif
     end if
     !
     if(mpiID==0) converged = check_convergence_local(dii,dmft_error,Nsuccess,nloop,id=0,file="error.err")
     !
#ifdef _MPI_INEQ     
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     call print_sc_out(converged)
     if(mpiID==0)call end_loop()
  enddo
#ifdef _MPI_INEQ
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)  
#endif
  !+******************************************************************+!
  !+******************************************************************+!
  !+******************************************************************+!

CONTAINS

  !+----------------------+!
  !+- AUXILIARY ROUTINES -+!
  !+----------------------+!
  ! print-out 
  subroutine print_sc_out(converged)
    integer                              :: i,j,is,row,col
    real(8)                              :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)              :: cdwii,rii,sii,zii
    real(8),dimension(Nrow,Ncol)         :: dij,nij,cij,pij
    real(8),dimension(Nrow)              :: grid_x
    real(8),dimension(Ncol)              :: grid_y
    real(8)                              :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)            :: data_covariance
    real(8),dimension(2,2)               :: covariance_nd
    real(8),dimension(2)                 :: data_mean,data_sdev
    logical                              :: converged
    complex(8),dimension(2,Lmats)        :: aGmats,aSmats
    complex(8),dimension(2,Lreal)        :: aGreal,aSreal
    character(len=4)                     :: loop
    integer,dimension(6)                 :: units,N_min
    real(8)                              :: Eint,Ekin,Epot,Eout(2)
    complex(8),dimension(1,1,Nlat,Lmats) :: Sigma_tmp
    complex(8),dimension(1,1,Nlat,Lmats) :: SigmaA_tmp

    if(mpiID==0)then
       write(loop,"(I4)")iloop
       !Get CDW "order parameter"
       do is=1,Nlat
          row=irow(is)
          col=icol(is)
          cdwii(is) = (-1.d0)**(row+col)*(nii(is)-1.d0)
       enddo
       nimp = sum(nii)/dble(Nlat)
       phi  = sum(pii)/dble(Nlat)
       docc = sum(dii)/dble(Nlat)
       ccdw = sum(cdwii)/dble(Nlat)
       print*,"<nimp>  =",nimp
       print*,"<phi>   =",phi
       print*,"<docc>  =",docc
       print*,"<ccdw>  =",ccdw
       call splot("nVSiloop.data",iloop,nimp,append=.true.)
       call splot("phiVSiloop.data",iloop,phi,append=.true.)
       call splot("doccVSiloop.data",iloop,docc,append=.true.)
       call splot("ccdwVSiloop.data",iloop,ccdw,append=.true.)
       call store_data("nVSisite.data",nii,(/(dble(i),i=1,Nlat)/))
       call store_data("phiVSisite.data",pii,(/(dble(i),i=1,Nlat)/))
       call store_data("doccVSisite.data",dii,(/(dble(i),i=1,Nlat)/))

       ! Plots at convergency
       if(converged)then
          ! GF & Sigma
          call store_data("LDelta_iw.data",Delta(1,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LGamma_iw.data",Delta(2,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LG_iw.data",Gmats(1,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LF_iw.data",Gmats(2,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LG_realw.data",Greal(1,1:Nlat,1,1,1,1,1:Lreal),wr(1:Lreal))
          call store_data("LF_realw.data",Greal(2,1:Nlat,1,1,1,1,1:Lreal),wr(1:Lreal))
          call store_data("LSigma_iw.data",Smats(1,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LSelf_iw.data",Smats(2,1:Nlat,1,1,1,1,1:Lmats),wm(1:Lmats))
          call store_data("LSigma_realw.data",Sreal(1,1:Nlat,1,1,1,1,1:Lreal),wr(1:Lreal))
          call store_data("LSelf_realw.data",Sreal(2,1:Nlat,1,1,1,1,1:Lreal),wr(1:Lreal))
          ! Observables
          do is=1,Nlat
             cdwii(is) = (-1.d0)**(is)*(nii(is)-1.d0)
             sii(is)   = dimag(Smats(1,is,1,1,1,1,1))-&
                  wm(1)*(dimag(Smats(1,is,1,1,1,1,2))-dimag(Smats(1,is,1,1,1,1,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(Gmats(1,is,1,1,1,1,1))-&
                  wm(1)*(dimag(Gmats(1,is,1,1,1,1,2))-dimag(Gmats(1,is,1,1,1,1,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(1,is,1,1,1,1,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)

          units = free_units(6)
          open(units(1),file='n_col.data')
          open(units(2),file='docc_col.data')
          open(units(3),file='phi_col.data')
          open(units(4),file='n_row.data')
          open(units(5),file='docc_row.data')
          open(units(6),file='phi_row.data')
          !
          do col=1,Ncol
             grid_y(col)=col
             do row=1,Nrow
                grid_x(row)  = row
                i            = ij2site(row,col)
                nij(row,col) = nii(i)
                dij(row,col) = dii(i)
                pij(row,col) = pii(i)
             enddo
          enddo
          !
          do row=1,Nrow
             write(units(1),'(100(f18.10))') dble(row),nij(row,:)
             write(units(2),'(100(f18.10))') dble(row),dij(row,:)
             write(units(3),'(100(f18.10))') dble(row),pij(row,:)
          end do
          !
          do col=1,Ncol
             write(units(4),'(100(f18.10))') dble(col),nij(:,col)
             write(units(5),'(100(f18.10))') dble(col),dij(:,col)
             write(units(6),'(100(f18.10))') dble(col),pij(:,col)
          end do
          !
          call store_data("cdwVSisite.data",cdwii,(/(dble(i),i=1,Nlat)/))
          call store_data("rhoVSisite.data",rii,(/(dble(i),i=1,Nlat)/))
          call store_data("sigmaVSisite.data",sii,(/(dble(i),i=1,Nlat)/))
          call store_data("zetaVSisite.data",zii,(/(dble(i),i=1,Nlat)/))
          call splot3d("3d_nVSij.data",grid_x,grid_y,nij)
          call splot3d("3d_doccVSij.data",grid_x,grid_y,dij)
          call splot3d("3d_phiVSij.data",grid_x,grid_y,pij)

       end if
    end if

    if(converged) then
       !
       Eint=0.d0
       Ekin=0.d0
       Epot=0.d0
       ! Sigma_tmp(1,1,:,:) = Smats(1,:,:)
       ! SigmaA_tmp(1,1,:,:) = Smats(2,:,:)
       Epot=sum(eii)/dble(Nlat)
       !Ekin=kinetic_energy_lattice(Hk,wt,Sigma_tmp,SigmaA_tmp)
       Eout = ed_kinetic_energy_lattice(Hk,Wt,Smats(1,:,:,:,:,:,:),Smats(2,:,:,:,:,:,:))
       !Ekin=ed_kinetic_energy_lattice(Hk,Wtk,Sigma,SigmaA)
       Ekin = Eout(1)
       Eint=Ekin+Epot
       unit=free_unit()
       if(mpiID==0) then
          open(unit,file='internal_energy.data')
          write(unit,'(10(F18.10))') Eint,Ekin,Epot
          close(unit)
       end if
       !
    end if

  end subroutine print_sc_out



  ! build Hk
  subroutine get_k_hamiltonian_stripe(Nrow,Ncol,pbc_col,pbc_row,k_grid)
    integer               :: Nrow
    integer               :: Ncol
    integer               :: Nsquare
    logical               :: pbc_row,pbc_col
    logical               :: symm
    logical               :: k_connect
    integer               :: i,jj,j,k,row,col,link(4),Lk,ik
    integer               :: unit,unitk
    real(8),optional      :: k_grid(:)
    real(8),allocatable   :: htmp(:,:),wk(:)    
    real(8),dimension(:,:),allocatable  :: tL,tR,H0
    if(Nlat /= Nrow*Ncol) stop "Nlat != Nrow*Ncol"
    if(present(k_grid)) then
       Lk=size(k_grid)
    else
       Lk=1
    end if
    ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk] 
    allocate(H0(Nlat,Nlat))
    allocate(Hk(Nlat,Nlat,Lk),tL(Nlat,Nlat),tR(Nlat,Nlat))

    !+-----------------------------------------------------+!
    !+- GET real-space block hamiltonian H0 (aka H-IRRED) -+!
    !+-----------------------------------------------------+!
    H0=0.d0
    unit=free_unit()
    if(mpiID==0) open(unit,file='rdmft_sites.lattice')
    allocate(icol(Nlat),irow(Nlat))
    allocate(ij2site(Nrow,Ncol))
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          irow(i)=row+1
          icol(i)=col+1
          ij2site(row+1,col+1)=i
          !
          if(mpiID==0) write(unit,*) dble(col+1),dble(row+1)

          ! right hop
          link(1)= i + 1     
          if((col+1)==Ncol) then
             if(pbc_col) then
                link(1)=1+row*Ncol  
             else
                link(1)=0  
             end if
          end if

          ! left  hop
          link(3)= i - 1    
          if((col-1)<0) then
             if(pbc_col) then
                link(3)=Ncol+row*Ncol
             else
                link(3)=0  
             end if
          end if

          ! up hop
          link(2)= i + Ncol 
          if((row+1)==Nrow) then
             if(pbc_row) then
                link(2)=col+1
             else
                link(2)=0  
             end if
          end if

          ! down  hop
          link(4)= i - Ncol 
          if((row-1)<0) then
             if(pbc_row) then
                link(4)=col+1+(Nrow-1)*Ncol
             else
                link(4)=0  
             end if
          end if
          !
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts 
          enddo
          !
       enddo
    enddo

    if(mpiID==0) close(unit)

    !+------------------------------------------------------------------------+!
    !+- GET real-space connection between blocks hamiltonian (aka H-CONNECT) -+!
    !+------------------------------------------------------------------------+
    Hk = 0.d0
    tL=0.d0;tR=0.d0
    !  tLeft 
    do row=0,Nrow-1
       col=0
       i=col+ 1 + row*Ncol
       link=0
       link(3)=Ncol+row*Ncol
       do jj=1,4
          if(link(jj)>0) then
             tL(i,link(jj))=-ts 
          end if
       enddo
    end do
    !+- build up tRight -+!
    do row=0,Nrow-1
       col=Ncol-1
       i=col+ 1 + row*Ncol
       link=0
       link(1)=1+row*Ncol  
       do jj=1,4
          if(link(jj)>0) then
             tR(i,link(jj))=-ts
          end if
       enddo
    end do
    !+- local fourier transform -+!
    do ik=1,Lk
       Hk(:,:,ik) = H0
       if(present(k_grid)) then
          Hk(:,:,ik) = Hk(:,:,ik) + tR(:,:)*exp(xi*k_grid(ik))+tL(:,:)*exp(-xi*k_grid(ik))
       end if
    end do

    unit=free_unit()
    open(unit,file="hk_symmetric_points")
    Hk_symm=.true.
    do ik=1,Lk
       do i=1,Nlat
          do j=1,Nlat
             if(Hk(i,j,ik) /= Hk(j,i,ik)) Hk_symm(ik)=.false.
          end do
       end do
       if(Hk_symm(ik)) write(unit,'(4(F18.10))') k_grid(ik)
    end do
  end subroutine get_k_hamiltonian_stripe





  function reflect(isite) result(jsite)
    integer             :: isite
    integer             :: row,col
    integer,allocatable :: jsite(:)
    integer             :: test_nsymm
    integer             :: rj,cj,isymm,itrans,ireflect
    allocate(jsite(Nsymm))
    !
    row=irow(isite)
    col=icol(isite)
    !
    !+- x-axis reflection -+!
    isymm=1
    rj=Nrow-(row-1)
    cj=col
    jsite(isymm) = ij2site(rj,cj)
  end function Reflect







end program ed_stripe

