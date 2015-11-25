!########################################################
!PURPOSE  :solve the attractive (A) Hubbard model
!          on a stripe geometry with periodic Hubbard moduation
!AUTHOR   :G Mazza & A Amaricci
!########################################################
program ed_stripe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE DMFT_VECTORS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none
  
  real(8),dimension(:),allocatable   :: wm,wr  
  !
  complex(8),allocatable             :: Hk(:,:,:)              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  complex(8),allocatable             :: Hloc(:,:,:,:,:)       ![Nlat][Nspin][Nspin][Norb][Norb]
  complex(8),allocatable             :: Hloc_(:,:,:,:,:)      ![Nindep][Nspin][Nspin][Norb][Norb]
  !
  complex(8),allocatable             :: Smats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Sreal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Gmats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Greal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Delta(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !
  complex(8),allocatable             :: Smats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Sreal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Gmats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Greal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Delta_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  !
  real(8),allocatable,dimension(:,:) :: bath,bath_old
  real(8),allocatable,dimension(:,:) :: bath_,bath_old_
  integer                            :: Nb
  !
  integer                            :: Nindep,Nsymm
  !
  real(8),allocatable,dimension(:)   :: nii,dii,pii,eii
  real(8),allocatable,dimension(:)   :: nii_,dii_,pii_,eii_
  !  
  real(8)                            :: ts
  real(8),allocatable,dimension(:,:) :: Usite,Usite_
  real(8),allocatable,dimension(:,:) :: Uij,Unodes

  logical                            :: converged
  real(8)                            :: Uamplitude
  real(8)                            :: wmixing
  real(8),allocatable,dimension(:)   :: wt
  logical,allocatable,dimension(:)   :: hk_symm,hk_symm_
  integer                            :: Uperiod,Nperiod
  integer                            :: i,iloop,ik
  integer                            :: Lk
  
  integer                            :: Nrow,Ncol
  integer                            :: row,col,ilat,i_ind
  
  integer                            :: unit
  logical                            :: pbc_row,pbc_col
  logical                            :: symmetry_flag,rdmft_phsym
  integer                            :: symmetry_type
  !
  integer                            :: Xperiod,Yperiod
  integer                            :: N_Xperiod,N_Yperiod
  logical                            :: Xpbc,Ypbc
  integer                            :: Nx,Ny,ix,iy
  !

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
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(Xperiod,"XPERIOD","inputRDMFT.in",default=10)
  call parse_input_variable(Yperiod,"YPERIOD","inputRDMFT.in",default=1)
  call parse_input_variable(N_Xperiod,"N_XPERIOD","inputRDMFT.in",default=10)
  call parse_input_variable(N_Yperiod,"N_YPERIOD","inputRDMFT.in",default=10)  
  call parse_input_variable(Xpbc,"Xpbc","inputRDMFT.in",default=.true.)
  call parse_input_variable(Ypbc,"Ypbc","inputRDMFT.in",default=.true.)
  call parse_input_variable(Uamplitude,"Uamplitude","inputRDMFT.in",default=0.d0)
  call parse_input_variable(symmetry_flag,"REFLECTION_SYMM","inputRDMFT.in",default=.true.)
  call parse_input_variable(rdmft_phsym,"RDMFT_PHSYM","inputRDMFT.in",default=.true.)
  !
  call ed_read_input("inputRDMFT.in")
  call set_store_size(1024)



  !+-----------------------------+!
  !+- BUILD LATTICE HAMILTONIAN -+!
  !+-----------------------------+!
  Nx=Xperiod*N_Xperiod;Ncol=Xperiod
  Ny=Yperiod*N_Yperiod;Nrow=Yperiod
  Lk = N_Xperiod*N_Yperiod
  !+- CONSTRAINTS ON INPUTS -+!
  if(Xperiod==1.and.N_Xperiod>1) Xpbc=.true.
  if(Yperiod==1.and.N_Yperiod>1) Ypbc=.true.
  !+-------------------------+!
  allocate(wt(Lk),hk_symm(Lk))
  call  get_Hk_2dsquare(Nx,Xpbc,Ny,Ypbc,Xperiod,Yperiod,Hk)
  wt=1.d0/dble(Lk)
  Nlat=size(Hk,1)
  
  allocate(icol(Nlat),irow(Nlat))
  allocate(ij2site(Nrow,Ncol))
  do row=0,Nrow-1
     do col=0,Ncol-1
        i=col+ 1 + row*Ncol
        irow(i)=row+1
        icol(i)=col+1
        ij2site(row+1,col+1)=i
     end do
  end do
  
  !+-------------------------+!
  !+- BUILD LATTICE DETAILS -+!
  !+-------------------------+!
  allocate(Usite(Nlat,Norb),Uij(Xperiod,Yperiod))
  Usite=Uloc(1)
  Uperiod = Xperiod
  if(mpiID==0) open(unit,file='Unodes.lattice')
  do iy=1,Yperiod
     do ix=1,Xperiod
        ilat = ix + (iy-1)*Xperiod 
        if(Uperiod.gt.1) Usite(ilat,:) = Usite(ilat,:) + Uamplitude*dsin(2.d0*pi*dble(ix-1)/dble(Uperiod))
        Uij(ix,iy) = Usite(ilat,1)
        write(*,*) Usite(ilat,1),ilat
        if(abs(dsin(2.d0*pi*dble(ix)/dble(Uperiod)))<1.d-10) then
           if(mpiID==0) write(unit,*) dble(ix),dble(iy)
        end if
     end do
  end do
  if(mpiID==0) call splot3d("Ustripe.ed",(/(dble(i),i=1,Xperiod)/),(/(dble(i),i=1,Yperiod)/),Uij)  
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
  allocate(bath(Nlat,Nb))
  allocate(bath_old(Nlat,Nb))
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
  !
  call  ed_init_solver_lattice(bath)
  Hloc=0.d0
  Hloc_=0.d0
  !
  if(Yperiod==1) symmetry_flag=.false.
  if(symmetry_flag) then     
     Nsymm=1
     call get_independent_sites(reflect,Nsymm,Nindep)
     !
     allocate(bath_(Nindep,Nb))
     allocate(bath_old_(Nindep,Nb))
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
     allocate(Usite_(Nindep,1))
     do i_ind=1,Nindep
        bath_(i_ind,:) = bath(indep_list(i_ind),:)
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
              call ph_symmetrize_bath(bath_(i_ind,:))
           enddo
        endif
        bath_old_=bath_
        call ed_solve_lattice(bath_,Hloc_,1,Uloc_ii=Usite_)
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
        call ed_get_gloc_lattice(Hk,Wt,Gmats,Greal,Smats,Sreal,1,hk_symm=hk_symm)
        !
        do i_ind=1,Nindep
           ilat=indep_list(i_ind)
           Gmats_(:,i_ind,:,:,:,:,:)=Gmats(:,ilat,:,:,:,:,:)
           Greal_(:,i_ind,:,:,:,:,:)=Greal(:,ilat,:,:,:,:,:)
        end do
        !
        call ed_get_weiss_lattice(Gmats_,Smats_,Delta_,Hloc_,1)                        
        call ed_chi2_fitgf_lattice(bath_,Delta_,Hloc_)
        !
        bath_=wmixing*bath_ + (1.d0-wmixing)*bath_old_
        !
     else
        if(rdmft_phsym)then
           do i=1,Nlat
              call ph_symmetrize_bath(bath(i,:))
           enddo
        endif
        bath_old=bath
        call ed_solve_lattice(bath,Hloc,1,Uloc_ii=Usite)
        nii = ed_get_dens_lattice(Nlat,1)
        dii = ed_get_docc_lattice(Nlat,1)
        eii = ed_get_epot_lattice(Nlat)
        pii = ed_get_phisc_lattice(Nlat,1)        
        call ed_get_sigma_matsubara_lattice(Smats(1,:,:,:,:,:,:),Nlat)
        call ed_get_sigma_real_lattice(Sreal(1,:,:,:,:,:,:),Nlat)
        call ed_get_self_matsubara_lattice(Smats(2,:,:,:,:,:,:),Nlat)
        call ed_get_self_real_lattice(Sreal(2,:,:,:,:,:,:),Nlat)
        !
        call ed_get_gloc_lattice(Hk,Wt,Gmats,Greal,Smats,Sreal,1,hk_symm=hk_symm)
        !
        call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc,1)        
        call ed_chi2_fitgf_lattice(bath,Delta,Hloc)
        !
        bath=wmixing*bath + (1.d0-wmixing)*bath_old
        if(rdmft_phsym)then
           do i=1,Nlat
              call ph_symmetrize_bath(bath(i,:))
           enddo
        endif
     end if
     !
     if(mpiID==0) converged = check_convergence_local(dii,dmft_error,Nsuccess,nloop,id=0,file="error.err")
#ifdef _MPI_INEQ     
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif
     if(mpiID==0)call print_sc_out(converged)
     if(mpiID==0)call end_loop()
  enddo
  call stripe_energy
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)  
  !+******************************************************************+!
  !+******************************************************************+!
  !+******************************************************************+!
CONTAINS
  !+----------------------+!
  !+- AUXILIARY ROUTINES -+!
  !+----------------------+!
  subroutine print_sc_out(converged)
    integer                              :: i,j,is,row,col,ilat,jlat
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

    !Get CDW "order parameter"
    do is = 1,Nlat
       ilat = mod(is-1,Nlat)+1
       jlat = (is-1)/Nlat+1
       cdwii(is) = (-1.d0)**(ilat+jlat)*(nii(is)-1.d0)
       write(*,'(3(I3),F18.10,I3)') is,ilat,jlat,cdwii(is),Nlat
    end do
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
    ! Plots at convergence
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
             ! i=col+ 1 + row*Ncol
             ! !
             ! irow(i)=row+1
             ! icol(i)=col+1
             ! ij2site(row+1,col+1)=i
             !i            = col + (row-1)*Ncol!ij2site(row,col)
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
    
  end subroutine print_sc_out


  subroutine stripe_energy
    integer,dimension(6)                 :: units,N_min
    real(8)                              :: Eint,Ekin,Epot,Eout(2)
    complex(8),dimension(1,1,Nlat,Lmats) :: Sigma_tmp
    complex(8),dimension(1,1,Nlat,Lmats) :: SigmaA_tmp
    !
    print*,""
    Eint=0.d0
    Ekin=0.d0
    Epot=0.d0
    Epot=sum(eii)/dble(Nlat)
    Eout = ed_kinetic_energy_lattice(Hk,Wt,Smats(1,:,:,:,:,:,:),Smats(2,:,:,:,:,:,:))
    Ekin = Eout(1)
    Eint=Ekin+Epot
    unit=free_unit()
    if(mpiID==0) then
       open(unit,file='internal_energy.data')
       write(unit,'(10(F18.10))') Eint,Ekin,Epot
       close(unit)
    end if
    !
  end subroutine stripe_energy
  


  subroutine  get_Hk_2dsquare(Nx,Xpbc,Ny,Ypbc,Lx,Ly,Hk_lat)
    implicit none
    integer                  :: Nx,Ny
    integer                  :: Lx,Ly
    logical                  :: Xpbc,Ypbc
    real(8),allocatable      :: kxgrid(:),kygrid(:),kxgrid_(:),kygrid_(:)    
    real(8),allocatable      :: epsik(:)
    real(8),allocatable      :: Hlat(:,:)
    integer                 ::ik,Nx_,Ny_,Nlat_,Nlat_sub,Lk_,Lk
    integer                  :: ilat,jlat,ilat_,jlat_
    integer                  :: ix_,iy_,ix_sub,iy_sub,ilat_sub,jlat_sub

    integer,allocatable      :: map_lattice_new2old_(:,:,:)
    integer,allocatable      :: map_lattice_new2old(:,:)
    integer,allocatable      :: map_lattice_old2new_sites(:),map_lattice_old2new_sublat(:)
    real(8),allocatable      :: thop_lattice(:,:,:,:)
    type(vect2D),allocatable :: Rlat_(:),kVect_(:)
    complex(8),allocatable   :: Hk_lat(:,:,:)
    complex(8)               :: psi_ki,psi_kj
    !
    if(mod(Nx,Lx)/=0) stop "wrong X-periodicity"
    if(mod(Ny,Ly)/=0) stop "wrong Y-periodicity"


    !+- STEP-1: build the lattice hamiltonian using FT of the energy-momentum dispersion with the right boundary conditions.
    allocate(kxgrid(Nx),kygrid(Ny))    
    if(Xpbc) then
       kxgrid = linspace(0.d0,2.d0*pi,Nx,istart=.true.,iend=.false.)
    else
       kxgrid = linspace(0.d0,pi,Nx,istart=.false.,iend=.false.)
    end if
    if(Ypbc) then
       kygrid = linspace(0.d0,2.d0*pi,Ny,istart=.true.,iend=.false.)
    else
       kygrid = linspace(0.d0,pi,Ny,istart=.false.,iend=.false.)  
    end if
    Nlat=Nx*Ny; Lk = Nx*Ny; allocate(epsik(Lk))
    ik=0
    do ix=1,Nx
       do iy=1,Ny
          ik = ik + 1
          epsik(ik) = -2*ts*( dcos(kxgrid(ix)) + dcos(kygrid(iy)))
       end do
    end do
    allocate(Hlat(Nlat,Nlat))
    call k2latticeFT(kxgrid,Xpbc,kygrid,Ypbc,epsik,Hlat)        
    !+- STEP-2: build the hopping matrices for a lattice of different periodicity -+!
    !           1 <= Lx =< Nx  ; 1 <= Ly =< Ny            !
    Nx_=Nx/Lx;Ny_=Ny/Ly    
    Nlat_=Nx_*Ny_;Nlat_sub=Lx*Ly    
    allocate(Rlat_(Nlat_))
    allocate(map_lattice_new2old(Nlat_,Nlat_sub),map_lattice_new2old_(Nlat_,Lx,Ly))
    allocate(map_lattice_old2new_sites(Nlat),map_lattice_old2new_sublat(Nlat))
    allocate(thop_lattice(Nlat_,Nlat_,Nlat_sub,Nlat_sub))    
    do ix_=1,Nx_
       do iy_=1,Ny_
          ilat_=square_lattice_stride(ix_,iy_,Nx_)        
          Rlat_(ilat_)%x= dble(ix_)*dble(Lx)
          Rlat_(ilat_)%y= dble(iy_)*dble(Ly)
          do ix_sub=1,Lx
             do iy_sub=1,Ly
                ilat_sub = square_lattice_stride(ix_sub,iy_sub,Lx)              
                ix = (ix_-1)*Lx+ix_sub
                iy = (iy_-1)*Ly+iy_sub
                ilat = square_lattice_stride(ix,iy,Nx)
                map_lattice_new2old_(ilat_,ix_sub,iy_sub) = ilat
                map_lattice_new2old(ilat_,ilat_sub) = ilat
                map_lattice_old2new_sites(ilat) = ilat_
                map_lattice_old2new_sublat(ilat) = ilat_sub
             end do
          end do
       end do
    end do
    thop_lattice=0.d0
    do ilat_=1,Nlat_
       do jlat_=1,Nlat_
          do ilat_sub=1,Nlat_sub
             do jlat_sub=1,Nlat_sub
                ilat=map_lattice_new2old(ilat_,ilat_sub)
                jlat=map_lattice_new2old(jlat_,jlat_sub)
                thop_lattice(ilat_,jlat_,ilat_sub,jlat_sub) = Hlat(ilat,jlat)
             end do
          end do
       end do
    end do
    !
    !+- STEP 3 -+! 
    !   build-up the new momentum grid corresponding to the periodicity (Lx,Ly)
    Lk_ =Nx_*Ny_;allocate(kxgrid_(Nx_),kygrid_(Ny_),kVect_(Lk_))
    if(Xpbc) then
       kxgrid_ = linspace(0.d0,2.d0*pi/dble(Lx),Nx_,istart=.true.,iend=.false.)
    else
       kxgrid_ = linspace(0.d0,pi/dble(Lx),Nx_,istart=.false.,iend=.false.)
    end if
    if(Ypbc) then
       kygrid_ = linspace(0.d0,2.d0*pi/dble(Ly),Ny_,istart=.true.,iend=.false.)
    else
       kygrid_ = linspace(0.d0,pi/dble(Ly),Ny_,istart=.false.,iend=.false.)
    end if
    ik=0
    do ix_=1,Nx_
       do iy_=1,Ny_
          ik=ik+1
          kVect_(ik)%x=kxgrid_(ix_)
          kVect_(ik)%y=kygrid_(iy_)
       end do
    end do
    !+- STEP 4 -+! 
    !   For each (ilat_sub,jlat_sub) transform back thop_lattice(:,:,ilat_sub,jlat_sub) to momentum space
    if(allocated(Hk_lat)) deallocate(Hk_lat)
    allocate(Hk_lat(Nlat_sub,Nlat_sub,Lk_))
    do ik=1,Lk_
       Hk_lat(:,:,ik) = 0.d0
       do ilat_=1,Nlat_
          do jlat_=1,Nlat_
             !
             psi_ki=1.d0
             psi_kj=1.d0             
             if(Xpbc) then
                psi_ki=psi_ki*sqrt(1.d0/dble(Nx_))*exp(-xi*kVect_(ik)%x*Rlat_(ilat_)%x)
                psi_kj=psi_kj*sqrt(1.d0/dble(Nx_))*exp(-xi*kVect_(ik)%x*Rlat_(jlat_)%x)
             else
                psi_ki=psi_ki*sqrt(2.d0/dble(Nx_+1))*sin(kVect_(ik)%x*Rlat_(ilat_)%x)
                psi_kj=psi_kj*sqrt(2.d0/dble(Nx_+1))*sin(kVect_(ik)%x*Rlat_(jlat_)%x)
             end if
             if(Ypbc) then
                psi_ki=psi_ki*sqrt(1.d0/dble(Ny_))*exp(-xi*kVect_(ik)%y*Rlat_(ilat_)%y)
                psi_kj=psi_kj*sqrt(1.d0/dble(Ny_))*exp(-xi*kVect_(ik)%y*Rlat_(jlat_)%y)
             else
                psi_ki=psi_ki*sqrt(2.d0/dble(Ny_+1))*sin(kVect_(ik)%y*Rlat_(ilat_)%y)
                psi_kj=psi_kj*sqrt(2.d0/dble(Ny_+1))*sin(kVect_(ik)%y*Rlat_(jlat_)%y)
             end if
             Hk_lat(:,:,ik) = Hk_lat(:,:,ik) + thop_lattice(ilat_,jlat_,:,:)*psi_kj*conjg(psi_ki)
          end do
       end do
    end do

    unit=free_unit()
    open(unit,file="hk_symmetric_points")
    Hk_symm=.true.
    do ik=1,Lk_
       do ilat_sub=1,Nlat_sub
          do jlat_sub=1,Nlat_sub
             if(Hk(ilat_sub,jlat_sub,ik) /= Hk(jlat_sub,ilat_sub,ik)) Hk_symm(ik)=.false.
          end do
       end do
       if(Hk_symm(ik)) write(unit,'(4(F18.10))') kVect_(ik)%x,kVect_(ik)%y
    end do
  end subroutine get_Hk_2dsquare



  subroutine k2latticeFT(kx_grid,Xpbc,ky_grid,Ypbc,epsik,Hlat) 
    real(8)                     :: kx_grid(:),ky_grid(:)
    logical                     :: Xpbc,Ypbc
    real(8)                     :: epsik(size(kx_grid)*size(ky_grid))
    real(8)                     :: Hlat(size(kx_grid)*size(ky_grid),size(kx_grid)*size(ky_grid))
    !
    type(vect2D),allocatable   :: kgrid(:),Rlat(:)
    integer                     :: Nx,Ny
    !    complex(8),dimension(Nx,Ny) :: Hk
    !type(vect2D),dimension(:),allocatable   :: 

    real(8) :: arg,test
    complex(8) :: psi_ki,psi_kj
    integer :: i,j,ilat,jlat,ill,ill_
    integer :: Lk,Nlat    
    !
    integer,allocatable :: stride_ill2ij(:)
    real(8),allocatable :: Hlat_tmp(:,:)
    !
    Nx=size(kx_grid);Ny=size(ky_grid);Nlat = Nx*Ny;Lk=Nlat
    allocate(Rlat(Nlat),kgrid(Lk))
    !+- build lattice vectors -+!
    ik=0
    do i=1,Nx
       do j=1,Ny
          ik=ik+1
          kgrid(ik)%x=kx_grid(i)
          kgrid(ik)%y=ky_grid(j)          
          ilat=(j-1)*Nx+i
          Rlat(ilat)%x = dble(i)
          Rlat(ilat)%y = dble(j)
       end do
    end do
    !+-------------------------+!
    allocate(Hlat_tmp(Nlat,Nlat))
    write(LOGfile,*) "Building full Real-space hamiltonian from Fourier transform."
    if(mpiID==0)call start_timer
    do ill=1+mpiID,Nlat*Nlat,mpiSIZE
       ilat=mod(ill-1,Nlat)+1
       jlat=(ill-1)/Nlat+1
       Hlat_tmp(ilat,jlat) = 0.d0
       do ik=1,Lk
          !
          psi_ki=1.d0
          psi_kj=1.d0             
          if(Xpbc) then
             psi_ki=psi_ki*sqrt(1.d0/dble(Nx))*exp(xi*kgrid(ik)%x*Rlat(ilat)%x)
             psi_kj=psi_kj*sqrt(1.d0/dble(Nx))*exp(xi*kgrid(ik)%x*Rlat(jlat)%x)
          else
             psi_ki=psi_ki*sqrt(2.d0/dble(Nx+1))*sin(kgrid(ik)%x*Rlat(ilat)%x)
             psi_kj=psi_kj*sqrt(2.d0/dble(Nx+1))*sin(kgrid(ik)%x*Rlat(jlat)%x)
          end if
          if(Ypbc) then
             psi_ki=psi_ki*sqrt(1.d0/dble(Ny))*exp(xi*kgrid(ik)%y*Rlat(ilat)%y)
             psi_kj=psi_kj*sqrt(1.d0/dble(Ny))*exp(xi*kgrid(ik)%y*Rlat(jlat)%y)
          else
             psi_ki=psi_ki*sqrt(2.d0/dble(Ny+1))*sin(kgrid(ik)%y*Rlat(ilat)%y)
             psi_kj=psi_kj*sqrt(2.d0/dble(Ny+1))*sin(kgrid(ik)%y*Rlat(jlat)%y)
          end if
          !
          Hlat_tmp(ilat,jlat) = Hlat_tmp(ilat,jlat) + epsik(ik)*psi_ki*conjg(psi_kj)
       end do
       !if(mpiID==0) write(*,*) ill,Nlat*Nlat
       if(mpiID==0)call eta(ill,Nlat*Nlat,unit=LOGfile)
    end do
    if(mpiId==0) call stop_timer
    call MPI_ALLREDUCE(Hlat_tmp,Hlat,Nlat*Nlat,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine k2latticeFT


  function square_lattice_stride(ix,iy,Ncol) result(ilat) 
    integer :: ix,iy,Ncol 
    integer :: ilat
    ilat = (iy-1)*Ncol + ix
    return
  end function square_lattice_stride





!!!!!!!!! OBSOLETE !!!!!!!!
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
    Hk_symm_=.true.
    do ik=1,Lk
       do i=1,Nlat
          do j=1,Nlat
             if(Hk(i,j,ik) /= Hk(j,i,ik)) Hk_symm_(ik)=.false.
          end do
       end do
       if(Hk_symm_(ik)) write(unit,'(4(F18.10))') k_grid(ik)
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





  !<DEBUG
  ! Ncol = Xperiod
  ! Nrow = Yperiod
  ! Nlat = Nrow*Ncol
  ! pbc_row = .false.
  ! pbc_col = .false.
  ! unit=free_unit()
  ! open(unit,file='k_grid.lattice')
  ! Lk=N_Xperiod
  ! allocate(wt_(Lk),epsik_(Lk),k_grid(Lk),hk_symm_(Lk))
  ! do ik=1,Lk
  !    k_grid(ik) = 2*pi/dble(Lk)*dble(ik-1)
  !    epsik_(ik) = 2.d0*cos(k_grid(ik))
  !    wt_(ik) = 1.d0/dble(Lk)
  !    write(unit,'(6(F18.10))') dble(ik),k_grid(ik),epsik_(ik),wt_(ik)
  ! end do
  ! close(unit)
  ! call get_k_hamiltonian_stripe(Nrow,Ncol,pbc_row,pbc_col,k_grid)  
  ! do ik=1,Lk
  !    write(401,*)
  !    write(401,*) ik
  !    write(401,*)
  !    write(*,*) ik,Lk,N_Xperiod,N_Yperiod
  !    do ilat=1,Nlat
  !       write(401,'(40(F7.2))') dreal(Hk(ilat,:,ik))
  !    end do
  ! end do
  !DEBUG>
