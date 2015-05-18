  !###############################################################!
  ! DRIVER FOR A ONE/TWO-BAND HUBBARD MODEL CALCULATION IN A SLAB GEOMETRY
  ! WITH APPLIED ELECTRIC FIELD.
  ! AUTHORS: G.Mazza & A.Amaricci @ SISSA - 2014
  !###############################################################!
!  include "SQUARE_LATTICE.f90"
  program ed_slab
    USE DMFT_ED
    USE SCIFOR
    USE DMFT_TOOLS
    !
    USE DMFT_VECTORS
 !   USE SQUARE_LATTICE
    !
#ifdef _MPI_INEQ
  USE MPI
#endif
    implicit none
    !
    complex(8),allocatable,dimension(:,:) :: Smats,Sreal !self_energies
    complex(8),allocatable,dimension(:,:) :: Smats_,Sreal_ !self_energies
    complex(8),allocatable,dimension(:,:) :: Gmats,Greal !local green's functions
    complex(8),allocatable,dimension(:,:) :: Gmats_,Greal_ !local green's functions
    complex(8),allocatable,dimension(:,:) :: Delta      
    complex(8),allocatable,dimension(:,:) :: Delta_      
    !
    complex(8),allocatable                :: Hk(:,:,:)
    complex(8),allocatable                :: Hloc(:,:,:,:,:)         ![Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),allocatable                :: Sigma_mats(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable                :: Sigma_real(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable                :: Gloc_mats(:,:,:,:,:,:)  ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable                :: Gloc_real(:,:,:,:,:,:)  ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable                :: Delta_bath(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]

    real(8),allocatable                   :: nsite(:,:),dsite(:,:)    ![Nlat][Norb]
    real(8),allocatable                   :: esite(:)                 ![Nlat]

    complex(8),allocatable                :: Hloc_(:,:,:,:,:)         ![Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),allocatable                :: Sigma_mats_(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable                :: Sigma_real_(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable                :: Gloc_mats_(:,:,:,:,:,:)  ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable                :: Gloc_real_(:,:,:,:,:,:)  ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable                :: Delta_bath_(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    real(8),allocatable                   :: nsite_(:,:),dsite_(:,:)  ![Nlat][Norb]
    real(8),allocatable                   :: esite_(:)                ![Nlat]

    !
    real(8),allocatable,dimension(:)      :: elocal
    real(8),allocatable,dimension(:)      :: elocal_
    real(8)                               :: crystal_field,Vhyb
    real(8),allocatable,dimension(:,:,:)  :: bath,bath_old
    real(8),allocatable,dimension(:,:,:)  :: bath_,bath_old_
    real(8),allocatable,dimension(:,:)    :: tmpBath
    !<DEBUG  comment
    real(8),allocatable,dimension(:,:)    :: nii_,dii_
    real(8),allocatable,dimension(:)      :: eii_  
    !>DEBUG
    logical                               :: converged
    real(8)                               :: DeltaV
    real(8)                               :: r,de,ts
    real(8)                               :: wmixing
    real(8),allocatable,dimension(:)      :: epsik,wt
    real(8),allocatable,dimension(:)      :: epsik_embedd
    logical,allocatable,dimension(:)      :: hk_symm
    real(8),allocatable,dimension(:)      :: docc_check,docc_check_
    real(8)                               :: ncheck,xmu_,kz
    integer                               :: i,is,js,iloop
    integer                               :: Nb(2),Nx,Lk,ik,ix,iy,Nkz,Nk
    integer                               :: ilat,i_ind,isymm,check_maps,iorb,jorb,ispin,iorb_,jorb_
    integer,dimension(2)                  :: unit
    integer                               :: tmp_unit,i1,i2
    integer                               :: eloc_unit
    character(len=5)                      :: tmp_suffix
    logical                               :: symmetry_flag,mu_ph,embedd
    type(vect2D),allocatable,dimension(:) :: kVect

    !<- supply the departure of SQUARE_LATTICE.f90
    type(vect2D)                          :: ai,aj,bi,bj
    integer,dimension(:),allocatable      :: ik2ix,ik2iy
    integer,dimension(:,:),allocatable    :: kindex
    real(8)                               :: peso,kx,ky
    real(8),dimension(:),allocatable      :: kxgrid,kygrid
    !>
    
    !<DEBUG 
    real(8),dimension(:),allocatable          :: wm,wr,eii
    real(8),dimension(:,:),allocatable        :: nii,dii ![Nlat][Norb]/[4]
    real(8)                                   :: Eint,Ekin,Epot,Eout(2)
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_tmp
    complex(8),dimension(:,:,:),allocatable   :: Hloc_tmp
    complex(8),dimension(:,:),allocatable     :: H_loc
    integer                                   :: unit_
    character(len=20)                         :: suffix 
    !>DEBUG


#ifdef _MPI_INEQ
    ! START MPI !
    call MPI_INIT(mpiERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
    write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

    !+--------------------+!
    ! READ INPUT VARIABLES !
    !+--------------------+!
    call parse_input_variable(Nx,"Nx","inputRDMFT.in",default=10)
    call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=1.d0)
    call parse_input_variable(DeltaV,"BIAS","inputRDMFT.in",default=0.d0)
    call parse_input_variable(crystal_field,"CFIELD","inputRDMFT.in",default=0.d0)
    call parse_input_variable(ts,"TS","inputRDMFT.in",default=0.5d0)
    call parse_input_variable(Vhyb,"Vhyb","inputRDMFT.in",default=0.d0)
    call parse_input_variable(symmetry_flag,"SYMMETRY_FLAG","inputRDMFT.in",default=.false.)
    call parse_input_variable(mu_ph,"MU_PH","inputRDMFT.in",default=.false.)
    call parse_input_variable(embedd,"EMBEDD","inputRDMFT.in",default=.false.)
    !
    call ed_read_input("inputRDMFT.in")
    call set_store_size(1024)

    !<TMP TEST: RESTORE PH SYMM IN THE FULLY POLARIZED INSULATING PHASE 
    if(mu_ph.and.nread==0.d0) xmu=-Uloc(1)-crystal_field*0.5d0
    !>

    !+-----------------------------+!
    !+- BUILD LATTICE HAMILTONIAN -+!
    !+-----------------------------+!
    Nlat = Nside
    !
    
    ! Nk = Nx/2+1
    ! Lk = Nk*(Nk+1)/2    
    ! ai=1.d0*Xver       ; aj=1.d0*Yver
    ! bi=(pi2/1.d0)*Xver ; bj=(pi2/1.d0)*Yver    
    ! allocate(epsik(Lk),wt(Lk),kVect(Lk),hk_symm(Lk))
    ! allocate(kindex(0:Nk,0:Nk))
    ! allocate(ik2ix(Lk),ik2iy(Lk))    
    ! ik=0
    ! do ix=0,Nx/2
    !    do iy=0,ix          
    !       ik=ik+1
    !       Kx=dble(ix)/dble(Nx)
    !       Ky=dble(iy)/dble(Nx)
    !       ik2ix(ik)=ix
    !       ik2iy(ik)=iy
    !       kVect(ik)=Kx*bi + Ky*bj - pi*Vone 
    !       kindex(ix,iy)=ik
    !       if (ix==0) then
    !          peso=1.d0       !center
    !       elseif(ix==Nx/2) then
    !          if (iy==0) then
    !             peso=2.d0    ! point (pi,0)
    !          elseif(iy==Nx/2) then 
    !             peso=1.d0    ! corner 
    !          else
    !             peso=4.d0    ! border
    !          endif
    !       else
    !          if (iy==ix) then
    !             peso=4.d0    ! diagonal
    !          elseif (iy==0) then
    !             peso=4.d0    ! x-axis
    !          else
    !             peso=8.d0    ! all other points
    !          endif
    !       endif
    !       wt(ik)=peso/dble(Nx**2)
    !    enddo
    ! enddo    
    ! do ik=1,Lk
    !    ix=ik2ix(ik)
    !    iy=ik2iy(ik)
    !    epsik(ik)=-2*ts*( dcos(kVect(ik)%x) + dcos(kVect(ik)%y) )
    ! enddo
    
    allocate(kxgrid(Nx),kygrid(Nx))
    kxgrid = linspace(0.d0,pi,Nx)
    kygrid = linspace(0.d0,pi,Nx)
    Lk = Nx*Nx
    allocate(epsik(Lk),wt(Lk),kVect(Lk),hk_symm(Lk))
    wt = 1.d0/dble(Lk)    
    ik=0
    do ix=1,Nx
       do iy=1,Nx
          ik = ik + 1
          epsik(ik) = -2*ts*( dcos(kxgrid(ix)) + dcos(kygrid(iy)) )
       end do
    end do
    
    call get_free_dos(epsik,wt)    
    Nkz=20
    allocate(epsik_embedd(Nkz))
    kz=0.d0
    do ik=1,Nkz
       kz = kz + pi*dble(ik)/dble(Nkz+1)
       epsik_embedd(ik) = -2.d0*ts*cos(kz)
    end do
    !
    call get_k_hamiltonian_ft(kVect)
    !

    !+-------------------------+!
    !+- BUILD LATTICE DETAILS -+!
    !+-------------------------+!
    allocate(elocal(Nlat))
    eloc_unit=free_unit()
    if(Nlat>1)then
       do i=1,Nlat
          elocal(i)= DeltaV*0.5d0 - DeltaV/dble(Nlat-1)*dble(i-1) 
       end do
    end if
    close(eloc_unit)
    call store_data("bias_used.data",elocal,(/(dble(i),i=1,Nlat)/))


    !+----------------------------------+!  
    !+- ALLOCATE GF & INITIALIZE BATHS -+!
    !+----------------------------------+!
    !Lattice sites

    ! Matsubara and Real freq
    allocate(wm(Lmats),wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
    ! Lattice sites baths
    Nb=get_bath_size()
    allocate(bath(Nlat,Nb(1),Nb(2)))
    allocate(bath_old(Nlat,Nb(1),Nb(2)))
    !
    allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
    allocate(Sigma_mats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(Sigma_real(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    allocate(Gloc_mats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(Gloc_real(Nlat,Nspin,Nspin,Norb,Norb,Lreal)) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    allocate(Delta_bath(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(nsite(Nlat,Norb),dsite(Nlat,Norb))            ![Nlat][Norb]
    allocate(esite(Nlat))                                  ![Nlat]
    allocate(docc_check(Nlat))
    !


    ! Observables
    allocate(nii(Nlat,Norb))
    allocate(dii(Nlat,Norb))
    allocate(eii(Nlat))
    ! Self energies
    allocate(Smats(Nlat,Lmats))
    allocate(Sreal(Nlat,Lreal))
    ! Green Function
    allocate(Gmats(Nlat,Lmats))
    allocate(Greal(Nlat,Lreal))
    ! Impurity-bath hybritizations
    allocate(Delta(Nlat,Lmats))
    ! 
    call  ed_init_solver_lattice(bath)

    !Independent sites
    call get_indep_sites


    allocate(Hloc_(Nindep,Nspin,Nspin,Norb,Norb))
    allocate(Sigma_mats_(Nindep,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(Sigma_real_(Nindep,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    allocate(Gloc_mats_(Nindep,Nspin,Nspin,Norb,Norb,Lmats))  ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(Gloc_real_(Nindep,Nspin,Nspin,Norb,Norb,Lreal))  ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    allocate(Delta_bath_(Nindep,Nspin,Nspin,Norb,Norb,Lmats)) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    allocate(nsite_(Nindep,Norb),dsite_(Nindep,Norb))         ![Nlat][Norb]
    allocate(esite_(Nindep))                                  ![Nlat]
    allocate(docc_check_(Nindep))


    !
    allocate(nii_(Nindep,Norb))
    allocate(dii_(Nindep,Norb))
    allocate(eii_(Nindep))
    !
    allocate(bath_(Nindep,Nb(1),Nb(2)))
    allocate(bath_old_(Nindep,Nb(1),Nb(2)))
    !
    allocate(Smats_(Nindep,Lmats))
    allocate(Sreal_(Nindep,Lreal))
    !
    allocate(Gmats_(Nindep,Lmats))
    allocate(Greal_(Nindep,Lreal))
    !
    allocate(Delta_(Nindep,Lmats))
    !
    allocate(elocal_(Nindep))
    !

    Hloc_=0.d0
    do i_ind=1,Nindep
       bath_(i_ind,:,:) = bath(indep_list(i_ind),:,:)
       elocal_(i_ind) = elocal(indep_list(i_ind))
       do ispin=1,Nspin
          do iorb=1,Norb
             Hloc_(i_ind,ispin,ispin,iorb,iorb) = elocal(indep_list(i_ind)) 
          end do
          if(Norb.eq.2) then
             Hloc_(i_ind,ispin,ispin,1,1) = Hloc_(i_ind,ispin,ispin,1,1) + crystal_field*0.5d0
             Hloc_(i_ind,ispin,ispin,2,2) = Hloc_(i_ind,ispin,ispin,2,2) - crystal_field*0.5d0
          end if
       end do
    end do

    Hloc=0.d0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             Hloc(ilat,ispin,ispin,iorb,iorb) = elocal(ilat) 
          end do
          if(Norb.eq.2) then
             Hloc(ilat,ispin,ispin,1,1) = Hloc(ilat,ispin,ispin,1,1) + crystal_field*0.5d0
             Hloc(ilat,ispin,ispin,2,2) = Hloc(ilat,ispin,ispin,2,2) - crystal_field*0.5d0
          end if
       end do
    end do


    !<DEBUG re-add the local part of the non-interacting Hamiltonian to H(k)
    !Hloc is defined as [Nlat][Nspin][Nspin][Norb][Norb]
    !H(k) is defiend as [Nlat*Norb][Nlat*Norb][Lk]
    !thus we have to remap Hloc to the size of H(k)
    allocate(H_loc(Nlat*Norb,Nlat*Norb))
    H_loc=zero
    do ilat=1,Nlat
       do iorb=1,Norb
          do jorb=1,Norb
             is = iorb + (ilat-1)*Norb
             js = jorb + (ilat-1)*Norb
             H_loc(is,js) = Hloc(ilat,1,1,iorb,jorb)
          enddo
       enddo
    enddo
    !RESTORE THE CORRECT H(k) INCLUDING THE LOCAL PART
    forall(ik=1:Lk)Hk(:,:,ik) = Hk(:,:,ik) + H_loc
    !>DEBUG



    !<DEBUG TO BE REMOVED TEST EKIN
    ! allocate(Sigma_tmp(Nlat,Norb,Norb,Lmats))
    ! allocate(Hloc_tmp(Nlat,Norb,Norb))
    ! Sigma_mats=zero
    ! esite=0d0
    ! Eint=0.d0
    ! Ekin=0.d0
    ! Epot=0.d0
    ! Epot=sum(esite)/dble(Nlat)
    ! !do ispin=1,Nspin
    ! Eout=ed_kinetic_energy_lattice(Hk,wt,Sigma_mats)
    ! !end do
    ! Ekin=Eout(1)!sum(Eout)
    ! Eint=Ekin+Epot
    ! unit_=free_unit()
    ! open(unit_,file='internal_energy.data')
    ! write(unit_,'(10(F18.10))') Eint,Ekin,Epot,Eint-xmu*sum(nii)/dble(Nlat)
    ! close(unit_)
    ! !if we substract again H_loc from Hk the kinetic energy 
    ! !as calculated below lack of the term Tr(Hloc.Gk)=Eloc (Eknot)
    ! !thus the two numbers have to be different.
    ! forall(ik=1:Lk)Hk(:,:,ik) = Hk(:,:,ik) - H_loc
    ! Ekin=0.d0
    ! Epot=0.d0
    ! Epot=sum(esite)/dble(Nlat)
    ! do ispin=1,Nspin
    !    Sigma_tmp=Sigma_mats(:,ispin,ispin,:,:,:)
    !    Hloc_tmp=Hloc(:,ispin,ispin,:,:)
    !    Ekin=Ekin+ kinetic_energy_lattice_OLD(Hk,Wt,Hloc_tmp,Sigma_tmp)
    ! end do
    ! Eint=Ekin+Epot
    ! unit_=free_unit()
    ! open(unit_,file='internal_energy2.data')
    ! write(unit_,'(10(F18.10))') Eint,Ekin,Epot,Eint-xmu*sum(nii)/dble(Nlat)
    ! close(unit_)
    !>DEBUG

    !+-------------+!
    !+- DMFT LOOP -+!
    !+-------------+!
    if(nread==0.d0.and.xmu/=0.d0) symmetry_flag=.false.
    iloop=0 ; converged=.false.
    do while(.not.converged.AND.iloop<nloop) 
       iloop=iloop+1
       if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   

       if(symmetry_flag.and.nread==0.d0) then
          ! Solve site-dependent impurities problems using lr ph symmetry !
          ! Save old baths
          bath_old_=bath_        
          ! Solve impurities
          !<DEBUG  comment
          call ed_solve_lattice(bath_,Hloc_)
          nsite_ = ed_get_dens_lattice(Nindep)
          dsite_ = ed_get_docc_lattice(Nindep)
          esite_ = ed_get_epot_lattice(Nindep)
          call ed_get_sigma_matsubara_lattice(Smats_,Nindep)
          call ed_get_sigma_real_lattice(Sreal_,Nindep)
          !>DEBUG
          ! Extend solution to full lattice
          do ilat=1,Nlat
             i_ind=map_lat2ind(ilat)
             if(abs(DeltaV).lt.1.d-8) then
                Sigma_mats(ilat,:,:,:,:,:) = Sigma_mats_(i_ind,:,:,:,:,:)
                Sigma_real(ilat,:,:,:,:,:) = Sigma_real_(i_ind,:,:,:,:,:)
             else
                if(i_ind /= ilat) then
                   select case(Norb)
                   case default
                      !Norb=1
                      Sigma_mats(ilat,:,:,:,:,:) = -conjg(Sigma_mats_(i_ind,:,:,:,:,:))  
                      Sigma_real(ilat,:,:,:,:,:) = -conjg(Sigma_real_(i_ind,:,:,:,:,:))  
                   case(2)
                      ! G_ab_ilat(iw) = -[G_{a}{b}_iind[iw]]*  {a} = not_a
                      do iorb=1,Norb
                         do jorb=1,Norb
                            iorb_=3-iorb
                            jorb_=3-jorb
                            Sigma_mats(ilat,:,:,iorb,jorb,:) = -conjg(Sigma_mats_(i_ind,:,:,iorb_,jorb_,:))  
                            Sigma_real(ilat,:,:,iorb,jorb,:) = -conjg(Sigma_real_(i_ind,:,:,iorb_,jorb_,:))  
                         end do
                      end do
                   end select
                else
                   Sigma_mats(ilat,:,:,:,:,:) = Sigma_mats_(i_ind,:,:,:,:,:)
                   Sigma_real(ilat,:,:,:,:,:) = Sigma_real_(i_ind,:,:,:,:,:)  
                end if
             end if
          end do
          ! Compute local GFs
          !+-(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)-+!
          !call rdmft_get_gloc_mb(Hk,wt,Hloc,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,hk_symm)
          call ed_get_gloc_lattice(Hk,Wt,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,1,hk_symm=hk_symm)
          ! Dump GLoc onto  independent sites GLoc
          do i_ind=1,Nindep
             ilat=indep_list(i_ind)
             Gloc_mats_(i_ind,:,:,:,:,:)=Gloc_mats(ilat,:,:,:,:,:)
             Gloc_real_(i_ind,:,:,:,:,:)=Gloc_real(ilat,:,:,:,:,:)
          end do
          ! Get Weiss fields
          !<DEBUG  comment
          !call rdmft_get_weiss_field_mb(Nindep,Gloc_mats_,Sigma_mats_,Delta_bath_,Hloc_)        
          !Nsites,Gloc,Smats,Weiss,Hloc
          call ed_get_weiss_lattice(Nindep,Gloc_mats_,Sigma_mats_,Delta_bath_,Hloc_)        
          !>DEBUG
          ! Fit baths
          !<DEBUG  comment
          ! call ed_fit_bath_sites_mb(bath_,Delta_bath_,Hloc_)        
          call ed_chi2_fitgf_lattice(bath_,Delta_bath_,Hloc_)
          !>DEBUG
          bath_=wmixing*bath_ + (1.d0-wmixing)*bath_old_
          ! constrain bath ph symmetry
          ! if(rdmft_phsym)then
          !    do i_ind=1,Nindep
          !       call ph_symmetrize_bath(bath_(i_ind,:,:))
          !    enddo
          ! endif
          ! check convergence
          docc_check_=sum(dsite_,2)
          if(mpiID==0) converged = check_convergence_local(docc_check_,dmft_error,Nsuccess,nloop,id=0,file="error.err")
#ifdef _MPI_INEQ
          call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
          ! Reshuffle baths to full slab and compute quantities at convergency
          if(converged) then            
             allocate(tmpBath(Nb(1),Nb(2)))
             do i_ind=1,Nindep              
                ilat=map_ind2lat(i_ind,1)
                bath(ilat,:,:) = bath_(i_ind,:,:)
                ilat=map_ind2lat(i_ind,2)
                bath(ilat,:,:) = bath_(i_ind,:,:)
                if(abs(DeltaV).lt.1.d-8) then
                   bath(ilat,:,:)=bath_(i_ind,:,:)
                else
                   tmpBath=bath_(i_ind,:,:)
                   call ph_trans_bath(tmpBath)
                   bath(ilat,:,:)=tmpBath
                end if
             end do
             deallocate(tmpBath)
             !<DEBUG  comment
             ! call ed_solve_impurity_sites_mb(bath,Hloc,Sigma_mats,Sigma_real,nsite,dsite,esite)
             call ed_solve_lattice(bath,Hloc)
             nsite = ed_get_dens_lattice(Nlat)
             dsite = ed_get_docc_lattice(Nlat)
             esite = ed_get_epot_lattice(Nlat)
             call ed_get_sigma_matsubara_lattice(Sigma_mats,Nlat)
             call ed_get_sigma_real_lattice(Sigma_real,Nlat)
             call ed_get_gloc_lattice(Hk,Wt,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,1,hk_symm=hk_symm)
             !call ed_get_gloc_lattice(Hk,Wt,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,hk_symm=hk_symm)    
             !>DEBUG
          end if



       else



          ! Solve site-dependent impurities problems without symmetries !        
          ! Save old baths
          bath_old=bath  
          xmu_=xmu
          ! Solve impurities
          !<DEBUG  comment
          call ed_solve_lattice(bath,Hloc)
          nsite = ed_get_dens_lattice(Nlat)
          dsite = ed_get_docc_lattice(Nlat)
          esite = ed_get_epot_lattice(Nlat)
          call ed_get_sigma_matsubara_lattice(Sigma_mats,Nlat)
          call ed_get_sigma_real_lattice(Sigma_real,Nlat)
          !>DEBUG
          !
          call kill_mixed_sigma(Sigma_mats,Sigma_real)
          ! Compute local GF
          ! if(embedd) then
          !    call rdmft_get_gloc_embedd(Hk,wt,Hloc,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,epsik_embedd,hk_symm=hk_symm)
          ! else
          !<DEBUG  comment
          !call ed_get_gloc_lattice(Hk,Wt,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,hk_symm=hk_symm)    
          call ed_get_gloc_lattice(Hk,Wt,Gloc_mats,Gloc_real,Sigma_mats,Sigma_real,1,hk_symm=hk_symm)
          !>DEBUG
          ! end if
          !<DEBUG  comment
          ! Compute Weiss fields (or hybridization functions)
          !call ed_get_weiss_lattice(Nindep,Gloc_mats_,Sigma_mats_,Delta_bath_,Hloc_)        
          call ed_get_weiss_lattice(Nlat,Gloc_mats,Sigma_mats,Delta_bath,Hloc)
          !>DEBUG
          !<DEBUG  comment
          ! fit baths and mix result with old baths
          ! call ed_fit_bath_sites_mb(bath,Delta_bath,Hloc)
          call ed_chi2_fitgf_lattice(bath,Delta_bath,Hloc)
          !>DEBUG
          bath=wmixing*bath + (1.d0-wmixing)*bath_old
          docc_check=sum(dsite,2)
          if(mpiID==0) converged = check_convergence_local(docc_check,dmft_error,Nsuccess,nloop,id=0,file="error.err")
          ncheck=0.d0
          do ilat=1,Nlat
             ncheck = ncheck + sum(nsite(ilat,:))/dble(Nlat)
          end do
          if(mpiID==0) then
             if(nread/=0.d0) call search_chemical_potential(xmu,ncheck,converged)
          end if
#ifdef _MPI_INEQ
          call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
          call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
          xmu=wmixing*xmu + (1.d0-wmixing)*xmu_
       end if
       call print_out_mb(converged)
       if(mpiID==0) call end_loop()
    enddo
  contains
    !
    subroutine kill_mixed_sigma(Smats,Sreal)
      complex(8),intent(inout) :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)  ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
      complex(8),intent(inout) :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)  ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
      integer :: iorb,jorb
      do iorb=1,Norb
         do jorb=1,Norb
            if(iorb/=jorb) then
               Smats(:,:,:,iorb,jorb,:) = zero
               Sreal(:,:,:,iorb,jorb,:) = zero
            end if
         end do
      end do
    end subroutine kill_mixed_sigma



    subroutine print_out_mb(converged)
      !
      logical                        :: converged
      !
      integer                        :: i,j,is,row,col,unit
      real(8),dimension(Norb)        :: nimp,docc
      real(8)                        :: docc_ave
      real(8) :: ccdw
      real(8),dimension(Nlat)        :: cdwii
      real(8),dimension(Nlat,Norb)   :: rii,sii,zii
      character(len=4)               :: loop
      real(8)                        :: Eint,Ekin,Epot,Eout(2)
      complex(8),dimension(Nlat,Norb,Norb,Lmats) :: Sigma_tmp
      complex(8),dimension(Nlat,Norb,Norb) :: Hloc_tmp
      integer                        :: iorb,jorb,ispin
      character(len=20)                              :: suffix


      if(mpiID==0)then


         if(symmetry_flag) then
            nimp=0.d0
            docc=0.d0          
            do ilat=1,Nindep
               nimp(:) = nimp(:) + nsite_(ilat,:)/dble(Nindep)
               docc(:) = docc(:) + dsite_(ilat,:)/dble(Nindep)
            end do
            docc_ave=0.d0
            do iorb=1,Norb
               suffix="_l"//reg(txtfy(iorb))
               call splot("n"//reg(suffix)//"VSiloop.data",iloop,nimp(iorb),append=.true.)
               call splot("docc"//reg(suffix)//"VSiloop.data",iloop,docc(iorb),append=.true.)
               call store_data("n"//reg(suffix)//"VSisite.data",nsite_(:,iorb),(/(dble(i),i=1,Nindep)/))
               call store_data("docc"//reg(suffix)//"VSisite.data",dsite_(:,iorb),(/(dble(i),i=1,Nindep)/))
               docc_ave = docc_ave + docc(iorb)/dble(Norb)
            end do
            call splot("docc_aveVSiloop.data",iloop,docc_ave,append=.true.)
            call store_data("eneLocVSisite.data",esite_,(/(dble(i),i=1,Nindep)/))
         else
            nimp=0.d0
            docc=0.d0
            do ilat=1,Nlat
               nimp(:) = nimp(:) + nsite(ilat,:)/dble(Nlat)
               docc(:) = docc(:) + dsite(ilat,:)/dble(Nlat)
            end do
            docc_ave=0.d0
            do iorb=1,Norb
               suffix="_l"//reg(txtfy(iorb))
               call splot("n"//reg(suffix)//"VSiloop.data",iloop,nimp(iorb),append=.true.)
               call splot("docc"//reg(suffix)//"VSiloop.data",iloop,docc(iorb),append=.true.)
               call store_data("n"//reg(suffix)//"VSisite.data",nsite(:,iorb),(/(dble(i),i=1,Nlat)/))
               call store_data("docc"//reg(suffix)//"VSisite.data",dsite(:,iorb),(/(dble(i),i=1,Nlat)/))
               docc_ave = docc_ave + docc(iorb)/dble(Norb)
            end do
            call splot("docc_aveVSiloop.data",iloop,docc_ave,append=.true.)
            call store_data("eneLocVSisite.data",esite_,(/(dble(i),i=1,Nlat)/))
         end if

         !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
         if(converged)then
            nimp=0.d0
            docc=0.d0
            do ilat=1,Nlat
               nimp(:) = nimp(:) + nsite(ilat,:)/dble(Nlat)
               docc(:) = docc(:) + dsite(ilat,:)/dble(Nlat)
            end do
            do iorb=1,Norb
               suffix="_l"//reg(txtfy(iorb))
               call store_data("n"//reg(suffix)//"VSisite.data",nsite(:,iorb),(/(dble(i),i=1,Nlat)/))
               call store_data("docc"//reg(suffix)//"VSisite.data",dsite(:,iorb),(/(dble(i),i=1,Nlat)/))
            end do
            call store_data("eneLocVSisite.data",esite,(/(dble(i),i=1,Nlat)/))
            !
            do iorb=1,Norb
               do jorb=1,Norb
                  suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
                  write(*,*) suffix
                  call store_data("LG_iw"//reg(suffix)//".data",Gloc_mats(1:Nlat,1,1,iorb,jorb,1:Lmats),wm(1:Lmats))
                  call store_data("LG_realw"//reg(suffix)//".data",Gloc_real(1:Nlat,1,1,iorb,jorb,1:Lreal),wr(1:Lreal))
                  call store_data("LSigma_iw"//reg(suffix)//".data",Sigma_mats(1:Nlat,1,1,iorb,iorb,1:Lmats),wm(1:Lmats))
                  call store_data("LSigma_real"//reg(suffix)//".data",Sigma_real(1:Nlat,1,1,iorb,jorb,1:Lreal),wr(1:Lreal))
               end do
            end do
            !
            do ilat=1,Nlat
               do iorb=1,Norb
                  zii(ilat,iorb)   = 1.d0/( 1.d0 + abs( dimag(Sigma_mats(ilat,1,1,iorb,iorb,1))/wm(1) ))
               end do
            enddo
            zii=abs(zii)
            !
            do iorb=1,Norb
               suffix="_l"//reg(txtfy(iorb))
               call store_data("zeta"//reg(suffix)//"VSisite.data",zii(:,iorb),(/(dble(i),i=1,Nlat)/))
            end do
            !
         end if
      end if


      if(converged) then       
         Eint=0.d0
         Ekin=0.d0
         Epot=0.d0
         Epot=sum(esite)/dble(Nlat)
         !do ispin=1,Nspin
         !Ekin=Ekin+ed_kinetic_energy_lattice(Hk,wt,Sigma_mats)
         Eout = ed_kinetic_energy_lattice(Hk,Wt,Sigma_mats)
         Ekin = Eout(1)
         !end do
         Eint=Ekin+Epot
         unit=free_unit()
         open(unit,file='internal_energy.data')
         write(unit,'(10(F18.10))') Eint,Ekin,Epot,Eint-xmu*sum(nii)/dble(Nlat)
         close(unit)

         !<DEBUG use the old routine from GIacomo which re-add Hloc by hand to H(k)
         !       as we now work with a compelte H(k) (which include local terms)
         !       we just set Hloc_tmp to zero
         ! Ekin=0.d0
         ! Epot=0.d0
         ! Epot=sum(esite)/dble(Nlat)
         ! do ispin=1,Nspin
         !    Sigma_tmp=Sigma_mats(:,ispin,ispin,:,:,:)
         !    Hloc_tmp=0d0!Hloc(:,ispin,ispin,:,:)
         !    Ekin=Ekin+ kinetic_energy_lattice_OLD(Hk,Wt,Hloc_tmp,Sigma_tmp)
         ! end do
         ! Eint=Ekin+Epot
         ! unit=free_unit()
         ! open(unit,file='internal_energy2.data')
         ! write(unit,'(10(F18.10))') Eint,Ekin,Epot,Eint-xmu*sum(nii)/dble(Nlat)
         ! close(unit)
         !>DEBUG
      end if

    end subroutine print_out_mb





    subroutine get_indep_sites
      !+- Number of independent sites -+!
      Nindep=Nlat/2
      ! allocate maps !
      allocate(indep_list(Nindep),map_lat2ind(Nlat),map_ind2lat(Nindep,2))
      unit=free_units(2)
      open(unit(1),file='independent_sites.lattice')  
      do ilat=1,Nindep
         indep_list(ilat)=ilat
         write(unit(1),'(2(F18.10))') dble(indep_list(ilat)),0.d0
         map_ind2lat(ilat,1)=ilat
         map_ind2lat(ilat,2)=Nlat+1-ilat
         write(tmp_suffix,'(I4.4)') ilat
         open(unit(2),file='equivalents'//trim(tmp_suffix)//'.lattice')
         write(unit(2),'(2(F18.10))') dble(map_ind2lat(ilat,1)),0.d0
         write(unit(2),'(2(F18.10))') dble(map_ind2lat(ilat,2)),0.d0
         close(unit(2))
         map_lat2ind(ilat)=ilat
         map_lat2ind(Nlat+1-ilat)=ilat
      end do
      close(unit(1))
      ! check maps !
      do i_ind=1,Nindep
         do isymm=1,2
            check_maps=map_ind2lat(i_ind,isymm)
            if(i_ind /= map_lat2ind(check_maps)) stop "WRONG MAPS"
         end do
      end do

    end subroutine get_indep_sites



    !+------------+!
    !+- build Hk -+!  
    !+------------+!
    subroutine get_k_hamiltonian_slab(kVect)
      type(vect2D)          :: kVect(:)
      integer               :: i,jj,j,k,row,col,link(4),Lk,ik,No
      integer               :: unit,unitk
      integer               :: iorb,jorb,io,jo
      real(8)               :: expk
      type(vect2D),dimension(4) :: Rhop
      real(8),dimension(:,:,:),allocatable  :: thop
      !
      Lk=size(kVect)
      No=Nlat*Norb
      allocate(H0(No,No))
      allocate(Hk(No,No,Lk),thop(No,No,4))


      !+-----------------------------------------------------+!
      !+- GET real-space block hamiltonian H0 (aka H-IRRED) -+!
      !+-----------------------------------------------------+!
      H0=0.d0
      do ilat=1,Nlat-1
         do iorb=1,Norb
            do jorb=1,Norb

               !+- hop from ilat to ilat+1 -+!
               io=(ilat-1)*Norb + iorb
               jo=(ilat+1-1)*Norb + jorb
               if(iorb.eq.jorb) H0(io,jo)=-ts
               !+- possibility of nn hybridization terms -+!

               !+- hop from ilat+1 to ilat -+!
               io=(ilat+1-1)*Norb + iorb
               jo=(ilat-1)*Norb + jorb
               if(iorb.eq.jorb) H0(io,jo)=-ts
               !+- possibility of nn hybridization terms -+!

            end do
         end do
      end do

      !+------------------------------------------------------------------------+!
      !+- GET real-space connection between blocks hamiltonian (aka H-CONNECT) -+!
      !+------------------------------------------------------------------------+!
      thop=0.d0
      do ilat=1,Nlat
         do iorb=1,Norb
            do jorb=1,Norb
               io=(ilat-1)*Norb + iorb
               jo=(ilat-1)*Norb + jorb
               if(iorb.eq.jorb) thop(io,jo,:) = -ts
               !+- possibility of nn hybridization terms -+!
            end do
         end do
      end do
      !
      Rhop(1)%x=-1.d0
      Rhop(1)%y= 0.d0
      !
      Rhop(2)%x= 1.d0
      Rhop(2)%y= 0.d0
      !
      Rhop(3)%x= 0.d0
      Rhop(3)%y= 1.d0
      !
      Rhop(4)%x= 0.d0
      Rhop(4)%y=-1.d0

      !+------------------------------------------------------------------------+!
      !+- Fourier transform of the real-space connections and get the full Hk  -+!
      !+------------------------------------------------------------------------+!
      Hk=zero
      do ik=1,Lk
         do j=1,4
            expk = Rhop(j).dot.kVect(ik)
            Hk(:,:,ik) =  Hk(:,:,ik) + thop(:,:,j)*exp(xi*expk)           
         end do
         Hk(:,:,ik) = Hk(:,:,ik) + H0
      end do

      do ik=1,Lk
         write(656,*)
         do io=1,No
            write(656,'(10(F18.10))') dreal(Hk(io,:,ik))
         end do
         write(656,*)
      end do

      !
      unit=free_unit()
      open(unit,file="hk_symmetric_points")
      Hk_symm=.true.
      do ik=1,Lk
         do io=1,No
            do jo=1,No
               if(Hk(io,jo,ik) /= Hk(jo,io,ik)) Hk_symm(ik)=.false.
            end do
         end do
         if(Hk_symm(ik)) write(unit,'(4(F18.10))') kVect(ik)%x,kVect(ik)%y
      end do
      close(unit)

    end subroutine get_k_hamiltonian_slab






    !+------------+!
    !+- build Hk -+!  
    !+------------+!
    subroutine get_k_hamiltonian_ft(kVect)
      type(vect2D)                            :: kVect(:)    
      integer                                 :: i,jj,j,k,row,col,link(4),ik,No
      integer                                 :: unit,unitk
      integer                                 :: iorb,jorb,io,jo,ilat,jlat
      integer                                 :: ikz,Nkz,Nk_tot,ik_
      real(8)                                 :: expk,kz,dk
      type(vect2D),dimension(4)               :: Rhop
      real(8),dimension(:,:,:),allocatable    :: thop
      complex(8),dimension(:,:),allocatable   :: psik
      complex(8),dimension(:,:,:),allocatable :: h_kspace
      real(8),dimension(:),allocatable        :: kz_grid,ekz
      !
      No=Nlat*Norb
      if(.not.allocated(Hk)) allocate(Hk(No,No,Lk))


      Nkz=Nlat
      allocate(kz_grid(Nkz),ekz(Nkz))
      dk = pi/dble(Nkz+1)
      kz = 0.d0
      do ikz = 1,Nkz
         kz = kz + dk
         kz_grid(ikz) = kz        
         ekz(ikz) = -2.d0*ts*cos(kz)
      end do

      allocate(psik(Nkz,Nlat))
      do ilat=1,Nlat
         do ikz=1,Nkz
            psik(ikz,ilat) = sqrt(2.d0/dble(Nkz+1))*sin(kz_grid(ikz)*dble(ilat))
         end do
      end do

      Nk_tot=Lk*Nkz
      allocate(h_kspace(Nk_tot,Norb,Norb))
      h_kspace=zero

      ik_=0
      do ik=1,Lk
         do ikz=1,Nkz
            ik_=ik_+1
            do iorb=1,Norb
               do jorb=1,Norb
                  if(iorb.eq.jorb) then
                     h_kspace(ik_,iorb,jorb) = epsik(ik) + ekz(ikz)
                  else
                     h_kspace(ik_,iorb,jorb) =  Vhyb*(dcos(kVect(ik)%x)-dcos(kVect(ik)%y))*dcos(kz_grid(ikz))
                  end if
               end do
            end do
         end do
      end do


      do ilat=1,Nlat
         do jlat=1,Nlat
            do iorb=1,Norb
               do jorb=1,Norb
                  io = iorb + (ilat-1)*Norb!lat_orb2lo(ilat,iorb)
                  jo = jorb + (jlat-1)*Norb!lat_orb2lo(jlat,jorb)
                  ik_=0
                  do ik=1,Lk
                     Hk(io,jo,ik)=zero
                     do ikz=1,Nkz
                        ik_=ik_+1
                        Hk(io,jo,ik) = Hk(io,jo,ik) + psik(ikz,ilat)*psik(ikz,jlat)*h_kspace(ik_,iorb,jorb)
                     end do
                  end do
               end do
            end do
         end do
      end do

      !
      unit=free_unit()
      open(unit,file="hk_symmetric_points")
      Hk_symm=.true.
      do ik=1,Lk
         do io=1,No
            do jo=1,No
               if(Hk(io,jo,ik) /= Hk(jo,io,ik)) Hk_symm(ik)=.false.
            end do
         end do
         if(Hk_symm(ik)) write(unit,'(4(F18.10))') kVect(ik)%x,kVect(ik)%y
      end do
      close(unit)

    end subroutine get_k_hamiltonian_ft




  end program ed_slab


