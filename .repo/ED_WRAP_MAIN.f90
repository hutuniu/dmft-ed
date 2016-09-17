module ED_WRAP_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_DIAG
  USE ED_BATH
  USE ED_MAIN
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: reg,store_data,txtfy
  USE SF_TIMER
  USE SF_LINALG
  implicit none
  private

  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara_lattice
     module procedure ed_get_sigma_matsubara_1
     module procedure ed_get_sigma_matsubara_2
     module procedure ed_get_sigma_matsubara_3
     module procedure ed_get_sigma_matsubara_11
     module procedure ed_get_sigma_matsubara_21
     module procedure ed_get_sigma_matsubara_31
  end interface ed_get_sigma_matsubara_lattice

  interface ed_get_self_matsubara_lattice
     module procedure ed_get_self_matsubara_1
     module procedure ed_get_self_matsubara_2
     module procedure ed_get_self_matsubara_3
     module procedure ed_get_self_matsubara_11
     module procedure ed_get_self_matsubara_21
     module procedure ed_get_self_matsubara_31
  end interface ed_get_self_matsubara_lattice

  interface ed_get_sigma_real_lattice
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
     module procedure ed_get_sigma_real_3
     module procedure ed_get_sigma_real_11
     module procedure ed_get_sigma_real_21
     module procedure ed_get_sigma_real_31
  end interface ed_get_sigma_real_lattice

  interface ed_get_self_real_lattice
     module procedure ed_get_self_real_1
     module procedure ed_get_self_real_2
     module procedure ed_get_self_real_3
     module procedure ed_get_self_real_11
     module procedure ed_get_self_real_21
     module procedure ed_get_self_real_31
  end interface ed_get_self_real_lattice


  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara_lattice
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
     module procedure ed_get_gimp_matsubara_3
     module procedure ed_get_gimp_matsubara_11
     module procedure ed_get_gimp_matsubara_21
     module procedure ed_get_gimp_matsubara_31
  end interface ed_get_gimp_matsubara_lattice

  interface ed_get_fimp_matsubara_lattice
     module procedure ed_get_fimp_matsubara_1
     module procedure ed_get_fimp_matsubara_2
     module procedure ed_get_fimp_matsubara_3
     module procedure ed_get_fimp_matsubara_11
     module procedure ed_get_fimp_matsubara_21
     module procedure ed_get_fimp_matsubara_31
  end interface ed_get_fimp_matsubara_lattice

  interface ed_get_gimp_real_lattice
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
     module procedure ed_get_gimp_real_3
     module procedure ed_get_gimp_real_11
     module procedure ed_get_gimp_real_21
     module procedure ed_get_gimp_real_31
  end interface ed_get_gimp_real_lattice

  interface ed_get_fimp_real_lattice
     module procedure ed_get_fimp_real_1
     module procedure ed_get_fimp_real_2
     module procedure ed_get_fimp_real_3
     module procedure ed_get_fimp_real_11
     module procedure ed_get_fimp_real_21
     module procedure ed_get_fimp_real_31
  end interface ed_get_fimp_real_lattice

  !Retrieve static common observables  
  interface ed_get_dens_lattice
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
  end interface ed_get_dens_lattice

  interface ed_get_mag_lattice
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
  end interface ed_get_mag_lattice

  interface ed_get_docc_lattice
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
  end interface ed_get_docc_lattice

  interface ed_get_phisc_lattice
     module procedure ed_get_phisc_1
     module procedure ed_get_phisc_2
  end interface ed_get_phisc_lattice


  !Overload lattice solver routines:
  interface ed_solve_lattice
     module procedure ed_solve_impurity_sites_eloc
     module procedure ed_solve_impurity_sites_hloc
  end interface ed_solve_lattice


  public :: ed_init_solver_lattice
  public :: ed_solve_lattice
  !
  public :: ed_get_sigma_matsubara_lattice
  public :: ed_get_self_matsubara_lattice
  public :: ed_get_sigma_real_lattice
  public :: ed_get_self_real_lattice
  !
  public :: ed_get_gimp_matsubara_lattice
  public :: ed_get_fimp_matsubara_lattice
  public :: ed_get_gimp_real_lattice
  public :: ed_get_fimp_real_lattice
  !
  public :: ed_get_dens_lattice
  public :: ed_get_mag_lattice
  public :: ed_get_docc_lattice
  public :: ed_get_phisc_lattice
  !
  public :: ed_get_eimp_lattice
  public :: ed_get_epot_lattice
  public :: ed_get_eint_lattice
  public :: ed_get_ehartree_lattice
  public :: ed_get_eknot_lattice
  !
  public :: ed_get_doubles_lattice
  public :: ed_get_dust_lattice
  public :: ed_get_dund_lattice
  public :: ed_get_dse_lattice
  public :: ed_get_dph_lattice
  !

  real(8),dimension(:,:),allocatable,save            :: nii,dii,mii,pii,ddii,eii ![Nlat][Norb/4]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii      ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: SAmatsii,SArealii    ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii      ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Fmatsii,Frealii      ![Nlat][Nspin][Nspin][Norb][Norb][L]
  real(8),dimension(:),allocatable :: wr,wm
  character(len=20)                :: suffix

contains


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize lattice baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:) :: bath ![Nlat][:]
    integer                :: ilat,Nineq
    logical                :: check_dim
    character(len=5)       :: tmp_suffix
    Nineq = size(bath,1)
    if(Nineq > Nlat)stop "init_lattice_bath error: size[bath,1] > Nlat"
    do ilat=1,Nineq
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       call ed_init_solver(bath(ilat,:))
    end do
#ifdef _MPI_INEQ
    call MPI_Barrier(MPI_COMM_WORLD,mpiERR)
#endif
  end subroutine ed_init_solver_lattice




  !-------------------------------------------------------------------------------------------
  !PURPOSE: solve the impurity problems for each independent
  ! lattice site using ED. 
  !-------------------------------------------------------------------------------------------
  subroutine ed_solve_impurity_sites_eloc(bath,Eloc,iprint,Uloc_ii,Ust_ii,Jh_ii)
    real(8)                  :: bath(:,:)                       ![Nlat][Nb]
    real(8),optional         :: Eloc(size(bath,1)*Norb*Nspin)   !A compact version of Hloc_ilat
    integer :: iprint
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI 
    complex(8)               :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Fmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Freal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: mii_tmp(size(bath,1),Norb)
    real(8)                  :: pii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1),4)
    real(8)                  :: ddii_tmp(size(bath,1),4)
    !
    integer                  :: ilat,i,Nsites,iorb,jorb,ispin,jspin
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    if(allocated(Fmatsii))deallocate(Fmatsii)
    if(allocated(Frealii))deallocate(Frealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Frealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero ; Smats_tmp  = zero
    Srealii  = zero ; Sreal_tmp  = zero
    SAmatsii = zero ; SAmats_tmp = zero
    SArealii = zero ; SAreal_tmp = zero
    Gmatsii  = zero ; Gmats_tmp  = zero
    Grealii  = zero ; Greal_tmp  = zero
    Fmatsii  = zero ; Fmats_tmp  = zero
    Frealii  = zero ; Freal_tmp  = zero
    nii      = 0d0  ; nii_tmp    = 0d0
    dii      = 0d0  ; dii_tmp    = 0d0
    mii      = 0d0  ; mii_tmp    = 0d0
    pii      = 0d0  ; pii_tmp    = 0d0
    eii      = 0d0  ; eii_tmp    = 0d0
    ddii     = 0d0  ; ddii_tmp   = 0d0
    !
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Set the local part of the Hamiltonian.
       Hloc=zero
       if(present(Eloc))then
          do ispin=1,Nspin
             do iorb=1,Norb
                i=iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                Hloc(ispin,ispin,iorb,iorb) = Eloc(i)
             enddo
          enddo
       endif
       call set_Hloc(Hloc)
       !
       !Solve the impurity problem for the ilat-th site
       call ed_solve(bath(ilat,:))
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Fmats_tmp(ilat,:,:,:,:,:)  = impFmats(:,:,:,:,:)
       Freal_tmp(ilat,:,:,:,:,:)  = impFreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)       = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)       = ed_docc(1:Norb)
       mii_tmp(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       pii_tmp(ilat,1:Norb)       = ed_phisc(1:Norb)
       eii_tmp(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii_tmp(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       !
    enddo
    if(mpiID==0)call stop_timer
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Fmats_tmp,Fmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Freal_tmp,Frealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Smatsii  =  Smats_tmp
    Srealii  =  Sreal_tmp
    SAmatsii = SAmats_tmp
    SArealii = SAreal_tmp
    Gmatsii  = Gmats_tmp
    Grealii  = Greal_tmp
    Fmatsii  = Fmats_tmp
    Frealii  = Freal_tmp
    nii      = nii_tmp
    dii      = dii_tmp
    mii      = mii_tmp
    pii      = pii_tmp
    eii      = eii_tmp
    ddii     = ddii_tmp
#endif
    if(ed_verbose>4)then
       if(mpiID==0)then
          if(allocated(wm))deallocate(wm)
          if(allocated(wr))deallocate(wr)
          allocate(wm(Lmats))
          allocate(wr(Lreal))
          wm = pi/beta*(2*arange(1,Lmats)-1)
          wr = linspace(wini,wfin,Lreal)
          select case(iprint)
          case (0)
             write(LOGfile,*)"Sigma not written on file."
          case(1)                  !print only diagonal elements
             write(LOGfile,*)"write spin-orbital diagonal elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,iorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,iorb,:),wr)
                   if(ed_mode=="superc")then
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,iorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,iorb,:),wr)
                   endif
                enddo
             enddo
          case(2)                  !print spin-diagonal, all orbitals 
             write(LOGfile,*)"write spin diagonal and all orbitals elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,jorb,:),wr)
                      if(ed_mode=="superc")then
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                         call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                         call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,jorb,:),wr)
                      endif
                   enddo
                enddo
             enddo
          case default                  !print all off-diagonals
             write(LOGfile,*)"write all elements:"
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                         call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                         call store_data("LSigma"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         if(ed_mode=="superc")then
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                            call store_data("LSelf"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                            call store_data("LSelf"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          end select
       endif
    endif
  end subroutine ed_solve_impurity_sites_eloc



  !-------------------------------- HLOC  -------------------------------------
  subroutine ed_solve_impurity_sites_hloc(bath,Hloc,iprint,Uloc_ii,Ust_ii,Jh_ii)
    !inputs
    real(8)                  :: bath(:,:) ![Nlat][Nb]
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                  :: iprint
    real(8),optional         :: Uloc_ii(size(bath,1),Norb)
    real(8),optional         :: Ust_ii(size(bath,1))
    real(8),optional         :: Jh_ii(size(bath,1))
    !MPI  auxiliary vars
    complex(8)               :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Fmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Freal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)                  :: nii_tmp(size(bath,1),Norb)
    real(8)                  :: dii_tmp(size(bath,1),Norb)
    real(8)                  :: mii_tmp(size(bath,1),Norb)
    real(8)                  :: pii_tmp(size(bath,1),Norb)
    real(8)                  :: eii_tmp(size(bath,1),4)
    real(8)                  :: ddii_tmp(size(bath,1),4)
    ! 
    integer                  :: ilat,iorb,jorb,ispin,jspin
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    if(allocated(Fmatsii))deallocate(Fmatsii)
    if(allocated(Frealii))deallocate(Frealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Frealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+mpiID,Nsites,mpiSIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero ; Smats_tmp  = zero
    Srealii  = zero ; Sreal_tmp  = zero
    SAmatsii = zero ; SAmats_tmp = zero
    SArealii = zero ; SAreal_tmp = zero
    Gmatsii  = zero ; Gmats_tmp  = zero
    Grealii  = zero ; Greal_tmp  = zero
    Fmatsii  = zero ; Fmats_tmp  = zero
    Frealii  = zero ; Freal_tmp  = zero
    nii      = 0d0  ; nii_tmp    = 0d0
    dii      = 0d0  ; dii_tmp    = 0d0
    mii      = 0d0  ; mii_tmp    = 0d0
    pii      = 0d0  ; pii_tmp    = 0d0
    eii      = 0d0  ; eii_tmp    = 0d0
    ddii     = 0d0  ; ddii_tmp   = 0d0
    !
    if(mpiID==0)call start_timer
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nsites,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Set the local part of the Hamiltonian.
       call set_Hloc(Hloc(ilat,:,:,:,:))
       ! 
       !Solve the impurity problem for the ilat-th site
       call ed_solve(bath(ilat,:))
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Fmats_tmp(ilat,:,:,:,:,:)  = impFmats(:,:,:,:,:)
       Freal_tmp(ilat,:,:,:,:,:)  = impFreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)       = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)       = ed_docc(1:Norb)
       mii_tmp(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       pii_tmp(ilat,1:Norb)       = ed_phisc(1:Norb)
       eii_tmp(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii_tmp(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
    enddo
    if(mpiID==0)call stop_timer
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Fmats_tmp,Fmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Freal_tmp,Frealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Smatsii  =  Smats_tmp
    Srealii  =  Sreal_tmp
    SAmatsii = SAmats_tmp
    SArealii = SAreal_tmp
    Gmatsii  = Gmats_tmp
    Grealii  = Greal_tmp
    Fmatsii  = Fmats_tmp
    Frealii  = Freal_tmp
    nii      = nii_tmp
    dii      = dii_tmp
    mii      = mii_tmp
    pii      = pii_tmp
    eii      = eii_tmp
    ddii     = ddii_tmp
#endif
    if(ed_verbose>4)then
       if(mpiID==0)then
          if(allocated(wm))deallocate(wm)
          if(allocated(wr))deallocate(wr)
          allocate(wm(Lmats))
          allocate(wr(Lreal))
          wm = pi/beta*(2*arange(1,Lmats)-1)
          wr = linspace(wini,wfin,Lreal)
          select case(iprint)
          case (0)
             write(LOGfile,*)"Sigma not written on file."
          case(1)                  !print only diagonal elements
             write(LOGfile,*)"write spin-orbital diagonal elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,iorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,iorb,:),wr)
                   if(ed_mode=="superc")then
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,iorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,iorb,:),wr)
                   endif
                enddo
             enddo
          case(2)                  !print spin-diagonal, all orbitals 
             write(LOGfile,*)"write spin diagonal and all orbitals elements:"
             do ispin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                      call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,ispin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                      call store_data("LSigma"//reg(suffix),Srealii(:,ispin,ispin,iorb,jorb,:),wr)
                      if(ed_mode=="superc")then
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                         call store_data("LSelf"//reg(suffix),SAmatsii(:,ispin,ispin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                         call store_data("LSelf"//reg(suffix),SArealii(:,ispin,ispin,iorb,jorb,:),wr)
                      endif
                   enddo
                enddo
             enddo
          case default                  !print all off-diagonals
             write(LOGfile,*)"write all elements:"
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                         call store_data("LSigma"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                         suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                         call store_data("LSigma"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         if(ed_mode=="superc")then
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                            call store_data("LSelf"//reg(suffix),Smatsii(:,ispin,jspin,iorb,jorb,:),wm)
                            suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                            call store_data("LSelf"//reg(suffix),Srealii(:,ispin,jspin,iorb,jorb,:),wr)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          end select
       endif
    endif
  end subroutine ed_solve_impurity_sites_hloc





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_sigma_matsubara_1(Smats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(1:Nsites,:,:,:,:,:) = Smatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_1
  !
  subroutine ed_get_sigma_matsubara_2(Smats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Smats(ilat,io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_matsubara_2
  !
  subroutine ed_get_sigma_matsubara_3(Smats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Smats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Smats(1:Nsites,:) = Smatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_matsubara_3
  !
  subroutine ed_get_sigma_matsubara_11(Smats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats(:,:,:,:,:) = Smatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_sigma_matsubara_11
  !
  subroutine ed_get_sigma_matsubara_21(Smats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Smats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Smats(io,jo,:) = Smatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_matsubara_21
  !
  subroutine ed_get_sigma_matsubara_31(Smats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Smats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Smats(:) = Smatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_matsubara_31



  !ANOMALous, MATSUBARA SELF-ENEGRGY
  subroutine ed_get_self_matsubara_1(SAmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
    SAmats(1:Nsites,:,:,:,:,:) = SAmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_self_matsubara_1
  !
  subroutine ed_get_self_matsubara_2(SAmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   SAmats(ilat,io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_matsubara_2
  !
  subroutine ed_get_self_matsubara_3(SAmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: SAmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAmats(1:Nsites,:) = SAmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_matsubara_3
  !
  subroutine ed_get_self_matsubara_11(SAmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
    SAmats(:,:,:,:,:) = SAmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_self_matsubara_11
  !
  subroutine ed_get_self_matsubara_21(SAmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAmats(io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_matsubara_21
  !
  subroutine ed_get_self_matsubara_31(SAmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: SAmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAmats(:) = SAmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_matsubara_31



  !NORMAL, REAL SELF-ENEGRGY
  subroutine ed_get_sigma_real_1(Sreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(1:Nsites,:,:,:,:,:) = Srealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_sigma_real_1
  !
  subroutine ed_get_sigma_real_2(Sreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Sreal(ilat,io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_real_2
  !
  subroutine ed_get_sigma_real_3(Sreal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Sreal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Sreal(1:Nsites,:) = Srealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_real_3
  !
  subroutine ed_get_sigma_real_11(Sreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal(:,:,:,:,:) = Srealii(ilat,:,:,:,:,:)
  end subroutine ed_get_sigma_real_11
  !
  subroutine ed_get_sigma_real_21(Sreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Sreal(io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_sigma_real_21
  !
  subroutine ed_get_sigma_real_31(Sreal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Sreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Sreal(:) = Srealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_sigma_real_31



  !ANOMALous, REAL SELF-ENERGY
  subroutine ed_get_self_real_1(SAreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
    SAreal(1:Nsites,:,:,:,:,:) = SArealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_self_real_1
  !
  subroutine ed_get_self_real_2(SAreal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   SAreal(ilat,io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_real_2
  !
  subroutine ed_get_self_real_3(SAreal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: SAreal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAreal(1:Nsites,:) = SArealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_real_3
  !
  subroutine ed_get_self_real_11(SAreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
    SAreal(:,:,:,:,:) = SArealii(ilat,:,:,:,:,:)
  end subroutine ed_get_self_real_11
  !
  subroutine ed_get_self_real_21(SAreal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                SAreal(io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_self_real_21
  !
  subroutine ed_get_self_real_31(SAreal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: SAreal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    SAreal(:) = SArealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_self_real_31














  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTION
  subroutine ed_get_gimp_matsubara_1(Gmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(1:Nsites,:,:,:,:,:) = Gmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_1
  !
  subroutine ed_get_gimp_matsubara_2(Gmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gmats(ilat,io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_matsubara_2
  !
  subroutine ed_get_gimp_matsubara_3(Gmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Gmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Gmats(1:Nsites,:) = Gmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_matsubara_3
  !
  subroutine ed_get_gimp_matsubara_11(Gmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats(:,:,:,:,:) = Gmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_gimp_matsubara_11
  !
  subroutine ed_get_gimp_matsubara_21(Gmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Gmats(io,jo,:) = Gmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_matsubara_21
  !
  subroutine ed_get_gimp_matsubara_31(Gmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Gmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Gmats(:) = Gmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_matsubara_31



  !ANOMALous, MATSUBARA GREEN'S FUNCTION
  subroutine ed_get_fimp_matsubara_1(Fmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
    Fmats(1:Nsites,:,:,:,:,:) = Fmatsii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_fimp_matsubara_1
  !
  subroutine ed_get_fimp_matsubara_2(Fmats,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Fmats(ilat,io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_matsubara_2
  !
  subroutine ed_get_fimp_matsubara_3(Fmats,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lmats),intent(inout) :: Fmats
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Fmats(1:Nsites,:) = Fmatsii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_matsubara_3
  !
  subroutine ed_get_fimp_matsubara_11(Fmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
    Fmats(:,:,:,:,:) = Fmatsii(ilat,:,:,:,:,:)
  end subroutine ed_get_fimp_matsubara_11
  !
  subroutine ed_get_fimp_matsubara_21(Fmats,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Fmats(io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_matsubara_21
  !
  subroutine ed_get_fimp_matsubara_31(Fmats,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lmats),intent(inout) :: Fmats
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Fmats(:) = Fmatsii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_matsubara_31



  !NORMAL, REAL GREEN'S FUNCTION
  subroutine ed_get_gimp_real_1(Greal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(1:Nsites,:,:,:,:,:) = Grealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_gimp_real_1
  !
  subroutine ed_get_gimp_real_2(Greal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Greal(ilat,io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_real_2
  !
  subroutine ed_get_gimp_real_3(Greal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Greal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Greal(1:Nsites,:) = Grealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_real_3
  !
  subroutine ed_get_gimp_real_11(Greal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal(:,:,:,:,:) = Grealii(ilat,:,:,:,:,:)
  end subroutine ed_get_gimp_real_11
  !
  subroutine ed_get_gimp_real_21(Greal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Greal(io,jo,:) = Grealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_gimp_real_21
  !
  subroutine ed_get_gimp_real_31(Greal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Greal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Greal(:) = Grealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_gimp_real_31



  !ANOMALous, REAL GREEN'S FUNCTION
  subroutine ed_get_fimp_real_1(Freal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
    Freal(1:Nsites,:,:,:,:,:) = Frealii(1:Nsites,:,:,:,:,:)
  end subroutine ed_get_fimp_real_1
  !
  subroutine ed_get_fimp_real_2(Freal,Nsites)
    integer                                                                :: Nsites
    complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
    integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nsites
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Freal(ilat,io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_real_2
  !
  subroutine ed_get_fimp_real_3(Freal,Nsites,ispin,jspin,iorb,jorb)
    integer                                          :: Nsites
    complex(8),dimension(Nsites,Lreal),intent(inout) :: Freal
    integer,optional                                 :: iorb,jorb,ispin,jspin
    integer                                          :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Freal(1:Nsites,:) = Frealii(1:Nsites,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_real_3
  !
  subroutine ed_get_fimp_real_11(Freal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
    Freal(:,:,:,:,:) = Frealii(ilat,:,:,:,:,:)
  end subroutine ed_get_fimp_real_11
  !
  subroutine ed_get_fimp_real_21(Freal,ilat)
    integer                                                         :: ilat
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
    integer                                                         :: io,jo,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Freal(io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end subroutine ed_get_fimp_real_21
  !
  subroutine ed_get_fimp_real_31(Freal,ilat,ispin,jspin,iorb,jorb)
    integer                                   :: ilat
    complex(8),dimension(Lreal),intent(inout) :: Freal
    integer,optional                          :: iorb,jorb,ispin,jspin
    integer                                   :: iorb_,jorb_,ispin_,jspin_
    ispin_=1 ; if(present(ispin))ispin_=ispin
    jspin_=1 ; if(present(jspin))jspin_=jspin
    iorb_=1  ; if(present(iorb))iorb_=iorb
    jorb_=1  ; if(present(jorb))jorb_=jorb
    Freal(:) = Frealii(ilat,ispin_,jspin_,iorb_,jorb_,:)
  end subroutine ed_get_fimp_real_31













  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of some local observables
  !+-----------------------------------------------------------------------------+!
  function ed_get_dens_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0    
    if(allocated(nii))then
       if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
       yii=nii
    end if
  end function ed_get_dens_1
  function ed_get_dens_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
    yii=0d0
    if(allocated(nii))then
       if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
       yii=nii(:,iorb)
    endif
  end function ed_get_dens_2


  function ed_get_docc_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(dii))then
       if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
       yii=dii
    endif
  end function ed_get_docc_1
  function ed_get_docc_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
    yii=0d0
    if(allocated(dii))then
       if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
       yii=dii(:,iorb)
    endif
  end function ed_get_docc_2



  function ed_get_mag_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(mii))then
       if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
       yii=mii
    endif
  end function ed_get_mag_1
  function ed_get_mag_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_mag error: orbital index > N_orbital"
    yii=0d0
    if(allocated(mii))then
       if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
       yii=mii(:,iorb)
    endif
  end function ed_get_mag_2



  function ed_get_phisc_1(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,Norb) :: yii
    yii=0d0
    if(allocated(pii))then
       if(Nlat>size(pii,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"   
       yii=pii
    endif
  end function ed_get_phisc_1
  function ed_get_phisc_2(Nlat,iorb) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    integer                 :: iorb
    if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
    yii=0d0
    if(allocated(pii))then
       if(Nlat>size(pii,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"
       yii=pii(:,iorb)
    endif
  end function ed_get_phisc_2




  function ed_get_eimp_lattice(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,4)    :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_eimp error: required N_sites > evaluated N_sites"
       yii=eii
    endif
  end function ed_get_eimp_lattice

  function ed_get_epot_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_epot error: required N_sites > evaluated N_sites"
       yii=eii(:,1)
    endif
  end function ed_get_epot_lattice

  function ed_get_eint_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_eint error: required N_sites > evaluated N_sites"
       yii=eii(:,2)
    endif
  end function ed_get_eint_lattice

  function ed_get_ehartree_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
       yii=eii(:,3)
    endif
  end function ed_get_ehartree_lattice

  function ed_get_eknot_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(eii))then
       if(Nlat>size(eii,1)) stop "ed_get_knot error: required N_sites > evaluated N_sites"
       yii=eii(:,4)
    endif
  end function ed_get_eknot_lattice





  function ed_get_doubles_lattice(Nlat) result(yii)
    integer                      :: Nlat
    real(8),dimension(Nlat,4)    :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
       yii=ddii(:,:)
    endif
  end function ed_get_doubles_lattice

  function ed_get_dust_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dust error: required N_sites > evaluated N_sites"
       yii=ddii(:,1)
    endif
  end function ed_get_dust_lattice

  function ed_get_dund_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dund error: required N_sites > evaluated N_sites"
       yii=ddii(:,2)
    endif
  end function ed_get_dund_lattice

  function ed_get_dse_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dse error: required N_sites > evaluated N_sites"
       yii=ddii(:,3)
    endif
  end function ed_get_dse_lattice

  function ed_get_dph_lattice(Nlat) result(yii)
    integer                 :: Nlat
    real(8),dimension(Nlat) :: yii
    yii=0d0
    if(allocated(ddii))then
       if(Nlat>size(ddii,1)) stop "ed_get_dph error: required N_sites > evaluated N_sites"
       yii=ddii(:,4)
    endif
  end function ed_get_dph_lattice



end module ED_WRAP_MAIN

