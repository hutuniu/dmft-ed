module ED_GLOC
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS,   only: reg,txtfy,splot,store_data
  USE SF_LINALG,    only: eye,inv,inv_sym,inv_tridiag,get_tridiag
  USE SF_ARRAYS,    only: linspace,arange
  USE SF_MISC,      only: assert_shape
  implicit none
  private
  !
  interface ed_get_gloc
     module procedure ed_get_gloc_normal_main
     module procedure ed_get_gloc_superc_main
     module procedure ed_get_gloc_normal_1b
     module procedure ed_get_gloc_normal_mb
     module procedure ed_get_gloc_superc_1b
     module procedure ed_get_gloc_superc_mb
  end interface ed_get_gloc

  interface ed_get_gloc_lattice
     module procedure ed_get_gloc_normal_lattice_main
     module procedure ed_get_gloc_superc_lattice_main
     module procedure ed_get_gloc_normal_lattice_1b
     module procedure ed_get_gloc_normal_lattice_mb
     module procedure ed_get_gloc_superc_lattice_1b
     module procedure ed_get_gloc_superc_lattice_mb
  end interface ed_get_gloc_lattice

  interface ed_get_gij_lattice
     module procedure ed_get_gij_normal_main
     module procedure ed_get_gij_superc_main
  end interface ed_get_gij_lattice


  public :: ed_get_gloc
  public :: ed_get_gloc_lattice
  public :: ed_get_gij_lattice

  real(8),dimension(:),allocatable :: wr,wm
  character(len=64)                :: suffix

contains

  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the NORMAL local Green's function given Hamiltonian and self-energy.
  ! Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! size [Nlat*Nspin*Norb]**2*[Nk]
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                            :: iprint
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    !local integers
    integer                                       :: Nspin,Norb,Nso,Lmats,Lreal,Lk
    integer                                       :: i,ik,iorb,jorb,ispin,jspin,io,jo
    type(effective_bath)                          :: dmft_bath_tmp


    integer                                       :: unit_ik
    !
    !Testing part:
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Nso   = Nspin*Norb
    Lmats = size(Smats,5)
    Lreal = size(Sreal,5)
    Lk    = size(Hk,3)
    call assert_shape(Hk,[Nso,Nso,Lk],"ed_get_gloc_normal_main","Hk")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_normal_main","Smats")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_normal_main","Sreal")
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_normal_main","Gmats")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_normal_main","Greal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(ED_MPI_MASTER)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_MASTER)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gkreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_mats(Nso,Nso,Lmats));zeta_mats=zero
    allocate(zeta_real(Nso,Nso,Lreal));zeta_real=zero
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo

    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_normal(zeta_mats,(Hk(:,:,ik)),hk_symm_(ik),Gkmats)      
       call add_to_gloc_normal(zeta_real,(Hk(:,:,ik)),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER.AND.ed_verbose<4)call print_Gloc(Gmats,Greal,"Gloc",iprint)
  end subroutine ed_get_gloc_normal_main

  subroutine ed_get_gloc_normal_lattice_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,&
       iprint,hk_symm,tridiag,Gamma_mats,Gamma_real,local_Gamma_mats,local_Gamma_real)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint    !
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    complex(8),dimension(:,:,:),optional            :: Gamma_mats![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),optional            :: Gamma_real![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:),optional          :: local_Gamma_mats![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:),optional          :: local_Gamma_real![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !local integers
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats,Lreal,Lk
    integer                                         :: i,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    !Testing part:
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"ed_get_gloc_normal_lattice_main","Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_normal_lattice_main","Smats")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_normal_lattice_main","Sreal")
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_normal_lattice_main","Gmats")
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_normal_lattice_main","Greal")
    if(present(Gamma_mats))&
         call assert_shape(Gamma_mats,[Nlso,Nlso,Lmats],"ed_get_gloc_normal_lattice_main","Gamma_mats")         
    if(present(Gamma_real))&
         call assert_shape(Gamma_real,[Nlso,Nlso,Lreal],"ed_get_gloc_normal_lattice_main","Gamma_real")
    if(present(local_Gamma_mats))&
         call assert_shape(local_Gamma_mats,[Nlat,Nso,Nso,Lmats],"ed_get_gloc_normal_lattice_main","local_Gamma_mats")         
    if(present(local_Gamma_real))&
         call assert_shape(local_Gamma_real,[Nlat,Nso,Nso,Lreal],"ed_get_gloc_normal_lattice_main","local_Gamma_real")         
    !
    if(present(Gamma_mats).AND.present(local_Gamma_mats))stop "ed_get_gloc_normal_lattice_main error: Gamma_mats & local_Gamma_mats present"
    if(present(Gamma_real).AND.present(local_Gamma_real))stop "ed_get_gloc_normal_lattice_main error: Gamma_real & local_Gamma_real present"
    if(tridiag_.AND.(present(Gamma_mats).OR.present(Gamma_real)))then
       if(ED_MPI_MASTER)write(LOGfile,"(A)")"ed_get_gloc_normal_lattice_main warning: called with tridiag=TRUE and Gamma Embded: Disreagarded."
    endif
    if(ED_MPI_MASTER)write(LOGfile,"(A)")"Get local GF (id=0):"
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    allocate(Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    select case(tridiag_)
    case default
       if(ED_MPI_MASTER)write(LOGfile,"(A)")"Direct Inversion:"
       do ik=1,Lk
          if(present(Gamma_mats))then
             call add_to_gloc_normal_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Gembed=Gamma_mats)
          else
             call add_to_gloc_normal_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          endif
          if(present(Gamma_real))then
             call add_to_gloc_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Gembed=Gamma_real)
          else
             call add_to_gloc_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          endif
          Gmats = Gmats + Gkmats*Wtk(ik)
          Greal = Greal + Gkreal*Wtk(ik)
          if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
       end do
    case(.true.)
       if(ED_MPI_MASTER)write(LOGfile,"(A)")"Tridiag Iterative Inversion:"
       do ik=1,Lk
          if(present(local_Gamma_mats))then
             call add_to_gloc_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Gembed=local_Gamma_mats)
          else
             call add_to_gloc_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          endif
          if(present(local_Gamma_real))then
             call add_to_gloc_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Gembed=local_Gamma_real)
          else
             call add_to_gloc_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          endif
          Gmats = Gmats + Gkmats*Wtk(ik)
          Greal = Greal + Gkreal*Wtk(ik)
          if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
       end do
    end select
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER)call print_gloc_lattice(Gmats,Greal,"LG",iprint)
  end subroutine ed_get_gloc_normal_lattice_main



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: evaluate the GF for a single k-point
  !+-----------------------------------------------------------------------------+!
  subroutine add_to_gloc_normal(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:),intent(in)        :: zeta    ![Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)          :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                            :: hk_symm                
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkout   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gktmp   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable         :: Gmatrix ![Nspin*Norb][Nspin*Norb]
    integer                                       :: Nspin,Norb,Nso,Lfreq
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    Nspin = size(Gkout,1)
    Norb  = size(Gkout,3)
    Lfreq = size(zeta,3)
    Nso   = Nspin*Norb
    call assert_shape(zeta,[Nso,Nso,Lfreq],"add_to_gloc_normal","zeta")
    call assert_shape(Hk,[Nso,Nso],"add_to_gloc_normal","Hk")
    call assert_shape(Gkout,[Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_normal","Gkout")
    allocate(Gktmp(Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nso,Nso))
    Gktmp=zero
    do i=1,Lfreq
       Gmatrix  = zeta(:,:,i) - Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
          !call inv1(100,Gmatrix,Nspin*Norb)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine add_to_gloc_normal

  subroutine add_to_gloc_normal_lattice(zeta,Hk,hk_symm,Gkout,Gembed)
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:),optional            :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"add_to_gloc_lattice_normal","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gloc_lattice_normal","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_lattice_normal","Gkout")
    if(present(Gembed))&
         call assert_shape(Gembed,[Nlso,Nlso,Lfreq],"add_to_gloc_lattice_normal","Gembed")   
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1+ED_MPI_ID,Lfreq,ED_MPI_SIZE
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i)) - Hk
       if(present(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
#ifdef _MPI
    call MPI_BARRIER(ED_MPI_COMM,ED_MPI_ERR)
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gloc_normal_lattice

  subroutine add_to_gloc_normal_tridiag(zeta,Hk,hk_symm,Gkout,Gembed)
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:),optional          :: Gembed  ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    !
    complex(8),dimension(:,:,:),allocatable         :: Diag
    complex(8),dimension(:,:,:),allocatable         :: Sub
    complex(8),dimension(:,:,:),allocatable         :: Over
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:),allocatable         :: Gmatrix ![Nlat][Nspin*Norb][Nspin*Norb]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"add_to_gloc_tridiag_normal","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gloc_tridiag_normal","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_tridiag_normal","Gkout")
    if(present(Gembed))&
         call assert_shape(Gembed,[Nlat,Nso,Nso,Lfreq],"add_to_gloc_tridiag_normal","Gembed")
    !
    allocate(Sub(Nlat-1,Nso,Nso))
    allocate(Diag(Nlat,Nso,Nso))
    allocate(Over(Nlat-1,Nso,Nso))
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlat,Nso,Nso))
    Gktmp=zero
    do i=1+ED_MPI_ID,Lfreq,ED_MPI_SIZE
       call get_tridiag(Nlat,Nso,Hk,Sub,Diag,Over)
       Diag = zeta(:,:,:,i) - Diag
       if(present(Gembed))Diag = Diag - Gembed(:,:,:,i)
       call inv_tridiag(Nlat,Nso,-Sub,Diag,-Over,Gmatrix)
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(ilat,io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
#ifdef _MPI
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gloc_normal_tridiag





  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  ADDITIONAL INTERFACES to the main procedures above 
  ! - _1b : one band case (input is independent of Nspin and Norb)
  ! - _mb : multi-bands case (input is independent of Nspin)
  !+-----------------------------------------------------------------------------+!
  subroutine ed_get_gloc_normal_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                     :: Smats(:)
    complex(8),intent(in)                     :: Sreal(:)
    complex(8),intent(inout)                  :: Gmats(size(Smats))
    complex(8),intent(inout)                  :: Greal(size(Sreal))
    logical,optional                          :: hk_symm(size(Hk,1))
    logical                                   :: hk_symm_(size(Hk,1))
    integer                                   :: iprint
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(1,1,1,1,size(Smats)) :: Gmats_
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Greal_
    complex(8),dimension(1,1,1,1,size(Smats)) :: Smats_
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Sreal_
    !
    Smats_(1,1,1,1,:) = Smats(:)
    Sreal_(1,1,1,1,:) = Sreal(:)
    Hk_(1,1,:)        = Hk(:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk_,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:) = Gmats_(1,1,1,1,:)
    Greal(:) = Greal_(1,1,1,1,:)
  end subroutine ed_get_gloc_normal_1b

  subroutine ed_get_gloc_normal_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)      :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                       :: Smats(:,:,:)
    complex(8),intent(in)                       :: Sreal(:,:,:)
    complex(8),intent(inout)                    :: Gmats(:,:,:)
    complex(8),intent(inout)                    :: Greal(:,:,:)
    logical,optional                            :: hk_symm(size(Hk,3))
    logical                                     :: hk_symm_(size(Hk,3))
    integer                                     :: iprint
    !
    complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats_
    complex(8),allocatable,dimension(:,:,:,:,:) :: Greal_
    complex(8),allocatable,dimension(:,:,:,:,:) :: Smats_
    complex(8),allocatable,dimension(:,:,:,:,:) :: Sreal_
    !
    integer :: Nspin,Norb,Nso,Lmats,Lreal
    !
    Nspin = 1
    Norb  = size(Smats,1)
    Lmats = size(Smats,3)
    Lreal = size(Sreal,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"ed_get_gloc_normal_mb","Hk")
    call assert_shape(Sreal,[Norb,Norb,Lreal],"ed_get_gloc_normal_mb","Sreal")
    call assert_shape(Gmats,[Norb,Norb,Lmats],"ed_get_gloc_normal_mb","Gmats")
    call assert_shape(Gmats,[Norb,Norb,Lreal],"ed_get_gloc_normal_mb","Greal")
    call assert_shape(Gmats,[Norb,Norb,Lmats],"ed_get_gloc_normal_mb","Gmats")
    !
    allocate(Smats_(1,1,Norb,Norb,Lmats))
    allocate(Gmats_(1,1,Norb,Norb,Lmats))
    allocate(Sreal_(1,1,Norb,Norb,Lreal))
    allocate(Greal_(1,1,Norb,Norb,Lreal))
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    Sreal_(1,1,:,:,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:,:) = Gmats_(1,1,:,:,:)
    Greal(:,:,:) = Greal_(1,1,:,:,:)
  end subroutine ed_get_gloc_normal_mb

  subroutine ed_get_gloc_normal_lattice_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)               :: Hk              ![Nk]
    real(8),intent(in)                                   :: Wtk(size(Hk,3)) ![Nk]
    complex(8),dimension(:,:),intent(in)                 :: Smats           ![Nlat][Lmats]
    complex(8),dimension(:,:),intent(in)                 :: Sreal           ![Nlat][Lreal]
    complex(8),dimension(:,:),intent(inout)              :: Gmats
    complex(8),dimension(:,:),intent(inout)              :: Greal
    logical,optional                                     :: hk_symm(size(Hk,1))
    logical                                              :: hk_symm_(size(Hk,1))
    integer                                              :: iprint
    !
    complex(8),dimension(size(Hk,1),1,1,1,1,size(Smats)) :: Gmats_
    complex(8),dimension(size(Hk,1),1,1,1,1,size(Sreal)) :: Greal_
    complex(8),dimension(size(Hk,1),1,1,1,1,size(Smats)) :: Smats_
    complex(8),dimension(size(Hk,1),1,1,1,1,size(Sreal)) :: Sreal_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    call ed_get_gloc_normal_lattice_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine ed_get_gloc_normal_lattice_1b

  subroutine ed_get_gloc_normal_lattice_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk              ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                         :: Smats(:,:,:,:)  ![Nlat][Norb][Norb][Lmats]
    complex(8),intent(in)                         :: Sreal(:,:,:,:)  ![Nlat][Norb][Norb][Lreal]
    complex(8),intent(inout)                      :: Gmats(:,:,:,:)
    complex(8),intent(inout)                      :: Greal(:,:,:,:)
    logical,optional                              :: hk_symm(size(Hk,3))
    logical                                       :: hk_symm_(size(Hk,3))
    integer                                       :: iprint
    !
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal_ ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_ ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    !
    integer                                       :: Nspin,Norb,Nso,Nlso,Lmats,Lreal
    !
    Nlat  = size(Smats,1)
    Nspin = 1    
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],"ed_get_gloc_normal_mb","Hk")
    call assert_shape(Sreal,[Nlat,Norb,Norb,Lreal],"ed_get_gloc_normal_mb","Sreal")
    call assert_shape(Gmats,[Nlat,Norb,Norb,Lmats],"ed_get_gloc_normal_mb","Gmats")
    call assert_shape(Greal,[Nlat,Norb,Norb,Lreal],"ed_get_gloc_normal_mb","Greal")
    call assert_shape(Gmats,[Nlat,Norb,Norb,Lmats],"ed_get_gloc_normal_mb","Gmats")
    !
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lmats))
    allocate(Gmats_(Nlat,1,1,Norb,Norb,Lmats))
    allocate(Sreal_(Nlat,1,1,Norb,Norb,Lreal))
    allocate(Greal_(Nlat,1,1,Norb,Norb,Lreal))
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_lattice_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine ed_get_gloc_normal_lattice_mb







  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the SUPERC local Green's function given Hamiltonian and self-energy.
  ! Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! size [Nlat*Nspin*Norb]**2*[Nk]
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_superc_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint
    logical,dimension(size(Hk,3)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    !local integers    
    integer                                         :: Nspin,Norb,Nso,Lmats,Lreal,Lk
    integer                                         :: i,ik,iorb,jorb,ispin,jspin,io,jo
    !Testing part:
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,Lk],"ed_get_gloc_superc_main","Hk")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_superc_main","Smats")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_superc_main","Sreal")
    call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_superc_main","Gmats")
    call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_superc_main","Greal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(ED_MPI_MASTER)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_MASTER)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    allocate(Gkmats(2,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                -        nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  -        nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  -        nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER.AND.ed_verbose<4)then
       call print_Gloc(Gmats(1,:,:,:,:,:),Greal(1,:,:,:,:,:),"Gloc",iprint)
       call print_Gloc(Gmats(2,:,:,:,:,:),Greal(2,:,:,:,:,:),"Floc",iprint)
    endif
  end subroutine ed_get_gloc_superc_main

  subroutine ed_get_gloc_superc_lattice_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !local integers
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lmats,Lreal,Lk
    integer                                           :: i,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    !Testing part:
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"ed_get_gloc_superc_lattice_main","Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_superc_lattice_main","Smats")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_superc_lattice_main","Sreal")
    call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gloc_superc_lattice_main","Gmats")
    call assert_shape(Greal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gloc_superc_lattice_main","Greal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(ED_MPI_MASTER)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_MASTER)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    allocate(Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    do ilat=1,Nlat
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do i=1,Lmats
          zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
       enddo
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]   
       do i=1,Lreal
          zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                -        nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(1,2,ilat,:,:,i) =                                                  -        nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,1,ilat,:,:,i) =                                                  -        nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
       enddo
    enddo
    !
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_superc_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER)then
       call print_gloc_lattice(Gmats(1,:,:,:,:,:,:),Greal(1,:,:,:,:,:,:),"LG",iprint)
       call print_gloc_lattice(Gmats(2,:,:,:,:,:,:),Greal(2,:,:,:,:,:,:),"LF",iprint)
    endif
  end subroutine ed_get_gloc_superc_lattice_main




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: evaluate the GF for a single k-point
  !+-----------------------------------------------------------------------------+!
  subroutine add_to_gloc_superc(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                              :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![2*Nspin*Norb][2*Nspin*Norb]
    integer                                         :: Nspin,Norb,Nso,Lfreq
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,5)
    Nso   = Nspin*Norb
    call assert_shape(zeta,[2,2,Nso,Nso,Lfreq],"add_to_gloc_superc","zeta")
    call assert_shape(Hk,[Nso,Nso],"add_to_gloc_superc","Hk")
    call assert_shape(Gkout,[2,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_superc","Gkout")
    !
    allocate(Gktmp(2,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nso,2*Nso))
    Gkout = zero
    Gktmp = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Nso,1:Nso)             = zeta(1,1,:,:,i) - Hk
       Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta(1,2,:,:,i)
       Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta(2,1,:,:,i)
       Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta(2,2,:,:,i) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(1,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   Gktmp(2,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nso+jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine add_to_gloc_superc

  subroutine add_to_gloc_superc_lattice(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![2*Nlat*Nspin*Norb][2*Nlat*Nspin*Norb]
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                           :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,3)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[2,2,Nlat,Nso,Nso,Lfreq],"add_to_gloc_lattice_superc","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gloc_lattice_superc","Hk")
    call assert_shape(Gkout,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_lattice_superc","Gkout")
    !
    allocate(Gktmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gkout = zero
    Gktmp = zero
    do i=1+ED_MPI_ID,Lfreq,ED_MPI_SIZE
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i)) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i)) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(1,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      Gktmp(2,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
#ifdef _MPI
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gloc_superc_lattice



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  ADDITIONAL INTERFACES to the main procedures above 
  ! - _1b : one band case (input is independent of Nspin and Norb)
  ! - _mb : multi-bands case (input is independent of Nspin)
  !+-----------------------------------------------------------------------------+!
  subroutine ed_get_gloc_superc_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:),intent(in) :: Hk              ![Nk]
    real(8),intent(in)                 :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)              :: Smats(:,:)
    complex(8),intent(in)              :: Sreal(:,:)
    complex(8),intent(inout)           :: Gmats(2,size(Smats,2))
    complex(8),intent(inout)           :: Greal(2,size(Sreal,2))
    logical,optional                   :: hk_symm(size(Hk,1))
    logical                            :: hk_symm_(size(Hk,1))
    integer                            :: iprint
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Gmats_
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Greal_
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Smats_
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Sreal_
    !
    call assert_shape(Smats,[2,size(Smats,2)],"ed_get_gloc_superc_1b","Smats")
    call assert_shape(Sreal,[2,size(Sreal,2)],"ed_get_gloc_superc_1b","Sreal")
    !
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    Hk_(1,1,:)          = Hk
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_superc_main(Hk_,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine ed_get_gloc_superc_1b

  subroutine ed_get_gloc_superc_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8),intent(in)                            :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                         :: Smats(:,:,:,:)
    complex(8),intent(in)                         :: Sreal(:,:,:,:)
    complex(8),intent(inout)                      :: Gmats(:,:,:,:)
    complex(8),intent(inout)                      :: Greal(:,:,:,:)
    logical,optional                              :: hk_symm(size(Hk,3))
    logical                                       :: hk_symm_(size(Hk,3))
    integer                                       :: iprint
    !
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats_
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal_
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats_
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_
    !
    integer                                       :: Nspin,Norb,Nso,Lmats,Lreal
    !
    Nspin = 1
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"ed_get_gloc_superc_mb","Hk")
    call assert_shape(Sreal,[2,Norb,Norb,Lreal],"ed_get_gloc_superc_mb","Sreal")
    call assert_shape(Gmats,[2,Norb,Norb,Lmats],"ed_get_gloc_superc_mb","Gmats")
    call assert_shape(Gmats,[2,Norb,Norb,Lreal],"ed_get_gloc_superc_mb","Greal")
    call assert_shape(Gmats,[2,Norb,Norb,Lmats],"ed_get_gloc_superc_mb","Gmats")
    !
    allocate(Smats_(2,1,1,Norb,Norb,Lmats))
    allocate(Gmats_(2,1,1,Norb,Norb,Lmats))
    allocate(Sreal_(2,1,1,Norb,Norb,Lreal))
    allocate(Greal_(2,1,1,Norb,Norb,Lreal))
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_superc_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine ed_get_gloc_superc_mb

  subroutine ed_get_gloc_superc_lattice_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                   :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                    :: Smats(:,:,:)
    complex(8),intent(in)                                    :: Sreal(:,:,:)
    complex(8),intent(inout)                                 :: Gmats(2,size(Hk,1),size(Smats,3))
    complex(8),intent(inout)                                 :: Greal(2,size(Hk,1),size(Sreal,3))
    logical,optional                                         :: hk_symm(size(Hk,1))
    logical                                                  :: hk_symm_(size(Hk,1))
    integer                                                  :: iprint
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Gmats_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Greal_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Smats_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Sreal_
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],"ed_get_gloc_superc_lattice_1b","Hk")
    call assert_shape(Smats,[2,size(Hk,1),size(Smats,3)],"ed_get_gloc_superc_lattice_1b","Smats")
    call assert_shape(Sreal,[2,size(Hk,1),size(Sreal,3)],"ed_get_gloc_superc_lattice_1b","Sreal")
    !
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    Sreal_(:,:,1,1,1,1,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_superc_lattice_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:,:) = Gmats_(:,:,1,1,1,1,:)
    Greal(:,:,:) = Greal_(:,:,1,1,1,1,:)
  end subroutine ed_get_gloc_superc_lattice_1b

  subroutine ed_get_gloc_superc_lattice_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk               ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                              :: Wtk(size(Hk,3))  ![Nk]
    complex(8),intent(in)                           :: Smats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    complex(8),intent(in)                           :: Sreal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                        :: Gmats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                        :: Greal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    logical,optional                                :: hk_symm(size(Hk,3))
    logical                                         :: hk_symm_(size(Hk,3))
    integer                                         :: iprint
    !
    complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats_           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Greal_           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Smats_           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Sreal_           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    !
    integer                                         :: Nspin,Norb,Nso,Nlso,Lmats,Lreal
    !
    Nlat  = size(Smats,2)
    Nspin = 1
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],"ed_get_gloc_superc_mb","Hk")
    call assert_shape(Sreal,[2,Nlat,Norb,Norb,Lreal],"ed_get_gloc_superc_mb","Sreal")
    call assert_shape(Gmats,[2,Nlat,Norb,Norb,Lmats],"ed_get_gloc_superc_mb","Gmats")
    call assert_shape(Gmats,[2,Nlat,Norb,Norb,Lreal],"ed_get_gloc_superc_mb","Greal")
    call assert_shape(Gmats,[2,Nlat,Norb,Norb,Lmats],"ed_get_gloc_superc_mb","Gmats")
    !
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lmats))
    allocate(Gmats_(2,Nlat,1,1,Norb,Norb,Lmats))
    allocate(Sreal_(2,Nlat,1,1,Norb,Norb,Lreal))
    allocate(Greal_(2,Nlat,1,1,Norb,Norb,Lreal))
    !
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    Sreal_(:,:,1,1,:,:,:) = Sreal(:,:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_superc_lattice_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,hk_symm_)
    Gmats(:,:,:,:,:) = Gmats_(:,:,1,1,:,:,:)
    Greal(:,:,:,:,:) = Greal_(:,:,1,1,:,:,:)
  end subroutine ed_get_gloc_superc_lattice_mb



  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the NORMALcomplete Green's function given Hamiltonian and self-energy.
  ! Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! size [Nlat*Nspin*Norb]**2*[Nk]
  !----------------------------------------------------------------------------------------!
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  ADDITIONAL INTERFACES to the main procedures intersite Gloc included
  !+-----------------------------------------------------------------------------+!
  subroutine ed_get_gij_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm,Gamma_mats,Gamma_real)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    complex(8),dimension(:,:,:),optional              :: Gamma_mats![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),optional              :: Gamma_real![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !local integers
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lmats,Lreal,Lk
    integer                                           :: i,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    !Testing part:
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"ed_get_gij_normal_main","Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gij_normal_main","Smats")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gij_normal_main","Sreal")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_gij_normal_main","Gmats")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"ed_get_gij_normal_main","Greal")
    if(present(Gamma_mats))&
         call assert_shape(Gamma_mats,[Nlso,Nlso,Lmats],"ed_get_gij_normal_main","Gamma_mats")         
    if(present(Gamma_real))&
         call assert_shape(Gamma_real,[Nlso,Nlso,Lreal],"ed_get_gij_normal_main","Gamma_real")         
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(ED_MPI_MASTER)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_MASTER)write(*,*)"Get local GF FULL (id=0):"
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       !if(ED_MPI_MASTER)write(*,*)ik
       if(present(Gamma_mats))then
          call add_to_gij_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Gembed=Gamma_mats)
       else
          call add_to_gij_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       endif
       if(present(Gamma_real))then
          call add_to_gij_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Gembed=Gamma_real)
       else
          call add_to_gij_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       endif
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER)call print_gij_lattice(Gmats,Greal,"Gij",iprint)
  end subroutine ed_get_gij_normal_main

  subroutine add_to_gij_normal(zeta,Hk,hk_symm,Gkout,Gembed)
    complex(8),dimension(:,:,:,:),intent(in)          :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:),optional              :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    !local integers
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                           :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"add_to_gij_normal","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gij_normal","Hk")
    call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gij_normal","Gkout")
    if(present(Gembed))&
         call assert_shape(Gembed,[Nlso,Nlso,Lfreq],"add_to_gij_normal","Gembed")   
    allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1+ED_MPI_ID,Lfreq,ED_MPI_SIZE
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i)) - Hk
       if(present(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
#ifdef _MPI
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gij_normal


  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the SUPERC complete Green's function given Hamiltonian and self-energy.
  ! Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! size [Nlat*Nspin*Norb]**2*[Nk]
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gij_superc_main(Hk,Wtk,Gmats,Fmats,Greal,Freal,Smats,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:)              :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                  :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)                 :: Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)                 :: Fmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)                 :: Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)                 :: Freal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)                 :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)                 :: Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                                  :: iprint
    !
    complex(8)                               :: zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                               :: zeta_real(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                               :: Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                               :: Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                               :: Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                               :: Fkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    logical,optional                         :: hk_symm(size(Hk,3))
    logical                                  :: hk_symm_(size(Hk,3))
    integer                                  :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    Lk=size(Hk,3)
    Nlso=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    if(ED_MPI_MASTER)write(LOGfile,*)"Get local GF (id=0):"
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    zeta_mats=zero
    zeta_real=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
             !G22(iw) = -[G11[iw]]*
             !G21(iw) =   G12[w]
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu !+ Eloc_(js)
             !
             !SYMMETRIES in real-frequencies   [assuming a real order parameter]
             !G22(w)  = -[G11[-w]]*
             !G21(w)  =   G12[w]             
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu !- Eloc_(js)
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )!- Eloc_(js) )
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
                   !
                   zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(ED_MPI_MASTER)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gij_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
       call add_to_gij_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Fmats = Fmats + Fkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       Freal = Freal + Fkreal*Wtk(ik)
       if(ED_MPI_MASTER)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_MASTER)call stop_timer
    if(ED_MPI_MASTER)then
       call print_gij_lattice(Gmats,Greal,"Gij",iprint)
       call print_gij_lattice(Fmats,Freal,"Fij",iprint)
    endif
  end subroutine ed_get_gij_superc_main


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: evaluate the GF for a single k-point
  !+-----------------------------------------------------------------------------+!
  subroutine add_to_gij_superc(zeta,Hk,hk_symm,Gkout,Fkout)
    complex(8)                               :: zeta(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)                               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                                  :: hk_symm                
    !output:
    complex(8),intent(inout)                 :: Gkout(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,6))
    complex(8),intent(inout)                 :: Fkout(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,6))
    complex(8)                               :: Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,6))
    complex(8)                               :: Fktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,6))
    !
    complex(8)                               :: Gmatrix(2*Nlat*Nspin*Norb , 2*Nlat*Nspin*Norb)
    integer                                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nlso
    if(size(zeta,1)/=2.OR.size(zeta,2)/=2)stop "add_to_gloc_superc error: zeta wrong size 1 or 2 = 2 (Nambu)"
    if(size(zeta,3)/=Nlat)stop "add_to_gloc_superc error: zeta wrong size 2 = Nlat"
    if(size(zeta,4)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta wrong size 4 = Nspin*Norb"
    if(size(zeta,5)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta wrong size 5 = Nspin*Norb"
    Lfreq = size(zeta,6)
    Nlso  = Nlat*Nspin*Norb
    Gkout = zero
    Gktmp = zero
    do i=1+ED_MPI_ID,Lfreq,ED_MPI_SIZE
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i)) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i)) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                         Fktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    enddo
#ifdef _MPI
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
    call MPI_ALLREDUCE(Fktmp,Fkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,ED_MPI_COMM,ED_MPI_ERR)
#else
    Gkout = Gktmp
    Fkout = Fktmp
#endif
  end subroutine add_to_gij_superc

  !+-----------------------------------------------------------------------------+!
  !PURPOSE: print a local GF according to iprint variable
  !+-----------------------------------------------------------------------------+!
  subroutine print_gloc(Gmats,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats,Greal
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    integer                                    :: Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nspin = size(Gmats,1)
    Norb  = size(Greal,3)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,size(Gmats,5)],"print_gloc","Gmats")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,size(Greal,5)],"print_gloc","Greal")
    !
    select case(iprint)
    case (0)
       if(ED_MPI_MASTER)write(LOGfile,*)"Gloc not written to file."
    case(1)                  !print only diagonal elements
       if(ED_MPI_MASTER)write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)//"_iw.ed"
             call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,iorb,:))
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
             call splot(reg(suffix),wr,-dimag(Greal(ispin,ispin,iorb,iorb,:))/pi,dreal(Greal(ispin,ispin,iorb,iorb,:)))
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       if(ED_MPI_MASTER)write(LOGfile,*)"write spin diagonal and all orbitals elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw"//reg(ed_file_suffix)//".ed"
                call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,jorb,:))
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
                call splot(reg(suffix),wr,-dimag(Greal(ispin,ispin,iorb,jorb,:))/pi,dreal(Greal(ispin,ispin,iorb,jorb,:)))
             enddo
          enddo
       enddo
    case default                  !print all off-diagonals
       if(ED_MPI_MASTER)write(LOGfile,*)"write all elements:"
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw"//reg(ed_file_suffix)//".ed"
                   call splot(reg(suffix),wm,Gmats(ispin,jspin,iorb,jorb,:))
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw"//reg(ed_file_suffix)//".ed"
                   call splot(reg(suffix),wr,-dimag(Greal(ispin,jspin,iorb,jorb,:))/pi,dreal(Greal(ispin,jspin,iorb,jorb,:)))
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine print_gloc

  subroutine print_gloc_lattice(Gmats,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats,Greal
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    integer                                      :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,2)
    Norb  = size(Greal,4)
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,6)],"print_gloc_lattice",reg(fname)//"_mats")
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,size(Greal,6)],"print_gloc_lattice",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(LOGfile,*)"Gloc not written to file."
    case(1)                  !print only diagonal elements
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call store_data(reg(suffix),Gmats(:,ispin,ispin,iorb,iorb,:),wm)
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call store_data(reg(suffix),Greal(:,ispin,ispin,iorb,iorb,:),wr)
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(LOGfile,*)"write spin diagonal and all orbitals elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data(reg(suffix),Gmats(:,ispin,ispin,iorb,jorb,:),wm)
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data(reg(suffix),Greal(:,ispin,ispin,iorb,jorb,:),wr)
             enddo
          enddo
       enddo
    case default                  !print all off-diagonals
       write(LOGfile,*)"write all elements:"
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                   call store_data(reg(suffix),Gmats(:,ispin,jspin,iorb,jorb,:),wm)
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                   call store_data(reg(suffix),Greal(:,ispin,jspin,iorb,jorb,:),wr)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine print_gloc_lattice

  subroutine print_gij_lattice(Gmats,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats,Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,3)
    Norb  = size(Greal,5)
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,7)],"print_gij_lattice",reg(fname)//"_mats")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Greal,7)],"print_gij_lattice",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(LOGfile,*)"Gloc not written on file."
    case(1)                  !print only diagonal elements
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call store_data(reg(suffix),Gmats(:,:,ispin,ispin,iorb,iorb,:),wm)
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call store_data(reg(suffix),Greal(:,:,ispin,ispin,iorb,iorb,:),wr)
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(LOGfile,*)"write spin diagonal and all orbitals elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data(reg(suffix),Gmats(:,:,ispin,ispin,iorb,jorb,:),wm)
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data(reg(suffix),Greal(:,:,ispin,ispin,iorb,jorb,:),wr)
             enddo
          enddo
       enddo
    case default              !print all off-diagonals
       write(LOGfile,*)"write all elements:"
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                   call store_data(reg(suffix),Gmats(:,:,ispin,jspin,iorb,jorb,:),wm)
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                   call store_data(reg(suffix),Greal(:,:,ispin,jspin,iorb,jorb,:),wr)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine print_gij_lattice















  !***************************************************************
  SUBROUTINE INV1(je,A,NP)
    !***************************************************************

    INTEGER,PARAMETER::NNP=100,nrhs=1
    complex (kind=8) :: A(NP,NP)
    COMPLEX*16 det(2)
    COMPLEX*16 AA(NP,NP),WORK(NP)
    INTEGER*4 IPVT(NP),LDA,NN,job,NP,JE,M,N
    REAL*8 RCOND

    NN=NP   ! dimension of the matrix AA
    LDA=NP ! dimension of the matrix to be inverted
    job=1   ! option flag
    DO M=1,NP
       DO N=1,NP
          AA(N,M)=DCMPLX(A(N,M))
       ENDDO
    ENDDO

    CALL ZGECO(AA,LDA,NN,IPVT,RCOND,WORK)

    IF (1.0D0+RCOND .EQ. 1.0D0) THEN
       WRITE(6,*) ' STOP FOR JE=',je
       WRITE(6,*) '*******************'
       WRITE(6,*) '*******************'
       WRITE(6,*) ' SINGULAR MATRIX !'
       WRITE(6,*) '*******************'
       WRITE(6,*) '*******************'
       DO N=1,NP
          WRITE(6,*) '-------------ROW  ',N
          WRITE(6,100) (REAL(AA(N,M)),AIMAG(AA(N,M)),M=1,NP)
       ENDDO
       STOP 
    ENDIF


    call zgedi(aa,lda,nn,ipvt,det,work,job)


    DO M=1,NP
       DO N=1,NP
          A(N,M)=AA(N,M)
       ENDDO
    ENDDO

100 FORMAT(30F7.3)	
    RETURN
  END subroutine INV1


  !****************************************************************
  subroutine zgeco(a,lda,n,ipvt,rcond,z)
    !****************************************************************
    integer lda,n,ipvt(1)
    complex*16 a(lda,1),z(1)
    double precision rcond
    !
    !     zgeco factors a complex*16 matrix by gaussian elimination
    !     and estimates the condition of the matrix.
    !
    !     if  rcond  is not needed, zgefa is slightly faster.
    !     to solve  a*x = b , follow zgeco by zgesl.
    !     to compute  inverse(a)*c , follow zgeco by zgesl.
    !     to compute  determinant(a) , follow zgeco by zgedi.
    !     to compute  inverse(a) , follow zgeco by zgedi.
    !
    !     on entry
    !
    !        a       complex*16(lda, n)
    !                the matrix to be factored.
    !
    !        lda     integer
    !                the leading dimension of the array  a .
    !
    !        n       integer
    !                the order of the matrix  a .
    !
    !     on return
    !
    !        a       an upper triangular matrix and the multipliers
    !                which were used to obtain it.
    !                the factorization can be written  a = l*u  where
    !                l  is a product of permutation and unit lower
    !                triangular matrices and  u  is upper triangular.
    !
    !        ipvt    integer(n)
    !                an integer vector of pivot indices.
    !
    !        rcond   double precision
    !                an estimate of the reciprocal condition of  a .
    !                for the system  a*x = b , relative perturbations
    !                in  a  and  b  of size  epsilon  may cause
    !                relative perturbations in  x  of size  epsilon/rcond .
    !                if  rcond  is so small that the logical expression
    !                           1.0 + rcond .eq. 1.0
    !                is true, then  a  may be singular to working
    !                precision.  in particular,  rcond  is zero  if
    !                exact singularity is detected or the estimate
    !                underflows.
    !
    !        z       complex*16(n)
    !                a work vector whose contents are usually unimportant.
    !                if  a  is close to a singular matrix, then  z  is
    !                an approximate null vector in the sense that
    !                norm(a*z) = rcond*norm(a)*norm(z) .
    !
    !     linpack. this version dated 08/14/78 .
    !     cleve moler, university of new mexico, argonne national lab.
    !
    !     subroutines and functions
    !
    !     linpack zgefa
    !     blas zaxpy,zdotc,zdscal,dzasum
    !     fortran dabs,dmax1,dcmplx,dconjg
    !
    !     internal variables
    !
    complex*16 ek,t,wk,wkm!zdotc
    double precision anorm,s,sm,ynorm!,dzasum
    integer info,j,k,kb,kp1,l
    !
    complex*16 zdum,zdum1,zdum2,csign1
    double precision cabs1
    double precision dreal,dimag
    complex*16 zdumr,zdumi
    dreal(zdumr) = zdumr
    dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
    cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
    csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
    !
    !     compute 1-norm of a
    !
    anorm = 0.0d0
    do 10 j = 1, n
       anorm = dmax1(anorm,dzasum(n,a(1,j),1))
10     continue
       !
       !     factor
       !
       call zgefa(a,lda,n,ipvt,info)
       !
       !     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
       !     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
       !     ctrans(a)  is the conjugate transpose of a .
       !     the components of  e  are chosen to cause maximum local
       !     growth in the elements of w  where  ctrans(u)*w = e .
       !     the vectors are frequently rescaled to avoid overflow.
       !
       !     solve ctrans(u)*w = e
       !
       ek = (1.0d0,0.0d0)
       do 20 j = 1, n
          z(j) = (0.0d0,0.0d0)
20        continue
          do 100 k = 1, n
             if (cabs1(z(k)) .ne. 0.0d0) ek = csign1(ek,-z(k))
             if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
             s = cabs1(a(k,k))/cabs1(ek-z(k))
             call zdscal(n,s,z,1)
             ek = dcmplx(s,0.0d0)*ek
30           continue
             wk = ek - z(k)
             wkm = -ek - z(k)
             s = cabs1(wk)
             sm = cabs1(wkm)
             if (cabs1(a(k,k)) .eq. 0.0d0) go to 40
             wk = wk/dconjg(a(k,k))
             wkm = wkm/dconjg(a(k,k))
             go to 50
40           continue
             wk = (1.0d0,0.0d0)
             wkm = (1.0d0,0.0d0)
50           continue
             kp1 = k + 1
             if (kp1 .gt. n) go to 90
             do 60 j = kp1, n
                sm = sm + cabs1(z(j)+wkm*dconjg(a(k,j)))
                z(j) = z(j) + wk*dconjg(a(k,j))
                s = s + cabs1(z(j))
60              continue
                if (s .ge. sm) go to 80
                t = wkm - wk
                wk = wkm
                do 70 j = kp1, n
                   z(j) = z(j) + t*dconjg(a(k,j))
70                 continue
80                 continue
90                 continue
                   z(k) = wk
100                continue
                   s = 1.0d0/dzasum(n,z,1)
                   call zdscal(n,s,z,1)
                   !
                   !     solve ctrans(l)*y = w
                   !
                   do 120 kb = 1, n
                      k = n + 1 - kb
                      if (k .lt. n) z(k) = z(k) + zdotc(n-k,a(k+1,k),1,z(k+1),1)
                      if (cabs1(z(k)) .le. 1.0d0) go to 110
                      s = 1.0d0/cabs1(z(k))
                      call zdscal(n,s,z,1)
110                   continue
                      l = ipvt(k)
                      t = z(l)
                      z(l) = z(k)
                      z(k) = t
120                   continue
                      s = 1.0d0/dzasum(n,z,1)
                      call zdscal(n,s,z,1)
                      !
                      ynorm = 1.0d0
                      !
                      !     solve l*v = y
                      !
                      do 140 k = 1, n
                         l = ipvt(k)
                         t = z(l)
                         z(l) = z(k)
                         z(k) = t
                         if (k .lt. n) call zaxpy(n-k,t,a(k+1,k),1,z(k+1),1)
                         if (cabs1(z(k)) .le. 1.0d0) go to 130
                         s = 1.0d0/cabs1(z(k))
                         call zdscal(n,s,z,1)
                         ynorm = s*ynorm
130                      continue
140                      continue
                         s = 1.0d0/dzasum(n,z,1)
                         call zdscal(n,s,z,1)
                         ynorm = s*ynorm
                         !
                         !     solve  u*z = v
                         !
                         do 160 kb = 1, n
                            k = n + 1 - kb
                            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
                            s = cabs1(a(k,k))/cabs1(z(k))
                            call zdscal(n,s,z,1)
                            ynorm = s*ynorm
150                         continue
                            if (cabs1(a(k,k)) .ne. 0.0d0) z(k) = z(k)/a(k,k)
                            if (cabs1(a(k,k)) .eq. 0.0d0) z(k) = (1.0d0,0.0d0)
                            t = -z(k)
                            call zaxpy(k-1,t,a(1,k),1,z(1),1)
160                         continue
                            !     make znorm = 1.0
                            s = 1.0d0/dzasum(n,z,1)
                            call zdscal(n,s,z,1)
                            ynorm = s*ynorm
                            !
                            if (anorm .ne. 0.0d0) rcond = ynorm/anorm
                            if (anorm .eq. 0.0d0) rcond = 0.0d0
                            return
                          end subroutine



                          !*******************************************************************
                          subroutine zgefa(a,lda,n,ipvt,info)
                            !*******************************************************************
                            integer lda,n,ipvt(1),info
                            complex*16 a(lda,1)
                            !
                            !     zgefa factors a complex*16 matrix by gaussian elimination.
                            !
                            !     zgefa is usually called by zgeco, but it can be called
                            !     directly with a saving in time if  rcond  is not needed.
                            !     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
                            !
                            !     on entry
                            !
                            !        a       complex*16(lda, n)
                            !                the matrix to be factored.
                            !
                            !        lda     integer
                            !                the leading dimension of the array  a .
                            !
                            !        n       integer
                            !                the order of the matrix  a .
                            !
                            !     on return
                            !
                            !        a       an upper triangular matrix and the multipliers
                            !                which were used to obtain it.
                            !                the factorization can be written  a = l*u  where
                            !                l  is a product of permutation and unit lower
                            !                triangular matrices and  u  is upper triangular.
                            !
                            !        ipvt    integer(n)
                            !                an integer vector of pivot indices.
                            !
                            !        info    integer
                            !                = 0  normal value.
                            !                = k  if  u(k,k) .eq. 0.0 .  this is not an error
                            !                     condition for this subroutine, but it does
                            !                     indicate that zgesl or zgedi will divide by zero
                            !                     if called.  use  rcond  in zgeco for a reliable
                            !                     indication of singularity.
                            !
                            !     linpack. this version dated 08/14/78 .
                            !     cleve moler, university of new mexico, argonne national lab.
                            !
                            !     subroutines and functions
                            !
                            !     blas zaxpy,zscal,izamax
                            !     fortran dabs
                            !
                            !     internal variables
                            !
                            complex*16 t
                            integer j,k,kp1,l,nm1!izamax
                            !
                            complex*16 zdum
                            double precision cabs1
                            double precision dreal,dimag
                            complex*16 zdumr,zdumi
                            dreal(zdumr) = zdumr
                            dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
                            cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
                            !
                            !     gaussian elimination with partial pivoting
                            !
                            info = 0
                            nm1 = n - 1
                            if (nm1 .lt. 1) go to 70
                            do 60 k = 1, nm1
                               kp1 = k + 1
                               !
                               !        find l = pivot index
                               !
                               l = izamax(n-k+1,a(k,k),1) + k - 1
                               ipvt(k) = l
                               !
                               !        zero pivot implies this column already triangularized
                               !
                               if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
                               !
                               !           interchange if necessary
                               !
                               if (l .eq. k) go to 10
                               t = a(l,k)
                               a(l,k) = a(k,k)
                               a(k,k) = t
10                             continue
                               !
                               !           compute multipliers
                               !
                               t = -(1.0d0,0.0d0)/a(k,k)
                               call zscal(n-k,t,a(k+1,k),1)
                               !
                               !           row elimination with column indexing
                               !
                               do 30 j = kp1, n
                                  t = a(l,j)
                                  if (l .eq. k) go to 20
                                  a(l,j) = a(k,j)
                                  a(k,j) = t
20                                continue
                                  call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30                                continue
                                  go to 50
40                                continue
                                  info = k
50                                continue
60                                continue
70                                continue
                                  ipvt(n) = n
                                  if (cabs1(a(n,n)) .eq. 0.0d0) info = n
                                  return
                                end subroutine

                                double complex function zdotc(n,zx,incx,zy,incy)
                                  !c
                                  !c     forms the dot product of a vector.
                                  !c     jack dongarra, 3/11/78.
                                  !c     modified 12/3/93, array(1) declarations changed to array(*)
                                  !c
                                  double complex zx(*),zy(*),ztemp
                                  integer i,incx,incy,ix,iy,n
                                  ztemp = (0.0d0,0.0d0)
                                  zdotc = (0.0d0,0.0d0)
                                  if(n.le.0)return
                                  if(incx.eq.1.and.incy.eq.1)go to 20
                                  !c
                                  !c        code for unequal increments or equal increments
                                  !c          not equal to 1
                                  !c
                                  ix = 1
                                  iy = 1
                                  if(incx.lt.0)ix = (-n+1)*incx + 1
                                  if(incy.lt.0)iy = (-n+1)*incy + 1
                                  do 10 i = 1,n
                                     ztemp = ztemp + dconjg(zx(ix))*zy(iy)
                                     ix = ix + incx
                                     iy = iy + incy
10                                   continue
                                     zdotc = ztemp
                                     return
                                     !c
                                     !c        code for both increments equal to 1
                                     !c
20                                   do 30 i = 1,n
                                        ztemp = ztemp + dconjg(zx(i))*zy(i)
30                                      continue
                                        zdotc = ztemp
                                        return
                                      end function

                                      subroutine  zdscal(n,da,zx,incx)
                                        !c
                                        !c     scales a vector by a constant.
                                        !c     jack dongarra, 3/11/78.
                                        !c     modified 3/93 to return if incx .le. 0.
                                        !c     modified 12/3/93, array(1) declarations changed to array(*)
                                        !c
                                        double complex zx(*)
                                        double precision da
                                        integer i,incx,ix,n
                                        !c
                                        if( n.le.0 .or. incx.le.0 )return
                                        if(incx.eq.1)go to 20
                                        !c
                                        !c        code for increment not equal to 1
                                        !c
                                        ix = 1
                                        do 10 i = 1,n
                                           zx(ix) = dcmplx(da,0.0d0)*zx(ix)
                                           ix = ix + incx
10                                         continue
                                           return
                                           !c
                                           !c        code for increment equal to 1
                                           !c
20                                         do 30 i = 1,n
                                              zx(i) = dcmplx(da,0.0d0)*zx(i)
30                                            continue
                                              return
                                            end subroutine

                                            subroutine zaxpy(n,za,zx,incx,zy,incy)
                                              !c
                                              !c     constant times a vector plus a vector.
                                              !c     jack dongarra, 3/11/78.
                                              !c     modified 12/3/93, array(1) declarations changed to array(*)
                                              !c
                                              double complex zx(*),zy(*),za
                                              integer i,incx,incy,ix,iy,n
                                              !double precision dcabs1
                                              if(n.le.0)return
                                              if (dcabs1(za) .eq. 0.0d0) return
                                              if (incx.eq.1.and.incy.eq.1)go to 20
                                              !c
                                              !c        code for unequal increments or equal increments
                                              !c          not equal to 1
                                              !c
                                              ix = 1
                                              iy = 1
                                              if(incx.lt.0)ix = (-n+1)*incx + 1
                                              if(incy.lt.0)iy = (-n+1)*incy + 1
                                              do 10 i = 1,n
                                                 zy(iy) = zy(iy) + za*zx(ix)
                                                 ix = ix + incx
                                                 iy = iy + incy
10                                               continue
                                                 return
                                                 !c
                                                 !c        code for both increments equal to 1
                                                 !c
20                                               do 30 i = 1,n
                                                    zy(i) = zy(i) + za*zx(i)
30                                                  continue
                                                    return
                                                  end subroutine


                                                  integer function izamax(n,zx,incx)
                                                    !c
                                                    !c     finds the index of element having max. absolute value.
                                                    !c     jack dongarra, 1/15/85.
                                                    !c     modified 3/93 to return if incx .le. 0.
                                                    !c     modified 12/3/93, array(1) declarations changed to array(*)
                                                    !c
                                                    double complex zx(*)
                                                    double precision smax
                                                    integer i,incx,ix,n
                                                    !double precision dcabs1
                                                    !c
                                                    izamax = 0
                                                    if( n.lt.1 .or. incx.le.0 )return
                                                    izamax = 1
                                                    if(n.eq.1)return
                                                    if(incx.eq.1)go to 20
                                                    !c
                                                    !c        code for increment not equal to 1
                                                    !c
                                                    ix = 1
                                                    smax = dcabs1(zx(1))
                                                    ix = ix + incx
                                                    do 10 i = 2,n
                                                       if(dcabs1(zx(ix)).le.smax) go to 5
                                                       izamax = i
                                                       smax = dcabs1(zx(ix))
5                                                      ix = ix + incx
10                                                     continue
                                                       return
                                                       !!c
                                                       !c        code for increment equal to 1
                                                       !c
20                                                     smax = dcabs1(zx(1))
                                                       do 30 i = 2,n
                                                          if(dcabs1(zx(i)).le.smax) go to 30
                                                          izamax = i
                                                          smax = dcabs1(zx(i))
30                                                        continue
                                                          return
                                                        end function

                                                        subroutine  zscal(n,za,zx,incx)
                                                          !c
                                                          !c     scales a vector by a constant.
                                                          !    jack dongarra, 3/11/78.
                                                          !c     modified 3/93 to return if incx .le. 0.
                                                          !c     modified 12/3/93, array(1) declarations changed to array(*)
                                                          !c
                                                          double complex za,zx(*)
                                                          integer i,incx,ix,n
                                                          !c
                                                          if( n.le.0 .or. incx.le.0 )return
                                                          if(incx.eq.1)go to 20
                                                          !c
                                                          !c        code for increment not equal to 1
                                                          !c
                                                          ix = 1
                                                          do 10 i = 1,n
                                                             zx(ix) = za*zx(ix)
                                                             ix = ix + incx
10                                                           continue
                                                             return
                                                             !c
                                                             !c        code for increment equal to 1
                                                             !c
20                                                           do 30 i = 1,n
                                                                zx(i) = za*zx(i)
30                                                              continue
                                                                return
                                                              end subroutine

                                                              double precision function dcabs1(z)
                                                                double complex z,zz
                                                                double precision t(2)
                                                                equivalence (zz,t(1))
                                                                zz = z
                                                                dcabs1 = dabs(t(1)) + dabs(t(2))
                                                                return
                                                              end function dcabs1

                                                              double precision function dzasum(n,zx,incx)
                                                                !c
                                                                !c     takes the sum of the absolute values.
                                                                !c     jack dongarra, 3/11/78.
                                                                !!c     modified 3/93 to return if incx .le. 0.
                                                                !c     modified 12/3/93, array(1) declarations changed to array(*)
                                                                !c
                                                                double complex zx(*)
                                                                double precision stemp!,dcabs1
                                                                integer i,incx,ix,n
                                                                !c
                                                                dzasum = 0.0d0
                                                                stemp = 0.0d0
                                                                if( n.le.0 .or. incx.le.0 )return
                                                                if(incx.eq.1)go to 20
                                                                !c
                                                                !c        code for increment not equal to 1
                                                                !c
                                                                ix = 1
                                                                do 10 i = 1,n
                                                                   stemp = stemp + dcabs1(zx(ix))
                                                                   ix = ix + incx
10                                                                 continue
                                                                   dzasum = stemp
                                                                   return
                                                                   !c
                                                                   !c        code for increment equal to 1
                                                                   !c
20                                                                 do 30 i = 1,n
                                                                      stemp = stemp + dcabs1(zx(i))
30                                                                    continue
                                                                      dzasum = stemp
                                                                      return
                                                                    end function

                                                                    !c***************************************************
                                                                    subroutine zgedi(a,lda,n,ipvt,det,work,job)
                                                                      !c***************************************************
                                                                      integer lda,ipvt(1),job
                                                                      integer*4 n
                                                                      complex*16 a(lda,1),det(2),work(1)
                                                                      !c
                                                                      !c     zgedi computes the determinant and inverse of a matrix
                                                                      !c     using the factors computed by zgeco or zgefa.
                                                                      !c
                                                                      !c     on entry
                                                                      !c
                                                                      !c        a       complex*16(lda, n)
                                                                      !c                the output from zgeco or zgefa.
                                                                      !c
                                                                      !c        lda     integer
                                                                      !c                the leading dimension of the array  a .
                                                                      !c
                                                                      !c        n       integer
                                                                      !c                the order of the matrix  a .
                                                                      !c
                                                                      !c        ipvt    integer(n)
                                                                      !c                the pivot vector from zgeco or zgefa.
                                                                      !c
                                                                      !c        work    complex*16(n)
                                                                      !c                work vector.  contents destroyed.
                                                                      !c
                                                                      !c        job     integer
                                                                      !c                = 11   both determinant and inverse.
                                                                      !c                = 01   inverse only.
                                                                      !c                = 10   determinant only.
                                                                      !c
                                                                      !c     on return
                                                                      !c
                                                                      !c        a       inverse of original matrix if requested.
                                                                      !c                otherwise unchanged.
                                                                      !c
                                                                      !c        det     complex*16(2)
                                                                      !c                determinant of original matrix if requested.
                                                                      !c                otherwise not referenced.
                                                                      !c                determinant = det(1) * 10.0**det(2)
                                                                      !c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
                                                                      !c                or  det(1) .eq. 0.0 .
                                                                      !c
                                                                      !c     error condition
                                                                      !c
                                                                      !c        a division by zero will occur if the input factor contains
                                                                      !c        a zero on the diagonal and the inverse is requested.
                                                                      !c        it will not occur if the subroutines are called correctly
                                                                      !c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
                                                                      !c        info .eq. 0 .
                                                                      !
                                                                      !c     linpack. this version dated 08/14/78 .
                                                                      !c     cleve moler, university of new mexico, argonne national lab.
                                                                      !c
                                                                      !c     subroutines and functions
                                                                      !c
                                                                      !c     blas zaxpy,zscal,zswap
                                                                      !c     fortran dabs,dcmplx,mod
                                                                      !c
                                                                      !c     internal variables
                                                                      !c
                                                                      complex*16 t
                                                                      double precision ten
                                                                      integer i,j,k,kb,kp1,l,nm1
                                                                      !c
                                                                      complex*16 zdum
                                                                      double precision cabs1
                                                                      double precision dreal,dimag
                                                                      complex*16 zdumr,zdumi
                                                                      dreal(zdumr) = zdumr
                                                                      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
                                                                      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
                                                                      !c
                                                                      !c     compute determinant
                                                                      !c
                                                                      if (job/10 .eq. 0) go to 70
                                                                      det(1) = (1.0d0,0.0d0)
                                                                      det(2) = (0.0d0,0.0d0)
                                                                      ten = 10.0d0
                                                                      do 50 i = 1, n
                                                                         if (ipvt(i) .ne. i) det(1) = -det(1)
                                                                         det(1) = a(i,i)*det(1)
                                                                         !c        ...exit
                                                                         if (cabs1(det(1)) .eq. 0.0d0) go to 60
10                                                                       if (cabs1(det(1)) .ge. 1.0d0) go to 20
                                                                         det(1) = dcmplx(ten,0.0d0)*det(1)
                                                                         det(2) = det(2) - (1.0d0,0.0d0)
                                                                         go to 10
20                                                                       continue
30                                                                       if (cabs1(det(1)) .lt. ten) go to 40
                                                                         det(1) = det(1)/dcmplx(ten,0.0d0)
                                                                         det(2) = det(2) + (1.0d0,0.0d0)
                                                                         go to 30
40                                                                       continue
50                                                                       continue
60                                                                       continue
70                                                                       continue
                                                                         !c
                                                                         !c     compute inverse(u)
                                                                         !c
                                                                         if (mod(job,10) .eq. 0) go to 150
                                                                         do 100 k = 1, n
                                                                            a(k,k) = (1.0d0,0.0d0)/a(k,k)
                                                                            t = -a(k,k)
                                                                            call zscal(k-1,t,a(1,k),1)
                                                                            kp1 = k + 1
                                                                            if (n .lt. kp1) go to 90
                                                                            do 80 j = kp1, n
                                                                               t = a(k,j)
                                                                               a(k,j) = (0.0d0,0.0d0)
                                                                               call zaxpy(k,t,a(1,k),1,a(1,j),1)
80                                                                             continue
90                                                                             continue
100                                                                            continue
                                                                               !c
                                                                               !c        form inverse(u)*inverse(l)
                                                                               !c
                                                                               nm1 = n - 1
                                                                               if (nm1 .lt. 1) go to 140
                                                                               do 130 kb = 1, nm1
                                                                                  k = n - kb
                                                                                  kp1 = k + 1
                                                                                  do 110 i = kp1, n
                                                                                     work(i) = a(i,k)
                                                                                     a(i,k) = (0.0d0,0.0d0)
110                                                                                  continue
                                                                                     do 120 j = kp1, n
                                                                                        t = work(j)
                                                                                        call zaxpy(n,t,a(1,j),1,a(1,k),1)
120                                                                                     continue
                                                                                        l = ipvt(k)
                                                                                        if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
130                                                                                     continue
140                                                                                     continue
150                                                                                     continue
                                                                                        return
                                                                                      end subroutine

                                                                                      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
                                                                                        !*     .. Scalar Arguments ..
                                                                                        INTEGER INCX,INCY,N
                                                                                        !*     ..
                                                                                        !*     .. Array Arguments ..
                                                                                        DOUBLE COMPLEX ZX(*),ZY(*)
                                                                                        !*     ..
                                                                                        !*
                                                                                        !*  Purpose
                                                                                        !*  =======
                                                                                        !*
                                                                                        !*     ZSWAP interchanges two vectors.
                                                                                        !*
                                                                                        !*  Further Details
                                                                                        !*  ===============
                                                                                        !*
                                                                                        !*     jack dongarra, 3/11/78.
                                                                                        !*     modified 12/3/93, array(1) declarations changed to array(*)
                                                                                        !*
                                                                                        !*  =====================================================================
                                                                                        !*
                                                                                        !*     .. Local Scalars ..
                                                                                        DOUBLE COMPLEX ZTEMP
                                                                                        INTEGER I,IX,IY
                                                                                        !*     ..
                                                                                        IF (N.LE.0) RETURN
                                                                                        IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
                                                                                           !*
                                                                                           !*       code for both increments equal to 1
                                                                                           DO I = 1,N
                                                                                              ZTEMP = ZX(I)
                                                                                              ZX(I) = ZY(I)
                                                                                              ZY(I) = ZTEMP
                                                                                           END DO
                                                                                        ELSE
                                                                                           !*
                                                                                           !*       code for unequal increments or equal increments not equal
                                                                                           !*         to 1
                                                                                           !*
                                                                                           IX = 1
                                                                                           IY = 1
                                                                                           IF (INCX.LT.0) IX = (-N+1)*INCX + 1
                                                                                           IF (INCY.LT.0) IY = (-N+1)*INCY + 1
                                                                                           DO I = 1,N
                                                                                              ZTEMP = ZX(IX)
                                                                                              ZX(IX) = ZY(IY)
                                                                                              ZY(IY) = ZTEMP
                                                                                              IX = IX + INCX
                                                                                              IY = IY + INCY
                                                                                           END DO
                                                                                        END IF
                                                                                        RETURN
                                                                                      END subroutine ZSWAP








                                                                                    end module ED_GLOC
