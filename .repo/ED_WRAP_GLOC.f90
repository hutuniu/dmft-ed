module ED_WRAP_GLOC
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS,   only: reg,txtfy,store_data
  USE SF_LINALG,    only:eye,inv,inv_sym
  USE SF_ARRAYS,    only:linspace,arange
  USE SF_MISC,      only:assert_shape
  implicit none
  private
  !
  interface ed_get_gloc_lattice
     module procedure ed_get_gloc_normal_lattice_main
     module procedure ed_get_gloc_superc_lattice_main
     module procedure ed_get_gloc_normal_lattice_1b
     module procedure ed_get_gloc_normal_lattice_mb
     module procedure ed_get_gloc_superc_lattice_1b
     module procedure ed_get_gloc_superc_lattice_mb
  end interface ed_get_gloc_lattice
  public :: ed_get_gloc_lattice


  interface ed_get_gij_lattice
     module procedure ed_get_gij_normal_main
     module procedure ed_get_gij_superc_main
  end interface ed_get_gij_lattice
  public :: ed_get_gij_lattice


  real(8),dimension(:),allocatable :: wr,wm
  character(len=20)                :: suffix



contains




  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the Normal local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_normal_lattice_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm,Gamma_mats,Gamma_real)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint    !
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    complex(8),dimension(:,:,:),optional            :: Gamma_mats(size(Hk,1),size(Hk,2),Lmats)
    complex(8),dimension(:,:,:),optional            :: Gamma_real(size(Hk,1),size(Hk,2),Lreal)
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
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
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
    ! zeta_mats=zero
    ! zeta_real=zero
    ! do ilat=1,Nlat
    !    do ispin=1,Nspin
    !       do iorb=1,Norb
    !          io = iorb + (ispin-1)*Norb
    !          js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
    !          zeta_mats(ilat,io,io,:) = xi*wm(:)       + xmu! - Eloc_(js)
    !          zeta_real(ilat,io,io,:) = wr(:) + xi*eps + xmu! - Eloc_(js)
    !       enddo
    !    enddo
    !    do ispin=1,Nspin
    !       do jspin=1,Nspin
    !          do iorb=1,Norb
    !             do jorb=1,Norb
    !                io = iorb + (ispin-1)*Norb
    !                jo = jorb + (jspin-1)*Norb
    !                zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Smats(ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_real(ilat,io,jo,:) = zeta_real(ilat,io,jo,:) - Sreal(ilat,ispin,jspin,iorb,jorb,:)
    !             enddo
    !          enddo
    !       enddo
    !    enddo
    ! enddo
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
    if(mpiID==0)call start_timer
    Gmats=zero
    Greal=zero
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
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)call print_gloc_lattice(Gmats,Greal,"LG",iprint)
  end subroutine ed_get_gloc_normal_lattice_main

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
    call assert_shape(zeta,[Nlat,Nlso,Nlso,Lfreq],"add_to_gloc_lattice_normal","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gloc_lattice_normal","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_lattice_normal","Gkout")
    if(present(Gembed))&
         call assert_shape(Gembed,[Nlso,Nlso,Lfreq],"add_to_gloc_lattice_normal","Gembed")   
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gloc_normal_lattice













  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the Nambu local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
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
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_ID==0)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
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
    ! zeta_mats=zero
    ! zeta_real=zero
    ! do ilat=1,Nlat
    !    do ispin=1,Nspin
    !       do iorb=1,Norb
    !          io = iorb + (ispin-1)*Norb
    !          js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
    !          !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
    !          !G22(iw) = -[G11[iw]]*
    !          !G21(iw) =   G12[w]
    !          zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu !- Eloc_(js)
    !          zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu !+ Eloc_(js)
    !          !
    !          !SYMMETRIES in real-frequencies   [assuming a real order parameter]
    !          !G22(w)  = -[G11[-w]]*
    !          !G21(w)  =   G12[w]             
    !          zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu !- Eloc_(js)
    !          zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )!- Eloc_(js) )
    !       enddo
    !    enddo
    !    do ispin=1,Nspin
    !       do jspin=1,Nspin
    !          do iorb=1,Norb
    !             do jorb=1,Norb
    !                io = iorb + (ispin-1)*Norb
    !                jo = jorb + (jspin-1)*Norb
    !                zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
    !                !
    !                zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
    !                zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
    !             enddo
    !          enddo
    !       enddo
    !    enddo
    ! enddo
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
    if(mpiID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_superc_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       call print_gloc_lattice(Gmats(1,:,:,:,:,:,:),Greal(1,:,:,:,:,:,:),"LG",iprint)
       call print_gloc_lattice(Gmats(2,:,:,:,:,:,:),Greal(2,:,:,:,:,:,:),"LF",iprint)
    endif
  end subroutine ed_get_gloc_superc_lattice_main

  subroutine add_to_gloc_superc_lattice(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm                
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
    call assert_shape(zeta,[2,2,Nlat,Nlso,Nlso,Lfreq],"add_to_gloc_lattice_superc","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"add_to_gloc_lattice_superc","Hk")
    call assert_shape(Gkout,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"add_to_gloc_lattice_superc","Gkout")
    !
    allocate(Gktmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gkout = zero
    Gktmp = zero
    do i=1+mpiID,Lfreq,mpiSIZE
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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gloc_superc_lattice




















  !+-----------------------------------------------------------------------------+!
  !PURPOSE: additional interfaces to the main procedures above to get Gloc 
  ! _1b : one band case (independent of Nspin and Norb)
  ! _mb : multi-bands case (independent of Nspin)
  !+-----------------------------------------------------------------------------+!
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
    call assert_shape(Gmats,[Nlat,Norb,Norb,Lreal],"ed_get_gloc_normal_mb","Greal")
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









  !>>>>> WARNING: THESE ARE NOT UPDATED <<<<<<
  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate all components of the Normal Green's function for a given 
  ! Hamiltonian matrix and self-energy functions. Hk is a big sparse matrix of the form 
  ! H(\k;R_i,R_j)_{ab}^{ss'} and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gij_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,hk_symm,Gamma_mats,Gamma_real)
    complex(8),dimension(:,:,:)              :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                  :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)                 :: Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)                 :: Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)                 :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)                 :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                                  :: iprint
    complex(8)                               :: zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                               :: zeta_real(Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                               :: Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                               :: Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    logical,optional                         :: hk_symm(size(Hk,3))
    logical                                  :: hk_symm_(size(Hk,3))
    complex(8),optional                      :: Gamma_mats(size(Hk,1),size(Hk,2),Lmats)
    complex(8),optional                      :: Gamma_real(size(Hk,1),size(Hk,2),Lreal)
    integer                                  :: i,ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    Lk=size(Hk,3)
    Nlso=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpiID==0)write(LOGfile,*)"Get full GF (id=0):"
    !here we create the "array" *zeta* of Nlat blocks, each of size (Nspin*Norb)
    !then we use a newly created function *blocks_to_matrix* to spread the blocks into
    !a matrix of rank 2 dimensions Nlso*Nlso
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    !
    zeta_mats=zero
    zeta_real=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(ilat,io,io,:) = xi*wm(:)       + xmu !- Eloc_(js)
             zeta_real(ilat,io,io,:) = wr(:) + xi*eps + xmu !- Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(ilat,io,jo,:) = zeta_mats(ilat,io,jo,:) - Smats(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(ilat,io,jo,:) = zeta_real(ilat,io,jo,:) - Sreal(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpiID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
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
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       call print_gij_lattice(Gmats,Greal,"Gij",iprint)
    endif
  end subroutine ed_get_gij_normal_main

  subroutine add_to_gij_normal(zeta,Hk,hk_symm,Gkout,Gembed)
    complex(8)                               :: zeta(:,:,:,:)              ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8)                               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                                  :: hk_symm
    !output:
    complex(8),intent(inout)                 :: Gkout(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,4))
    complex(8)                               :: Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta,4))
    !
    complex(8)                               :: Gmatrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex(8),optional                      :: Gembed(Nlat*Nspin*Norb,Nlat*Nspin*Norb,size(zeta,4)) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    integer                                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    if(size(zeta,1)/=Nlat)stop "add_to_gij_normal error: zeta wrong size 1 = Nlat"
    if(size(zeta,2)/=Nspin*Norb)stop "add_to_gij_normal error: zeta wrong size 2 = Nspin*Norb"
    if(size(zeta,3)/=Nspin*Norb)stop "add_to_gij_normal error: zeta wrong size 3 = Nspin*Norb"
    Lfreq = size(zeta,4)
    Gktmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
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
       !
    enddo
    Gkout=zero
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gij_normal





  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate all the components of the Nambu Green's function for a given 
  ! Hamiltonian matrix and self-energy functions. Hk is a big sparse matrix of the form 
  ! H(k;R_i,R_j)_{ab}^{ss'} and size [Nk]*[Nlat*Nspin*Norb]**2
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
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
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
    if(mpiID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gij_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
       call add_to_gij_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Fmats = Fmats + Fkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       Freal = Freal + Fkreal*Wtk(ik)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       call print_gij_lattice(Gmats,Greal,"Gij",iprint)
       call print_gij_lattice(Fmats,Freal,"Fij",iprint)
    endif
  end subroutine ed_get_gij_superc_main

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
    do i=1+mpiID,Lfreq,mpiSIZE
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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Gktmp,Gkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Fktmp,Fkout,size(Gkout),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Gkout = Gktmp
    Fkout = Fktmp
#endif
  end subroutine add_to_gij_superc








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: print a local GF according to iprint variable
  !+-----------------------------------------------------------------------------+!
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

end module ED_WRAP_GLOC
