!###############################################################
! PROGRAM  : RDMFT_FUNX
! PURPOSE  : Compute Local GFunction for generel real-space scheme
!###############################################################
module ED_WRAP_GLOC
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS,   only:reg,sread,free_unit
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym,matrix_diagonalize
  USE SF_ARRAYS,    only:linspace,arange
  implicit none
  private
  !
  interface ed_get_gloc_lattice
     module procedure ed_get_gloc_normal,ed_get_gloc_superc
  end interface ed_get_gloc_lattice
  public :: ed_get_gloc_lattice

  ! !OBSOLETE
  ! interface rdmft_get_gloc
  !    module procedure rdmft_get_gloc_normal,rdmft_get_gloc_superc
  ! end interface rdmft_get_gloc
  ! public :: rdmft_get_gloc
  ! public :: rdmft_get_gloc_mb
  ! public :: rdmft_get_gloc_embedd
  ! ! interface ed_get_gloc_superc
  ! !    module procedure rdmft_get_gloc_superc
  ! ! end interface ed_get_gloc_superc


  real(8),dimension(:),allocatable        :: wr,wm,tau

contains


  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the Normal local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_normal(Hk,Wtk,Gmats,Greal,Smats,Sreal,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: ik,Lk,Nside,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    Lk=size(Hk,3)
    Nside=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nside.OR.size(Hk,2)/=Nside) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    Gmats=zero
    Greal=zero
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
    !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
    !then we use a newly created function *blocks_to_matrix* to spread the blocks into
    !a matrix of rank 2 dimensions Nside*Nside
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
             zeta_mats(ilat,io,io,:) = xi*wm(:)       + xmu - Eloc_(js)
             zeta_real(ilat,io,io,:) = wr(:) + xi*eps + xmu - Eloc_(js)
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
    !one you have *zeta_site* array you pass it to the routines that invert 
    !G for each k-point
    if(mpiID==0)call start_timer
    do ik=1,Lk
       call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Gmats)
       call add_to_gloc_normal(zeta_real,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Greal)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
  end subroutine ed_get_gloc_normal

  subroutine add_to_gloc_normal(zeta_site,Hk,Wtk,hk_symm,Gloc)
    complex(8)               :: zeta_site(:,:,:,:)              ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    real(8)                  :: Wtk                    
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gloc(Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
    !
    complex(8)               :: Gmatrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex(8)               :: Gk_tmp(Nlat,Nspin*Norb,Nspin*Norb,size(zeta_site,4))
    complex(8)               :: Gk(Nlat,Nspin*Norb,Nspin*Norb,size(zeta_site,4))
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    if(size(zeta_site,1)/=Nlat)stop "get_gloc_kpoint error: zeta_site wrong size 1 = Nlat"
    if(size(zeta_site,2)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 2 = Nspin*Norb"
    if(size(zeta_site,3)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 3 = Nspin*Norb"
    Lfreq = size(zeta_site,4)
    Gk    =zero
    Gk_tmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
       Gmatrix  = zero
       Gmatrix  = blocks_to_matrix(zeta_site(:,:,:,i)) - Hk
       if(hk_symm) then
          call matrix_inverse_sym(Gmatrix)
       else
          call matrix_inverse(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       Gk_tmp(:,:,:,i) = matrix_to_blocks(Gmatrix)
    enddo
    call MPI_ALLREDUCE(Gk_tmp,Gk,size(Gk),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !add to k-summation
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gloc(ilat,ispin,jspin,iorb,jorb,:) = Gloc(ilat,ispin,jspin,iorb,jorb,:) + Gk(ilat,io,jo,:)*Wtk
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine add_to_gloc_normal







  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the Nambu local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_superc(Hk,Wtk,Gmats,Greal,Smats,Sreal,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    Lk=size(Hk,3)
    Nlso=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    Gmats=zero
    Greal=zero
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
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:)          + xmu - Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:)          - xmu + Eloc_(js)
             !
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu - Eloc_(js)
             zeta_real(2,2,ilat,io,io,:) = -conjg(dcmplx(wr(Lreal:1:-1),eps)) + xmu - Eloc_(js)
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
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - conjg(Sreal(2,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1))
                   zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) - Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(mpiID==0)call start_timer
    do ik=1,Lk
       call add_to_gloc_superc(zeta_mats,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Gmats)
       call add_to_gloc_superc(zeta_real,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Greal)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
  end subroutine ed_get_gloc_superc

  subroutine add_to_gloc_superc(zeta_site,Hk,Wtk,hk_symm,Gloc)
    complex(8)               :: zeta_site(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    real(8)                  :: Wtk                    
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gloc(2,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,6))
    !
    complex(8)               :: Gmatrix(2*Nlat*Nspin*Norb , 2*Nlat*Nspin*Norb)
    complex(8)               :: Gk_tmp(2,Nlat,Nspin*Norb,Nspin*Norb,size(zeta_site,6))
    complex(8)               :: Gk(2,Nlat,Nspin*Norb,Nspin*Norb,size(zeta_site,6))
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nlso
    if(size(zeta_site,1)/=2.OR.size(zeta_site,2)/=2)stop "add_to_gloc_superc error: zeta_site wrong size 1 or 2 = 2 (Nambu)"
    if(size(zeta_site,3)/=Nlat)stop "add_to_gloc_superc error: zeta_site wrong size 2 = Nlat"
    if(size(zeta_site,4)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 4 = Nspin*Norb"
    if(size(zeta_site,5)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 5 = Nspin*Norb"
    Lfreq = size(zeta_site,6)
    Nlso  = Nlat*Nspin*Norb
    Gk    =zero
    Gk_tmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta_site(1,1,:,:,:,i)) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_site(1,2,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta_site(2,1,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_site(2,2,:,:,:,i)) + Hk
       if(hk_symm) then
          call matrix_inverse_sym(Gmatrix)
       else
          call matrix_inverse(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       Gk_tmp(1,:,:,:,i) = matrix_to_blocks(Gmatrix(1:Nlso,1:Nlso))         !block 11 of size Nlso*Nlso
       Gk_tmp(2,:,:,:,i) = matrix_to_blocks(Gmatrix(1:Nlso,Nlso+1:2*Nlso))  !block 12 of size Nlso*Nlso
    enddo
    call MPI_ALLREDUCE(Gk_tmp,Gk,size(Gk),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !add to k-summation
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gloc(1,ilat,ispin,jspin,iorb,jorb,:) = Gloc(1,ilat,ispin,jspin,iorb,jorb,:) + Gk(1,ilat,io,jo,:)*Wtk
                   Gloc(2,ilat,ispin,jspin,iorb,jorb,:) = Gloc(2,ilat,ispin,jspin,iorb,jorb,:) + Gk(2,ilat,io,jo,:)*Wtk
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine add_to_gloc_superc



















  ! !###################################################################################################
  ! !###################################################################################################
  ! !               POSSIBLY OBSOLETE ROUTINES NOW SUPERSEDED BY TOP ROUTINE HERE
  ! !               GETTING GLOC FOR ANY Norb*Nspin*Nlat*Nk NORMAL GREEN'S FUNCTION
  ! !                (these routines are left temporarily for back-compatibility)
  ! !                                     ( to be removed)
  ! !###################################################################################################
  ! !###################################################################################################

  ! !+------------------------+!
  ! !+- SINGLE BAND ROUTINES -+!
  ! !+------------------------+!
  ! subroutine rdmft_get_gloc_normal(Hk,wt,Gmats,Greal,Smats,Sreal,Eloc,hk_symm)
  !   complex(8),dimension(:,:,:)            :: Hk
  !   real(8),dimension(size(Hk,3))          :: wt
  !   real(8),optional                       :: Eloc(Nlat)
  !   real(8)                                :: Eloc_(Nlat)
  !   complex(8),intent(inout)               :: Gmats(Nlat,Lmats)
  !   complex(8),intent(inout)               :: Greal(Nlat,Lreal)
  !   complex(8),intent(inout)               :: Smats(Nlat,Lmats)
  !   complex(8),intent(inout)               :: Sreal(Nlat,Lreal)
  !   logical,dimension(size(Hk,3)),optional :: hk_symm
  !   logical,dimension(size(Hk,3))          :: hk_symm_
  !   integer                                :: ik,Lk
  !   if(allocated(wm))deallocate(wm)
  !   if(allocated(wr))deallocate(wr)
  !   allocate(wm(Lmats))
  !   allocate(wr(Lreal))
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !   Lk=size(Hk,3)
  !   if(size(Hk,1)/=Nlat.or.size(Hk,2)/=Nlat) stop "get_gloc_k_lattice :Wrong dimensions Hk"
  !   eloc_=0.d0;if(present(eloc))eloc_=eloc
  !   hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !   Gmats=zero
  !   Greal=zero
  !   if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
  !   if(mpiID==0)call start_timer
  !   do ik=1,Lk
  !      call get_gloc_mats(eloc_,Smats,Gmats,Hk(:,:,ik),wt(ik),hk_symm_(ik))
  !      call get_gloc_real(eloc_,Sreal,Greal,Hk(:,:,ik),wt(ik),hk_symm_(ik))
  !      if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !   end do
  !   if(mpiID==0)call stop_timer
  ! end subroutine rdmft_get_gloc_normal

  ! !mats freq
  ! subroutine get_gloc_mats(elocal,sigma,fg,Hk,wtk,hk_symm)
  !   real(8)                  :: elocal(Nlat)
  !   complex(8)               :: Hk(Nlat,Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Lmats),sigma(Nlat,Lmats)
  !   real(8)                  :: ek,wtk
  !   logical                  :: hk_symm
  !   complex(8)               :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,Lmats),fg_k(Nlat,Lmats)
  !   integer                  :: i,is
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lmats,mpiSIZE
  !      Gloc  = zero
  !      zeta  = xi*wm(i) + xmu
  !      Gloc  = -Hk
  !      do is=1,Nlat
  !         Gloc(is,is) = Gloc(is,is) + zeta  - elocal(is) - sigma(is,i)
  !      enddo
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
  !      end if
  !      forall(is=1:Nlat)
  !         gf_tmp(is,i) = Gloc(is,is)
  !      end forall
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_mats

  ! ! real freq
  ! subroutine get_gloc_real(elocal,sigma,fg,Hk,wtk,hk_symm)
  !   real(8)                  :: elocal(Nlat)
  !   complex(8)               :: Hk(Nlat,Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Lreal),sigma(Nlat,Lreal)
  !   real(8)                  :: ek,wtk
  !   logical                  :: hk_symm
  !   complex(8)               :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,Lreal),fg_k(Nlat,Lreal)
  !   integer                  :: i,is
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lreal,mpiSIZE
  !      Gloc  = zero
  !      zeta  = cmplx(wr(i),eps,8) + xmu
  !      Gloc  = - Hk
  !      do is=1,Nlat
  !         Gloc(is,is) = Gloc(is,is) +  zeta  - elocal(is) - sigma(is,i)
  !      enddo
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)
  !      end if
  !      forall(is=1:Nlat)
  !         gf_tmp(is,i) = Gloc(is,is)
  !      end forall
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_real


  ! !+----------------------+!
  ! !+- MULTIBAND ROUTINES -+!   FOR THE MOMENT ONLY REPULSIVE CASE IS CONSIDERED
  ! !+----------------------+!
  ! subroutine rdmft_get_gloc_mb(Hk,wt,Hloc,Gmats,Greal,Smats,Sreal,hk_symm)
  !   complex(8)                       :: Hk(:,:,:)          ![Norb*Nlat][Norb*Nlat][Lk]
  !   real(8)                          :: wt(:)              ![Lk]
  !   complex(8)                       :: Hloc(:,:,:,:,:)    ![Nlat][Nspin][Nspin][Norb][Norb]
  !   complex(8),intent(inout)         :: Gmats(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !   complex(8),intent(inout)         :: Greal(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  !   complex(8),intent(inout)         :: Smats(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !   complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  !   logical,dimension(:),optional    :: hk_symm
  !   logical,dimension(:),allocatable :: hk_symm_
  !   integer                          :: ik,Lk,ispin
  !   integer                          :: No
  !   if(allocated(wm))deallocate(wm)
  !   if(allocated(wr))deallocate(wr)
  !   allocate(wm(Lmats))
  !   allocate(wr(Lreal))
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !   !checks!
  !   No = Norb*Nlat
  !   Lk=size(Hk,3)
  !   ! check Hk !
  !   if(size(Hk,1)/=No) stop "size Hk(:,*,*) /= No"
  !   if(size(Hk,2)/=size(Hk,1)) stop "size Hk(*,:,*) /= Hk(:,*,*) "
  !   ! check wt !
  !   if(size(wt)/=Lk) stop "size wt(*) /= Lk"
  !   ! check Hloc !
  !   if(size(Hloc,1)/=Nlat) stop "size Hloc(:,*,*,*,*) /= Nlat"
  !   if(size(Hloc,2)/=size(Hloc,3)) stop "size Hloc(*,:,*,*,*) /= Hloc(*,*,:,*,*) "
  !   if(size(Hloc,2)/=Nspin) stop "size Hloc(*,:,*,*,*) /= Nspin"
  !   if(size(Hloc,4)/=size(Hloc,5)) stop "size Hloc(*,*,*,:,*) /= Hloc(*,*,*,*,:) "
  !   if(size(Hloc,4)/=Norb) stop "size Hloc(*,*,*,:,*) /= Norb"
  !   ! check Smats !
  !   if(size(Smats,1)/=Nlat) stop "size Smats(:,*,*,*,*,*) /= Nlat"
  !   if(size(Smats,2)/=size(Smats,3)) stop "size Smats(*,:,*,*,*,*) /= Smats(*,*,:,*,*,*) "
  !   if(size(Smats,2)/=Nspin) stop "size Smats(*,:,*,*,*,*) /= Nspin"
  !   if(size(Smats,4)/=size(Smats,5)) stop "size Smats(*,*,*,:,*,*) /= Smats(*,*,*,*,:,*) "
  !   if(size(Smats,4)/=Norb) stop "size Smats(*,*,*,*,*,:) /= Norb"
  !   if(size(Smats,6)/=Lmats) stop "size Smats(*,*,*,*,*,:) /= Lmats"
  !   ! check Sreal !
  !   if(size(Sreal,1)/=Nlat) stop "size Sreal(:,*,*,*,*,*) /= Nlat"
  !   if(size(Sreal,2)/=size(Sreal,3)) stop "size Sreal(*,:,*,*,*,*) /= Sreal(*,*,:,*,*,*) "
  !   if(size(Sreal,2)/=Nspin) stop "size Sreal(*,:,*,*,*,*) /= Nspin"
  !   if(size(Sreal,4)/=size(Sreal,5)) stop "size Sreal(*,*,*,:,*,*) /= Sreal(*,*,*,*,:,*) "
  !   if(size(Sreal,4)/=Norb) stop "size Sreal(*,*,*,*,*,:) /= Norb"
  !   if(size(Sreal,6)/=Lmats) stop "size Sreal(*,*,*,*,*,:) /= Lmats"
  !   ! check Gmats !
  !   if(size(Gmats,1)/=Nlat) stop "size Gmats(:,*,*,*,*,*) /= Nlat"
  !   if(size(Gmats,2)/=size(Gmats,3)) stop "size Gmats(*,:,*,*,*,*) /= Gmats(*,*,:,*,*,*) "
  !   if(size(Gmats,2)/=Nspin) stop "size Gmats(*,:,*,*,*,*) /= Nspin"
  !   if(size(Gmats,4)/=size(Gmats,5)) stop "size Gmats(*,*,*,:,*,*) /= Gmats(*,*,*,*,:,*) "
  !   if(size(Gmats,4)/=Norb) stop "size Gmats(*,*,*,*,*,:) /= Norb"
  !   if(size(Gmats,6)/=Lmats) stop "size Gmats(*,*,*,*,*,:) /= Lmats"
  !   ! check Greal !
  !   if(size(Greal,1)/=Nlat) stop "size Greal(:,*,*,*,*,*) /= Nlat"
  !   if(size(Greal,2)/=size(Greal,3)) stop "size Greal(*,:,*,*,*,*) /= Greal(*,*,:,*,*,*) "
  !   if(size(Greal,2)/=Nspin) stop "size Greal(*,:,*,*,*,*) /= Nspin"
  !   if(size(Greal,4)/=size(Greal,5)) stop "size Greal(*,*,*,:,*,*) /= Greal(*,*,*,*,:,*) "
  !   if(size(Greal,4)/=Norb) stop "size Greal(*,*,*,*,*,:) /= Norb"
  !   if(size(Greal,6)/=Lreal) stop "size Greal(*,*,*,*,*,:) /= Lreal"
  !   ! Check Hk_symm !
  !   if(present(hk_symm)) then
  !      if(size(hk_symm)/=Lk) stop "size hk_symm(:) /= Lk"
  !   end if
  !   allocate(hk_symm_(Lk))
  !   hk_symm_=.false.
  !   if(present(hk_symm)) hk_symm_=hk_symm
  !   No=Norb*Nlat
  !   Gmats=zero
  !   Greal=zero
  !   if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
  !   if(mpiID==0)call start_timer
  !   do ik=1,Lk
  !      !+- NOTE: here I pass only SPIN-diagonal quantities -+!
  !      do ispin=1,Nspin
  !         call get_gloc_mats_mb(Hk(:,:,ik),wt(ik),Hloc(:,ispin,ispin,:,:),Smats(:,ispin,ispin,:,:,:),Gmats(:,ispin,ispin,:,:,:),hk_symm_(ik))
  !         call get_gloc_real_mb(Hk(:,:,ik),wt(ik),Hloc(:,ispin,ispin,:,:),Sreal(:,ispin,ispin,:,:,:),Greal(:,ispin,ispin,:,:,:),hk_symm_(ik))
  !      end do
  !      if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !   end do
  !   if(mpiID==0)call stop_timer
  ! end subroutine rdmft_get_gloc_mb
  ! !+- NOTES: here I pass only SPIN-diagonal quantities -+!
  ! subroutine get_gloc_mats_mb(Hk,wtk,Hloc,Smats,fg,hk_symm)
  !   complex(8)               :: Hloc(Nlat,Norb,Norb)
  !   complex(8)               :: Hk(Norb*Nlat,Norb*Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Norb,Norb,Lmats) 
  !   complex(8),intent(inout) :: Smats(Nlat,Norb,Norb,Lmats)
  !   logical                  :: hk_symm
  !   real(8)                  :: ek,wtk
  !   complex(8)               :: Gloc(Nlat*Norb,Nlat*Norb)
  !   real(8)                  :: Id(Nlat*Norb,Nlat*Norb)
  !   complex(8)               :: fg_k(Nlat,Norb,Norb,Lmats)
  !   complex(8)               :: gf_tmp(Nlat,Norb,Norb,Lmats)
  !   complex(8)               :: zeta
  !   integer                  :: i,is,io,jo,ilat,iorb,jorb
  !   integer                  :: No
  !   !
  !   No = Norb*Nlat
  !   Id=0.d0
  !   do io=1,No
  !      Id(io,io) = 1.d0
  !   end do
  !   !
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lmats,mpiSIZE
  !      Gloc(:,:) = (xi*wm(i)+xmu)*Id(:,:) - Hk(:,:)
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               Gloc(io,jo) = Gloc(io,jo) - Smats(ilat,iorb,jorb,i) 
  !               Gloc(io,jo) = Gloc(io,jo) - Hloc(ilat,iorb,jorb) 
  !            end do
  !         end do
  !      end do
  !      !
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)  
  !      else
  !         call matrix_inverse(Gloc)  
  !      end if
  !      !
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb             
  !               ! dump the local GF onto gf_tmp !
  !               gf_tmp(ilat,iorb,jorb,i) = Gloc(io,jo)
  !            end do
  !         end do
  !      end do
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_mats_mb
  ! ! real freq
  ! subroutine get_gloc_real_mb(Hk,wtk,Hloc,Sreal,fg,hk_symm)
  !   complex(8)                  :: Hloc(Nlat,Norb,Norb)
  !   complex(8)               :: Hk(Norb*Nlat,Norb*Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Norb,Norb,Lreal) 
  !   complex(8),intent(inout) :: Sreal(Nlat,Norb,Norb,Lreal) 
  !   real(8)                  :: ek,wtk
  !   logical                  :: hk_symm
  !   complex(8)               :: Gloc(Nlat*Norb,Nlat*Norb)    
  !   real(8)                  :: Id(Nlat*Norb,Nlat*Norb)
  !   complex(8)               :: fg_k(Nlat,Norb,Norb,Lreal)
  !   complex(8)               :: gf_tmp(Nlat,Norb,Norb,Lreal)
  !   complex(8)               :: zeta
  !   integer                  :: i,is,io,jo,iorb,jorb,ilat
  !   integer                  :: No
  !   !
  !   No = Norb*Nlat
  !   Id=0.d0
  !   do io=1,No
  !      Id(io,io) = 1.d0
  !   end do
  !   !
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lreal,mpiSIZE
  !      Gloc(:,:) = (cmplx(wr(i),eps,8) + xmu)*Id(:,:) - Hk(:,:)
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               Gloc(io,jo) = Gloc(io,jo) - Sreal(ilat,iorb,jorb,i) 
  !               Gloc(io,jo) = Gloc(io,jo) - Hloc(ilat,iorb,jorb) 
  !            end do
  !         end do
  !      end do
  !      !
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)
  !      end if
  !      !
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb             
  !               ! dump local GF onto gf_tmp !
  !               gf_tmp(ilat,iorb,jorb,i) = Gloc(io,jo)
  !            end do
  !         end do
  !      end do
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_real_mb


  ! ! superconducting case !
  ! subroutine rdmft_get_gloc_superc(Hk,wt,Gmats,Greal,Smats,Sreal,Eloc,hk_symm)
  !   complex(8),dimension(:,:,:)            :: Hk
  !   real(8),dimension(size(Hk,3))          :: wt
  !   real(8),optional                       :: Eloc(Nlat)
  !   real(8)                                :: Eloc_(Nlat)
  !   complex(8),intent(inout)               :: Gmats(2,Nlat,Lmats)
  !   complex(8),intent(inout)               :: Greal(2,Nlat,Lreal)
  !   complex(8),intent(inout)               :: Smats(2,Nlat,Lmats)
  !   complex(8),intent(inout)               :: Sreal(2,Nlat,Lreal)
  !   logical,dimension(size(Hk,3)),optional :: hk_symm
  !   logical,dimension(size(Hk,3))          :: hk_symm_    
  !   integer                                :: ik,Lk
  !   Lk=size(Hk,3)
  !   if(size(Hk,1)/=Nlat.or.size(Hk,2)/=Nlat) stop "Wrong dimensions Hk"
  !   !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
  !   !then we use a newly created function *blocks_to_matrix* to spread the blocks into
  !   !a matrix of rank 2 dimensions Nside*Nside
  !   if(allocated(wm))deallocate(wm)
  !   if(allocated(wr))deallocate(wr)
  !   allocate(wm(Lmats))
  !   allocate(wr(Lreal))
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !   eloc_=0.d0;if(present(eloc))eloc_=eloc
  !   hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !   Gmats=zero
  !   Greal=zero
  !   if(mpiID==0)call start_timer
  !   if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
  !   do ik=1,Lk
  !      call get_sc_gloc_mats(eloc_,Smats,Gmats,Hk(:,:,ik),wt(ik),hk_symm_(ik))
  !      call get_sc_gloc_real(eloc_,Sreal,Greal,Hk(:,:,ik),wt(ik),hk_symm_(ik))
  !      if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !   end do
  !   if(mpiID==0)call stop_timer
  ! end subroutine rdmft_get_gloc_superc

  ! ! mats freq
  ! subroutine get_sc_gloc_mats(elocal,sigma,fg,Hk,wtk,hk_symm)
  !   real(8)    :: elocal(Nlat) 
  !   complex(8) :: Hk(Nlat,Nlat),Hk_dag(Nlat,Nlat)
  !   complex(8),intent(inout) :: fg(2,Nlat,Lmats),sigma(2,Nlat,Lmats)
  !   real(8)    :: ek,wtk
  !   logical    :: hk_symm
  !   complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lmats),fg_k(2,Nlat,Lmats)
  !   integer    :: i,j,is,js
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lmats,mpiSIZE
  !      Gloc=zero
  !      Gloc(1:Nlat,1:Nlat)              =  -Hk
  !      Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat)=   Hk !!! REMEMBER CONJG(H(-k)) = H(K) for a bipartite lattice       
  !      do is=1,Nlat
  !         Gloc(is,is)           =  Gloc(is,is) + xi*wm(i) - sigma(1,is,i)        - elocal(is) + xmu
  !         Gloc(Nlat+is,Nlat+is) =  Gloc(Nlat+is,Nlat+is) + xi*wm(i) + conjg(sigma(1,is,i)) + elocal(is) - xmu !==-conjg(Gloc(is,is))
  !         Gloc(is,Nlat+is)      =  Gloc(is,Nlat+is) - sigma(2,is,i)
  !         Gloc(Nlat+is,is)      =  Gloc(Nlat+is,is) - sigma(2,is,i)                                         
  !      enddo
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)
  !      end if
  !      forall(is=1:Nlat)
  !         gf_tmp(1,is,i) = Gloc(is,is)
  !         !#ACHTUNG
  !         gf_tmp(2,is,i) = dreal(Gloc(is,Nlat+is))
  !      end forall
  !   enddo
  !   ! reduce gf_tmp on fg_k 
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   fg= fg + fg_k*wtk
  ! end subroutine get_sc_gloc_mats

  ! ! real freq
  ! subroutine get_sc_gloc_real(elocal,sigma,fg,Hk,wtk,hk_symm)
  !   real(8)    :: elocal(Nlat)
  !   complex(8) :: Hk(Nlat,Nlat)
  !   complex(8) :: fg(2,Nlat,Lreal),sigma(2,Nlat,Lreal)
  !   real(8)    :: ek,wtk
  !   logical    :: hk_symm
  !   complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lreal),zeta1,zeta2,fg_k(2,Nlat,Lreal)
  !   integer    :: i,is,js
  !   fg_k=zero; 
  !   gf_tmp=zero
  !   do i=1+mpiID,Lreal,mpiSIZE
  !      Gloc=zero
  !      Gloc(1:Nlat,1:Nlat)               = - Hk
  !      Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat) =   Hk  !!REMEMEBER CONJG(H(-k))=H(k) for a bipartite lattice
  !      do is=1,Nlat
  !         ! diagonal parts of Gk^(-1)
  !         zeta1=        cmplx(wr(i),eps,8) + xmu - sigma(1,is,i)       - elocal(is)
  !         zeta2= -conjg(cmplx(wr(Lreal+1-i),eps,8) + xmu - sigma(1,is,Lreal+1-i)) + elocal(is)
  !         !
  !         Gloc(is,is)           =  Gloc(is,is) + zeta1
  !         Gloc(Nlat+is,Nlat+is) =  Gloc(Nlat+is,Nlat+is) + zeta2
  !         Gloc(is,Nlat+is)      =  Gloc(is,Nlat+is) - sigma(2,is,i)
  !         Gloc(Nlat+is,is)      =  Gloc(Nlat+is,is) - conjg(sigma(2,is,Lreal+1-i)) 
  !      enddo
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)
  !      end if
  !      forall(is=1:Nlat)
  !         gf_tmp(1,is,i) = Gloc(is,is)
  !         gf_tmp(2,is,i) = Gloc(is,Nlat+is)
  !      end forall
  !   enddo
  !   ! reduce gf_tmp
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   ! k sum
  !   fg = fg + fg_k*wtk
  ! end subroutine get_sc_gloc_real



  ! !+- EMBEDD LAYERS -+!
  ! subroutine rdmft_get_gloc_embedd(Hk,wt,Hloc,Gmats,Greal,Smats,Sreal,epsik,hk_symm)
  !   complex(8)                       :: Hk(:,:,:)          ![Norb*Nlat][Norb*Nlat][Lk]
  !   real(8)                          :: wt(:)              ![Lk]
  !   complex(8)                       :: Hloc(:,:,:,:,:)    ![Nlat][Nspin][Nspin][Norb][Norb]
  !   complex(8),intent(inout)         :: Gmats(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !   complex(8),intent(inout)         :: Greal(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  !   complex(8),intent(inout)         :: Smats(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !   complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:) ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  !   logical,dimension(:),optional    :: hk_symm
  !   real(8),dimension(:)             :: epsik
  !   logical,dimension(:),allocatable :: hk_symm_
  !   integer                          :: ik,Lk,ispin
  !   integer                          :: No
  !   if(allocated(wm))deallocate(wm)
  !   if(allocated(wr))deallocate(wr)
  !   allocate(wm(Lmats))
  !   allocate(wr(Lreal))
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !   !checks!
  !   No = Norb*Nlat
  !   Lk=size(Hk,3)
  !   ! check Hk !
  !   if(size(Hk,1)/=No) stop "size Hk(:,*,*) /= No"
  !   if(size(Hk,2)/=size(Hk,1)) stop "size Hk(*,:,*) /= Hk(:,*,*) "
  !   ! check wt !
  !   if(size(wt)/=Lk) stop "size wt(*) /= Lk"
  !   ! check Hloc !
  !   if(size(Hloc,1)/=Nlat) stop "size Hloc(:,*,*,*,*) /= Nlat"
  !   if(size(Hloc,2)/=size(Hloc,3)) stop "size Hloc(*,:,*,*,*) /= Hloc(*,*,:,*,*) "
  !   if(size(Hloc,2)/=Nspin) stop "size Hloc(*,:,*,*,*) /= Nspin"
  !   if(size(Hloc,4)/=size(Hloc,5)) stop "size Hloc(*,*,*,:,*) /= Hloc(*,*,*,*,:) "
  !   if(size(Hloc,4)/=Norb) stop "size Hloc(*,*,*,:,*) /= Norb"
  !   ! check Smats !
  !   if(size(Smats,1)/=Nlat) stop "size Smats(:,*,*,*,*,*) /= Nlat"
  !   if(size(Smats,2)/=size(Smats,3)) stop "size Smats(*,:,*,*,*,*) /= Smats(*,*,:,*,*,*) "
  !   if(size(Smats,2)/=Nspin) stop "size Smats(*,:,*,*,*,*) /= Nspin"
  !   if(size(Smats,4)/=size(Smats,5)) stop "size Smats(*,*,*,:,*,*) /= Smats(*,*,*,*,:,*) "
  !   if(size(Smats,4)/=Norb) stop "size Smats(*,*,*,*,*,:) /= Norb"
  !   if(size(Smats,6)/=Lmats) stop "size Smats(*,*,*,*,*,:) /= Lmats"
  !   ! check Sreal !
  !   if(size(Sreal,1)/=Nlat) stop "size Sreal(:,*,*,*,*,*) /= Nlat"
  !   if(size(Sreal,2)/=size(Sreal,3)) stop "size Sreal(*,:,*,*,*,*) /= Sreal(*,*,:,*,*,*) "
  !   if(size(Sreal,2)/=Nspin) stop "size Sreal(*,:,*,*,*,*) /= Nspin"
  !   if(size(Sreal,4)/=size(Sreal,5)) stop "size Sreal(*,*,*,:,*,*) /= Sreal(*,*,*,*,:,*) "
  !   if(size(Sreal,4)/=Norb) stop "size Sreal(*,*,*,*,*,:) /= Norb"
  !   if(size(Sreal,6)/=Lmats) stop "size Sreal(*,*,*,*,*,:) /= Lmats"
  !   ! check Gmats !
  !   if(size(Gmats,1)/=Nlat) stop "size Gmats(:,*,*,*,*,*) /= Nlat"
  !   if(size(Gmats,2)/=size(Gmats,3)) stop "size Gmats(*,:,*,*,*,*) /= Gmats(*,*,:,*,*,*) "
  !   if(size(Gmats,2)/=Nspin) stop "size Gmats(*,:,*,*,*,*) /= Nspin"
  !   if(size(Gmats,4)/=size(Gmats,5)) stop "size Gmats(*,*,*,:,*,*) /= Gmats(*,*,*,*,:,*) "
  !   if(size(Gmats,4)/=Norb) stop "size Gmats(*,*,*,*,*,:) /= Norb"
  !   if(size(Gmats,6)/=Lmats) stop "size Gmats(*,*,*,*,*,:) /= Lmats"
  !   ! check Greal !
  !   if(size(Greal,1)/=Nlat) stop "size Greal(:,*,*,*,*,*) /= Nlat"
  !   if(size(Greal,2)/=size(Greal,3)) stop "size Greal(*,:,*,*,*,*) /= Greal(*,*,:,*,*,*) "
  !   if(size(Greal,2)/=Nspin) stop "size Greal(*,:,*,*,*,*) /= Nspin"
  !   if(size(Greal,4)/=size(Greal,5)) stop "size Greal(*,*,*,:,*,*) /= Greal(*,*,*,*,:,*) "
  !   if(size(Greal,4)/=Norb) stop "size Greal(*,*,*,*,*,:) /= Norb"
  !   if(size(Greal,6)/=Lreal) stop "size Greal(*,*,*,*,*,:) /= Lreal"
  !   ! Check Hk_symm !
  !   if(present(hk_symm)) then
  !      if(size(hk_symm)/=Lk) stop "size hk_symm(:) /= Lk"
  !   end if
  !   allocate(hk_symm_(Lk))
  !   hk_symm_=.false.
  !   if(present(hk_symm)) hk_symm_=hk_symm

  !   No=Norb*Nlat
  !   Gmats=zero
  !   Greal=zero
  !   if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
  !   if(mpiID==0)call start_timer
  !   do ik=1,Lk
  !      !+- NOTE: here I pass only SPIN-diagonal quantities -+!
  !      do ispin=1,Nspin
  !         call get_gloc_mats_embedd(Hk(:,:,ik),wt(ik),Hloc(:,ispin,ispin,:,:),Smats(:,ispin,ispin,:,:,:),Gmats(:,ispin,ispin,:,:,:),hk_symm_(ik),epsik)
  !         call get_gloc_real_embedd(Hk(:,:,ik),wt(ik),Hloc(:,ispin,ispin,:,:),Sreal(:,ispin,ispin,:,:,:),Greal(:,ispin,ispin,:,:,:),hk_symm_(ik),epsik)
  !      end do
  !      if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
  !   end do
  !   if(mpiID==0)call stop_timer
  ! end subroutine rdmft_get_gloc_embedd
  ! !+- NOTES: here I pass only SPIN-diagonal quantities -+!
  ! subroutine get_gloc_mats_embedd(Hk,wtk,Hloc,Smats,fg,hk_symm,epsik)
  !   complex(8)               :: Hloc(Nlat,Norb,Norb)
  !   complex(8)               :: Hk(Norb*Nlat,Norb*Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Norb,Norb,Lmats) 
  !   complex(8),intent(inout) :: Smats(Nlat,Norb,Norb,Lmats)
  !   real(8),dimension(:)     :: epsik
  !   complex(8)               :: sL(Norb,Norb,Lmats),sR(Norb,Norb,Lmats)
  !   complex(8)               :: gL(Norb,Norb,Lmats),gR(Norb,Norb,Lmats)
  !   complex(8)               :: tL(Norb,Norb),tR(Norb,Norb)
  !   complex(8)               :: tmpG(Norb,Norb)
  !   logical                  :: hk_symm
  !   real(8)                  :: ek,wtk
  !   complex(8)               :: Gloc(Nlat*Norb,Nlat*Norb)
  !   real(8)                  :: Id(Nlat*Norb,Nlat*Norb)
  !   complex(8)               :: fg_k(Nlat,Norb,Norb,Lmats)
  !   complex(8)               :: gf_tmp(Nlat,Norb,Norb,Lmats)
  !   complex(8)               :: zeta
  !   integer                  :: i,is,io,jo,ilat,iorb,jorb
  !   integer                  :: No,ikz,Nkz
  !   No = Norb*Nlat
  !   Id=0.d0
  !   do io=1,No
  !      Id(io,io) = 1.d0
  !   end do
  !   !
  !   !+- compute embedded potential -+!
  !   !
  !   Nkz=size(epsik)    
  !   tL = zero; tR=zero
  !   do iorb=1,Norb
  !      tL(iorb,iorb) = -ts
  !      tR(iorb,iorb) = -ts
  !   end do
  !   gL=zero
  !   gR=zero
  !   do i=1,Lmats
  !      do ikz=1,Nkz
  !         !left
  !         tmpG=zero
  !         do iorb=1,Norb          
  !            tmpG(iorb,iorb) = xi*wm(i) + xmu
  !         end do
  !         ilat=1
  !         tmpG = tmpG - Hloc(ilat,:,:) - epsik(ikz) !- Smats(1,:,:,i)
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               tmpG(iorb,jorb) = tmpG(iorb,jorb) - Hk(io,jo)
  !            end do
  !         end do
  !         call matrix_inverse_sym(tmpG)
  !         gL(:,:,i) = gL(:,:,i) + tmpG/dble(Nkz)
  !         !right
  !         tmpG=zero
  !         do iorb=1,Norb          
  !            tmpG(iorb,iorb) = xi*wm(i) + xmu
  !         end do
  !         ilat=Nlat
  !         tmpG = tmpG - Hloc(ilat,:,:) - epsik(ikz) !- Smats(Nlat,:,:,i)
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               tmpG(iorb,jorb) = tmpG(iorb,jorb) - Hk(io,jo)
  !            end do
  !         end do
  !         call matrix_inverse_sym(tmpG)
  !         gR(:,:,i) = gR(:,:,i) + tmpG/dble(Nkz)
  !      end do
  !      sL(:,:,i) = matmul(tL,gL(:,:,i)); sL(:,:,i) = matmul(sL(:,:,i),tL)
  !      sR(:,:,i) = matmul(tR,gR(:,:,i)); sR(:,:,i) = matmul(sR(:,:,i),tR)
  !   end do
  !   !
  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lmats,mpiSIZE
  !      Gloc(:,:) = (xi*wm(i)+xmu)*Id(:,:) - Hk(:,:)
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               Gloc(io,jo) = Gloc(io,jo) - Smats(ilat,iorb,jorb,i) 
  !               Gloc(io,jo) = Gloc(io,jo) - Hloc(ilat,iorb,jorb) 
  !               if(ilat.eq.1) Gloc(io,jo) = Gloc(io,jo) - sL(iorb,jorb,i)
  !               if(ilat.eq.Nlat) Gloc(io,jo) = Gloc(io,jo) - sR(iorb,jorb,i)
  !            end do
  !         end do
  !      end do
  !      !
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)  
  !      else
  !         call matrix_inverse(Gloc)  
  !      end if
  !      !
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb             
  !               ! dump the local GF onto gf_tmp !
  !               gf_tmp(ilat,iorb,jorb,i) = Gloc(io,jo)
  !            end do
  !         end do
  !      end do
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_mats_embedd
  ! ! real freq
  ! subroutine get_gloc_real_embedd(Hk,wtk,Hloc,Sreal,fg,hk_symm,epsik)
  !   complex(8)                  :: Hloc(Nlat,Norb,Norb)
  !   complex(8)               :: Hk(Norb*Nlat,Norb*Nlat)
  !   complex(8),intent(inout) :: fg(Nlat,Norb,Norb,Lreal) 
  !   complex(8),intent(inout) :: Sreal(Nlat,Norb,Norb,Lreal) 
  !   real(8)                  :: ek,wtk
  !   logical                  :: hk_symm
  !   real(8),dimension(:)     :: epsik
  !   complex(8)               :: sL(Norb,Norb,Lreal),sR(Norb,Norb,Lreal)
  !   complex(8)               :: gL(Norb,Norb,Lreal),gR(Norb,Norb,Lreal)
  !   complex(8)               :: tL(Norb,Norb),tR(Norb,Norb)
  !   complex(8)               :: tmpG(Norb,Norb)
  !   complex(8)               :: Gloc(Nlat*Norb,Nlat*Norb)    
  !   real(8)                  :: Id(Nlat*Norb,Nlat*Norb)
  !   complex(8)               :: fg_k(Nlat,Norb,Norb,Lreal)
  !   complex(8)               :: gf_tmp(Nlat,Norb,Norb,Lreal)
  !   complex(8)               :: zeta
  !   integer                  :: i,is,io,jo,iorb,jorb,ilat
  !   integer                  :: No,Nkz,ikz
  !   !
  !   No = Norb*Nlat
  !   Id=0.d0
  !   do io=1,No
  !      Id(io,io) = 1.d0
  !   end do
  !   !
  !   Nkz=size(epsik)    
  !   tL = zero; tR=zero
  !   do iorb=1,Norb
  !      tL(iorb,iorb) = -ts
  !      tR(iorb,iorb) = -ts
  !   end do
  !   gL=zero
  !   gR=zero
  !   do i=1,Lreal
  !      do ikz=1,Nkz
  !         !left
  !         tmpG=zero
  !         do iorb=1,Norb          
  !            tmpG(iorb,iorb) = cmplx(wr(i),eps,8) + xmu
  !         end do
  !         tmpG = tmpG - Hloc(1,:,:) - epsik(ikz) - Sreal(1,:,:,i)
  !         call matrix_inverse_sym(tmpG)
  !         gL(:,:,i) = gL(:,:,i) + tmpG/dble(Nkz)
  !         !right
  !         tmpG=zero
  !         do iorb=1,Norb          
  !            tmpG(iorb,iorb) = xi*wm(i) + xmu
  !         end do
  !         tmpG = tmpG - Hloc(Nlat,:,:) - epsik(ikz) - Sreal(Nlat,:,:,i)
  !         call matrix_inverse_sym(tmpG)
  !         gR(:,:,i) = gR(:,:,i) + tmpG/dble(Nkz)
  !      end do
  !      sL(:,:,i) = matmul(tL,gL(:,:,i)); sL(:,:,i) = matmul(sL(:,:,i),tL)
  !      sR(:,:,i) = matmul(tR,gR(:,:,i)); sR(:,:,i) = matmul(sR(:,:,i),tR)
  !   end do

  !   fg_k=zero
  !   gf_tmp=zero
  !   do i=1+mpiID,Lreal,mpiSIZE
  !      Gloc(:,:) = (cmplx(wr(i),eps,8) + xmu)*Id(:,:) - Hk(:,:)
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb
  !               Gloc(io,jo) = Gloc(io,jo) - Sreal(ilat,iorb,jorb,i) 
  !               Gloc(io,jo) = Gloc(io,jo) - Hloc(ilat,iorb,jorb) 
  !               sL(:,:,i) = matmul(tL,gL(:,:,i)); sL(:,:,i) = matmul(sL(:,:,i),tL)
  !               sR(:,:,i) = matmul(tR,gR(:,:,i)); sR(:,:,i) = matmul(sR(:,:,i),tR)
  !            end do
  !         end do
  !      end do
  !      !
  !      if(hk_symm) then
  !         call matrix_inverse_sym(Gloc)
  !      else
  !         call matrix_inverse(Gloc)
  !      end if
  !      !
  !      do ilat=1,Nlat
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               io = (ilat-1)*Norb + iorb
  !               jo = (ilat-1)*Norb + jorb             
  !               ! dump local GF onto gf_tmp !
  !               gf_tmp(ilat,iorb,jorb,i) = Gloc(io,jo)
  !            end do
  !         end do
  !      end do
  !   enddo
  !   !+- reduce gf_tmp on fg_k -+!
  !   call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   !+- k-sums -+!
  !   fg = fg + fg_k*wtk
  ! end subroutine get_gloc_real_embedd


end module ED_WRAP_GLOC
