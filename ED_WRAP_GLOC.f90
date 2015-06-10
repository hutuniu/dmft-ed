module ED_WRAP_GLOC
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_WRAP_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS,   only: reg,txtfy,store_data
  USE SF_LINALG,    only:inv,inv_sym,matrix_inverse,matrix_inverse_sym,matrix_diagonalize
  USE SF_ARRAYS,    only:linspace,arange
  implicit none
  private
  !
  interface ed_get_gloc_lattice
     module procedure ed_get_gloc_normal_main
     module procedure ed_get_gloc_superc_main
     module procedure ed_get_gloc_normal_1b
     module procedure ed_get_gloc_normal_mb
     module procedure ed_get_gloc_superc_1b
     module procedure ed_get_gloc_superc_mb
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
  subroutine ed_get_gloc_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm,Gamma_mats,Gamma_real)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                     :: iprint
    complex(8)                  :: zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    complex(8),optional         :: Gamma_mats(size(Hk,1),size(Hk,2),Lmats)
    complex(8),optional         :: Gamma_real(size(Hk,1),size(Hk,2),Lreal)
    integer                     :: i,ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    Lk=size(Hk,3)
    Nlso=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    !
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
    !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
    !then we use a newly created function *blocks_to_matrix* to spread the blocks into
    !a matrix of rank 2 dimensions Nlso*Nlso
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
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpiID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       if(present(Gamma_mats))then
          call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Gembed=Gamma_mats)
       else
          call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       endif
       if(present(Gamma_real))then
          call add_to_gloc_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Gembed=Gamma_real)
       else
          call add_to_gloc_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       endif
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       select case(iprint)
       case (0)
          write(LOGfile,*)"Gloc not written on file."
       case(1)                  !print only diagonal elements
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data("LG"//reg(suffix),Gmats(:,ispin,ispin,iorb,iorb,:),wm)
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data("LG"//reg(suffix),Greal(:,ispin,ispin,iorb,iorb,:),wr)
             enddo
          enddo
       case(2)                  !print spin-diagonal, all orbitals 
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("LG"//reg(suffix),Gmats(:,ispin,ispin,iorb,jorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("LG"//reg(suffix),Greal(:,ispin,ispin,iorb,jorb,:),wr)
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
                      call store_data("LG"//reg(suffix),Gmats(:,ispin,jspin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                      call store_data("LG"//reg(suffix),Greal(:,ispin,jspin,iorb,jorb,:),wr)
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_gloc_normal_main

  subroutine add_to_gloc_normal(zeta_site,Hk,hk_symm,Gkout,Gembed)
    complex(8)               :: zeta_site(:,:,:,:)              ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                  :: hk_symm
    !output:
    complex(8),intent(inout) :: Gkout(Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
    complex(8)               :: Gktmp(Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
    !
    complex(8)               :: Gmatrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex(8),optional      :: Gembed(Nlat*Nspin*Norb,Nlat*Nspin*Norb,size(zeta_site,4)) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    if(size(zeta_site,1)/=Nlat)stop "get_gloc_kpoint error: zeta_site wrong size 1 = Nlat"
    if(size(zeta_site,2)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 2 = Nspin*Norb"
    if(size(zeta_site,3)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 3 = Nspin*Norb"
    Lfreq = size(zeta_site,4)
    Gktmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
       Gmatrix  = blocks_to_matrix(zeta_site(:,:,:,i)) - Hk
       if(present(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call matrix_inverse_sym(Gmatrix)
       else
          call matrix_inverse(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
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
  end subroutine add_to_gloc_normal











  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the Nambu local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k;R_i,R_j)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gloc_superc_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
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
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu - Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu + Eloc_(js)
             !
             !SYMMETRIES in real-frequencies   [assuming a real order parameter]
             !G22(w)  = -[G11[-w]]*
             !G21(w)  =   G12[w]             
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu - Eloc_(js)
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu - Eloc_(js) )
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
       call add_to_gloc_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       select case(iprint)
       case (0)
          write(LOGfile,*)"Gloc not written on file."
       case(1)                  !print only diagonal elements
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data("LG"//reg(suffix),Gmats(1,:,ispin,ispin,iorb,iorb,:),wm)
                call store_data("LF"//reg(suffix),Gmats(2,:,ispin,ispin,iorb,iorb,:),wm)
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data("LG"//reg(suffix),Greal(1,:,ispin,ispin,iorb,iorb,:),wr)
                call store_data("LF"//reg(suffix),Greal(2,:,ispin,ispin,iorb,iorb,:),wr)
             enddo
          enddo
       case(2)                  !print spin-diagonal, all orbitals 
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("LG"//reg(suffix),Gmats(1,:,ispin,ispin,iorb,jorb,:),wm)
                   call store_data("LF"//reg(suffix),Gmats(2,:,ispin,ispin,iorb,jorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("LG"//reg(suffix),Greal(1,:,ispin,ispin,iorb,jorb,:),wr)
                   call store_data("LF"//reg(suffix),Greal(2,:,ispin,ispin,iorb,jorb,:),wr)
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
                      call store_data("LG"//reg(suffix),Gmats(1,:,ispin,jspin,iorb,jorb,:),wm)
                      call store_data("LF"//reg(suffix),Gmats(2,:,ispin,jspin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                      call store_data("LG"//reg(suffix),Greal(1,:,ispin,jspin,iorb,jorb,:),wr)
                      call store_data("LF"//reg(suffix),Greal(2,:,ispin,jspin,iorb,jorb,:),wr)
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_gloc_superc_main

  subroutine add_to_gloc_superc(zeta_site,Hk,hk_symm,Gkout)
    complex(8)               :: zeta_site(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gkout(2,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,6))
    complex(8)               :: Gktmp(2,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,6))
    !
    complex(8)               :: Gmatrix(2*Nlat*Nspin*Norb , 2*Nlat*Nspin*Norb)
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nlso
    if(size(zeta_site,1)/=2.OR.size(zeta_site,2)/=2)stop "add_to_gloc_superc error: zeta_site wrong size 1 or 2 = 2 (Nambu)"
    if(size(zeta_site,3)/=Nlat)stop "add_to_gloc_superc error: zeta_site wrong size 2 = Nlat"
    if(size(zeta_site,4)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 4 = Nspin*Norb"
    if(size(zeta_site,5)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 5 = Nspin*Norb"
    Lfreq = size(zeta_site,6)
    Nlso  = Nlat*Nspin*Norb
    Gkout = zero
    Gktmp = zero
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
  end subroutine add_to_gloc_superc










  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate all components of the Normal Green's function for a given 
  ! Hamiltonian matrix and self-energy functions. Hk is a big sparse matrix of the form 
  ! H(\k;R_i,R_j)_{ab}^{ss'} and size [Nk]*[Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine ed_get_gij_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm,Gamma_mats,Gamma_real)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                     :: iprint
    complex(8)                  :: zeta_mats(Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    complex(8),optional         :: Gamma_mats(size(Hk,1),size(Hk,2),Lmats)
    complex(8),optional         :: Gamma_real(size(Hk,1),size(Hk,2),Lreal)
    integer                     :: i,ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    Lk=size(Hk,3)
    Nlso=Nlat*Norb*Nspin
    if(size(Hk,1)/=Nlso.OR.size(Hk,2)/=Nlso) stop "rdmft_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    !
    if(mpiID==0)write(LOGfile,*)"Get full GF (id=0):"
    !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
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
       select case(iprint)
       case (0)
          write(LOGfile,*)"Gloc not written on file."
       case(1)                  !print only diagonal elements
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data("Gij"//reg(suffix),Gmats(:,:,ispin,ispin,iorb,iorb,:),wm)
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data("Gij"//reg(suffix),Greal(:,:,ispin,ispin,iorb,iorb,:),wr)
             enddo
          enddo
       case(2)                  !print spin-diagonal, all orbitals 
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("Gij"//reg(suffix),Gmats(:,:,ispin,ispin,iorb,jorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("Gij"//reg(suffix),Greal(:,:,ispin,ispin,iorb,jorb,:),wr)
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
                      call store_data("Gij"//reg(suffix),Gmats(:,:,ispin,jspin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                      call store_data("Gij"//reg(suffix),Greal(:,:,ispin,jspin,iorb,jorb,:),wr)
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_gij_normal_main

  subroutine add_to_gij_normal(zeta_site,Hk,hk_symm,Gkout,Gembed)
    complex(8)               :: zeta_site(:,:,:,:)              ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                  :: hk_symm
    !output:
    complex(8),intent(inout) :: Gkout(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
    complex(8)               :: Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,4))
    !
    complex(8)               :: Gmatrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex(8),optional      :: Gembed(Nlat*Nspin*Norb,Nlat*Nspin*Norb,size(zeta_site,4)) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    if(size(zeta_site,1)/=Nlat)stop "add_to_gij_normal error: zeta_site wrong size 1 = Nlat"
    if(size(zeta_site,2)/=Nspin*Norb)stop "add_to_gij_normal error: zeta_site wrong size 2 = Nspin*Norb"
    if(size(zeta_site,3)/=Nspin*Norb)stop "add_to_gij_normal error: zeta_site wrong size 3 = Nspin*Norb"
    Lfreq = size(zeta_site,4)
    Gktmp=zero
    do i=1+mpiID,Lfreq,mpiSIZE
       Gmatrix  = blocks_to_matrix(zeta_site(:,:,:,i)) - Hk
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
  subroutine ed_get_gij_superc_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: zeta_mats(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(2,2,Nlat,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
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
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu - Eloc_(js)
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu + Eloc_(js)
             !
             !SYMMETRIES in real-frequencies   [assuming a real order parameter]
             !G22(w)  = -[G11[-w]]*
             !G21(w)  =   G12[w]             
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu - Eloc_(js)
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu - Eloc_(js) )
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
       call add_to_gloc_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(mpiID==0)call stop_timer
    if(mpiID==0)then
       select case(iprint)
       case (0)
          write(LOGfile,*)"Gloc not written on file."
       case(1)                  !print only diagonal elements
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data("Gij"//reg(suffix),Gmats(1,:,:,ispin,ispin,iorb,iorb,:),wm)
                call store_data("Fij"//reg(suffix),Gmats(2,:,:,ispin,ispin,iorb,iorb,:),wm)
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                call store_data("Gij"//reg(suffix),Greal(1,:,:,ispin,ispin,iorb,iorb,:),wr)
                call store_data("Fij"//reg(suffix),Greal(2,:,:,ispin,ispin,iorb,iorb,:),wr)
             enddo
          enddo
       case(2)                  !print spin-diagonal, all orbitals 
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   call store_data("Gij"//reg(suffix),Gmats(1,:,:,ispin,ispin,iorb,jorb,:),wm)
                   call store_data("Fij"//reg(suffix),Gmats(2,:,:,ispin,ispin,iorb,jorb,:),wm)
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
                   call store_data("Gij"//reg(suffix),Greal(1,:,:,ispin,ispin,iorb,jorb,:),wr)
                   call store_data("Fij"//reg(suffix),Greal(2,:,:,ispin,ispin,iorb,jorb,:),wr)
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
                      call store_data("Gij"//reg(suffix),Gmats(1,:,:,ispin,jspin,iorb,jorb,:),wm)
                      call store_data("Fij"//reg(suffix),Gmats(2,:,:,ispin,jspin,iorb,jorb,:),wm)
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                      call store_data("Gij"//reg(suffix),Greal(1,:,:,ispin,jspin,iorb,jorb,:),wr)
                      call store_data("Fij"//reg(suffix),Greal(2,:,:,ispin,jspin,iorb,jorb,:),wr)
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_gij_superc_main

  subroutine add_to_gij_superc(zeta_site,Hk,hk_symm,Gkout)
    complex(8)               :: zeta_site(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)               :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gkout(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,6))
    complex(8)               :: Gktmp(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta_site,6))
    !
    complex(8)               :: Gmatrix(2*Nlat*Nspin*Norb , 2*Nlat*Nspin*Norb)
    integer                  :: i,is,Lfreq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nlso
    if(size(zeta_site,1)/=2.OR.size(zeta_site,2)/=2)stop "add_to_gloc_superc error: zeta_site wrong size 1 or 2 = 2 (Nambu)"
    if(size(zeta_site,3)/=Nlat)stop "add_to_gloc_superc error: zeta_site wrong size 2 = Nlat"
    if(size(zeta_site,4)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 4 = Nspin*Norb"
    if(size(zeta_site,5)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 5 = Nspin*Norb"
    Lfreq = size(zeta_site,6)
    Nlso  = Nlat*Nspin*Norb
    Gkout = zero
    Gktmp = zero
    do i=1+mpiID,Lfreq,mpiSIZE
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta_site(1,1,:,:,:,i)) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_site(1,2,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta_site(2,1,:,:,:,i))
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta_site(2,2,:,:,:,i)) + Hk
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
                         Gktmp(1,ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                         Gktmp(2,ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
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
#else
    Gkout = Gktmp
#endif
  end subroutine add_to_gij_superc








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: additional interfaces to the main procedures above to get Gloc 
  ! _1b : one band case (independent of Nspin and Norb)
  ! _mb : multi-bands case (independent of Nspin)
  !+-----------------------------------------------------------------------------+!
  subroutine ed_get_gloc_normal_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nlat,Lmats)
    complex(8),intent(inout)    :: Greal(Nlat,Lreal)
    complex(8),intent(inout)    :: Smats(Nlat,Lmats)
    complex(8),intent(inout)    :: Sreal(Nlat,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: Gmats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    if(Norb>1)stop "ed_get_gloc_normal_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_gloc_normal_1b error: Nspin > 1 in 1-band routine" 
    Gmats_(:,1,1,1,1,:) = Gmats(:,:)
    Greal_(:,1,1,1,1,:) = Greal(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Sreal(:,:) = Sreal_(:,1,1,1,1,:)
  end subroutine ed_get_gloc_normal_1b

  subroutine ed_get_gloc_normal_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nlat,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nlat,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nlat,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nlat,Norb,Norb,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: Gmats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    if(Nspin>1)stop "ed_get_gloc_normal_mb error: Nspin > 1 in M-band routine" 
    Gmats_(:,1,1,:,:,:) = Gmats(:,:,:,:)
    Greal_(:,1,1,:,:,:) = Greal(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Sreal(:,:,:,:) = Sreal_(:,1,1,:,:,:)
  end subroutine ed_get_gloc_normal_mb

  subroutine ed_get_gloc_superc_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nlat,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nlat,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nlat,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nlat,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: Gmats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    if(Norb>1)stop "ed_get_gloc_superc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_gloc_superc_1b error: Nspin > 1 in 1-band routine" 
    Gmats_(:,:,1,1,1,1,:) = Gmats(:,:,:)
    Greal_(:,:,1,1,1,1,:) = Greal(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    Sreal_(:,:,1,1,1,1,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_gloc_superc_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:,:) = Gmats_(:,:,1,1,1,1,:)
    Greal(:,:,:) = Greal_(:,:,1,1,1,1,:)
    Smats(:,:,:) = Smats_(:,:,1,1,1,1,:)
    Sreal(:,:,:) = Sreal_(:,:,1,1,1,1,:)
  end subroutine ed_get_gloc_superc_1b

  subroutine ed_get_gloc_superc_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nlat,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nlat,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nlat,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nlat,Norb,Norb,Lreal)
    integer                     :: iprint
    !
    complex(8)                  :: Gmats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Nlat*Norb*Nspin)
    real(8)                     :: Eloc_(Nlat*Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    if(Nspin>1)stop "ed_get_gloc_superc_mb error: Nspin > 1 in M-band routine" 
    Gmats_(:,:,1,1,:,:,:) = Gmats(:,:,:,:,:)
    Greal_(:,:,1,1,:,:,:) = Greal(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    Sreal_(:,:,1,1,:,:,:) = Sreal(:,:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_gloc_superc_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:,:,:,:) = Gmats_(:,:,1,1,:,:,:)
    Greal(:,:,:,:,:) = Greal_(:,:,1,1,:,:,:)
    Smats(:,:,:,:,:) = Smats_(:,:,1,1,:,:,:)
    Sreal(:,:,:,:,:) = Sreal_(:,:,1,1,:,:,:)
  end subroutine ed_get_gloc_superc_mb

end module ED_WRAP_GLOC
