!###############################################################
! purpose  : Compute Local GFunction 
!###############################################################
module ED_GLOC
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_SPECIAL,   only: gfbethe,gfbether
  USE SF_IOTOOLS,   only: reg,txtfy,splot
  USE SF_LINALG,    only: matrix_inverse,matrix_inverse_sym
  USE SF_ARRAYS,    only: linspace,arange
  implicit none
  private
  !
  interface ed_get_gloc
     module procedure ed_get_gloc_normal_main
     module procedure ed_get_gloc_normal_1b
     module procedure ed_get_gloc_normal_mb
     !
     module procedure ed_get_gloc_superc_main
     module procedure ed_get_gloc_superc_1b
     module procedure ed_get_gloc_superc_mb
     !
     module procedure ed_get_gloc_normal_bethe
     module procedure ed_get_gloc_normal_bethe_1b
     module procedure ed_get_gloc_normal_bethe_mb
  end interface ed_get_gloc

  public :: ed_get_gloc

  real(8),dimension(:),allocatable :: wr,wm
  character(len=64)                :: suffix


contains


  !----------------------------------------------------------------------------------------!
  ! PURPOSE: evaluate the local Green's function for a given Hamiltonian matrix and
  ! self-energy functions. Hk is a big sparse matrix of the form H(k)_{ab}^{ss'}
  ! and size [Nk]*[Nlat*Nspin*Norb]**2
  !   NORMAL:_normal_main
  !   SUPERC:_superc_main
  !----------------------------------------------------------------------------------------!
  !NORMAL
  subroutine ed_get_gloc_normal_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: zeta_mats(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(Nspin,Nspin,Norb,Norb,Lreal)
    real(8),optional            :: Eloc(Norb*Nspin)
    real(8)                     :: Eloc_(Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: iprint
    integer                     :: i,ik,Lk,Nso,iorb,jorb,ispin,jspin,io,jo,js
    Lk=size(Hk,3)
    Nso=Norb*Nspin
    if(size(Hk,1)/=Nso.OR.size(Hk,2)/=Nso) stop "ed_get_gloc_normal error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    if(ED_MPI_ID==0)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_ID==0)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    !here we create the "array" *zeta_site* of Nlat blocks, each of size (Nspin*Norb)
    !then we use a newly created function *blocks_to_matrix* to spread the blocks into
    !a matrix of rank 2 dimensions Nso*Nso
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    !
    zeta_mats=zero
    zeta_real=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_mats(io,io,:) = xi*wm(:)       + xmu - Eloc_(io)
          zeta_real(io,io,:) = wr(:) + xi*eps + xmu - Eloc_(io)
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_mats(io,jo,:) = zeta_mats(io,jo,:) - Smats(ispin,jspin,iorb,jorb,:)
                zeta_real(io,jo,:) = zeta_real(io,jo,:) - Sreal(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(ED_MPI_ID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
       call add_to_gloc_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_ID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_ID==0)call stop_timer
    !
    !
    !PRINT:
    call print_gloc_normal(Gmats,Greal,iprint)
    !
  end subroutine ed_get_gloc_normal_main

  !SUPERCONDUCTING:
  subroutine ed_get_gloc_superc_main(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: zeta_mats(2,2,Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(2,2,Nspin*Norb,Nspin*Norb,Lreal)
    complex(8)                  :: Gkmats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal)
    real(8),optional            :: Eloc(Norb*Nspin)
    real(8)                     :: Eloc_(Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: ik,Lk,Nso,iorb,jorb,ispin,jspin,io,jo,js
    integer                     :: iprint
    Lk=size(Hk,3)
    Nso=Norb*Nspin
    !if(size(Hk,1)/=Nso.OR.size(Hk,2)/=Nso) stop "rdmft_get_gloc_superc error: wrong dimensions of Hk"
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc    
    if(ED_MPI_ID==0)write(LOGfile,*)"Get local GF (id=0):"
    if(ED_MPI_ID==0)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    !
    zeta_mats=zero
    zeta_real=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          !
          !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
          !G22(iw) = -[G11[iw]]*
          !G21(iw) =   G12[w]
          zeta_mats(1,1,io,io,:) = xi*wm(:) + xmu - Eloc_(io)
          zeta_mats(2,2,io,io,:) = xi*wm(:) - xmu + Eloc_(io)
          !
          !SYMMETRIES in real-frequencies   [assuming a real order parameter]
          !G22(w)  = -[G11[-w]]*
          !G21(w)  =   G12[w]
          zeta_real(1,1,io,io,:) =  dcmplx(wr(:),eps)                  + xmu - Eloc_(io)    
          zeta_real(2,2,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu - Eloc_(io) )
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                !
                zeta_mats(1,1,io,jo,:) = zeta_mats(1,1,io,jo,:) - Smats(1,ispin,jspin,iorb,jorb,:)
                zeta_mats(1,2,io,jo,:) = zeta_mats(1,2,io,jo,:) - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_mats(2,1,io,jo,:) = zeta_mats(2,1,io,jo,:) - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_mats(2,2,io,jo,:) = zeta_mats(2,2,io,jo,:) + conjg(Smats(1,ispin,jspin,iorb,jorb,:))
                !
                zeta_real(1,1,io,jo,:) = zeta_real(1,1,io,jo,:) - Sreal(1,ispin,jspin,iorb,jorb,:)
                zeta_real(1,2,io,jo,:) = zeta_real(1,2,io,jo,:) - Sreal(2,ispin,jspin,iorb,jorb,:)
                zeta_real(2,1,io,jo,:) = zeta_real(2,1,io,jo,:) - Sreal(2,ispin,jspin,iorb,jorb,:)
                zeta_real(2,2,io,jo,:) = zeta_real(2,2,io,jo,:) + conjg( Sreal(1,ispin,jspin,iorb,jorb,Lreal:1:-1) )
                !
             enddo
          enddo
       enddo
    enddo
    !
    if(ED_MPI_ID==0)call start_timer
    Gmats=zero
    Greal=zero
    do ik=1,Lk
       call add_to_gloc_superc(zeta_mats,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Gkmats)
       call add_to_gloc_superc(zeta_real,Hk(:,:,ik),Wtk(ik),hk_symm_(ik),Gkreal)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Greal = Greal + Gkreal*Wtk(ik)
       if(ED_MPI_ID==0)call eta(ik,Lk,unit=LOGfile)
    end do
    if(ED_MPI_ID==0)call stop_timer
    !
    !
    !PRINT:
    call print_gloc_superc(Gmats,Greal,iprint)
    !
  end subroutine ed_get_gloc_superc_main


  !NORMAL BETHE
  subroutine ed_get_gloc_normal_bethe(Gmats,Greal,Smats,Sreal,Hloc,iprint,wband)
    complex(8),intent(inout)    :: Gmats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                     :: iprint
    real(8),optional            :: wband
    real(8)                     :: wband_
    complex(8)                  :: zeta_mats(Nspin*Norb,Lmats)
    complex(8)                  :: zeta_real(Nspin*Norb,Lreal)
    integer                     :: i,ik,Lk,Nso,iorb,jorb,ispin,jspin,io,jo,js
    !
    wband_=1d0;if(present(wband))wband_=wband
    Nso=Norb*Nspin
    if(ED_MPI_ID==0)write(LOGfile,*)"Get Bethe GF (id=0):"
    if(ED_MPI_ID==0)write(LOGfile,*)"print in mode "//reg(txtfy(iprint))
    if(bath_type/="normal")stop "ed_get_gloc_normal_bethe error: this routine works only for diagonal problems. set bath_type=normal"
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    allocate(wm(Lmats))
    allocate(wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    !
    if(ED_MPI_ID==0)call start_timer
    !only Diagonal component are allowed.
    Gmats=zero
    Greal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_mats(io,:) = xi*wm(:)          + xmu - Hloc(ispin,ispin,iorb,iorb) - Smats(ispin,ispin,iorb,iorb,:)
          zeta_real(io,:) = dcmplx(wr(:),eps) + xmu - Hloc(ispin,ispin,iorb,iorb) - Sreal(ispin,ispin,iorb,iorb,:)
          do i=1,Lmats
             Gmats(ispin,ispin,iorb,iorb,i) = gfbethe(wm(i),zeta_mats(io,i),wband_)
          enddo
          do i=1,Lreal
             Greal(ispin,ispin,iorb,iorb,i) = gfbether(wr(i),zeta_real(io,i),wband_)
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0)call stop_timer
    !
    !
    !PRINT:
    call print_gloc_normal(Gmats,Greal,iprint)
    !
  end subroutine ed_get_gloc_normal_bethe


















  !----------------------------------------------------------------------------------------!
  ! PURPOSE: routines that actually evaluate the G_k(iw//w+) 
  !----------------------------------------------------------------------------------------!
  subroutine add_to_gloc_normal(zeta_site,Hk,hk_symm,Gkout)
    complex(8)               :: zeta_site(:,:,:)              ![Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)               :: Hk(Nspin*Norb,Nspin*Norb) 
    real(8)                  :: Wtk                    
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gkout(Nspin,Nspin,Norb,Norb,size(zeta_site,3))
    complex(8)               :: Gktmp(Nspin,Nspin,Norb,Norb,size(zeta_site,3))
    !
    complex(8)               :: Gmatrix(Nspin*Norb,Nspin*Norb)
    integer                  :: i,is,Lfreq,iorb,jorb,ispin,jspin,io,jo
    if(size(zeta_site,1)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 2 = Nspin*Norb"
    if(size(zeta_site,2)/=Nspin*Norb)stop "get_gloc_kpoint error: zeta_site wrong size 3 = Nspin*Norb"
    Lfreq = size(zeta_site,3)
    Gktmp=zero
    do i=1,Lfreq
       Gmatrix  = zeta_site(:,:,i) - Hk
       if(Nspin*Norb==1)then
          Gmatrix = one/Gmatrix
       else
          if(hk_symm) then
             call matrix_inverse_sym(Gmatrix)
          else
             call matrix_inverse(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
          endif
       endif
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

  subroutine add_to_gloc_superc(zeta_site,Hk,Wtk,hk_symm,Gkout)
    complex(8)               :: zeta_site(:,:,:,:,:)              ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)               :: Hk(Nspin*Norb,Nspin*Norb) 
    real(8)                  :: Wtk                    
    logical                  :: hk_symm                
    !output:
    complex(8),intent(inout) :: Gkout(2,Nspin,Nspin,Norb,Norb,size(zeta_site,5))
    complex(8)               :: Gktmp(2,Nspin,Nspin,Norb,Norb,size(zeta_site,5))
    !
    complex(8)               :: Gmatrix(2*Nspin*Norb ,2*Nspin*Norb)
    integer                  :: i,is,Lfreq,iorb,jorb,ispin,jspin,io,jo,Nso
    if(size(zeta_site,1)/=2.OR.size(zeta_site,2)/=2)stop "add_to_gloc_superc error: zeta_site wrong size 1 or 2 = 2 (Nambu)"
    if(size(zeta_site,3)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 4 = Nspin*Norb"
    if(size(zeta_site,4)/=Nspin*Norb)stop "add_to_gloc_superc error: zeta_site wrong size 5 = Nspin*Norb"
    Lfreq = size(zeta_site,5)
    Nso   = Nspin*Norb
    Gkout = zero
    Gktmp = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Nso,1:Nso)             = zeta_site(1,1,:,:,i) - Hk
       Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta_site(1,2,:,:,i)
       Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta_site(2,1,:,:,i)
       Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta_site(2,2,:,:,i) + Hk
       if(hk_symm) then
          call matrix_inverse_sym(Gmatrix)
       else
          call matrix_inverse(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
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









  !----------------------------------------------------------------------------------------!
  ! PURPOSE: write local GF to file
  !----------------------------------------------------------------------------------------!
  subroutine print_gloc_normal(Gmats,Greal,iprint)
    complex(8),intent(in)    :: Gmats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(in)    :: Greal(Nspin,Nspin,Norb,Norb,Lreal)
    integer                  :: iprint
    integer                  :: ik,Lk,Nso,iorb,jorb,ispin,jspin,io,jo,js
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_iw"//reg(ed_file_suffix)//".ed"
                call splot("Gloc"//reg(suffix),wm,Gmats(ispin,ispin,iorb,iorb,:))
                suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
                call splot("Gloc"//reg(suffix),wr,-dimag(Greal(ispin,ispin,iorb,iorb,:))/pi,dreal(Greal(ispin,ispin,iorb,iorb,:)))
             enddo
          enddo
       case(2)                  !print spin-diagonal, all orbitals 
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_iw"//reg(ed_file_suffix)//".ed"
                   call splot("Gloc"//reg(suffix),wm,Gmats(ispin,ispin,iorb,jorb,:))
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
                   call splot("Gloc"//reg(suffix),wr,-dimag(Greal(ispin,ispin,iorb,jorb,:))/pi,dreal(Greal(ispin,ispin,iorb,jorb,:)))
                enddo
             enddo
          enddo
       case default                  !print all off-diagonals
          write(LOGfile,*)"write all elements:"
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw"//reg(ed_file_suffix)//".ed"
                      call splot("Gloc"//reg(suffix),wm,Gmats(ispin,jspin,iorb,jorb,:))
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw"//reg(ed_file_suffix)//".ed"
                      call splot("Gloc"//reg(suffix),wr,-dimag(Greal(ispin,jspin,iorb,jorb,:))/pi,dreal(Greal(ispin,jspin,iorb,jorb,:)))
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine print_gloc_normal

  subroutine print_gloc_superc(Gmats,Greal,iprint)
    complex(8),intent(in)    :: Gmats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8),intent(in)    :: Greal(2,Nspin,Nspin,Norb,Norb,Lreal)
    integer                  :: iprint
    integer                  :: ik,Lk,Nso,iorb,jorb,ispin,jspin,io,jo,js
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       select case(iprint)
       case(1)
          write(LOGfile,*)"write spin-orbital diagonal elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_iw"//reg(ed_file_suffix)//".ed"
                call splot("Gloc"//reg(suffix),wm,Gmats(1,ispin,ispin,iorb,iorb,:))
                call splot("Floc"//reg(suffix),wm,Gmats(2,ispin,ispin,iorb,iorb,:))
                suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
                call splot("Gloc"//reg(suffix),wr,-dimag(Greal(1,ispin,ispin,iorb,iorb,:))/pi,dreal(Greal(1,ispin,ispin,iorb,iorb,:)))
                call splot("Floc"//reg(suffix),wr,-dimag(Greal(2,ispin,ispin,iorb,iorb,:))/pi,dreal(Greal(2,ispin,ispin,iorb,iorb,:)))
             enddo
          enddo
       case(2)
          write(LOGfile,*)"write spin diagonal and all orbitals elements:"
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_iw"//reg(ed_file_suffix)//".ed"
                   call splot("Gloc"//reg(suffix),wm,Gmats(1,ispin,ispin,iorb,jorb,:))
                   call splot("Floc"//reg(suffix),wm,Gmats(2,ispin,ispin,iorb,jorb,:))
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(ispin))//"_realw"//reg(ed_file_suffix)//".ed"
                   call splot("Gloc"//reg(suffix),wr,-dimag(Greal(1,ispin,ispin,iorb,jorb,:))/pi,dreal(Greal(1,ispin,ispin,iorb,jorb,:)))
                   call splot("Floc"//reg(suffix),wr,-dimag(Greal(2,ispin,ispin,iorb,jorb,:))/pi,dreal(Greal(2,ispin,ispin,iorb,jorb,:)))
                enddo
             enddo
          enddo
       case(3)
          write(LOGfile,*)"write all elements:"
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw"//reg(ed_file_suffix)//".ed"
                      call splot("Gloc"//reg(suffix),wm,Gmats(1,ispin,jspin,iorb,jorb,:))
                      call splot("Floc"//reg(suffix),wm,Gmats(2,ispin,jspin,iorb,jorb,:))
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw"//reg(ed_file_suffix)//".ed"
                      call splot("Gloc"//reg(suffix),wr,-dimag(Greal(1,ispin,jspin,iorb,jorb,:))/pi,dreal(Greal(1,ispin,jspin,iorb,jorb,:)))
                      call splot("Floc"//reg(suffix),wr,-dimag(Greal(2,ispin,jspin,iorb,jorb,:))/pi,dreal(Greal(2,ispin,jspin,iorb,jorb,:)))
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine print_gloc_superc










  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  additional interfaces to the main procedures above for the 
  ! normal and superc case:
  ! - _1b : one band case (input is independent of Nspin and Norb)
  ! - _mb : multi-bands case (input is independent of Nspin)
  !+-----------------------------------------------------------------------------+!
  subroutine ed_get_gloc_normal_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:)            :: Hk              ![Nk]
    real(8)                            :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(inout)           :: Gmats(Lmats)
    complex(8),intent(inout)           :: Greal(Lreal)
    complex(8),intent(inout)           :: Smats(Lmats)
    complex(8),intent(inout)           :: Sreal(Lreal)
    !
    complex(8),dimension(1,1,size(Hk)) :: Hk_              ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8)                         :: Gmats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                         :: Greal_(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                         :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                         :: Sreal_(Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional                   :: Eloc(Norb*Nspin)
    real(8)                            :: Eloc_(Norb*Nspin)
    logical,optional                   :: hk_symm(size(Hk,1))
    logical                            :: hk_symm_(size(Hk,1))
    integer                            :: iprint
    if(Norb>1)stop "ed_get_gloc_normal_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop"ed_get_gloc_normal_1b error: Nspin > 1 in 1-band routine" 
    Gmats_(1,1,1,1,:) = Gmats(:)
    Greal_(1,1,1,1,:) = Greal(:)
    Smats_(1,1,1,1,:) = Smats(:)
    Sreal_(1,1,1,1,:) = Sreal(:)
    Hk_(1,1,:)        = Hk(:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk_,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:) = Gmats_(1,1,1,1,:)
    Greal(:) = Greal_(1,1,1,1,:)
    Smats(:) = Smats_(1,1,1,1,:)
    Sreal(:) = Sreal_(1,1,1,1,:)
  end subroutine ed_get_gloc_normal_1b

  subroutine ed_get_gloc_normal_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(Norb,Norb,Lreal)
    !
    complex(8)                  :: Gmats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Norb*Nspin)
    real(8)                     :: Eloc_(Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: iprint
    if(Nspin>1)stop"ed_get_gloc_normal_mb error: Nspin > 1 in M-band routine" 
    Gmats_(1,1,:,:,:) = Gmats(:,:,:)
    Greal_(1,1,:,:,:) = Greal(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    Sreal_(1,1,:,:,:) = Sreal(:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call ed_get_gloc_normal_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:,:) = Gmats_(1,1,:,:,:)
    Greal(:,:,:) = Greal_(1,1,:,:,:)
    Smats(:,:,:) = Smats_(1,1,:,:,:)
    Sreal(:,:,:) = Sreal_(1,1,:,:,:)
  end subroutine ed_get_gloc_normal_mb

  subroutine ed_get_gloc_superc_1b(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:)     :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Lmats)
    complex(8),intent(inout)    :: Greal(2,Lreal)
    complex(8),intent(inout)    :: Smats(2,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Lreal)
    !
    complex(8)                  :: Hk_(1,1,size(Hk))              ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8)                  :: Gmats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(2,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(2,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Norb*Nspin)
    real(8)                     :: Eloc_(Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,1))
    logical                     :: hk_symm_(size(Hk,1))
    integer                     :: iprint
    if(Norb>1)stop"ed_get_gloc_superc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop"ed_get_gloc_superc_1b error: Nspin > 1 in 1-band routine" 
    Gmats_(:,1,1,1,1,:) = Gmats(:,:)
    Greal_(:,1,1,1,1,:) = Greal(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    Hk_(1,1,:)          = Hk
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_gloc_superc_main(Hk_,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Sreal(:,:) = Sreal_(:,1,1,1,1,:)
  end subroutine ed_get_gloc_superc_1b

  subroutine ed_get_gloc_superc_mb(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint,Eloc,hk_symm)
    complex(8),dimension(:,:,:) :: Hk              ![Norb*Nspin][Norb*Nspin][Nk]
    real(8)                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout)    :: Gmats(2,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Greal(2,Norb,Norb,Lreal)
    complex(8),intent(inout)    :: Smats(2,Norb,Norb,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Norb,Norb,Lreal)
    !
    complex(8)                  :: Gmats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Greal_(2,Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                  :: Sreal_(2,Nspin,Nspin,Norb,Norb,Lreal)
    !
    real(8),optional            :: Eloc(Norb*Nspin)
    real(8)                     :: Eloc_(Norb*Nspin)
    logical,optional            :: hk_symm(size(Hk,3))
    logical                     :: hk_symm_(size(Hk,3))
    integer                     :: iprint
    if(Nspin>1)stop "ed_get_gloc_superc_mb error: Nspin > 1 in M-band routine" 
    Gmats_(:,1,1,:,:,:) = Gmats(:,:,:,:)
    Greal_(:,1,1,:,:,:) = Greal(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_gloc_superc_main(Hk,Wtk,Gmats_,Greal_,Smats_,Sreal_,iprint,Eloc_,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Sreal(:,:,:,:) = Sreal_(:,1,1,:,:,:)
  end subroutine ed_get_gloc_superc_mb

  subroutine ed_get_gloc_normal_bethe_1b(Gmats,Greal,Smats,Sreal,Hloc,iprint,wband)
    complex(8),intent(inout) :: Gmats(Lmats)
    complex(8),intent(inout) :: Greal(Lreal)
    complex(8),intent(inout) :: Smats(Lmats)
    complex(8),intent(inout) :: Sreal(Lreal)
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                  :: iprint
    real(8),optional         :: wband
    real(8)                  :: wband_
    complex(8)               :: Gmats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Greal_(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(Nspin,Nspin,Norb,Norb,Lreal)
    !
    if(Norb>1)stop "ed_get_gloc_normal_bethe_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop"ed_get_gloc_normal_bethe_1b error: Nspin > 1 in 1-band routine" 
    !
    Gmats_(1,1,1,1,:) = Gmats(:)
    Greal_(1,1,1,1,:) = Greal(:)
    Smats_(1,1,1,1,:) = Smats(:)
    Sreal_(1,1,1,1,:) = Sreal(:)
    wband_=1d0       ;if(present(wband))wband_=wband
    call ed_get_gloc_normal_bethe(Gmats_,Greal_,Smats_,Sreal_,Hloc,iprint,wband_)
    Gmats(:) = Gmats_(1,1,1,1,:)
    Greal(:) = Greal_(1,1,1,1,:)
    Smats(:) = Smats_(1,1,1,1,:)
    Sreal(:) = Sreal_(1,1,1,1,:)
  end subroutine ed_get_gloc_normal_bethe_1b

  subroutine ed_get_gloc_normal_bethe_mb(Gmats,Greal,Smats,Sreal,Hloc,iprint,wband)
    complex(8),intent(inout) :: Gmats(Norb,Norb,Lmats)
    complex(8),intent(inout) :: Greal(Norb,Norb,Lreal)
    complex(8),intent(inout) :: Smats(Norb,Norb,Lmats)
    complex(8),intent(inout) :: Sreal(Norb,Norb,Lreal)
    complex(8)               :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                  :: iprint
    real(8),optional         :: wband
    real(8)                  :: wband_
    complex(8)               :: Gmats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Greal_(Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)               :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Sreal_(Nspin,Nspin,Norb,Norb,Lreal)
    if(Nspin>1)stop"ed_get_gloc_normal_bethe_mb error: Nspin > 1 in M-band routine" 
    Gmats_(1,1,:,:,:) = Gmats(:,:,:)
    Greal_(1,1,:,:,:) = Greal(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    Sreal_(1,1,:,:,:) = Sreal(:,:,:)
    wband_=1d0       ;if(present(wband))wband_=wband
    call ed_get_gloc_normal_bethe(Gmats_,Greal_,Smats_,Sreal_,Hloc,iprint,wband_)
    Gmats(:,:,:) = Gmats_(1,1,:,:,:)
    Greal(:,:,:) = Greal_(1,1,:,:,:)
    Smats(:,:,:) = Smats_(1,1,:,:,:)
    Sreal(:,:,:) = Sreal_(1,1,:,:,:)
  end subroutine ed_get_gloc_normal_bethe_mb

end module ED_GLOC
