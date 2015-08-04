module ED_WEISS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_IOTOOLS,   only:reg,txtfy,splot,store_data
  USE SF_TIMER
  USE SF_ARRAYS,    only:arange
  USE SF_LINALG,    only:eye,inv,inv_sym
  USE SF_MISC,      only:assert_shape
  implicit none
  private


  interface ed_get_weiss
     module procedure ed_get_weiss_field_normal_main
     module procedure ed_get_weiss_field_superc_main
     module procedure ed_get_weiss_field_normal_1b
     module procedure ed_get_weiss_field_superc_1b
     module procedure ed_get_weiss_field_normal_mb
     module procedure ed_get_weiss_field_superc_mb
  end interface ed_get_weiss


  interface ed_get_weiss_lattice
     module procedure ed_get_weiss_field_normal_lattice
     module procedure ed_get_weiss_field_superc_lattice
     module procedure ed_get_weiss_field_normal_lattice_1b
     module procedure ed_get_weiss_field_superc_lattice_1b
     module procedure ed_get_weiss_field_normal_lattice_mb
     module procedure ed_get_weiss_field_superc_lattice_mb
  end interface ed_get_weiss_lattice


  public :: ed_get_weiss
  public :: ed_get_weiss_lattice



  real(8),dimension(:),allocatable :: wm
  character(len=20)                :: suffix

contains




  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! NORMAL PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_normal_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)   :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)   :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    integer                                      :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable      :: zeta_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable      :: Smats_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable      :: invGloc_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable      :: calG0_site![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                      :: Nspin,Norb,Nso,Lmats
    integer                                      :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Testing part:
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lmats = size(Gloc,5)
    Nso   = Nspin*Norb
    call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_main","Gloc")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_main","Smats")
    call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_main","Weiss")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"ed_get_weiss_field_normal_main","Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix and create the zeta_site
    do i=1,Lmats
       zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc,Nspin,Norb) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
       invGloc_site(:,:,i) = nn2so_reshape(Gloc(:,:,:,:,i),Nspin,Norb)
       Smats_site(:,:,i)   = nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    if(cg_scheme=="weiss")then
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
    else
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
    endif
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nspsin,Nspin,Norb,Norb] matrix
    do i=1,Lmats
       Weiss(:,:,:,:,i) = so2nn_reshape(calG0_site(:,:,i),Nspin,Norb)
    enddo
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       call print_weiss(Weiss,iprint)
    endif
  end subroutine ed_get_weiss_field_normal_main

  subroutine ed_get_weiss_field_normal_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:),intent(in)             :: Gloc
    complex(8),dimension(size(Gloc)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc)),intent(inout) :: Weiss
    complex(8),dimension(1,1,1,1)                  :: Hloc
    integer                                        :: iprint
    !aux
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Gloc_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Smats_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Weiss_
    Gloc_(1,1,1,1,:) = Gloc(:)
    Smats_(1,1,1,1,:) = Smats(:)
    call ed_get_weiss_field_normal_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:) = Weiss_(1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_1b

  subroutine ed_get_weiss_field_normal_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)      :: Gloc  ![Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:),intent(in)    :: Hloc  ![Nspin][Nspin][Norb][Norb]
    integer                                     :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:),allocatable :: Gloc_ ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:),allocatable :: Weiss_!
    integer                                     :: Nspin,Norb,Lfreq
    Nspin = 1
    Norb  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    call assert_shape(Gloc,[Norb,Norb,Lfreq],"ed_get_weiss_field_normal_mb","Gloc")
    call assert_shape(Smats,[Norb,Norb,Lfreq],"ed_get_weiss_field_normal_mb","Smats")
    call assert_shape(Weiss,[Norb,Norb,Lfreq],"ed_get_weiss_field_normal_mb","Weiss")
    call assert_shape(Hloc,[1,1,Norb,Norb],"ed_get_weiss_field_normal_mb","Hloc")
    allocate(Gloc_(1,1,Norb,Norb,Lfreq))
    allocate(Smats_(1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(1,1,Norb,Norb,Lfreq))
    Gloc_(1,1,:,:,:) = Gloc(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    call ed_get_weiss_field_normal_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:) = Weiss_(1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_mb








  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! SUPERCONDUCTING PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_superc_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        !
    complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    integer                                         :: Nspin,Norb,Nso,Nso2,Lmats
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Testing part:
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_main","Gloc")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_main","Smats")
    call assert_shape(Weiss,[2,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_main","Weiss")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"ed_get_weiss_field_superc_main","Hloc")
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_site(io,io,:)         = xi*wm(:) + xmu 
          zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu 
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Hloc(ispin,jspin,iorb,jorb) - Smats(1,ispin,jspin,iorb,jorb,:)
                zeta_site(io,jo+Nso,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + Hloc(ispin,jspin,iorb,jorb) + conjg(Smats(1,ispin,jspin,iorb,jorb,:))
                !
                invGloc_site(io,jo,:)        = Gloc(1,ispin,jspin,iorb,jorb,:)
                invGloc_site(io,jo+Nso,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ispin,jspin,iorb,jorb,:))
                !
                Smats_site(io,jo,:)          = Smats(1,ispin,jspin,iorb,jorb,:)
                Smats_site(io,jo+Nso,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ispin,jspin,iorb,jorb,:))
             enddo
          enddo
       enddo
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    if(cg_scheme=="weiss")then
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
    else
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
    endif
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Weiss(1,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                Weiss(2,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
             enddo
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       call print_weiss(Weiss(1,:,:,:,:,:),iprint,"_normal")
       call print_weiss(Weiss(2,:,:,:,:,:),iprint,"_anomal")
    endif
  end subroutine ed_get_weiss_field_superc_main

  subroutine ed_get_weiss_field_superc_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)               :: Gloc
    complex(8),dimension(2,size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2)),intent(inout) :: Weiss
    complex(8),dimension(1,1,1,1)                      :: Hloc
    integer                                            :: iprint
    !aux
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Weiss_
    call assert_shape(Gloc,[2,size(Gloc,2)],"ed_get_weiss_field_superc_1b","Gloc")
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call ed_get_weiss_field_superc_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_1b

  subroutine ed_get_weiss_field_superc_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc   ![2][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats  !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Weiss  !
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc   ![Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_  ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_ !
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Weiss_ !
    integer                                       :: Nspin,Norb,Lfreq
    !
    Nspin=1
    Norb =size(Gloc,2)
    Lfreq=size(Gloc,4)
    call assert_shape(Gloc,[2,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_mb","Gloc")
    call assert_shape(Smats,[2,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_mb","Smats")
    call assert_shape(Weiss,[2,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_mb","Weiss")
    call assert_shape(Hloc,[1,1,Norb,Norb],"ed_get_weiss_field_superc_mb","Hloc")
    allocate(Gloc_(2,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(2,1,1,Norb,Norb,Lfreq))
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call ed_get_weiss_field_superc_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_mb










  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! NORMAL PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_normal_lattice(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)   :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)   :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                        :: iprint       
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable  :: Weiss_tmp    ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable        :: zeta_site    ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable        :: Smats_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable        :: invGloc_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable        :: calG0_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                        :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                        :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_lattice","Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_lattice","Smats")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_normal_lattice","Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"ed_get_weiss_field_normal_lattice","Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss_tmp = zero
    Weiss     = zero
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
          invGloc_site(:,:,i) = nn2so_reshape(Gloc(ilat,:,:,:,:,i),Nspin,Norb)
          Smats_site(:,:,i)   = nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
          calG0_site(:,:,:) = invGloc_site(:,:,:) + Smats_site(:,:,:)
          do i=1,Lmats
             call inv(calG0_site(:,:,i))
          enddo
       else
          ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          calG0_site(:,:,:) = zeta_site(:,:,:) - invGloc_site(:,:,:)
       endif
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss_tmp(ilat,ispin,jspin,iorb,jorb,1:Lmats) = calG0_site(io,jo,1:Lmats)
                enddo
             enddo
          enddo
       enddo
    end do mpi_site_loop
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_normal_lattice

  subroutine ed_get_weiss_field_normal_lattice_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)                          :: Gloc
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,1),1,1,1,1)                    :: Hloc
    integer                                                       :: iprint
    !aux
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Weiss_
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call ed_get_weiss_field_normal_lattice(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_lattice_1b

  subroutine ed_get_weiss_field_normal_lattice_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc  ![Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                       :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,1)
    Nspin = 1
    Norb  = size(Gloc,2)
    Lfreq = size(Gloc,4)
    call assert_shape(Gloc,[Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_normal_lattice_mb","Gloc")
    call assert_shape(Smats,[Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_normal_lattice_mb","Smats")
    call assert_shape(Weiss,[Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_normal_lattice_mb","Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],"ed_get_weiss_field_normal_lattice_mb","Hloc")
    allocate(Gloc_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(Nlat,1,1,Norb,Norb,Lfreq))
    !
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call ed_get_weiss_field_normal_lattice(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_lattice_mb











  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! SUPERCONDUCTING PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_superc_lattice(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)   :: Gloc         ! [2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)   :: Smats        ! 
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                          :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable  :: Weiss_tmp    ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable          :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable          :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable          :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable          :: calG0_site   !
    integer                                          :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                          :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    !
    !Testing part:
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lmats = size(Gloc,7)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_lattice","Gloc")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_lattice","Smats")
    call assert_shape(Weiss,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"ed_get_weiss_field_superc_lattice","Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"ed_get_weiss_field_superc_lattice","Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Weiss_tmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss_tmp   = zero
    Weiss       = zero
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   !
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(1,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
          calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
          do i=1,Lmats
             call inv(calG0_site(:,:,i))
          enddo
       else
          ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
       endif
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss_tmp(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   Weiss_tmp(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do mpi_site_loop
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_superc_lattice

  subroutine ed_get_weiss_field_superc_lattice_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)                          :: Gloc ![2][Nlat][Lmats]
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,2),1,1,1,1)                      :: Hloc
    integer                                                         :: iprint
    !aux
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Gloc_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Smats_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Weiss_
    call assert_shape(Gloc,[2,size(Gloc,2),size(Gloc,3)],"ed_get_weiss_field_superc_lattice_1b","Gloc")
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call ed_get_weiss_field_superc_lattice(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:) = Weiss_(:,:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_lattice_1b

  subroutine ed_get_weiss_field_superc_lattice_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc  ![2][Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gloc_ ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                         :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,2)
    Nspin = 1
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)
    call assert_shape(Gloc,[2,Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_lattice_mb","Gloc")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_lattice_mb","Smats")
    call assert_shape(Weiss,[2,Nlat,Norb,Norb,Lfreq],"ed_get_weiss_field_superc_lattice_mb","Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],"ed_get_weiss_field_superc_lattice_mb","Hloc")
    allocate(Gloc_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(2,Nlat,1,1,Norb,Norb,Lfreq))
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call ed_get_weiss_field_superc_lattice(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:,:) = Weiss_(:,:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_lattice_mb










  !+-----------------------------------------------------------------------------+!
  !PURPOSE: print a local GF according to iprint variable
  !+-----------------------------------------------------------------------------+!
  subroutine print_weiss(Weiss,iprint,sfx)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Weiss ![Nspin][Nspin][Norb][Norb][Lfreq]
    character(len=*),optional                  :: sfx
    character(len=16)                          :: fname
    integer,intent(in)                         :: iprint
    integer                                    :: Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nspin = size(Weiss,1)
    Norb  = size(Weiss,3)
    call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,size(Weiss,5)],"print_weiss","Weiss")
    !
    if(cg_scheme=="weiss")then
       fname="WeissField"
    else
       fname="Delta"
    endif
    if(present(sfx))fname=reg(fname)//reg(sfx)
    !
    select case(iprint)
    case (0)
       write(LOGfile,*)"Delta//WeissField not written to file."
    case(1)
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call splot(reg(suffix),wm,Weiss(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
    case(2)
       write(LOGfile,*)"write spin diagonal and all orbitals elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call splot(reg(suffix),wm,Weiss(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
    case default
       write(LOGfile,*)"write all elements:"
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                   call splot(reg(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine print_weiss

  subroutine print_weiss_lattice(Weiss,iprint,sfx)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Weiss ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    character(len=*),optional                    :: sfx
    character(len=16)                            :: fname
    integer,intent(in)                           :: iprint
    integer                                      :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nlat  = size(Weiss,1)
    Nspin = size(Weiss,2)
    Norb  = size(Weiss,3)
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,size(Weiss,6)],"print_weiss_lattice","Weiss")
    !
    if(cg_scheme=="weiss")then
       fname="WeissField"
    else
       fname="Delta"
    endif
    if(present(sfx))fname=reg(fname)//reg(sfx)
    !
    select case(iprint)
    case (0)
       write(LOGfile,*)"Delta//WeissField not written to file."
    case(1)
       write(LOGfile,*)"write spin-orbital diagonal elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call store_data(reg(suffix),Weiss(:,ispin,ispin,iorb,iorb,:),wm)
          enddo
       enddo
    case(2)
       write(LOGfile,*)"write spin diagonal and all orbitals elements:"
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                call store_data(reg(suffix),Weiss(:,ispin,ispin,iorb,jorb,:),wm)
             enddo
          enddo
       enddo
    case default
       write(LOGfile,*)"write all elements:"
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//"_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                   call store_data(reg(suffix),Weiss(:,ispin,jspin,iorb,jorb,:),wm)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine print_weiss_lattice



end module ED_WEISS
