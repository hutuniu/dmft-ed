!###############################################################
! PROGRAM  : RDMFT_FUNX
! PURPOSE  : Compute Local GFunction for generel real-space scheme
!###############################################################
module ED_WRAP_WEISS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_ARRAYS,    only: arange
  USE SF_TIMER
  USE SF_IOTOOLS,   only:reg,sread,free_unit
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym,matrix_diagonalize,matrix_inverse_gj
  implicit none
  private


  interface ed_get_weiss_normal
     module procedure ed_get_weiss_field_normal_eloc,ed_get_weiss_field_normal_hloc
  end interface ed_get_weiss_normal

  interface ed_get_weiss_superc
     module procedure ed_get_weiss_field_superc_eloc,ed_get_weiss_field_superc_hloc
  end interface ed_get_weiss_superc
  public :: ed_get_weiss_normal
  public :: ed_get_weiss_superc



  !OBSOLETE:
  interface rdmft_get_weiss_field
     module procedure rdmft_get_weiss_field_normal,rdmft_get_weiss_field_superc
  end interface rdmft_get_weiss_field
  public :: rdmft_get_weiss_field
  public :: rdmft_get_weiss_field_mb


  real(8),dimension(:),allocatable        :: wm

contains


  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! NORMAL PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_normal_eloc(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    !aux
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    complex(8)                                  :: Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: zeta_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: Hloc(Nspin*Norb,Nspin*Norb)
    complex(8)                                  :: Smats_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(Nspin*Norb,Nspin*Norb,Lmats)
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
       !and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_site(io,io,:) = xi*wm(:) + xmu - Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_site(io,jo,:)    = zeta_site(io,jo,:) - Smats(ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo,:) = Gloc(ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo,:)   = Smats(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
       !
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call matrix_inverse(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          !if calG0 is required get it as:
          ![calG0]_ilat^-1 = [Gloc]_ilat^-1 + [Smats]_ilat
          ![calG0]_ilat = [[calG0]_ilat^-1]^-1
          calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
          do i=1,Lmats
             call matrix_inverse(calG0_site(:,:,i))
          enddo
       else
          !else if Delta is required get is as:
          ! [Delta]_ilat = [Zeta]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          !              = [Zeta]_ilat - [Eloc]_ilat - [Gloc]_ilat^-1
          !              = [iw+mu-Eloc]_ilat         - [Gloc]_ilat^-1
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
                   Weiss_tmp(ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                enddo
             enddo
          enddo
       enddo
    end do mpi_site_loop
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_get_weiss_field_normal_eloc

  subroutine ed_get_weiss_field_normal_hloc(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    !aux
    complex(8)                                  :: Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: zeta_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(Nspin*Norb,Nspin*Norb,Lmats)
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
       !and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_site(io,io,:) = xi*wm(:) + xmu
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_site(io,jo,:)    = zeta_site(io,jo,:) - Hloc(ilat,ispin,jspin,iorb,jorb) - Smats(ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo,:) = Gloc(ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo,:)   = Smats(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
       !
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call matrix_inverse(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          !if calG0 is required get it as:
          ![calG0]_ilat^-1 = [Gloc]_ilat^-1 + [Smats]_ilat
          ![calG0]_ilat = [[calG0]_ilat^-1]^-1
          calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
          do i=1,Lmats
             call matrix_inverse(calG0_site(:,:,i))
          enddo
       else
          !else if Delta is required get is as:
          ! [Delta]_ilat = [Zeta]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          !              = [Zeta]_ilat - [Eloc]_ilat - [Gloc]_ilat^-1
          !              = [iw+mu-Eloc]_ilat         - [Gloc]_ilat^-1
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
                   Weiss_tmp(ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                enddo
             enddo
          enddo
       enddo
    end do mpi_site_loop
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_get_weiss_field_normal_hloc



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! SUPERCONDUCTING PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_superc_eloc(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    !aux
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    complex(8)                                  :: Weiss_tmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: zeta_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    integer                                     :: Nlso,Nso
    !
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Nso =Nspin*Norb
    Nlso=Nlat*Nspin*Norb
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
       !and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Eloc_(js)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Eloc_(js)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
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
                   Smats_site(io+Nso,jo+Nso,:)  = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call matrix_inverse(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          !if calG0 is required get it as:
          ![calG0]_ilat^-1 = [Gloc]_ilat^-1 + [Smats]_ilat
          ![calG0]_ilat = [[calG0]_ilat^-1]^-1
          calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
          do i=1,Lmats
             call matrix_inverse(calG0_site(:,:,i))
          enddo
       else
          !else if Delta is required get is as:
          ! [Delta]_ilat = [Zeta]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          !              = [Zeta]_ilat - [Eloc]_ilat - [Gloc]_ilat^-1
          !              = [iw+mu-Eloc]_ilat         - [Gloc]_ilat^-1
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
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_get_weiss_field_superc_eloc



  subroutine ed_get_weiss_field_superc_hloc(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    !aux
    complex(8)                                  :: Weiss_tmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: zeta_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    integer                                     :: Nlso,Nso
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Nso =Nspin*Norb
    Nlso=Nlat*Nspin*Norb
    mpi_site_loop: do ilat=1+mpiID,Nlat,mpiSIZE
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
       !and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,jspin,iorb,jorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,jspin,iorb,jorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
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
                   Smats_site(io+Nso,jo+Nso,:)  = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call matrix_inverse(invGloc_site(:,:,i))
       enddo
       !
       if(cg_scheme=="weiss")then
          !if calG0 is required get it as:
          ![calG0]_ilat^-1 = [Gloc]_ilat^-1 + [Smats]_ilat
          ![calG0]_ilat = [[calG0]_ilat^-1]^-1
          calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
          do i=1,Lmats
             call matrix_inverse(calG0_site(:,:,i))
          enddo
       else
          !else if Delta is required get is as:
          ! [Delta]_ilat = [Zeta]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
          !              = [Zeta]_ilat - [Eloc]_ilat - [Gloc]_ilat^-1
          !              = [iw+mu-Eloc]_ilat         - [Gloc]_ilat^-1
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
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_get_weiss_field_superc_hloc











  !###################################################################################################
  !###################################################################################################
  !               POSSIBLY OBSOLETE ROUTINES NOW SUPERSEDED BY TOP ROUTINE HERE
  !               GETTING G0 FOR ANY Norb*Nspin*Nlat*Nk NORMAL GREEN'S FUNCTION
  !                (these routines are left temporarily for back-compatibility)
  !                                     ( to be removed)
  !###################################################################################################
  !###################################################################################################
  subroutine rdmft_get_weiss_field_normal(Nsites,Gmats,Smats,Delta,Eloc)
    integer                  :: Nsites
    complex(8),intent(inout) :: Delta(Nsites,Lmats)
    complex(8),intent(in)    :: Gmats(Nsites,Lmats)
    complex(8),intent(in)    :: Smats(Nsites,Lmats)    
    real(8),optional      :: Eloc(Nsites)
    !
    complex(8)               :: Delta_tmp(Nsites,Lmats)
    complex(8),allocatable   :: Hloc(:,:,:,:)
    integer                  :: ilat,i
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta=zero
    Delta_tmp=zero
    allocate(Hloc(1,1,1,1))
    do ilat=1+mpiID,Nsites,mpiSIZE
       ! set Hloc
       Hloc(1,1,1,1)=zero
       if(present(eloc))Hloc(1,1,1,1)=eloc(ilat)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             Delta_tmp(ilat,i)  =   one/(one/Gmats(ilat,i) + Smats(ilat,i))
          else
             ! same as for the weiss case
             Delta_tmp(ilat,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(ilat,i) - one/Gmats(ilat,i)
          endif
       end do
    end do
    deallocate(Hloc)
    call MPI_ALLREDUCE(Delta_tmp(:,1:Lmats),Delta,Nsites*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine rdmft_get_weiss_field_normal

  subroutine rdmft_get_weiss_field_superc(Nsites,Gmats,Smats,Delta,Eloc)
    integer                  :: Nsites
    real(8),optional         :: eloc(Nsites)
    complex(8),intent(inout) :: Delta(2,Nsites,Lmats)
    complex(8),intent(inout) :: Gmats(2,Nsites,Lmats)
    complex(8),intent(inout) :: Smats(2,Nsites,Lmats)
    complex(8)               :: Delta_tmp(2,Nsites,Lmats)
    complex(8)               :: calG(2,Lmats),cdet
    integer                  :: ilat,i
    complex(8),allocatable   :: Hloc(:,:,:,:)
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta=zero
    Delta_tmp=zero
    allocate(Hloc(1,1,1,1))
    !+- GET INDEPENDENT SITES HYBRIDIZATION FUNCTION AND FIT THE BATHS -+!
    do ilat=1+mpiID,Nsites,mpiSIZE
       Hloc(1,1,1,1)=zero
       if(present(eloc))Hloc(1,1,1,1)=eloc(ilat)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             cdet                = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             calG(1,i)           = conjg(Gmats(1,ilat,i))/cdet + Smats(1,ilat,i)
             calG(2,i)           =  Gmats(2,ilat,i)/cdet + Smats(2,ilat,i) 
             cdet                =  abs(calG(1,i))**2 + (calG(2,i))**2
             Delta_tmp(1,ilat,i) =  conjg(calG(1,i))/cdet
             Delta_tmp(2,ilat,i) =  calG(2,i)/cdet
          else
             cdet                = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             Delta_tmp(1,ilat,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(1,ilat,i) - conjg(Gmats(1,ilat,i))/cdet 
             Delta_tmp(2,ilat,i) = -(Gmats(2,ilat,i)/cdet + Smats(2,ilat,i))
          endif
       end do
    end do
    deallocate(Hloc)
    call MPI_ALLREDUCE(Delta_tmp,Delta,2*Nsites*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine rdmft_get_weiss_field_superc

  subroutine rdmft_get_weiss_field_mb(Nsites,Gmats,Smats,Delta,Hloc)
    ! inputs
    integer                  :: Nsites
    complex(8),intent(inout) :: Delta(Nsites,Nspin,Nspin,Norb,Norb,Lmats)     ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(in)    :: Gmats(Nsites,Nspin,Nspin,Norb,Norb,Lmats)     ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(in)    :: Smats(Nsites,Nspin,Nspin,Norb,Norb,Lmats) ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8)               :: Hloc(Nsites,Nspin,Nspin,Norb,Norb)  ! [Nlat][Nspin][Nspin][Norb][Norb]
    ! auxiliary vars
    complex(8)               :: Delta_tmp(Nsites,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Gloc_tmp(Norb,Norb)
    real(8),allocatable      :: Id(:,:)    
    integer                  :: ilat,i,iorb,ispin,jorb
    logical                  :: check
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    allocate(Id(Norb,Norb))
    !
    Id=0.d0
    do iorb=1,Norb
       Id(iorb,iorb)=1.d0
    end do
    check=.true.
    !
    Delta=zero
    Delta_tmp=zero
    do ilat=1+mpiID,Nsites,mpiSIZE
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             ! (ONLY DIAGONAL SPIN CALCULATIONS!!!!!)
             do ispin=1,Nspin
                if(Norb.gt.1) then
                   ! select case(bath_type)
                   ! case default
                   !    do iorb=1,Norb
                   !       Delta_tmp(ilat,ispin,ispin,iorb,iorb,i) = one/Gmats(ilat,ispin,ispin,iorb,iorb,i) + Smats(ilat,ispin,ispin,iorb,iorb,i)
                   !       Delta_tmp(ilat,ispin,ispin,iorb,iorb,i) = one/Delta_tmp(ilat,ispin,ispin,iorb,iorb,i)                         
                   !    end do
                   ! case ('hybrid')
                   !    call matrix_inverse(Gmats(ilat,ispin,ispin,:,:,i)) 
                   !    Delta_tmp(ilat,ispin,ispin,:,:,i) = Gmats(ilat,ispin,ispin,:,:,i) + Smats(ilat,ispin,ispin,:,:,i)
                   !    call matrix_inverse(Delta_tmp(ilat,ispin,ispin,:,:,i))
                   ! end select 
                   Gloc_tmp(:,:)=Gmats(ilat,ispin,ispin,:,:,i)
                   call matrix_inverse(Gloc_tmp)
                   Delta_tmp(ilat,ispin,ispin,:,:,i) =  Gloc_tmp + Smats(ilat,ispin,ispin,:,:,i)                    
                   call  matrix_inverse(Delta_tmp(ilat,ispin,ispin,:,:,i))
                else 
                   Delta_tmp(ilat,ispin,ispin,1,1,i) = one/Gmats(ilat,ispin,ispin,1,1,i) +  Smats(ilat,ispin,ispin,1,1,i)                    
                   Delta_tmp(ilat,ispin,ispin,1,1,i) =    one/Delta_tmp(ilat,ispin,ispin,1,1,i) 
                end if
             end do
          else 
             do ispin=1,Nspin 
                if(Norb.gt.1) then
                   ! select case(bath_type)
                   ! case default
                   !    do iorb=1,Norb
                   !       Delta_tmp(ilat,ispin,ispin,iorb,iorb,i) = Id(iorb,iorb)*(xi*wm(i) + xmu) - Hloc(ilat,ispin,ispin,iorb,iorb) &
                   !            - Smats(ilat,ispin,ispin,iorb,iorb,i) - one/Gmats(ilat,ispin,ispin,iorb,iorb,i)
                   !    end do
                   ! case ('hybrid')
                   !    call matrix_inverse(Gmats(ilat,ispin,ispin,:,:,i)) 
                   !    Delta_tmp(ilat,ispin,ispin,:,:,i) = Id(:,:)*(xi*wm(i) + xmu) - Hloc(ilat,ispin,ispin,:,:) &
                   !         - Smats(ilat,ispin,ispin,:,:,i) - Gmats(ilat,ispin,ispin,:,:,i)                   
                   ! end select
                   Gloc_tmp(:,:)=Gmats(ilat,ispin,ispin,:,:,i)
                   call matrix_inverse(Gloc_tmp)
                   Delta_tmp(ilat,ispin,ispin,:,:,i) = Id(:,:)*(xi*wm(i) + xmu) - Hloc(ilat,ispin,ispin,:,:) &
                        - Smats(ilat,ispin,ispin,:,:,i) - GLoc_tmp
                else
                   Delta_tmp(ilat,ispin,ispin,1,1,i) = Id(1,1)*(xi*wm(i) + xmu) - Hloc(ilat,ispin,ispin,1,1) &
                        - Smats(ilat,ispin,ispin,1,1,i) - one/Gmats(ilat,ispin,ispin,1,1,i)
                end if
             end do
          endif
       end do
    end do

    call MPI_ALLREDUCE(Delta_tmp(:,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats),Delta,Nsites*Lmats*Nspin*Nspin*Norb*Norb, &
         MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)

  end subroutine rdmft_get_weiss_field_mb














end module ED_WRAP_WEISS
