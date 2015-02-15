module ED_WRAP_WEISS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_ARRAYS,    only: arange
  USE SF_TIMER
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym,matrix_diagonalize,matrix_inverse_gj
  implicit none
  private


  interface ed_get_weiss_lattice
     module procedure &
          ed_get_weiss_field_normal_eloc,   &
          ed_get_weiss_field_normal_eloc_1b,&
          ed_get_weiss_field_normal_eloc_mb,&
          ed_get_weiss_field_normal_hloc,   &
          ed_get_weiss_field_normal_hloc_1b,&
          ed_get_weiss_field_normal_hloc_mb,&
          ed_get_weiss_field_superc_eloc,   &
          ed_get_weiss_field_superc_eloc_1b,&
          ed_get_weiss_field_superc_eloc_mb,&
          ed_get_weiss_field_superc_hloc,   &
          ed_get_weiss_field_superc_hloc_1b,&
          ed_get_weiss_field_superc_hloc_mb
  end interface ed_get_weiss_lattice
  public :: ed_get_weiss_lattice


  real(8),dimension(:),allocatable        :: wm

contains


  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! NORMAL PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_normal_eloc_1b(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(Nlat,Lmats)
    complex(8)                                  :: Smats(Nlat,Lmats)
    complex(8)                                  :: Weiss(Nlat,Lmats)
    !
    complex(8)                                  :: Gloc_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    if(Norb>1)stop "ed_get_weiss_field_normal_eloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_normal_eloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_normal_eloc(Gloc_,Smats_,Weiss_,Eloc_)
    Gloc(:,:) = Gloc_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_eloc_1b

  subroutine ed_get_weiss_field_normal_eloc_mb(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nlat,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Gloc_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    if(Nspin>1)stop "ed_get_weiss_field_normal_eloc_1m error: Nspin > 1 in M-band routine" 
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_normal_eloc(Gloc_,Smats_,Weiss_,Eloc_)
    Gloc(:,:,:,:) = Gloc_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_eloc_mb

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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_normal_eloc


  subroutine ed_get_weiss_field_normal_hloc_1b(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(Nlat,Lmats)
    complex(8)                                  :: Smats(Nlat,Lmats)
    complex(8)                                  :: Weiss(Nlat,Lmats)
    !
    complex(8)                                  :: Gloc_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    if(Norb>1)stop "ed_get_weiss_field_normal_hloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_normal_hloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call ed_get_weiss_field_normal_hloc(Gloc_,Smats_,Weiss_,Hloc)
    Gloc(:,:) = Gloc_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_hloc_1b

  subroutine ed_get_weiss_field_normal_hloc_mb(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nlat,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Gloc_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    if(Nspin>1)stop "ed_get_weiss_field_normal_hloc_mb error: Nspin > 1 in 1-band routine" 
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call ed_get_weiss_field_normal_hloc(Gloc_,Smats_,Weiss_,Hloc)
    Gloc(:,:,:,:) = Gloc_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_hloc_mb

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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_normal_hloc



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! SUPERCONDUCTING PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_superc_eloc_1b(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(2,Nlat,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Lmats)
    !
    complex(8)                                  :: Gloc_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    if(Norb>1)stop "ed_get_weiss_field_superc_eloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_superc_eloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_superc_eloc(Gloc_,Smats_,Weiss_,Eloc_)
    Gloc(:,:,:) = Gloc_(:,:,1,1,1,1,:)
    Smats(:,:,:) = Smats_(:,:,1,1,1,1,:)
    Weiss(:,:,:) = Weiss_(:,:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_eloc_1b

  subroutine ed_get_weiss_field_superc_eloc_mb(Gloc,Smats,Weiss,Eloc)
    complex(8)                                  :: Gloc(2,Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Gloc_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Nlat*Norb*Nspin)
    real(8)                                     :: Eloc_(Nlat*Norb*Nspin)
    if(Nspin>1)stop "ed_get_weiss_field_superc_eloc_Mb error: Nspin > 1 in M-band routine" 
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_superc_eloc(Gloc_,Smats_,Weiss_,Eloc_)
    Gloc(:,:,:,:,:) = Gloc_(:,:,1,1,:,:,:)
    Smats(:,:,:,:,:) = Smats_(:,:,1,1,:,:,:)
    Weiss(:,:,:,:,:) = Weiss_(:,:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_eloc_mb

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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_superc_eloc


  subroutine ed_get_weiss_field_superc_hloc_1b(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(2,Nlat,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Lmats)
    !
    complex(8)                                  :: Gloc_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    if(Norb>1)stop "ed_get_weiss_field_superc_hloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_superc_hloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call ed_get_weiss_field_superc_hloc(Gloc_,Smats_,Weiss_,Hloc)
    Gloc(:,:,:) = Gloc_(:,:,1,1,1,1,:)
    Smats(:,:,:) = Smats_(:,:,1,1,1,1,:)
    Weiss(:,:,:) = Weiss_(:,:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_hloc_1b

  subroutine ed_get_weiss_field_superc_hloc_mb(Gloc,Smats,Weiss,Hloc)
    complex(8)                                  :: Gloc(2,Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nlat,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nlat,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Gloc_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats)
    !
    complex(8)                                  :: Hloc(Nlat,Nspin,Nspin,Norb,Norb)
    if(Nspin>1)stop "ed_get_weiss_field_superc_hloc_mb error: Nspin > 1 in M-band routine" 
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call ed_get_weiss_field_superc_hloc(Gloc_,Smats_,Weiss_,Hloc)
    Gloc(:,:,:,:,:) = Gloc_(:,:,1,1,:,:,:)
    Smats(:,:,:,:,:) = Smats_(:,:,1,1,:,:,:)
    Weiss(:,:,:,:,:) = Weiss_(:,:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_hloc_mb

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
#ifdef _MPI_INEQ
    call MPI_ALLREDUCE(Weiss_tmp,Weiss,size(Weiss),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
#else
    Weiss = Weiss_tmp
#endif
  end subroutine ed_get_weiss_field_superc_hloc


end module ED_WRAP_WEISS
