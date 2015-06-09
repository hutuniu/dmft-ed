module ED_WEISS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_IOTOOLS,   only:reg,txtfy,splot
  USE SF_TIMER
  USE SF_ARRAYS,    only: arange
  USE SF_LINALG,    only:matrix_inverse,matrix_inverse_sym
  implicit none
  private

  interface ed_get_weiss
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
  end interface ed_get_weiss

  public :: ed_get_weiss


  real(8),dimension(:),allocatable :: wm
  character(len=20)                :: suffix

contains


  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! NORMAL PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_normal_eloc_1b(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(Lmats)
    complex(8)                                  :: Smats(Lmats)
    complex(8)                                  :: Weiss(Lmats)
    !
    complex(8)                                  :: Gloc_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nspin,Nspin,Norb,Norb,Lmats)
    !
    integer                                     :: iprint
    real(8),optional                            :: Eloc(Norb*Nspin)
    real(8)                                     :: Eloc_(Norb*Nspin)
    if(Norb>1)stop "ed_get_weiss_field_normal_eloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_normal_eloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(1,1,1,1,:) = Gloc(:)
    Smats_(1,1,1,1,:) = Smats(:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_normal_eloc(Gloc_,Smats_,Weiss_,iprint,Eloc_)
    Gloc(:) = Gloc_(1,1,1,1,:)
    Smats(:) = Smats_(1,1,1,1,:)
    Weiss(:) = Weiss_(1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_eloc_1b

  subroutine ed_get_weiss_field_normal_eloc_mb(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Norb,Norb,Lmats)
    !
    complex(8)                                  :: Gloc_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nspin,Nspin,Norb,Norb,Lmats)
    !
    integer                                     :: iprint
    real(8),optional                            :: Eloc(Norb*Nspin)
    real(8)                                     :: Eloc_(Norb*Nspin)
    if(Nspin>1)stop "ed_get_weiss_field_normal_eloc_1m error: Nspin > 1 in M-band routine" 
    Gloc_(1,1,:,:,:) = Gloc(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_normal_eloc(Gloc_,Smats_,Weiss_,iprint,Eloc_)
    Gloc(:,:,:) = Gloc_(1,1,:,:,:)
    Smats(:,:,:) = Smats_(1,1,:,:,:)
    Weiss(:,:,:) = Weiss_(1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_eloc_mb

  subroutine ed_get_weiss_field_normal_eloc(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nspin,Nspin,Norb,Norb,Lmats)
    integer                                     :: iprint
    real(8),optional                            :: Eloc(Norb*Nspin)
    !aux
    real(8)                                     :: Eloc_(Norb*Nspin)
    complex(8)                                  :: zeta_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(Nspin*Norb,Nspin*Norb,Lmats)
    integer                                     :: i,iorb,jorb,ispin,jspin,io,jo,js
    !
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix
    !and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          js = iorb + (ispin-1)*Norb
          zeta_site(io,io,:) = xi*wm(:) + xmu - Eloc_(io)
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)    = zeta_site(io,jo,:) - Smats(ispin,jspin,iorb,jorb,:)
                invGloc_site(io,jo,:) = Gloc(ispin,jspin,iorb,jorb,:)
                Smats_site(io,jo,:)   = Smats(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    !
    !Invert the [Norb*Nspin]**2 Gloc block matrix 
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
                Weiss(ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
             enddo
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       select case(iprint)
       case(1)
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                if(cg_scheme=="weiss")then
                   call splot("WeissField"//reg(suffix),wm,Weiss(ispin,ispin,iorb,iorb,:))
                else
                   call splot("Delta"//reg(suffix),wm,Weiss(ispin,ispin,iorb,iorb,:))
                endif
             enddo
          enddo
       case(2)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   if(cg_scheme=="weiss")then
                      call splot("WeissField"//reg(suffix),wm,Weiss(ispin,ispin,iorb,jorb,:))
                   else
                      call splot("Delta"//reg(suffix),wm,Weiss(ispin,ispin,iorb,jorb,:))
                   endif
                enddo
             enddo
          enddo
       case default
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                      if(cg_scheme=="weiss")then
                         call splot("WeissField"//reg(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                      else
                         call splot("Delta"//reg(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                      endif
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_weiss_field_normal_eloc


  subroutine ed_get_weiss_field_normal_hloc_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(Lmats)
    complex(8)                                  :: Smats(Lmats)
    complex(8)                                  :: Weiss(Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nspin,Nspin,Norb,Norb,Lmats)
    if(Norb>1)stop "ed_get_weiss_field_normal_hloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_normal_hloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(1,1,1,1,:) = Gloc(:)
    Smats_(1,1,1,1,:) = Smats(:)
    call ed_get_weiss_field_normal_hloc(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Gloc(:) = Gloc_(1,1,1,1,:)
    Smats(:) = Smats_(1,1,1,1,:)
    Weiss(:) = Weiss_(1,1,1,1,:)
  end subroutine ed_get_weiss_field_normal_hloc_1b

  subroutine ed_get_weiss_field_normal_hloc_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(Nspin,Nspin,Norb,Norb,Lmats)
    if(Nspin>1)stop "ed_get_weiss_field_normal_hloc_mb error: Nspin > 1 in 1-band routine" 
    Gloc_(1,1,:,:,:) = Gloc(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    call ed_get_weiss_field_normal_hloc(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Gloc(:,:,:) = Gloc_(1,1,:,:,:)
    Smats(:,:,:) = Smats_(1,1,:,:,:)
    Weiss(:,:,:) = Weiss_(1,1,:,:,:)
  end subroutine ed_get_weiss_field_normal_hloc_mb

  subroutine ed_get_weiss_field_normal_hloc(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !aux
    complex(8)                                  :: zeta_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(Nspin*Norb,Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(Nspin*Norb,Nspin*Norb,Lmats)
    integer                                     :: i,iorb,jorb,ispin,jspin,io,jo
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix
    !and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_site(io,io,:) = xi*wm(:) + xmu
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)    = zeta_site(io,jo,:) - Hloc(ispin,jspin,iorb,jorb) - Smats(ispin,jspin,iorb,jorb,:)
                invGloc_site(io,jo,:) = Gloc(ispin,jspin,iorb,jorb,:)
                Smats_site(io,jo,:)   = Smats(ispin,jspin,iorb,jorb,:)
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
    !output structure of [Nspsin,Nspin,Norb,Norb] matrix
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Weiss(ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
             enddo
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       select case(iprint)
       case(1)
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                if(cg_scheme=="weiss")then
                   call splot("WeissField"//reg(suffix),wm,Weiss(ispin,ispin,iorb,iorb,:))
                else
                   call splot("Delta"//reg(suffix),wm,Weiss(ispin,ispin,iorb,iorb,:))
                endif
             enddo
          enddo
       case(2)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   if(cg_scheme=="weiss")then
                      call splot("WeissField"//reg(suffix),wm,Weiss(ispin,ispin,iorb,jorb,:))
                   else
                      call splot("Delta"//reg(suffix),wm,Weiss(ispin,ispin,iorb,jorb,:))
                   endif
                enddo
             enddo
          enddo
       case default
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                      if(cg_scheme=="weiss")then
                         call splot("WeissField"//reg(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                      else
                         call splot("Delta"//reg(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                      endif
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_weiss_field_normal_hloc



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Get the local Weiss Field calG0 or Hybridization function \Delta using 
  ! self-consistency equations and given G_loc and Sigma.
  ! SUPERCONDUCTING PHASE
  !-------------------------------------------------------------------------------------------
  subroutine ed_get_weiss_field_superc_eloc_1b(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(2,Lmats)
    complex(8)                                  :: Smats(2,Lmats)
    complex(8)                                  :: Weiss(2,Lmats)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Norb*Nspin)
    real(8)                                     :: Eloc_(Norb*Nspin)
    if(Norb>1)stop "ed_get_weiss_field_superc_eloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_superc_eloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_superc_eloc(Gloc_,Smats_,Weiss_,iprint,Eloc_)
    Gloc(:,:) = Gloc_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_eloc_1b

  subroutine ed_get_weiss_field_superc_eloc_mb(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(2,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Norb,Norb,Lmats)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nspin,Nspin,Norb,Norb,Lmats)
    !
    real(8),optional                            :: Eloc(Norb*Nspin)
    real(8)                                     :: Eloc_(Norb*Nspin)
    if(Nspin>1)stop "ed_get_weiss_field_superc_eloc_Mb error: Nspin > 1 in M-band routine" 
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    call ed_get_weiss_field_superc_eloc(Gloc_,Smats_,Weiss_,iprint,Eloc_)
    Gloc(:,:,:,:) = Gloc_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_eloc_mb

  subroutine ed_get_weiss_field_superc_eloc(Gloc,Smats,Weiss,iprint,Eloc)
    complex(8)                                  :: Gloc(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nspin,Nspin,Norb,Norb,Lmats)
    integer                                     :: iprint
    real(8),optional                            :: Eloc(Norb*Nspin)
    !aux
    real(8)                                     :: Eloc_(Norb*Nspin)
    complex(8)                                  :: zeta_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    integer                                     :: i,iorb,jorb,ispin,jspin,io,jo,js
    integer                                     :: Nso
    !
    Eloc_=0d0       ;if(present(Eloc))Eloc_=Eloc
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Nso =Nspin*Norb
    !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
    !and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          js = iorb + (ispin-1)*Norb
          zeta_site(io,io,:)         = xi*wm(:) + xmu - Eloc_(io)
          zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Eloc_(io)
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ispin,jspin,iorb,jorb,:)
                zeta_site(io,jo+Nso,:)       =-Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo,:)       =-Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ispin,jspin,iorb,jorb,:))
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
    !output structure of [Nspsin,Nspin,Norb,Norb] matrix
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
       select case(iprint)
       case(1)
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                if(cg_scheme=="weiss")then
                   call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,iorb,:))
                   call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,iorb,:))
                else
                   call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,iorb,:))
                   call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,iorb,:))
                endif
             enddo
          enddo
       case(2)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   if(cg_scheme=="weiss")then
                      call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,jorb,:))
                      call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,jorb,:))
                   else
                      call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,jorb,:))
                      call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,jorb,:))
                   endif
                enddo
             enddo
          enddo
       case default
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                      if(cg_scheme=="weiss")then
                         call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,jspin,iorb,jorb,:))
                         call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,jspin,iorb,jorb,:))
                      else
                         call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,jspin,iorb,jorb,:))
                         call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,jspin,iorb,jorb,:))
                      endif
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_weiss_field_superc_eloc


  subroutine ed_get_weiss_field_superc_hloc_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(2,Lmats)
    complex(8)                                  :: Smats(2,Lmats)
    complex(8)                                  :: Weiss(2,Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nspin,Nspin,Norb,Norb,Lmats)
    !
    if(Norb>1)stop "ed_get_weiss_field_superc_hloc_1b error: Norb > 1 in 1-band routine" 
    if(Nspin>1)stop "ed_get_weiss_field_superc_hloc_1b error: Nspin > 1 in 1-band routine" 
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call ed_get_weiss_field_superc_hloc(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Gloc(:,:) = Gloc_(:,1,1,1,1,:)
    Smats(:,:) = Smats_(:,1,1,1,1,:)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine ed_get_weiss_field_superc_hloc_1b

  subroutine ed_get_weiss_field_superc_hloc_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(2,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !
    complex(8)                                  :: Gloc_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats_(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss_(2,Nspin,Nspin,Norb,Norb,Lmats)
    !
    if(Nspin>1)stop "ed_get_weiss_field_superc_hloc_mb error: Nspin > 1 in M-band routine" 
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call ed_get_weiss_field_superc_hloc(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Gloc(:,:,:,:) = Gloc_(:,1,1,:,:,:)
    Smats(:,:,:,:) = Smats_(:,1,1,:,:,:)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine ed_get_weiss_field_superc_hloc_mb

  subroutine ed_get_weiss_field_superc_hloc(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8)                                  :: Gloc(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Smats(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Weiss(2,Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)                                  :: Hloc(Nspin,Nspin,Norb,Norb)
    integer                                     :: iprint
    !aux
    complex(8)                                  :: zeta_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: Smats_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: invGloc_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    complex(8)                                  :: calG0_site(2*Nspin*Norb,2*Nspin*Norb,Lmats)
    integer                                     :: i,iorb,jorb,ispin,jspin,io,jo
    integer                                     :: Nso
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Nso =Nspin*Norb
    !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix
    !and create the zeta_site
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
                Weiss(1,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                Weiss(2,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
             enddo
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0.AND.ed_verbose<4)then
       select case(iprint)
       case(1)
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                if(cg_scheme=="weiss")then
                   call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,iorb,:))
                   call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,iorb,:))
                else
                   call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,iorb,:))
                   call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,iorb,:))
                endif
             enddo
          enddo
       case(2)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
                   if(cg_scheme=="weiss")then
                      call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,jorb,:))
                      call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,jorb,:))
                   else
                      call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,ispin,iorb,jorb,:))
                      call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,ispin,iorb,jorb,:))
                   endif
                enddo
             enddo
          enddo
       case default
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                      if(cg_scheme=="weiss")then
                         call splot("WeissField_normal"//reg(suffix),wm,Weiss(1,ispin,jspin,iorb,jorb,:))
                         call splot("WeissField_anomal"//reg(suffix),wm,Weiss(2,ispin,jspin,iorb,jorb,:))
                      else
                         call splot("Delta_normal"//reg(suffix),wm,Weiss(1,ispin,jspin,iorb,jorb,:))
                         call splot("Delta_anomal"//reg(suffix),wm,Weiss(2,ispin,jspin,iorb,jorb,:))
                      endif
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine ed_get_weiss_field_superc_hloc



end module ED_WEISS
