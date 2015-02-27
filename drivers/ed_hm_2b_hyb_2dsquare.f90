!                             MODEL Hamiltonian is:
!
!  | Mh - e0*[cos(kx)+cos(ky)]  ,       - lambda*[cos(kx)+cos(ky)]  |
!  | - lambda*[cos(kx)+cos(ky)] ,       -Mh - e0*[cos(kx)+cos(ky)] |
!
!
program ed_hm_2bands_nohyb_2Dsquare
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS


  implicit none
  integer                :: iloop
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  complex(8),allocatable :: Hloc(:,:)
  !variables for the model:
  real(8)                :: e0,Mh,lambda
  integer                :: Lk,Nkstep,Nkx,Nky
  real(8)                :: wmixing
  character(len=32)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym
  
  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_HM2B.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(Nkx,"NKX",finput,default=20)
  call parse_input_variable(Nky,"NKY",finput,default=20)
  call parse_input_variable(Nkstep,"NKSTEP",finput,default=100)
  call parse_input_variable(e0,"E0",finput,default=1.d0)
  call parse_input_variable(lambda ,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(Mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call ed_read_input(trim(finput))

  if(Norb/=2)stop "Wrong setup from input file: Norb=2"
  if(Nspin>1)stop "Wrong setup from input file: Nspin>1"

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))


  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))


  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb(1),Nb(2)))
  allocate(Bath_(Nb(1),Nb(2)))
  call ed_init_solver(bath)	!!!!!!!!!! F90 is NOT case-sensitive !!!!!!!!!!
  call set_Hloc(Hloc,1)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)	!!!!!!!!!! F90 is NOT case-sensitive !!!!!!!!!!

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta
     !Fit the new bath, starting from the old bath + the supplied delta
     call ed_chi2_fitgf(delta(:,:,:),bath,ispin=1)	!!!!!!!!!! F90 is NOT case-sensitive !!!!!!!!!!

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  !Get Kinetic Energy too
  call ed_kinetic_energy(Hk,dos_wt,impSmats(1,1,:,:,:))

contains




  !---------------------------------------------------------------------
  !PURPOSE: GET DELTA FUNCTION
  !---------------------------------------------------------------------
  subroutine get_delta
    integer                                     :: i,j,ik,iorb,jorb,ispin,jspin,iso,unit
    complex(8),dimension(Norb,Norb)             :: zeta,fg,gdelta
    complex(8),dimension(:,:,:),allocatable     :: Smats
    complex(8),dimension(:,:,:),allocatable     :: gloc
    complex(8)                                  :: iw
    real(8)                                     :: wm(Lmats),wr(Lreal)
    character(len=20)                           :: suffix
    !
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta = zero
    !
    !MATSUBARA AXIS		!!!!!!!!!! OBTAIN G_LOC FOR MATSUBARA FREQUENCIES !!!!!!!!!!
    print*,"Get Gloc_iw:"
    allocate(gloc(Norb,Norb,Lmats))
    do i=1,Lmats		!!!!!!!!!! BEGIN LOOP ALL OVER MATSUBARA FREQUENCIES !!!!!!!!!!
       iw = xi*wm(i)
       forall(iorb=1:Norb)zeta(iorb,iorb)=iw+xmu     ! diagonal elements of z are iw_n + \mu
       zeta(:,:) = zeta(:,:) - impSmats(1,1,:,:,i)   ! z is (iw_n + \mu)
       fg=zero		!!!!! dummy variable for gloc at each mats.freq. !!!!!!!!!!
       do ik=1,Lk       !!!!! sum over k of the inverse of [z-h_k] renormalized by number of k-points !!!!!!!!!!
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo		
       gloc(:,:,i) = fg	!!!!! the value stored in the dummy variable fg is assigned to gloc for the current mats.freq. !!!!!!!!!!
       !
       !Get Delta=\Delta or G_0 i.e. hybridization function or weiss field !!!!!!!!!!
       if(lambda==0.d0)then	!!!!!!!!!! diagonal case !!!!!!!!!!
          forall(iorb=1:Norb)fg(iorb,iorb)=one/fg(iorb,iorb)	!!!!!!!!!! both in G_0 and Delta gloc enters as its inverse so fg becomes the container of (G_loc)^-1 !!!!!!!!!!
          if(cg_scheme=='weiss')then !!!!!!!!!! case: gdelta = G_0 !!!!!!!!!!
             gdelta=zero
             forall(iorb=1:Norb)gdelta(iorb,iorb) = one/(fg(iorb,iorb) + impSmats(1,1,iorb,iorb,i))
          else !!!!!!!!!! case: gdelta = Delta (hybridization function)
             gdelta=zero
             forall(iorb=1:Norb)gdelta(iorb,iorb) = zeta(iorb,iorb) - Hloc(iorb,iorb) - fg(iorb,iorb)
          endif
       else	!!!!!!!!!! off-diagonal case
          call matrix_inverse(fg)
          if(cg_scheme=='weiss')then	!!!!!!!!!! case: gdelta = G_0 !!!!!!!!!!
             gdelta = fg + impSmats(1,1,:,:,i)
             call matrix_inverse(gdelta)
          else !!!!!!!!!! case: gdelta = Delta (hybridization function) !!!!!!!!!!
             gdelta = zeta - Hloc - fg
          endif
       endif
       !
       delta(:,:,i) = gdelta	!!!!!!!!!! gdelta is either G_0 or Delta depending on the variable cg_scheme
       !
    enddo			!!!!!!!!!! END LOOP ALL OVER MATSUBARA FREQUENCIES !!!!!!!!!!
    do iorb=1,Norb		!!!!!!!!!! print for all mats.freq. G_loc and G_0 (or Delta) !!!!!!!!!!
       suffix="_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed"	!!!!!!!!!! suffixes l1m1 ecc for the two bands !!!!!!!!!!
       call splot("Delta"//reg(suffix),wm,delta(iorb,iorb,:))	!!!!!!!!!! Delta_loc_l iorb m iorb VS mast.freq. !!!!!!!!!!
       call splot("Gloc"//reg(suffix),wm,gloc(iorb,iorb,:))	!!!!!!!!!! G_loc_l iorb m iorb VS mast.freq.  !!!!!!!!!!
    enddo        
    deallocate(gloc)	!!!!!!!!!! deallocate gloc because it'll be used for store gloc in real freq. !!!!!!!!!!
    !
    !
    !REAL AXIS		!!!!!!!!!! OBTAIN G_LOC FOR REAL FREQUENCIES !!!!!!!!!!
    allocate(gloc(Norb,Norb,Lreal))		!!!!!!!!!! same steps as for FOR MATSUBARA FREQUENCIES !!!!!!!!!!
    print*,"Get Gloc_realw:"
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       forall(iorb=1:Norb)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - impSreal(1,1,:,:,i)
       fg=zero
       do ik=1,Lk
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
    enddo
    do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed"
          call splot("Gloc"//reg(suffix),wr,-dimag(gloc(iorb,iorb,:))/pi,dreal(gloc(iorb,iorb,:)))   
    enddo
  deallocate(gloc)
  end subroutine get_delta






  !---------------------------------------------------------------------
      !PURPOSE: CONSTRUCT MODEL HAMILTONIAN (from the NonInteracting code)
      !---------------------------------------------------------------------
      subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky 
    integer                             :: iorb,jorb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Norb,Norb,Lmats) :: Gmats
    complex(8),dimension(Norb,Norb,Lreal) :: Greal
    real(8)                               :: wm(Lmats),wr(Lreal),dw,n0(Norb)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)    !!!!!!!! taking odd pi/beta multiples , real(x,8) is a cast with precision 8 arange is like linspace but on complex numbers  !!!!!!!!
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    !
    write(LOGfile,*)"Build model H(k) for HM 2bands without hybr. in 2D square lattice"
    !
    Lk=Nkx*Nky
    !
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Norb,Norb,Lk))	!!!!!!!! Hk is 2x2 x #_kpoints !!!!!!!!
    allocate(dos_wt(Lk))
    if(present(file))then	!!!!!!!!  ???????? !!!!!!!!
       unit=free_unit()
       open(unit,file=file)
       Gmats=zero			!!!!!!!! set to zero Gmats which is 2x2 x #_mats.freq.  !!!!!!!!
       Greal=zero			!!!!!!!! set to zero Greal which is 2x2 x #_realfreq.  !!!!!!!!
    endif
    ik=0				!!!!!!!! Hk is not stored respect to (kx,ky) but respect to ik which runs from 1 to Lk=Nkx*Nky  !!!!!!!!
    call start_timer
    do ix=1,Nkx
       kx = -pi + 2.d0*pi*(ix-1)/Nkx	!!!!!!!! kx runs from -pi to +pi  !!!!!!!!
       do iy=1,Nky
          ky = -pi + 2.d0*pi*(iy-1)/Nky	!!!!!!!! ky runs from -pi to +pi  !!!!!!!!
            !
             ik=ik+1			!!!!!!!! increase ik  !!!!!!!!
             Hk(:,:,ik) = Hk_model(kx,ky)	!!!!!!!! assign ik-th value to the four components of Hk  !!!!!!!!
             !
             if(present(file))then		!!!!!!!! Only the first processor and ???????? !!!!!!!!
                write(unit,"(3(F10.7,1x))")kx,ky
                do iorb=1,Norb
                   write(unit,"(100(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Norb)		!!!!!!!! ???????? !!!!!!!!
                enddo
                !		!!!!!!!! Compute Greal and Gmats non interacting (Sigma=0) to store it because they are the initial ones !!!!!!!!
                do i=1,Lreal	!!!!!!!! it's like in get_delta but here the inverse of z+h is summed directly to G for each kpoint  !!!!!!!!
                   Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk	!!!!!!!! look at inverse_g0k !!!!!!!!
                enddo
                do i=1,Lmats
                   Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(xi*wm(i)+xmu,Hk(:,:,ik))/Lk
                enddo
             endif
    
             call eta(ik,Lk)
       enddo
    enddo
    call stop_timer	!!!!!!!! ???????? !!!!!!!!
    if(present(file))then	!!!!!!!! ???????? !!!!!!!!
       write(unit,*)""			!!!!!!!! ???????? !!!!!!!!
       close(unit)				!!!!!!!! ???????? !!!!!!!!
    endif
    !
    dos_wt=1.d0/Lk
    !
    if(present(file))then	!!!!!!!!  ???????? !!!!!!!!
       do iorb=1,Norb
          n0(iorb)=fft_get_density(Gmats(iorb,iorb,:),beta)	!!!!!!!! Get dos for both orbitals !!!!!!!!
          call splot("U0_Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))	!!!!!!!! write dos's in files !!!!!!!!
          call splot("U0_Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,-dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
	  !call splot("DOS_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,-dimag(Greal(iorb,iorb,:))/pi)
       enddo
       unit=free_unit()
       open(unit,file="U0_observables.ed")
       write(unit,"(24F20.12)")e0,lambda,Mh,xmu,(n0(i),i=1,Norb),sum(n0)
       close(unit)
       write(LOGfile,"(24F20.12)")e0,lambda,Mh,xmu,(n0(i),i=1,Norb),sum(n0)
    endif
    !
    allocate(Hloc(Norb,Norb))
    Hloc = sum(Hk(:,:,:),dim=3)/Lk	!!!!!!!! each of the four element H_loc_(i,j)=(#_kpoints)^{-1}sum_k H_k(i,j) !!!!!!!!
    where(abs(dreal(Hloc))<1.d-9)Hloc=0.d0
    !
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of bands        :",Norb
      end subroutine build_hk
    
    
    
  !---------------------------------------------------------------------
  !PURPOSE: GET THE model HAMILTONIAN ALONG THE Gamma-X-M-Gamma path
  !---------------------------------------------------------------------
  subroutine build_hk_GXMG()
    integer                                     :: i,j,ik=0
    integer                                     :: ix,iy
    integer                                     :: iorb,jorb
    integer                                     :: isporb,jsporb
    integer                                     :: ispin,jspin
    integer                                     :: iso,unit
    real(8)                                     :: foo
    real(8)                                     :: kx,ky    
    complex(8),dimension(Norb,Norb)               :: zeta,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal,Smats
    complex(8),dimension(:,:,:,:),allocatable   :: gk,gfoo,ReSmat
    complex(8),dimension(:,:,:),allocatable     :: Hktilde
    real(8)                                     :: eig(Norb)
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is used in get_delta with getak=T
    if(ed_mpi_id==0)write(LOGfile,*)"Build H(k) model along the path GXMG:"
    unit=free_unit() 
    if(ed_mpi_id==0)open(unit,file="Eigenbands.dat")
    Lk=3*Nkstep
    ik = 0
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Norb,Norb,Lk))
    !From \Gamma=(0,0) to X=(pi,0): Nkstep
    do ix=1,Nkstep
       ik=ik+1
       kx = 0.d0 + pi*real(ix-1,8)/dble(Nkstep)
       ky = 0.d0
       Hk(:,:,ik)=hk_model(kx,ky)
       eig = Eigk(hk_model(kx,ky))
       if(ed_mpi_id==0)write(unit,"(I3,16F25.12)")ik,(eig(i),i=1,Norb)
    enddo
    !From X=(pi,0) to M=(pi,pi): Nk steps
    do iy=1,Nkstep
       ik=ik+1
       kx = pi
       ky = 0.d0 + pi*real(iy-1,8)/dble(Nkstep)
       Hk(:,:,ik)=hk_model(kx,ky)
       eig = Eigk(hk_model(kx,ky))
       if(ed_mpi_id==0)write(unit,"(I3,16F25.12)")ik,(eig(i),i=1,Norb)
    enddo
    !From M=(pi,pi) to \Gamma=(0,0): Nk steps
    do ix=1,Nkstep
       ik=ik+1
       iy=ix
       kx = pi - pi*real(ix-1,8)/dble(Nkstep)
       ky = pi - pi*real(iy-1,8)/dble(Nkstep)
       Hk(:,:,ik)=hk_model(kx,ky)
       eig = Eigk(hk_model(kx,ky))
       if(ed_mpi_id==0)write(unit,"(I3,16F25.12)")ik,(eig(i),i=1,Norb)
    enddo
    if(ed_mpi_id==0)close(unit)
  end subroutine build_hk_GXMG


  function Eigk(hk) result(eig)
    complex(8),dimension(2,2) :: hk
    real(8),dimension(2)      :: eig
    !call matrix_diagonalize(hk,eig)
    eig(1)=hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig(2)=hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig = eig/2.d0
  end function Eigk
    
    
    
    
      !--------------------------------------------------------------------!
      !MODEL HAMILTONIAN H(k)
      !--------------------------------------------------------------------!
      function Hk_model(kx,ky) result(Hk)
    real(8)                         :: kx,ky
    real(8)                         :: epsik,vpsik
    complex(8),dimension(Norb,Norb) :: Hk
    epsik = cos(kx)+cos(ky)
    vpsik = cos(kx)-cos(ky)
    Hk(1,1) = Mh  - e0*epsik
    Hk(2,2) = -Mh - e0*epsik
    Hk(1,2) = - lambda*vpsik
    Hk(2,1) = - lambda*vpsik
      end function Hk_model



  !--------------------------------------------------------------------!
  !ADDITIONAL ROUTINES:
  !--------------------------------------------------------------------!
  function inverse_gk(zeta,hk) result(gk)
    complex(8),dimension(Norb,Norb) :: zeta,hk
    complex(8),dimension(Norb,Norb) :: gk
    complex(8)                      :: h11,h22,h12
    integer                         :: iorb
    select case(lambda==0.d0)
    case default                !non-diagonal Hamiltonian
       h11 = zeta(1,1) - hk(1,1)
       h22 = zeta(2,2) - hk(2,2)
       !h12 = zeta(1,2) - hk(1,2)   !we should assume that z(1,2)==0.0
       h12 =           - hk(1,2)
       gk(1,1) =  one/(h11 - abs(h12)**2/h22)
       gk(2,2) =  one/(h22 - abs(h12)**2/h11)
       gk(1,2) = -h12/(h11*h22 - abs(h12)**2)
       gk(2,1) = conjg(gk(1,2))
    case (.true.)               !diagonal Hamiltonian
       gk = zero
       forall(iorb=1:Norb)gk(iorb,iorb)=one/(zeta(iorb,iorb)-hk(iorb,iorb))
    end select
  end function inverse_gk
  !
  function inverse_g0k(zeta,hk) result(gk)	!!!!!!!!!! the difference with inverse_gk is that I'm taking zeta as diagonal !!!!!!!!!!
    complex(8)                      :: zeta
    complex(8),dimension(Norb,Norb) :: hk
    complex(8),dimension(Norb,Norb) :: gk
    complex(8)                      :: h11,h22,h12
    integer                         :: iorb
    select case(lambda==0.d0)
    case default                !non-diagonal Hamiltonian
       h11 = zeta - hk(1,1)
       h22 = zeta - hk(2,2)
       h12 =      - hk(1,2)
       gk(1,1) =  one/(h11 - abs(h12)**2/h22)
       gk(2,2) =  one/(h22 - abs(h12)**2/h11)
       gk(1,2) = -h12/(h11*h22 - abs(h12)**2)
       gk(2,1) = conjg(gk(1,2))
    case (.true.)               !diagonal Hamiltonian
       gk = zero
       h11 = zeta - hk(1,1)
       h22 = zeta - hk(2,2)
       gk(1,1) =  one/h11
       gk(2,2) =  one/h22
    end select
  end function inverse_g0k




  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do iorb=1,Norb
       unit=free_unit()
       if(L==Lreal)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw.ed"
       elseif(L==Lmats)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed"
       endif
       write(*,*)"read from file=","impSigma"//reg(suffix)
       open(unit,file="impSigma"//reg(suffix),status='old')
       do i=1,L
          read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
          forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
       enddo
       close(unit)
    enddo
  end subroutine read_sigma

end program ed_hm_2bands_nohyb_2Dsquare



