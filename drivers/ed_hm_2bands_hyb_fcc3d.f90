!                             MODEL Hamiltonian is:
!
!  | Mh - e0*[cos(kx)+cos(ky)+cos(kz)]  ,        lambda*(cos(kx)-cos(ky))*cos(kz)  |
!  | lambda*(cos(kx)-cos(ky))*cos(kz)   ,       -Mh - e0*[cos(kx)+cos(ky)+cos(kz)] |
!
!
program ed_hm2bhyb_fcc
  USE DMFT_ED
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif

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
  real(8)                :: e0,lambda,Mh
  integer                :: Lk,Nk,Nkstep,Nkx,Nky,Nkz
  real(8)                :: wmixing
  character(len=32)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym

#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_HM2B.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(Nkx,"NKX",finput,default=20)
  call parse_input_variable(Nky,"NKY",finput,default=20)
  call parse_input_variable(Nkz,"NKZ",finput,default=20)
  call parse_input_variable(Nkstep,"NKSTEP",finput,default=100)
  call parse_input_variable(e0,"E0",finput,default=1.d0)
  call parse_input_variable(lambda ,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(Mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  !
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
  call init_ed_solver(bath)
  call set_Hloc(Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta
     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:,:,:),bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath

     if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0)call end_loop
  enddo

  !Get Kinetic Energy too
  if(ED_MPI_ID==0)call ed_kinetic_energy(impSmats(1,1,:,:,:),Hk,dos_wt)

#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif
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
    !MATSUBARA AXIS
    if(ED_MPI_ID==0)print*,"Get Gloc_iw:"
    allocate(gloc(Norb,Norb,Lmats))
    do i=1,Lmats
       iw = xi*wm(i)
       forall(iorb=1:Norb)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - impSmats(1,1,:,:,i)
       fg=zero
       do ik=1,Lk
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
       !
       !Get Delta=\Delta or G_0
       if(lambda==0.d0)then
          forall(iorb=1:Norb)fg(iorb,iorb)=one/fg(iorb,iorb)
          if(cg_scheme=='weiss')then
             gdelta=zero
             forall(iorb=1:Norb)gdelta(iorb,iorb) = one/(fg(iorb,iorb) + impSmats(1,1,iorb,iorb,i))
          else
             gdelta=zero
             forall(iorb=1:Norb)gdelta(iorb,iorb) = zeta(iorb,iorb) - Hloc(iorb,iorb) - fg(iorb,iorb)
          endif
       else
          call matrix_inverse(fg)
          if(cg_scheme=='weiss')then
             gdelta = fg + impSmats(1,1,:,:,i)
             call matrix_inverse(gdelta)
          else
             gdelta = zeta - Hloc - fg
          endif
       endif
       !
       delta(:,:,i) = gdelta
       !
    enddo
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed"
          call splot("Delta"//reg(suffix),wm,delta(iorb,iorb,:))
          call splot("Gloc"//reg(suffix),wm,gloc(iorb,iorb,:))
       enddo
    endif
    deallocate(gloc)
    !
    !
    !REAL AXIS
    allocate(gloc(Norb,Norb,Lreal))
    if(ED_MPI_ID==0)print*,"Get Gloc_realw:"
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
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed"
          call splot("Gloc"//reg(suffix),wr,-dimag(gloc(iorb,iorb,:))/pi,dreal(gloc(iorb,iorb,:)))
       enddo
    endif
    deallocate(gloc)
  end subroutine get_delta









  !---------------------------------------------------------------------
  !PURPOSE: CONSTRUCT MODEL HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy,iz
    real(8)                             :: kx,ky,kz 
    integer                             :: iorb,jorb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Norb,Norb,Lmats) :: Gmats
    complex(8),dimension(Norb,Norb,Lreal) :: Greal
    real(8)                               :: wm(Lmats),wr(Lreal),dw,n0(Norb)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    !
    if(ED_MPI_ID==0)write(LOGfile,*)"Build model H(k) for HM 2bands w/ hyb in 3D fcc lattice:"
    !
    Lk=Nkx*Nky*Nkz
    !
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Norb,Norb,Lk))
    allocate(dos_wt(Lk))
    if(ED_MPI_ID==0.AND.present(file))then
       unit=free_unit()
       open(unit,file=file)
       Gmats=zero
       Greal=zero
       if(ED_MPI_ID==0)call start_progress
    endif
    ik=0
    if(ED_MPI_ID==0)call start_progress
    do ix=1,Nkx
       kx = -pi + 2.d0*pi*(ix-1)/Nkx
       do iy=1,Nky
          ky = -pi + 2.d0*pi*(iy-1)/Nky
          do iz=1,Nkz
             kz = -pi + 2.d0*pi*(iz-1)/Nkz
             !
             ik=ik+1
             Hk(:,:,ik) = hk_model(kx,ky,kz)
             !
             if(ED_MPI_ID==0.AND.present(file))then
                write(unit,"(3(F10.7,1x))")kx,ky,kz
                do iorb=1,Norb
                   write(unit,"(100(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Norb)
                enddo
                !
                do i=1,Lreal
                   Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
                enddo
                do i=1,Lmats
                   Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(xi*wm(i)+xmu,Hk(:,:,ik))/Lk
                enddo
             endif

             if(ED_MPI_ID==0)call progress_bar(ik,Lk)
          enddo
       enddo
    enddo
    if(ED_MPI_ID==0)call stop_progress
    if(ED_MPI_ID==0.AND.present(file))then
       write(unit,*)""
       close(unit)
    endif
    !
    dos_wt=1.d0/Lk
    !
    if(ED_MPI_ID==0.AND.present(file))then
       do iorb=1,Norb
          n0(iorb)=get_density_fromFFT(Gmats(iorb,iorb,:),beta)
          call splot("U0_Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
          call splot("U0_Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,-dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
       enddo
       unit=free_unit()
       open(unit,file="U0_observables.ed")
       write(unit,"(24F20.12)")e0,lambda,Mh,xmu,(n0(i),i=1,Norb),sum(n0)
       close(unit)
       write(LOGfile,"(24F20.12)")e0,lambda,Mh,xmu,(n0(i),i=1,Norb),sum(n0)
    endif
    !
    allocate(Hloc(Norb,Norb))
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-9)Hloc=0.d0
    !
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of bands        :",Norb
  end subroutine build_hk








  !---------------------------------------------------------------------
  !PURPOSE: GET THE MODEL HAMILTONIAN ALONG A PATH IN THE 3D BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path()
    integer                                     :: i,j,ik,iorb,Lkpath
    integer                                     :: ix,iy,iz
    real(8)                                     :: kx,ky,kz
    integer                                     :: unit
    real(8)                                     :: Eval(Norb)
    real(8),allocatable                         :: Kpath(:,:)
    complex(8)                                  :: pHk(Norb,Norb)
    if(ed_mpi_id==0)write(LOGfile,*)"Build the model H(k) along a path in the 3D BZ:"
    !
    Lkpath=6*Nkstep
    allocate(kpath(Lk,3))
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Norb,Norb,Lk))
    ik = 0
    !1: (0,0,0)-->(0,0,pi)
    do iz=1,Nkstep
       ik = ik+1
       kx = 0.d0
       ky = 0.d0
       kz = 0.d0 + pi*(iz-1)/Nkstep
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !2: (0,0,pi)-->(0,pi,pi)
    do iy=1,Nkstep
       ik = ik+1
       kx = 0.d0
       ky = 0.d0 + pi*(iy-1)/Nkstep
       kz = pi
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !3: (0,pi,pi)-->(0,pi,0)
    do iz=1,Nkstep
       ik=ik+1
       kx = 0.d0
       ky = pi
       kz = pi - pi*(iz-1)/Nkstep
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !4: (0,pi,0)-->(pi,pi,0)
    do ix=1,Nkstep
       ik  =ik+1
       kx = 0.d0 + pi*(ix-1)/Nkstep
       ky = pi
       kz = 0.d0
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !5: (pi,pi,0)-->(pi,pi,pi)
    do iz=1,Nkstep
       ik = ik+1
       kx = pi
       ky = pi
       kz = 0.d0 + pi*(iz-1)/Nkstep
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !6: (pi,pi,pi)-->(0,0,0)
    do ix=1,Nkstep
       ik = ik+1
       kx = pi - pi*(ix-1)/Nkstep
       ky = pi - pi*(ix-1)/Nkstep
       kz = pi - pi*(ix-1)/Nkstep
       kpath(ik,:)=[kx,ky,kz]
    enddo
    !
    unit=free_unit() 
    if(ED_MPI_ID==0)open(unit,file="Eigenbands.dat")
    do ik=1,Lk
       Hk(:,:,ik) = Hk_model(kpath(ik,1),kpath(ik,2),kpath(ik,3))
       pHk = Hk(:,:,ik)
       call Eigensolve(pHk,Eval)
       if(ED_MPI_ID==0)write(unit,'(I,10F18.12)')ik,(Eval(iorb),iorb=1,Norb)
    enddo
    if(ED_MPI_ID==0)close(unit)
    !
  end subroutine build_hk_path













  !--------------------------------------------------------------------!
  !MODEL HAMILTONIAN H(k)
  !--------------------------------------------------------------------!
  function Hk_model(kx,ky,kz) result(Hk)
    real(8)                         :: kx,ky,kz
    real(8)                         :: epsik,vpsik
    complex(8),dimension(Norb,Norb) :: Hk
    epsik = cos(kx)+cos(ky)+cos(kz)
    vpsik = (cos(kx)-cos(ky))*cos(kz)
    Hk(1,1) = mh - e0*epsik
    Hk(2,2) =-mh - e0*epsik
    Hk(1,2) = lambda*vpsik
    Hk(2,1) = lambda*vpsik
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
  function inverse_g0k(zeta,hk) result(gk)
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

  subroutine eigensolve(hk,Evals)
    complex(8),dimension(Norb,Norb) :: hk
    real(8),dimension(Norb)         :: evals
    real(8)                         :: eplus,eminus
    complex(8),dimension(2,2)       :: uk
    complex(8)                      :: delta00,ek,u,v
    !Evals
    eplus   = hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eplus   = eplus/2.d0
    eminus  = hk(1,1)+hk(2,2) -sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus  = eminus/2.d0
    !Eigenvectors
    delta00  = -(hk(1,1)-hk(2,2))
    ek      = sqrt( delta00**2 + 4.d0*hk(1,2)*hk(2,1) )
    u       = sqrt(0.5d0*(1.d0 - delta00/ek))
    v       = sqrt(0.5d0*(1.d0 + delta00/ek))
    uk(1,1) = u
    uk(2,2) = u
    uk(1,2) = v
    uk(2,1) =-v
    !
    Hk = Uk
    Evals(1)=Eplus;Evals(2)=Eminus
  end subroutine Eigensolve


  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT








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

end program ed_hm2bhyb_fcc



