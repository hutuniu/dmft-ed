!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: iloop,Nb(2),Ne,ie,mu_loop,Nk
  logical                :: converged,dos_3d
  real(8),allocatable    :: wm(:),wr(:)
  real(8)                :: wband,ts,de
  real(8),allocatable    :: dens(:),docc(:) 
  !Bath:
  real(8),allocatable    :: bath(:,:),bath_old(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  real(8),allocatable    :: epsik(:),wt(:)
  real(8)                :: crystal_field,var,wmix
  integer                :: unit,iorb
  
  complex(8),allocatable,dimension(:,:,:) :: Hk


#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  call parse_input_variable(wband,"WBAND","inputED.in",default=1.d0)
  call parse_input_variable(Ne,"NE","inputED.in",default=2000)
  call parse_input_variable(Nk,"Nk","inputED.in",default=10)
  call parse_input_variable(dos_3d,"DOS","inputED.in",default=.false.)
  call parse_input_variable(crystal_field,"CRYSTAL_FIELD","inputED.in",default=0.d0)
  call parse_input_variable(wmix,"WMIX","inputED.in",default=1.d0)
  call ed_read_input("inputED.in")
  !+- allocate frequency arrays -+!
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr = linspace(wini,wfin,Lreal)
  !+- allocate Weiss Field +-!
  allocate(delta(Norb,Norb,Lmats))
  !+- allocate observables -+!
  allocate(dens(Norb),docc(Norb))
  !+- allocate k-grids -+! 
  ts=wband/6.d0
  if(dos_3d) then
     call get_cubic_dos(100000000)
  else
     call get_cubic_k(Nk)
  end if
  !+- setup ED_solver -+!
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  allocate(bath_old(Nb(1),Nb(2)))
  call init_ed_solver(bath)
  bath_old=bath
  Hloc(1,1,1,1) =  crystal_field*0.5d0
  Hloc(1,1,2,2) = -crystal_field*0.5d0
  if(ED_MPI_ID==0)then
     write(LOGfile,*)"Updated Hloc:"
     call print_Hloc(Hloc)
  endif
  call build_Hk
  !+- DMFT loop -+!
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0) call start_loop(iloop,nloop,"DMFT-loop")
     bath_old=bath
     !+- Solve the EFFECTIVE 2BAND-IMPURITY PROBLEM -+!
     call ed_solver(bath) 
     !+- Get the WEISS FILED/HYBRIDIZATION function to be fitted (user defined) -+!
     call get_delta  
     !+- Perform the SELF-CONSISTENCY ----> fitting the new bath -+!
     call chi2_fitgf(delta,bath,ispin=1)
     bath=wmix*bath + (1.d0-wmix)*bath_old
     !+- check CONVERGENCE and adjust CHEMICAL POTENTIAL (if required) -+!
     if(ED_MPI_ID==0)  converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0) call search_chemical_potential(xmu,ed_dens(1)+ed_dens(2),converged)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
     if(ED_MPI_ID==0) call end_loop
  enddo
  !Get Kinetic Energy too
  if(ED_MPI_ID==0) call ed_kinetic_energy(impSmats(1,1,:,:,:),Hk,wt)

#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif

contains

  subroutine get_delta
    integer                     :: i,j,ie,iorb,jorb
    complex(8)                  :: iw,zita,g0and,g0loc,gg
    complex(8),dimension(Lmats) :: self
    complex(8),dimension(Lreal) :: selfr,grloc
    !
    complex(8),allocatable      :: GLoc_mats(:,:,:),GLoc_real(:,:,:)
    complex(8),allocatable      :: tmp_gloc(:,:)
    !
    allocate(GLoc_mats(Norb,Norb,Lmats),GLoc_real(Norb,Norb,Lreal))
    allocate(tmp_gloc(Norb,Norb))    
    !+- MATSUBARA FREQ -+!
     if(ED_MPI_ID==0) print*,"Get Gloc_iw:"
    delta = zero
    do i=1,Lmats
       iw = xi*wm(i)+xmu
       !+- compute loacl greens function -+!
       GLoc_mats(:,:,i)=zero
       do ie=1,Ne
          GLoc_mats(:,:,i) = GLoc_mats(:,:,i) +  gk(iw,impSmats(1,1,:,:,i),Hk(:,:,ie))*wt(ie)
       end do
       if(cg_scheme=='weiss') then
          tmp_gloc=GLoc_mats(:,:,i)
          call matrix_inverse(tmp_gloc)
          delta(:,:,i) = tmp_gloc + impSmats(1,1,:,:,i)
          call matrix_inverse(delta(:,:,i))
       else
          tmp_gloc=GLoc_mats(:,:,i)
          call matrix_inverse(tmp_gloc)
          delta(:,:,i) = inverse_g0imp(iw) - tmp_gloc - impSmats(1,1,:,:,i)
       end if
    end do
    !+- REAL FREQ -+!
    if(ED_MPI_ID==0) print*,"Get Gloc_real:"
    do i=1,Lreal
       iw=cmplx(wr(i),eps)+xmu
       GLoc_real(:,:,i) = zero       
       do ie=1,Ne
          GLoc_real(:,:,i) = GLoc_real(:,:,i) + gk(iw,impSreal(1,1,:,:,i),Hk(:,:,ie))*wt(ie)
       enddo
    enddo
    !
    if(ED_MPI_ID==0) then
       if(bath_type.eq.'normal') then
          do iorb=1,Norb
             call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed",wm,GLoc_mats(iorb,iorb,:))
             call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_real.ed",wr,GLoc_real(iorb,iorb,:))
             call splot("DOS_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_real.ed",wr,-dimag(GLoc_real(iorb,iorb,:))/pi)
             call splot("Delta_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,iorb,:))
          end do
       else
          do iorb=1,Norb
             do jorb=1,Norb          
                call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed",wm,GLoc_mats(iorb,jorb,:))
                call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_real.ed",wr,GLoc_real(iorb,jorb,:))
                call splot("Delta_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
             end do
          end do
       end if
    end if
    deallocate(GLoc_mats,GLoc_real)

 
  end subroutine get_delta

  !+- 2BAND GREENS FUNCTION -+!  
  function g0k(iw,hk) 
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k,g0k_
    complex(8)                  :: det,ginv11,ginv22
    complex(8)                  :: vmix12,vmix21
    g0k =zero
    g0k_=zero
    !+- build G0^(-1) -+!
    g0k_(1,1) = iw 
    g0k_(2,2) = iw 
    g0k_=g0k_-hk
    det = g0k_(1,1)*g0k_(2,2) - g0k_(1,2)*g0k_(2,1)
    !+- compute inverse -+!
    g0k(1,1) =  g0k_(2,2)/det
    g0k(2,2) =  g0k_(1,1)/det
    g0k(1,2) = -g0k_(1,2)/det
    g0k(2,1) = -g0k_(2,1)/det
  end function g0k
  !
  function gk(iw,sigma,hk) 
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: sigma
    complex(8),dimension(2,2)   :: hk
    complex(8),dimension(2,2)   :: gk,gk_
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: det
    gk  = zero
    gk_ = zero
    g0k = zero
    g0k = inverse_g0k(iw,hk)
    !+- build G^(-1) -+!
    gk_ = g0k-sigma
    !+- compute inverse -+!
    det = gk_(1,1)*gk_(2,2) - gk_(1,2)*gk_(2,1)
    gk(1,1) =  gk_(2,2)/det
    gk(2,2) =  gk_(1,1)/det
    gk(1,2) = -gk_(1,2)/det
    gk(2,1) = -gk_(2,1)/det
  end function gk

  !+- INVERSE 2BAND GREENS FUNCTION -+!
  function inverse_g0k(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: ginv11,ginv22
    complex(8)                  :: vmix12,vmix21
    g0k=zero
    !+- build G0^(-1) -+!
    g0k(1,1) = iw
    g0k(2,2) = iw
    g0k = g0k - hk
  end function inverse_g0k
  !
  function inverse_g0imp(iw) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: ginv11,ginv22
    complex(8)                  :: vmix12,vmix21
    g0k=zero
    !+- build G0_imp^(-1) -+!
    g0k(1,1) = iw
    g0k(2,2) = iw    
    g0k = g0k - Hloc(1,1,:,:)
  end function inverse_g0imp
  !
  function inverse_gk(iw,sigma,hk) result(gk)
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: sigma
    complex(8),dimension(2,2)   :: hk
    complex(8),dimension(2,2)   :: gk,gk_
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: det
    gk  = zero
    g0k = zero
    g0k = inverse_g0k(iw,hk)
    !+- build G^(-1) -+!
    gk = g0k-sigma
  end function inverse_gk
  !

  !+- HK routines -+!
  subroutine build_hk
    integer  :: ik
    allocate(Hk(Norb,Norb,Ne))
    Hk=0.d0
    do ik=1,Ne
       call get_hk(Hk(:,:,ik),ik)
    end do
  end subroutine build_hk

  subroutine get_hk(hk,ik)
    complex(8),dimension(Norb,Norb),intent(inout) :: hk
    integer  :: ik
    hk=zero
    do iorb=1,Norb
       hk(iorb,iorb) = epsik(ik) 
    end do
    hk = hk + Hloc(1,1,:,:)
  end subroutine get_hk
  
  
  
  subroutine get_cubic_dos(Nrnd)
    integer :: Nrnd
    integer :: i,ik,unit
    integer(8),allocatable,dimension(:) :: count_
    real(8) :: kx,ky,kz,ek
    allocate(wt(Ne),epsik(Ne))
    allocate(count_(Ne))
    epsik=linspace(-wband,wband,Ne,mesh=de)
    count_=0
    call random_seed(put=[1234567])

    if(ED_MPI_ID==0)write(LOGfile,*)"Building 3d DOS"
    do i=1,Nrnd
       call random_number(kx)
       kx=(kx-1.d0)*pi
       call random_number(ky)
       ky=(ky-1.d0)*pi
       call random_number(kz)    
       kz=(kz-1.d0)*pi              
       ek=-2.d0*ts*(cos(kx)+cos(ky)+cos(kz))
       ik=1+ek/de + wband/de
       count_(ik)=count_(ik)+1   
    end do
    !+- symmetryze -+!
    do i=1,Ne/2
       count_(i) = count_(Ne-(i-1))
    end do
    wt=dble(count_)/dble(Nrnd) 
    wt=wt/trapz(de,wt)
    wt=wt*de
    if(ED_MPI_ID==0) then
       unit=free_unit()
       open(unit,file='DOS.3d')
       do i=1,Ne
          write(unit,*) epsik(i),wt(i)
       end do
       close(unit)
    end if
  end subroutine get_cubic_dos

  subroutine get_cubic_k(Nk)
    integer :: Nk
    integer :: i,ik
    integer :: ix,iy,iz
    integer(8),allocatable,dimension(:) :: count_
    real(8),allocatable,dimension(:) :: kx
    real(8) :: tmp
    
    allocate(kx(Nk))
    kx=linspace(0.d0,pi,Nk)
    Ne=8.d0*Nk*Nk*Nk
    allocate(wt(Ne),epsik(Ne))    
    ik=0
    do ix=1,Nk
       do iy=1,Nk
          do iz=1,Nk
             ik=ik+1
             epsik(ik)=-2.d0*ts*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz)))
             wt(ik) = 8.d0/dble(Ne)
          end do
       end do
    end do
    call get_free_dos(epsik,wt,file='DOS_free.kgrid')
  end subroutine get_cubic_k


end program lancED





  !+-
  ! tmp_bath=bath
  ! call allocate_bath(tmp_dmft_bath)
  ! call set_bath(tmp_bath,tmp_dmft_bath)
  ! tmp_dmft_bath%e(:,1,:)=tmp_dmft_bath%e(:,1,:) + crystal_field!*0.5d0
  ! tmp_dmft_bath%e(:,2,:)=tmp_dmft_bath%e(:,2,:) - crystal_field!*0.5d0
  ! call copy_bath(tmp_dmft_bath,tmp_bath)
  ! call write_bath(tmp_dmft_bath,890)
  ! call deallocate_bath(tmp_dmft_bath)
  ! bath=tmp_bath
  !+-
  ! real(8),allocatable    :: tmp_bath(:,:)
  ! type(effective_bath)   :: tmp_dmft_bath



  ! subroutine get_internal_energy
  !   integer :: ie,i,unit
  !   complex(8),dimension(Lmats) :: G0_loc,G_loc,kin_G_loc,kin_G0_loc
  !   complex(8) :: zita0,zita
  !   real(8) :: kin_ene,pot_ene
  !   real(8) :: kin_ene0,kin_ene0_gl

  !   kin_ene0=0.d0
  !   do ie=1,Ne
  !      kin_ene0 = kin_ene0 + 2.d0*epsik(ie)*wt(ie)*fermi(epsik(ie),beta)
  !   end do


  !   do i=1,Lmats
  !      zita0 = xi*wm(i) 
  !      G0_loc(i)=zero
  !      kin_G0_loc(i)=zero
  !      do ie=1,Ne
  !         G0_loc(i) = G0_loc(i) + wt(ie)/(zita0-epsik(ie))
  !         kin_G0_loc(i) = kin_G0_loc(i) + wt(ie)/(zita0-epsik(ie))*epsik(ie)
  !      end do
  !   end do

  !   unit=free_unit()
  !   open(unit,file='kin_Gloc')
  !   do i=1,Lmats
  !      zita    = xi*wm(i) - impSmats(1,1,1,1,i)
  !      G_loc(i) = zero
  !      kin_G_loc(i)=zero
  !      do ie=1,Ne
  !         G_loc(i)=G_loc(i)+wt(ie)/(zita-epsik(ie))
  !         kin_G_loc(i)=kin_G_loc(i)+wt(ie)/(zita-epsik(ie))*epsik(ie)
  !      enddo
  !      write(unit,'(6(F18.10))') wm(i),kin_G_loc(i)
  !   enddo

  !   kin_ene = 0.d0
  !   do i=1,Lmats
  !      zita0 = xi*wm(i)
  !      kin_ene = kin_ene + 2.d0*dreal(zita0*(G_loc(i)-G0_loc(i)))
  !      kin_ene = kin_ene - 2.d0*dreal(impSmats(1,1,1,1,i)*G_loc(i))
  !   end do

  !   kin_ene = 2.d0/beta*kin_ene

  !   pot_ene = Uloc(1)*ed_docc(1)
  !   unit = free_unit()
  !   open(unit,file='internal_energy.data')
  !   write(unit,'(6(F18.10))') Uloc(1),kin_ene+pot_ene,kin_ene,pot_ene
  !   close(unit)

  !   kin_ene=0.d0
  !   kin_ene0_gl=0.d0
  !   do i=1,Lmats
  !      zita0 = xi*wm(i)
  !      kin_ene = kin_ene + 2.d0*dreal(kin_G_loc(i))
  !      kin_ene0_gl = kin_ene0_gl + 2.d0*dreal(kin_G0_loc(i))
  !   end do
  !   kin_ene = 2.d0/beta*kin_ene 
  !   kin_ene0_gl = 2.d0/beta*kin_ene0_gl 

  !   unit = free_unit()
  !   open(unit,file='internal_energy_.data')
  !   write(unit,'(6(F18.10))') Uloc(1),kin_ene+pot_ene,kin_ene,pot_ene
  !   close(unit)


  !   unit = free_unit()
  !   open(unit,file='internal_energy_0.data')
  !   write(unit,'(6(F18.10))') kin_ene0,kin_ene0_gl
  !   close(unit)

  ! end subroutine get_internal_energy


