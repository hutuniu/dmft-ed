program ed_ah
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                :: iloop,Nb(2),Lk,Nx
  logical                :: converged
  real(8)                :: wband,ts,wmixing
  !Bath:
  real(8),allocatable    :: Bath(:,:),BathOld(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:)
  complex(8),allocatable :: Gmats(:,:),Greal(:,:)
  complex(8),allocatable :: Smats(:,:),Sreal(:,:)
  character(len=16)      :: finput,fhloc
  logical                :: phsym,normal_bath
  real(8),allocatable    :: Hk(:),wt(:),kxgrid(:),kygrid(:),wm(:),wr(:)

  ! call MPI_INIT(mpiERR)
  ! call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  ! call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  ! write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  ! call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.5d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(phsym,"phsym",finput,default=.false.,comment="Flag to enforce p-h symmetry of the bath.")
  call parse_input_variable(normal_bath,"normal",finput,default=.false.,comment="Flag to enforce no symmetry braking in the bath.")
  !
  call ed_read_input(trim(finput))

  !Allocate Weiss Field:
  allocate(delta(2,Norb,Norb,Lmats))
  allocate(Gmats(2,Lmats),Greal(2,Lreal))
  allocate(Smats(2,Lmats),Sreal(2,Lreal))

  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  allocate(bathold(Nb(1),Nb(2)))
  call ed_init_solver(bath)


  Lk = Nx*Nx
  allocate(Hk(Lk),Wt(Lk))
  allocate(kxgrid(Nx),kygrid(Nx))
  write(*,*) "Using Nk_total="//txtfy(Lk)
  kxgrid = kgrid(Nx)
  kygrid = kgrid(Nx)
  Hk     = build_hk_model(hk_model,kxgrid,kygrid,[0d0])
  Wt     = 1d0/Lk
  call write_hk_w90("Hk2d.dat",1,1,0,1,dcmplx(Hk,0d0),kxgrid,kygrid,[0d0])
  call get_free_dos(Hk,Wt)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 
     call ed_get_sigma_matsubara(Smats(1,:))
     call ed_get_sigma_real(Sreal(1,:))
     call ed_get_self_matsubara(Smats(2,:))
     call ed_get_self_real(Sreal(2,:))

     !Get the Weiss field/Delta function to be fitted (user defined)
     !call get_delta
     call ed_get_gloc(dcmplx(Hk,0d0),Wt,Gmats,Greal,Smats,Sreal)
     call ed_get_weiss(Gmats,Smats,Delta(:,1,1,:))
     call splot("Gloc_iw.ed",wm,Gmats(1,:))
     call splot("Floc_iw.ed",wm,Gmats(2,:))
     call splot("Delta_iw.ed",wm,delta(1,1,1,:),delta(2,1,1,:))
     call splot("Gloc_realw.ed",wr,-dimag(Greal(1,:))/pi,dreal(Greal(1,:)))
     call splot("Floc_realw.ed",wr,-dimag(Greal(2,:))/pi,dreal(Greal(2,:)))


     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(delta,bath,ispin=1)
     if(phsym)call ph_symmetrize_bath(bath)
     if(normal_bath)call enforce_normal_bath(bath)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
     BathOld=Bath
     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,:)+delta(2,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu,converged)
     call end_loop
  enddo

  !call get_sc_optical_conductivity
  !call get_sc_internal_energy(Lmats)

  call  ed_kinetic_energy(Hk,wt,impSmats(1,1,1,1,:),impSAmats(1,1,1,1,:))

  ! call MPI_FINALIZE(mpiERR)


contains


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)           :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -one*2d0*ts*(cos(kx)+cos(ky))
  end function hk_model









  ! !+----------------------------------------+
  ! subroutine get_delta
  !   integer                    :: i,j,iorb,ik
  !   complex(8)                 :: iw,zita,g0loc,cdet,zita1,zita2
  !   complex(8),dimension(Lreal)   :: zeta
  !   complex(8),dimension(2,Lmats) :: gloc,calG
  !   complex(8),dimension(2,Lreal) :: grloc
  !   real(8)                    :: wm(Lmats),wr(Lreal),tau(0:Lmats)
  !   complex(8),dimension(Norb,Norb) :: Hloc
  !   wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  !   wr = linspace(wini,wfin,Lreal)
  !   call get_Hloc(Hloc,1)
  !   delta=zero
  !   do i=1,Lmats
  !      iw = xi*wm(i)
  !      zita    = iw + xmu - Hloc(1,1) - impSmats(1,1,1,1,i) 
  !      gloc(:,i)=zero
  !      do ik=1,Lk
  !         cdet = abs(zita-Hk(ik))**2 + impSAmats(1,1,1,1,i)**2
  !         gloc(1,i)=gloc(1,i) + wt(ik)*(conjg(zita)-Hk(ik))/cdet
  !         gloc(2,i)=gloc(2,i) - wt(ik)*impSAmats(1,1,1,1,i)/cdet
  !      enddo
  !      if(cg_scheme=='weiss')then
  !         !Get G0^{-1} matrix components:
  !         cdet      =  abs(gloc(1,i))**2 + (gloc(2,i))**2
  !         calG(1,i) =  conjg(gloc(1,i))/cdet + impSmats(1,1,1,1,i)
  !         calG(2,i) =  gloc(2,i)/cdet        + impSAmats(1,1,1,1,i) 
  !         !Get Weiss field G0 components:
  !         cdet            =  abs(calG(1,i))**2 + (calG(2,i))**2
  !         delta(1,1,1,i)  =  conjg(calG(1,i))/cdet
  !         delta(2,1,1,i)  =  calG(2,i)/cdet
  !      else
  !         cdet            = abs(gloc(1,i))**2 + (gloc(2,i))**2
  !         delta(1,1,1,i)  = zita  - conjg(gloc(1,i))/cdet 
  !         delta(2,1,1,i)  = - impSAmats(1,1,1,1,i) - gloc(2,i)/cdet 
  !      endif
  !   enddo
  !   !
  !   zeta(:) = cmplx(wr(:),eps,8) + xmu - impSreal(1,1,1,1,:)
  !   do i=1,Lreal
  !      zita1 = zeta(i)
  !      zita2 = conjg(zeta(Lreal+1-i))
  !      grloc(:,i) = zero
  !      do ik=1,Lk
  !         cdet = (zita1-Hk(ik))*(zita2-Hk(ik)) + impSAreal(1,1,1,1,i)*impSAreal(1,1,1,1,i)
  !         grloc(1,i) = grloc(1,i) + wt(ik)*(zita2-Hk(ik))/cdet
  !         grloc(2,i) = grloc(2,i) + wt(ik)*impSAreal(1,1,1,1,i)/cdet
  !      enddo
  !   enddo
  !   call splot("Gloc_iw.ed",wm,gloc(1,:))
  !   call splot("Floc_iw.ed",wm,gloc(2,:))
  !   call splot("Gloc_realw.ed",wr,grloc(1,:))
  !   call splot("Floc_realw.ed",wr,grloc(2,:))
  !   call splot("DOS.ed",wr,-dimag(grloc(1,:))/pi)
  !   call splot("Delta_iw.ed",wm,delta(1,1,1,:),delta(2,1,1,:))
  ! end subroutine get_delta
  ! !+----------------------------------------+




  ! subroutine  get_sc_optical_conductivity()  
  !   integer                :: i,ik,iv,iw  
  !   real(8),allocatable    :: oc(:), Ak(:,:,:),cDOS(:,:),ocw(:)
  !   complex(8)            :: zeta(Lreal)
  !   integer               :: Nw
  !   real(8)               :: wr(Lreal),dw
  !   complex(8)            :: cdet,zita1,zita2,fg(2)
  !   real(8),dimension(Lk) :: epsik,wt
  !   real(8)                :: vel2,dos,Dfermi,A2,B2,ock


  !   call bethe_lattice(wt,epsik,Lk,wband)

  !   print*,"Get OC with:",Lreal,"freq."
  !   Nw=Lreal/2

  !   allocate(Ak(2,Lk,Lreal))
  !   wr = linspace(wini,wfin,Lreal,mesh=dw)

  !   allocate(cDOS(2,Lreal))
  !   cDOS=0.d0
  !   zeta(:) = cmplx(wr(:),eps,8) + xmu - impSreal(1,1,1,1,:)
  !   do i=1,Lreal
  !      zita1 = zeta(i)
  !      zita2 = conjg(zeta(Lreal+1-i))
  !      do ik=1,Lk
  !         cdet = (zita1-epsik(ik))*(zita2-epsik(ik)) + impSAreal(1,1,1,1,i)*impSAreal(1,1,1,1,i)
  !         fg(1)=(zita2-epsik(ik))/cdet
  !         fg(2)=impSAreal(1,1,1,1,i)/cdet
  !         Ak(1,ik,i)=-dimag(fg(1))/pi
  !         Ak(2,ik,i)=-dimag(fg(2))/pi
  !         cDOS(1,i)=cDOS(1,i)+Ak(1,ik,i)*wt(ik)
  !         cDOS(2,i)=cDOS(2,i)+Ak(2,ik,i)*wt(ik)
  !      enddo
  !   enddo
  !   call splot("ocDOS.last",wr,cDOS(1,:),cDOS(2,:))
  !   deallocate(cDOS)



  !   !Changing the loop order does not affect the calculation.
  !   allocate(oc(Nw),ocw(Lreal))
  !   oc=0.d0
  !   call start_progress
  !   do iv=1,Nw
  !      ocw   = 0.d0
  !      do iw=1,Lreal-iv
  !         !Dfermi  = istep(wr(iw)) - istep(wr(iw+iv))
  !         Dfermi  = fermi(wr(iw),beta) - fermi(wr(iw+iv),beta)
  !         ock=0.d0
  !         do ik=1,Lk
  !            A2 = Ak(1,ik,iw)*Ak(1,ik,iw+iv)
  !            B2 = Ak(2,ik,iw)*Ak(2,ik,iw+iv)
  !            vel2= (wband**2-epsik(ik)**2)/3.d0
  !            ock = ock + vel2*(A2-B2)*wt(ik)
  !         enddo
  !         ocw(iw) = Dfermi*ock
  !      enddo
  !      oc(iv)=trapz(dw,ocw(:Lreal-iv))/wr(Nw+iv)
  !      call progress(iv,Nw)
  !   enddo
  !   call stop_progress
  !   call splot("OC_realw.ed",wr(Nw+1:Lreal),oc(:))
  !   call splot("OC_integral.ed",uloc(1),trapz(dw,oc))

  ! end subroutine get_sc_optical_conductivity





  ! subroutine get_sc_internal_energy(L)
  !   integer                       :: L
  !   real(8)                       :: wm(Lmats)
  !   complex(8)                    :: fg(2,L),sigma(2,L)
  !   real(8)                       :: matssum,fmatssum,checkP,checkdens,vertex,Dssum
  !   complex(8)                    :: iw,gkw,fkw,g0kw,f0kw
  !   real(8)                       :: Epot,Etot,Eint,kin,kinsim,Ds,docc
  !   real(8)                       :: Sigma_infty,S_infty,det,det_infty,csi,Ei,thermal_factor
  !   real(8)                       :: free(Lk),Ffree(Lk),n_k(Lk),n,delta,u,ts
  !   integer                       :: i,j,iorb,ik
  !   complex(8)                    :: zita,g0loc,cdet,zita1,zita2
  !   complex(8),dimension(Lmats)   :: zeta
  !   real(8),dimension(Lk)         :: epsik,wt

  !   wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  !   call bethe_lattice(wt,epsik,Lk,1.d0)
  !   ts=0.5d0
  !   sigma(1,:)=impSmats(1,1,1,1,:)
  !   sigma(2,:)=impSAmats(1,1,1,1,:)!-conjg(impSAmats(1,1,1,1,:))
  !   fg=zero
  !   do i=1,L
  !      iw   = xi*wm(i)
  !      zita = iw + xmu - sigma(1,i)
  !      do ik=1,Lk
  !         cdet = abs(zita-epsik(ik))**2 + sigma(2,i)**2
  !         fg(1,i)=fg(1,i) + wt(ik)*(conjg(zita)-epsik(ik))/cdet
  !         fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
  !      enddo
  !   enddo
  !   !fg(2,:)=-conjg(fg(2,:))
  !   u    = uloc(1)
  !   n    = ed_dens(1)/2.d0
  !   delta= ed_phisc(1)*u

  !   !Get asymptotic self-energies
  !   Sigma_infty =   dreal(sigma(1,L))
  !   S_infty     =   dreal(sigma(2,L))

  !   checkP=0.d0 ; checkdens=0.d0 ;          ! test variables

  !   kin=0.d0                      ! kinetic energy (generic)
  !   Ds=0.d0                       ! superfluid stiffness (Bethe)
  !   do ik=1,Lk
  !      csi            = epsik(ik)-(xmu-Sigma_infty)
  !      Ei             = dsqrt(csi**2 + S_infty**2)
  !      thermal_factor = dtanh(0.5d0*beta*Ei)
  !      free(ik)        = 0.5d0*(1.d0 - csi/Ei)*thermal_factor
  !      Ffree(ik)       =-(0.5d0*S_infty)/Ei*thermal_factor
  !      fmatssum= 0.d0
  !      matssum = 0.d0
  !      Dssum   = 0.d0
  !      vertex=(4.d0*ts**2-epsik(ik)**2)/3.d0
  !      do i=1,L
  !         iw       = xi*wm(i)
  !         det      = abs(iw+xmu-epsik(ik)-sigma(1,i))**2 + dreal(sigma(2,i))**2
  !         det_infty= wm(i)**2 + (epsik(ik)-(xmu-Sigma_infty))**2 + S_infty**2
  !         gkw = (-iw+xmu - epsik(ik) - conjg(sigma(1,i)) )/det
  !         fkw = -sigma(2,i)/det
  !         g0kw= (-iw - (epsik(ik)-(xmu-Sigma_infty)))/det_infty
  !         f0kw=-S_infty/det_infty
  !         matssum =  matssum +  dreal(gkw)-dreal(g0kw)
  !         fmatssum= fmatssum +  dreal(fkw)-dreal(f0kw)
  !         Dssum   = Dssum    +  fkw*fkw
  !      enddo
  !      n_k(ik)   = 4.d0/beta*matssum + 2.d0*free(ik)
  !      checkP    = checkP    - wt(ik)*(2.d0/Beta*fmatssum+Ffree(ik))
  !      checkdens = checkdens + wt(ik)*n_k(ik)
  !      kin    = kin    + wt(ik)*n_k(ik)*epsik(ik)
  !      Ds=Ds + 8.d0/beta* wt(ik)*vertex*Dssum
  !   enddo

  !   kinsim=0.d0
  !   kinsim = sum(fg(1,:)*fg(1,:)+conjg(fg(1,:)*fg(1,:))-2.d0*fg(2,:)*fg(2,:))*2.d0*ts**2/beta

  !   Epot=zero
  !   Epot = sum(fg(1,:)*sigma(1,:) + fg(2,:)*sigma(2,:))/beta*2.d0

  !   docc = 0.5d0*n**2
  !   if(u > 0.01d0)docc=-Epot/u + n - 0.25d0

  !   Eint=kin+Epot

  !   Ds=zero
  !   Ds = sum(fg(2,:)*fg(2,:))/beta*2.d0

  !   write(*,*)"Asymptotic Self-Energies",Sigma_infty, S_infty
  !   write(*,*)"n,delta",n,delta
  !   write(*,*)"Dn% ,Ddelta%",(n-0.5d0*checkdens)/n,(delta + u*checkP)/delta ! u is positive
  !   write(*,*)'========================================='
  !   write(*,*)"Kinetic energy",kin
  !   write(*,*)'========================================='
  !   write(*,*)"double occupancy   =",docc
  !   write(*,*)'========================================='
  !   write(*,*) 'Kinetic Energy TEST (simple formula)'
  !   write(*,*) '###ACTHUNG: FOR BETHE ONLY####',kinsim
  !   write(*,*) 'Dkin%',(kin-kinsim)/kin
  !   write(*,*)'========================================='
  !   write(*,*) 'Superfluid stiffness',Ds
  !   write(*,*) 'Potential Energy U(n_up-1/2)(n_do-1/2)',Epot
  !   write(*,*) 'Internal Energy',Eint
  !   write(*,*)'========================================='
  !   call splot("nk_distribution.ed",epsik,n_k/2.d0,free)
  !   open(100,file="columns.ed")
  !   write(100,"(11A21)")"1vbias","2u","3beta","4n","5kin","6docc","7Ds","8Epot","9Eint"
  !   close(100)
  !   open(200,file="thermodynamics.ed")
  !   write(200,"(11F21.12)")0.d0,u,beta,n,kinsim,docc,Ds,Epot,Eint
  !   close(200)
  !   open(300,file="energy_bethe.ed")
  !   write(300,"(11F21.12)") kin,Epot,Eint
  !   close(300)
  !   return 
  ! end subroutine get_sc_internal_energy

  ! elemental function istep(x) result(out)
  !   real(8),intent(in) :: x
  !   real(8)            :: out
  !   if(x < 0.d0) then
  !      out = 1.0d0
  !   elseif(x==0.d0)then
  !      out = 0.50d0
  !   else
  !      out = 0.0d0
  !   endif
  ! end function istep



end program ed_ah



