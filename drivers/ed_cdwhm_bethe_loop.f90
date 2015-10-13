!#####################################################################
!PURPOSE : Solution of DMFT problem for EHM model with CO long-range order.
!          Mixed use of single-site and lattice ("wrap") functions.           
!AUTHORS : K.J. Kapcia
!#####################################################################

program ed_cdw_loop

  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR

  implicit none
  
  integer                                       :: iloop,ineq
  logical                                       :: converged
  real(8)                                       :: wband,W0
  character(len=16)                             :: finput
  real(8)                                       :: wmixing
  integer                                       :: Nk
  real(8)                                       :: n0(2),nQ
  real(8)					:: Ekin(3),Epot(3),Docc(2),Emu,EWint,Etot,OmegaTot
  real(8)					:: Eout(2)
  real(8)					:: temporaryreal, error
  real(8)					:: xmubar, xmu0, xmusub(2)
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath(:,:),Bath_prev(:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta ![Nlat][Nspin][Nspin][Norb][Norb][Nfreq]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_loc
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_loc

  !Hamiltonian
  real(8),allocatable                           :: BetheDOS(:)
  complex(8),allocatable                        :: BetheEkImp(:,:,:)
  complex(8),allocatable                        :: BetheEk(:,:,:) 
  complex(8),allocatable                        :: Hloc(:)
  character(len=20)                             :: tmp_suffix
  
  ! tag for printing, tag for distinguish files,
      
  real(8)  		:: loopincrement,loopfinish,loopstart
  real(8)   		:: loop_min,loop_max,looplocvalue
  integer  		:: itertemp,loopnumberofsteps,looptag 
  integer 		:: printingtag, numfiles 
  character(len=80) 	:: temporaryfilename
  real(8) 		:: sumtest
  integer		:: i  

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  !
  call parse_input_variable(printingtag,"printingtag",finput,default=0)
  call parse_input_variable(numfiles,"numfiles",finput,default=1)  
  call parse_input_variable(loopstart,"LoopStart",finput,default=-0.5d0)
  call parse_input_variable(loopfinish,"LoopFinish",finput,default=0d0) 
  call parse_input_variable(loopnumberofsteps,"LoopNumberofSteps",finput,default=0)
  call parse_input_variable(looptag,"LoopTag",finput,default=1)
  !
  call parse_input_variable(Nk,"Nk",finput,default=500)
  call parse_input_variable(n0,"n0",finput,default=[0.1d0,1.9d0])
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(W0,"W0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0) 
  !
  call ed_read_input(trim(finput))
  
  ! starting three main files for results
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_it_',numfiles,'.dat'
  open(35,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')      
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_par_',numfiles,'.dat'
  open(36,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')
  WRITE(temporaryfilename,'(a,i1.1,a)')'energy_vs_par_',numfiles,'.dat'
  open(37,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')

  xmubar = xmu ! from starting file xmubar = 0 if n = 1

  Nlat=2      ! the two inequivalent sublattices A,B
  if(Nlat/=2)stop "This drivers requires Nlat==2"
  if(Norb/=1)stop "This drivers requires Norb==1"
  if(Nspin/=1)stop "This drivers requires Nspin==1"

  allocate(BetheDOS(Nk),BetheEk(Nlat,Nlat,Nk),BetheEkImp(Nlat,Nlat,Nk),Hloc(Nlat))
  
  !building the Hamiltonian
  call build_Hbethe_lat_2x2(n0) 
  !call build_Hloc(n0)
  !call build_Hbethe
  
  !check of the density of states
  !call splot("testDOSbethe.dat",dreal(BetheEk(1,2,:)),BetheDOS(:))
  !sumtest=0d0
  !do i=1,Nk
  !  sumtest = sumtest + BetheDOS(i)
  !enddo
  !print*,"sumtest = ",sumtest

  !Allocate Weiss Field and other:
  allocate(Delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_loc(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_loc(Nlat,Nspin,Nspin,Norb,Norb,Lreal)) 

  !get bath size, setup solver
  Nb=get_bath_size()
  print*,"Nb = ",Nb
  allocate(bath(Nlat,Nb),bath_prev(Nlat,Nb))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! set loop parameters  
! if loopnumberofsteps==0 only one eveluation
      if(loopnumberofsteps.eq.0) then
       loopincrement = 10000000
      else
       loopincrement=(loopfinish-loopstart)/loopnumberofsteps
      endif  
! loop with changing mu/W0/U/...
if (looptag.eq.1) then
       WRITE(temporaryfilename,'(a)')'mu -- LOOP'
      else if  (looptag.eq.2) then
       WRITE(temporaryfilename,'(a)')'W0 -- LOOP'
      else if  (looptag.eq.3) then
       WRITE(temporaryfilename,'(a)')'U -- LOOP'
      else
       WRITE(temporaryfilename,'(a)')'WRONG TAG !!!'
      end if      
! preparing for loop      
! determination of max and min values in looop - range of the loop
!      
      if (loopincrement > 0 ) then
         loop_max =   loopfinish
         loop_min =   loopstart
      else
         loop_max =   loopstart
         loop_min =   loopfinish
      end if
! looplocvalue is current value of the loop parameter
! we set a proper start 
      looplocvalue = loopstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  !do the mu/W0/U/... loop when current value is between min and max
  itertemp = 0

  do while ((looplocvalue - &
      0.5*abs(loopincrement) < loop_max)  &
      .and. (looplocvalue +  &
      0.5*abs(loopincrement) > loop_min) &
      .and. (itertemp.le.loopnumberofsteps ) )

  if (looptag.eq.1) then        
     	xmubar = looplocvalue
  else if  (looptag.eq.2) then
  	W0 = looplocvalue
      else if  (looptag.eq.3) then
           Uloc = looplocvalue
      else
          print*,'---- WRONG LOOPTAG ---- '
      end if
      !set proper value of mu to use in solvers      
      xmu0 = xmubar + W0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  !init solvers
  do ineq=1,Nlat
     write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
     ed_file_suffix="_path"//trim(tmp_suffix)
     call ed_init_solver(Bath(ineq,:))
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     
     iloop=iloop+1
     
     call start_loop(iloop,nloop,"DMFT-loop")
     
     bath_prev=bath

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     !Solve the impurities on each inequivalent site:
     do ineq=1,Nlat
        xmu=xmu0-W0*n0(3-ineq)
        
        write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
        ed_file_suffix="_path"//trim(tmp_suffix)

        call ed_solve(bath(ineq,:))
        n0(ineq) = ed_dens(1)

        call ed_get_sigma_matsubara(Smats(ineq,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(ineq,:,:,:,:,:))
        call ed_get_gimp_matsubara(Gmats(ineq,:,:,:,:,:))
        call ed_get_gimp_real(Greal(ineq,:,:,:,:,:))
     enddo
     
     !calculating proper local Green function (needed if not bethe !)
     call build_Hbethe_lat_2x2(n0) 
     !call build_Hloc(n0)
     
     xmu = xmu0 
     call ed_get_gloc_lattice(BetheEkImp,BetheDOS,Gmats_loc,Greal_loc,Smats,Sreal,iprint=1)
     call get_gloc_CDW(BetheEKimp(1,2,:),BetheDOS,Gmats_loc,Greal_loc,Smats,Sreal,iprint=1)
     
     print*,""
!     print'(A4,F12.10,A4,F12.10)',"mi1= ",xmu0-W0*n0(2)," mi2=",xmu0-W0*n0(1)
     print'(A4,F12.10,A4,F12.10)',"n1 = ",n0(1)," n2= ",n0(2)
     print*,""

     do ineq=1,Nlat
        xmu=xmu0-W0*n0(3-ineq)
        write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
        ed_file_suffix="_path"//trim(tmp_suffix)
        !!! call ed_get_weiss(Gmats_loc(ineq,:,:,:,:,:),Smats(ineq,:,:,:,:,:),Delta(3-ineq,:,:,:,:,:),iprint=1)
        !!! call ed_get_weiss(Gmats(ineq,1,1,1,1,:),Smats(ineq,1,1,1,1,:),Delta(ineq,1,1,1,1,:),iprint=1)
        !only in such way, G_loc needed if not bethe !:
        if (cg_scheme=='delta') then
          Delta(ineq,1,1,1,1,:) = (wband**2)/4d0*Gmats_loc(3-ineq,1,1,1,1,:)        
        else
          Delta(ineq,1,1,1,1,:) = one/( xi*pi/beta*(2*arange(1,Lmats)-1) + xmu - (wband**2)/4d0*Gmats_loc(3-ineq,1,1,1,1,:) )
        endif
        call ed_chi2_fitgf(Delta(ineq,1,1,:,:,:),bath(ineq,:),ispin=1)
     enddo
     
      ! G_loc needed if not bethe !
      !xmu = xmu0 
      
      !call ed_get_weiss_lattice(Nlat,Gmats(:,1,1,1,1,:),Smats(:,1,1,1,1,:),Delta(:,1,1,1,1,:),Eloc=Hloc) ! always for delta nad weiss scheme
      
      !call ed_chi2_fitgf_lattice(bath(:,:),Delta(:,1,1,1,1,:),ispin=1) !OK only for weiss scheme
      !call ed_chi2_fitgf_lattice(bath(:,:),Delta(:,1,1,1,1,:),Eloc=Hloc,ispin=1) !for delta and weiss scheme 

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,1,1,:),dmft_error,nsuccess,nloop)

     !if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu0,converged)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !'results  of the interation - to file'  
      
      write(35,'(2(I4.1),8(F15.10))')itertemp,iloop, &
      Uloc(1),W0,xmu0-W0,1.d0/beta, &
      n0(1),n0(2),0.5d0*(n0(1)+n0(2)),abs(0.5d0*(n0(1)-n0(2)))
      call flush(35) 
      call end_loop
  enddo !finish of DMFT loop
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Printing output, after finish dmft loop
  !in this funtion Hloc is important and Hk must include Hloc !
    xmu = xmu0 
    Eout = ed_kinetic_energy_lattice(BetheEkImp(:,:,:),BetheDOS,Smats(:,1,1,1,1,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !WRITE(temporaryfilename,'(a9)')'error.err'
  !open(40,file = temporaryfilename, form ='formatted',status = 'old')
  
  !'the most important results of the single loop'
  !file 1 - order parameters
  write(36,'(I3.1,8(F15.10),I4.1)')itertemp, &
      Uloc(1),W0,xmubar,1.d0/beta, &
      n0(1),n0(2),0.5d0*(n0(1)+n0(2)),abs(0.5d0*(n0(1)-n0(2))), &
      iloop
  call flush(36)
  
  !get informations from files
  !Ekin
  !do ineq=1,2
  !  WRITE(temporaryfilename,'(a17,i1.1,a5,i1.1,a3)')'kinetic_last_path',numfiles,'_site',ineq,'.ed'     
  !  open(41,file = temporaryfilename, form ='formatted',status = 'old')
  !  read(41,*)Ekin(ineq)
  !  close(41)
  !  print*,"Ekin =", temporaryfilename, Ekin(ineq)
  !enddo
  
  !Epot
  do ineq=1,Nlat
    WRITE(temporaryfilename,'(a16,i1.1,a5,i1.1,a3)')'energy_last_path',numfiles,'_site',ineq,'.ed'     
    open(41,file = temporaryfilename, form ='formatted',status = 'old')
    read(41,*)temporaryreal,Epot(ineq)
    close(41)
  !  print*,"Epot =", temporaryfilename, Epot(ineq)
  enddo

  !Ekin(3) = 0.5*(Ekin(1) + Ekin(2))
  Ekin(3) = Eout(1)
  Epot(3) = 0.5*(Epot(1) + Epot(2))
  EWint =  W0* (n0(1) * n0(2) ) * 0.5
  Emu = - 0.5*(n0(1)+n0(2))*(xmubar+W0+Uloc(1)/float(2))
  Etot = Ekin(3) + Epot(3) + EWint
  OmegaTot = Etot + Emu
 
  !docc
  !do ineq=1,2
  !  WRITE(temporaryfilename,'(a21,i1.1,a5,i1.1,a3)')'observables_last_path',numfiles,'_site',ineq,'.ed'     
  !  open(41,file = temporaryfilename, form ='formatted',status = 'old')
  !  read(41,*)temporaryreal,Docc(ineq)
  !  close(41)
  !  print*,"Docc =", temporaryfilename, Docc(ineq)
  !enddo

  !file 2 - energy
  write(37,'(I3.1,10(F15.10))')itertemp, &
      Uloc(1),W0,xmubar,1.d0/beta, &
      Ekin(3), &!Ekin(2),Ekin(1),&
      Epot(3), &!Epot(2),Epot(1),&
      EWint, Emu, Etot, OmegaTot
  call flush(37)
  
  !change of mu/W0/U/... loop parameters
  itertemp = itertemp+1
  looplocvalue = looplocvalue + loopincrement

  enddo !finish of mu/W0/U/... loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  close(35)
  close(36)
  close(37)
      
  ! checking file
  WRITE(temporaryfilename,'(a,i1.1,a)') 'AA_finish_',numfiles,'.dat'
  open(44,file=temporaryfilename,form='formatted',status='unknown')
  write(44,'(1(I3.3))')itertemp 
  write(44,'(1(A30))')'Finish of the CDW-EHM program'  
  close(44) 
   
  write(*,*)"Finish of the CDW-EHM program"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine build_Hbethe_lat_2x2 (dens)
    real(8),dimension(Nk) :: ek,wtk
    real(8),dimension(2),intent(in)  :: dens
    real(8)               :: de
    integer               :: i
    ek  = linspace(-Wband,Wband,Nk,mesh=de)
    Wtk = dens_bethe(ek,wband)*de
    BetheDOS        = zero
    BetheEk         = 0d0
    BetheDOS(:)     = Wtk(:)
    BetheEkImp(1,2,:)  = ek
    BetheEkImp(2,1,:)  = ek
    do i=1,Nk
      BetheEkImp(1,1,i)  = W0*dens(2) 
      BetheEkImp(2,2,i)  = W0*dens(1)
    enddo
  end subroutine build_Hbethe_lat_2x2
  
  subroutine build_Hbethe()
    real(8),dimension(Nk) :: ek,wtk
    real(8)               :: de
    integer               :: i
    ek  = linspace(-Wband,Wband,Nk,mesh=de)
    Wtk = dens_bethe(ek,wband)*de
    BetheDOS        = zero
    BetheEk         = 0d0
    BetheDOS(:)     = Wtk(:)
    BetheEk(1,2,:)  = ek
    BetheEk(2,1,:)  = ek
  end subroutine build_Hbethe

  subroutine build_Hloc (dens)
    real(8),dimension(2),intent(in)  :: dens
    Hloc = 0d0
    Hloc(1) = W0*dens(2)
    Hloc(2) = W0*dens(1)
  end subroutine build_Hloc

  ! calculate Gloc for Nlat=2 case 'by hand' (for check)
  subroutine get_gloc_CDW (Hk,Wtk,Gmats_loc,Greal_loc,Smats,Sreal,iprint)
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                     :: Wtk(size(Hk,1)) ![Nk]
    complex(8),intent(inout)    :: Gmats_loc(2,Lmats)
    complex(8),intent(inout)    :: Greal_loc(2,Lreal)
    complex(8),intent(inout)    :: Smats(2,Lmats)
    complex(8),intent(inout)    :: Sreal(2,Lreal)
    complex(8)                  :: zeta_mats(2,Lmats)
    complex(8)                  :: zeta_real(2,Lreal)
    complex(8)                  :: Gkmats(Lmats)
    complex(8)                  :: Gkreal(Lreal)
    real(8) :: wm(Lmats), wr(Lreal)
    integer :: i,iprint,Nk
    Nk=size(Hk)
    do ineq=1,2
        xmusub(ineq) = xmu0 - W0*n0(3-ineq)
    end do
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)

    zeta_mats(1,:) = xi*wm(:) + xmusub(1) - Smats(1,:)
    zeta_mats(2,:) = xi*wm(:) + xmusub(2) - Smats(2,:)
    
    zeta_real(1,:) = wr(:) + xi*eps + xmusub(1) - Sreal(1,:) 
    zeta_real(2,:) = wr(:) + xi*eps + xmusub(2) - Sreal(2,:)
    
    Gkmats(:) = zero
    Gkreal(:) = zero
    do i=1,Nk
       Gkmats(:) = Gkmats(:) + Wtk(i)/(zeta_mats(1,:) * zeta_mats(2,:) - Hk(i)**2) 
       Gkreal(:) = Gkreal(:) + Wtk(i)/(zeta_real(1,:) * zeta_real(2,:) - Hk(i)**2) 
    enddo
	
    Gmats_loc(1,:) = zeta_mats(2,:) * Gkmats(:)
    Gmats_loc(2,:) = zeta_mats(1,:) * Gkmats(:)
    
    Greal_loc(1,:) = zeta_real(2,:) * Gkreal(:)
    Greal_loc(2,:) = zeta_real(1,:) * Gkreal(:)

    if (iprint .eq. 1) then
       call splot("Gloc_iw_site1.ed",wm,Gmats_loc(1,:))
       call splot("Gloc_iw_site2.ed",wm,Gmats_loc(2,:))
       !call splot("Gloc_realw_site1.ed",wr,-dimag(Greal_loc(1,:))/pi,dreal(Greal_loc(1,:)))
       !call splot("Gloc_realw_site2.ed",wr,-dimag(Greal_loc(2,:))/pi,dreal(Greal_loc(2,:)))
       call splot("Gloc_realw_site1.ed",wr,dimag(Greal_loc(1,:)),dreal(Greal_loc(1,:)))
       call splot("Gloc_realw_site2.ed",wr,dimag(Greal_loc(2,:)),dreal(Greal_loc(2,:)))
    
    endif

  end subroutine get_gloc_CDW

end program ed_cdw_loop



