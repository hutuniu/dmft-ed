program ed_ehm_loop

  USE DMFT_ED
  USE DMFT_TOOLS
  USE SCIFOR
  implicit none
  
  integer                                       :: iloop
  logical                                       :: converged
  real(8)                                       :: wband,W0
  character(len=16)                             :: finput
  real(8)                                       :: wmixing
  integer                                       :: Nk
  real(8)                                       :: dens,dens_old,n0(2)
  real(8)					:: Ekin,Epot,Emu,EWint,Etot,OmegaTot
  real(8)					:: temporaryreal, error
  real(8)					:: xmubar, xmu0
  !Bath:
  integer                                       :: Nb(2)
  real(8),allocatable                           :: Bath(:,:,:),Bath_prev(:,:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Delta ![Nlat][Nspin][Nspin][Norb][Norb][Nfreq]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_loc
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_loc
  !Hamiltonian
  real(8),allocatable                           :: BetheDOS(:)
  real(8),allocatable                           :: BetheEk(:,:,:)
  complex(8),allocatable                        :: Hloc(:,:,:,:,:)
  character(len=20)                             :: tmp_suffix
  
  ! tag for printing, tag for distinguish files,
      
  real(8)  		:: loopincrement,loopfinish,loopstart
  real(8)   		:: loop_min,loop_max,looplocvalue
  integer  		:: itertemp,loopnumberofsteps,looptag 
  integer 		:: printingtag, numfiles 
  !real(8) 		:: printingtag_r, numfiles_r
  character(len=80) 	:: temporaryfilename
  
  !WRITE(temporaryfilename,'(a8,i1.1,a3)')'inputED_',numfiles,'.in'
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  !
  call parse_input_variable(printingtag,"printingtag",finput,default=0)
  call parse_input_variable(numfiles,"numfiles",finput,default=5)  
  call parse_input_variable(loopstart,"LoopStart",finput,default=0.0d0)
  call parse_input_variable(loopfinish,"LoopFinish",finput,default=0d0) 
  call parse_input_variable(loopnumberofsteps,"LoopNumberofSteps",finput,default=0)
  call parse_input_variable(looptag,"LoopTag",finput,default=1)
  !
  call parse_input_variable(Nk,"Nk",finput,default=500)
  call parse_input_variable(dens,"dens",finput,default=1.0d0)
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(W0,"W0",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0) 
  !
  call ed_read_input(trim(finput))
  
  ! starting two main files
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_it_',numfiles,'.dat'
  open(35,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')      
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_par_',numfiles,'.dat'
  open(36,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')
  WRITE(temporaryfilename,'(a,i1.1,a)')'energy_vs_par_',numfiles,'.dat'
  open(37,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')

  xmubar = xmu ! from starting file

  Nlat=1   
  if(Norb/=1)stop "This drivers requires Norb==1"
  if(Nspin/=1)stop "This drivers requires Nspin==1"

  allocate(BetheDOS(Nk),BetheEk(Nlat,Nlat,Nk))
  call build_Hbethe2x2 
  !call splot("testDOSbethe.dat",BetheEk(1,1,:),BetheDOS(:))


  !Allocate Weiss Field:
  allocate(Delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_loc(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_loc(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))  

  !setup solver
  Nb=get_bath_size()
  print*,"Nb = ",Nb(1)," ",Nb(2)
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_prev(Nlat,Nb(1),Nb(2)))
 
     write(tmp_suffix,'(I1.1)') numfiles
     ed_file_suffix="_path"//trim(tmp_suffix)
     call ed_init_solver(Bath(1,:,:))

!     'if loopnumberofsteps ==0 only one eveluation
      if(loopnumberofsteps.eq.0) then
       loopincrement = 10000000
      else
       loopincrement=(loopfinish-loopstart)/loopnumberofsteps
      endif  

if (looptag.eq.1) then
       WRITE(temporaryfilename,'(a)')'mu -- LOOP'
      else if  (looptag.eq.2) then
       WRITE(temporaryfilename,'(a)')'W0 -- LOOP'
      else if  (looptag.eq.3) then
       WRITE(temporaryfilename,'(a)')'U -- LOOP'
      else
       WRITE(temporaryfilename,'(a)')'WRONG TAG !!!'
      end if 
!      
! preparing for loop      
! 'determination of max and min values in looop - range of the loop'
!      
      if (loopincrement > 0 ) then
         loop_max =   loopfinish
         loop_min =   loopstart
      else
         loop_max =   loopstart
         loop_min =   loopfinish
      end if
! 'looplocvalue is current value of the loop parameter'
! 'we set a proper start 
      looplocvalue = loopstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !'do the loop when current value is between min and max'
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
      
      xmu0 = xmubar + W0
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     
     call start_loop(iloop,nloop,"DMFT-loop")
     bath_prev=bath
     dens_old=dens

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     ! solve the impurities on each inequivalent site:
        xmu=xmu0-W0*dens
        
        write(tmp_suffix,'(I1.1)') numfiles
        ed_file_suffix="_path"//trim(tmp_suffix)

        call ed_solve(bath(1,:,:))
        dens = ed_dens(1)
        call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(1,:,:,:,:,:))
        call ed_get_gimp_matsubara(Gmats(1,:,:,:,:,:))
        call ed_get_gimp_real(Greal(1,:,:,:,:,:))

!        xmu=xmu0-W0*dens
        
        !if printing for every step in par loop       
        if (printingtag .eq. 1) then
          ED_VERBOSE = 0
          write(tmp_suffix,'(I1.1,A5,I3.3)')numfiles,"_iter",itertemp
          ed_file_suffix="_path"//trim(tmp_suffix)
        endif

	  call ed_get_gloc(one*BetheEk(1,1,:),BetheDOS,&
             Gmats_loc(1,1,1,1,1,:),&
             Greal_loc(1,1,1,1,1,:),&
             Smats(1,1,1,1,1,:),&
             Sreal(1,1,1,1,1,:),iprint=1)
          
     ED_VERBOSE = 6

         
     
     print*,""
     print'(A4,F12.10,A6,F12.8)',"n = ",dens," mu = ", xmu
     print*,""

!        xmu=xmu0-W0*dens

	!Get Delta and fit it
        write(tmp_suffix,'(I1.1)') numfiles
        ed_file_suffix="_path"//trim(tmp_suffix)

        !Delta(1,:,:,:,:,:) = (wband**2d0)/4d0*Gmats(1,:,:,:,:,:)!Gmats_loc(3-ineq,:,:,:,:,:)
        call ed_get_weiss(Gmats_loc(1,:,:,:,:,:),Smats(1,:,:,:,:,:),Delta(1,:,:,:,:,:),iprint=1)
        call ed_chi2_fitgf(Delta(1,1,1,:,:,:),bath(1,:,:),ispin=1)


     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     !if(iloop>1)dens = wmixing*dens + (1.d0-wmixing)*dens_old

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,1,1,:),dmft_error,nsuccess,nloop)

     !if(nread/=0.d0)call search_chemical_potential(ed_dens(1),xmu0,converged)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !'results  of the interation - to file'  
      n0(1) = dens
      n0(2) = dens

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
 
!     xmu = xmu0 - W0*dens

     write(tmp_suffix,'(I1.1)') numfiles
     ed_file_suffix="_path"//trim(tmp_suffix)

     call ed_kinetic_energy(one*BetheEk(1,1,:),BetheDOS,Smats(1,1,1,1,1,:))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !WRITE(temporaryfilename,'(a9)')'error.err'
  !open(40,file = temporaryfilename, form ='formatted',status = 'old')
  
  !'the most important results of the single loop'
  !file 1 - order parameters
  write(36,'(I3.1,8(F15.10),1(I4.1))')itertemp, &
      Uloc(1),W0,xmubar,1.d0/beta, &
      n0(1),n0(2),0.5d0*(n0(1)+n0(2)),abs(0.5d0*(n0(1)-n0(2))), &
      iloop
  call flush(36)
  
  !get informations from files
  !Ekin

    WRITE(temporaryfilename,'(a17,i1.1,a3)')'kinetic_last_path',numfiles,'.ed'     
    open(41,file = temporaryfilename, form ='formatted',status = 'old')
    read(41,*)Ekin
    close(41)
  !  print*,"Ekin =", temporaryfilename, Ekin(ineq)

  
  !Epot

    WRITE(temporaryfilename,'(a16,i1.1,a3)')'energy_last_path',numfiles,'.ed'     
    open(41,file = temporaryfilename, form ='formatted',status = 'old')
    read(41,*)temporaryreal,Epot
    close(41)
  !  print*,"Epot =", temporaryfilename, Epot(ineq)


  EWint =  W0* (n0(1) * n0(2) ) * 0.5
  Emu = - 0.5*(n0(1)+n0(2))*(xmubar+W0+Uloc(1)/float(2))
  Etot = Ekin + Epot + EWint
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
      Ekin, &!Ekin(2),Ekin(1),&
      Epot, &!Epot(2),Epot(1),&
      EWint, Emu, Etot, OmegaTot
  call flush(37)
  
  itertemp = itertemp+1
  looplocvalue = looplocvalue + loopincrement

  enddo !finish of mu/W/U/T loop

  close(35)
  close(36)
  close(37)
      
  !'checking file'
  WRITE(temporaryfilename,'(a,i1.1,a)') 'AA_finish_',numfiles,'.dat'
  open(44,file=temporaryfilename,form='formatted',status='unknown')
  write(44,'(1(I3.3))')itertemp    
  close(44) 
   
  write(*,*)"Finish of the program"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine build_Hbethe2x2
    real(8),dimension(Nk) :: ek,wtk
    real(8)               :: de
    ek  = linspace(-Wband,Wband,Nk,mesh=de)
    Wtk = dens_bethe(ek,wband)*de
    BetheDOS        = zero
    BetheEk         = 0d0
    BetheDOS(:)     = Wtk(:)
    BetheEk(1,1,:)  = ek
!    BetheEk(2,2,:)  = ek
  end subroutine build_Hbethe2x2

end program ed_ehm_loop



