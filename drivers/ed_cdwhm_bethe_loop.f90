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
  
  ! starting two main files
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_it_',numfiles,'.dat'
  open(35,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')      
  WRITE(temporaryfilename,'(a,i1.1,a)')'dens_vs_par_',numfiles,'.dat'
  open(36,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')
  WRITE(temporaryfilename,'(a,i1.1,a)')'energy_vs_par_',numfiles,'.dat'
  open(37,file=temporaryfilename,form='formatted',status='unknown', Access = 'append')

  xmubar = xmu ! from starting file

  Nlat=2      !the two ineq. sublattices A,B
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
  do ineq=1,2
     !write(tmp_suffix,'(I4.4)') ineq
     !ed_file_suffix="_site"//trim(tmp_suffix)
     write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
     ed_file_suffix="_path"//trim(tmp_suffix)
     call ed_init_solver(Bath(ineq,:,:))
  enddo

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

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     ! solve the impurities on each inequivalent site:
     do ineq=1,2
        xmu=xmu0-W0*n0(3-ineq)
        
        !write(tmp_suffix,'(I4.4)') ineq
        !ed_file_suffix="_site"//trim(tmp_suffix)
        write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
        ed_file_suffix="_path"//trim(tmp_suffix)

        call ed_solve(bath(ineq,:,:))
        n0(ineq) = ed_dens(1)
        call ed_get_sigma_matsubara(Smats(ineq,:,:,:,:,:))
        call ed_get_sigma_real(Sreal(ineq,:,:,:,:,:))
        call ed_get_gimp_matsubara(Gmats(ineq,:,:,:,:,:))
        call ed_get_gimp_real(Greal(ineq,:,:,:,:,:))
     enddo
     
     do ineq=1,2
        xmu=xmu0-W0*n0(3-ineq)

        !write(tmp_suffix,'(I4.4)') ineq
        !ed_file_suffix="_site"//trim(tmp_suffix)
        !if printing for every step in par loop       
        if (printingtag .eq. 1) then
          ED_VERBOSE = 0
          write(tmp_suffix,'(I1.1,A5,I1.1,A5,I3.3)')numfiles,"_site",ineq,"_iter",itertemp
          ed_file_suffix="_path"//trim(tmp_suffix)
        endif
        
	  call ed_get_gloc(one*BetheEk(1,1,:),BetheDOS,&
             Gmats_loc(ineq,1,1,1,1,:),&
             Greal_loc(ineq,1,1,1,1,:),&
             Smats(ineq,1,1,1,1,:),&
             Sreal(ineq,1,1,1,1,:),iprint=1)
          ED_VERBOSE = 6

     enddo
     
     print*,""
     print'(A4,F8.6,A4,F8.6)',"n1 = ",n0(1)," n2= ",n0(2)
     print*,""

     do ineq=1,2
        xmu=xmu0-W0*n0(3-ineq)
        write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
        ed_file_suffix="_path"//trim(tmp_suffix)
        call ed_get_weiss(Gmats_loc(ineq,:,:,:,:,:),Smats(ineq,:,:,:,:,:),Delta(ineq,:,:,:,:,:),iprint=1)
        call ed_chi2_fitgf(Delta(ineq,1,1,:,:,:),bath(ineq,:,:),ispin=1)
     enddo

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

  do ineq=1,2
     xmu = xmu0 - W0*n0(3-ineq)
     !write(tmp_suffix,'(I1.1,A5,I3.3)')ineq,"_iter",itertemp
     !write(tmp_suffix,'(I1.1)') ineq
     !ed_file_suffix="_site"//trim(tmp_suffix)

     write(tmp_suffix,'(I1.1,A5,I1.1)') numfiles,"_site",ineq
     ed_file_suffix="_path"//trim(tmp_suffix)

     call ed_kinetic_energy(one*BetheEk(1,1,:),BetheDOS,Smats(ineq,1,1,1,1,:))
  enddo

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
  do ineq=1,2
    WRITE(temporaryfilename,'(a17,i1.1,a5,i1.1,a3)')'kinetic_last_path',numfiles,'_site',ineq,'.ed'     
    open(41,file = temporaryfilename, form ='formatted',status = 'old')
    read(41,*)Ekin(ineq)
    close(41)
  !  print*,"Ekin =", temporaryfilename, Ekin(ineq)
  enddo
  
  !Epot
  do ineq=1,2
    WRITE(temporaryfilename,'(a16,i1.1,a5,i1.1,a3)')'energy_last_path',numfiles,'_site',ineq,'.ed'     
    open(41,file = temporaryfilename, form ='formatted',status = 'old')
    read(41,*)temporaryreal,Epot(ineq)
    close(41)
  !  print*,"Epot =", temporaryfilename, Epot(ineq)
  enddo

  Ekin(3) = 0.5*(Ekin(1) + Ekin(2))
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
    BetheEk(2,2,:)  = ek
  end subroutine build_Hbethe2x2

end program ed_cdw_loop



