!MULTI-PURPOSE DRIVER FOR DMFT CALCULATIONS
!READ multi-band H(k)  from file or generate 
!a Bethe lattice multi-band H(e) using input 
!variables.
program ed_lda
  USE DMFT_ED
  USE SCIFOR
  implicit none
  !GLOBAL VARIABLE USED HERE:
  integer                :: iloop
  integer                :: Nkpts
  integer                :: Norb_d
  integer                :: Norb_p
  integer                :: Norb_tot
  integer                :: N_ineq
  integer                :: Ntype
  character(len=32)      :: Finput
  character(len=32)      :: hkfile
  character(len=10)      :: hktype
  logical                :: converged
  !HAMILTONIAN 
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !HOPPING PARAMETERS TO BE READ FROM FILE
  real(8),allocatable    :: ts(:)
  !FIX DENSITY VARIABLE
  real(8)                :: nobj
  !BATH:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !THE LOCAL HYBRIDIZATION FUNCTION:
  complex(8),allocatable :: Delta(:,:,:)
  complex(8),allocatable :: Smats(:,:,:)


  !Parse general input file and variables
  call parse_cmd_variable(finput,'FINPUT',default='inputED.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default='hkfile.in')
  call parse_input_variable(hktype,"HKTYPE",finput,default='lda')
  call parse_input_variable(Nkpts,"NKPTS",finput,default=1000)
  call parse_input_variable(Norb_p,"NORB_P",finput,default=0)
  call parse_input_variable(N_ineq,"N_INEQ",finput,default=0)
  call parse_input_variable(ntype,"NTYPE",finput,default=0)
  !
  call ed_read_input(trim(finput))

  !SETUP SOLVER
  !READ HLOC
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)


  !read Hkfile. Invoked AFTER init_ed_solver, because 
  !needs to read Hloc from file if exist first
  call read_hkfile(trim(hkfile))


  !Allocate Weiss Field:
  allocate(delta(Norb_tot,Norb_tot,Lmats))



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(xmu,nobj,converged)
     call end_loop
  enddo


  allocate(Smats(Norb_tot,Norb_tot,Lmats))
  Smats=zero
  Smats(:Norb_d,:Norb_d,:Lmats) = impSmats(1,1,:Norb_d,:Norb_d,:Lmats)
  call ed_kinetic_energy(Smats,Hk,dos_wt)
  deallocate(Smats)


contains


  !-------------------------------------------------------------------------------
  ! READ THE MODEL HAMILTONIAN OR GENERATE A BETHE LATTICE EQUIVALENT
  !-------------------------------------------------------------------------------
  subroutine read_hkfile(file)
    character(len=*) :: file
    logical          :: ioexist,fbethe
    character(len=1) :: cchar
    integer          :: Norb_d_,Norb_p_,N_ineq_
    integer          :: i,j,ik,iorb,jorb
    real(8)          :: de,e,kx,ky,kz,foo
    Norb_d = Norb
    inquire(file=trim(file),exist=ioexist)
    if(.not.ioexist)stop "HKfile does not exist or can not be read."
    open(50,file=file,status='old')
    !READ HEADER OF THE H(k) FILE:
    read(50,*)cchar          !first character should be a comment line:
    if(cchar/='#')stop "read_hkfile error: file does not start with a comment #"
    backspace(50)
    read(50,*)cchar,Nkpts,Norb_d_,Norb_p_,N_ineq_
    write(*,*)"Nkpts    =",Nkpts
    write(*,*)"Norb_d   =",Norb_d_
    write(*,*)"Norb_p   =",Norb_p_
    write(*,*)"N_ineq   =",N_ineq_
    if(Norb_d_/=Norb_d)stop "read_hkfile error: read  Norb_d != Norb_d from input file"
    if(Norb_p_/=Norb_p)stop "read_hkfile error: read  Norb_p != Norb_p from input file"
    if(N_ineq_/=N_ineq)stop "read_hkfile error: read  N_ineq != N_ineq from input file"
    Norb_tot=Norb_d + Norb_p
    write(*,*)"Norb_tot =",Norb_tot
    if(Ntype > Norb_tot)stop "read_hkfile error: required Ntype > Norb_tot"
    !
    allocate(Hk(Norb_tot,Norb_tot,Nkpts))
    allocate(dos_wt(Nkpts))

    if(trim(hktype)=='lda')then
       !
       write(*,*)"Reading LDA H(k) from file:"
       do ik=1,Nkpts
          read(50,"(3(F10.7,1x))")kx,ky,kz
          do iorb=1,Norb_tot
             read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Norb_tot)
          enddo
       enddo
       close(50)
       Hloc(1,1,:Norb,:Norb) = sum(Hk(:Norb,:Norb,:),dim=3)/dble(Nkpts)
       dos_wt=1.d0/dble(Nkpts)
       !
    elseif(trim(hktype)=='bethe')then
       !
       write(*,*)"Reading t_hopping for the Bethe lattice:"
       allocate(ts(Norb_tot))
       do iorb=1,Norb_tot
          read(50,*)ts(iorb)
          write(*,*)ts(iorb)
       enddo
       close(50)
       !
       de=2.d0/dble(Nkpts)
       do ik=1,Nkpts
          e = -1.d0 + dble(ik-1)*de
          Hk(1:Norb_d,1:Norb_d,ik) = Hloc(1,1,1:Norb,1:Norb) !set local part
          do iorb=1,Norb_tot
             Hk(iorb,iorb,ik) = Hk(iorb,iorb,ik) - 2.d0*ts(iorb)*e
          enddo
          dos_wt(ik)=dens_bethe(e,1.d0)*de
       enddo
    else
       stop "read_hkfile error: this value of hktype is not supported [lda,bethe]"
    endif
    print*,"Hloc:"
    call print_Hloc(Hloc)
  end subroutine read_hkfile






  !-------------------------------------------------------------------------------
  ! UPDATE THE WEISS/HYBRIDIZATION FUNCTION. PERFORM SELF-CONSISTENCY
  !-------------------------------------------------------------------------------
  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8),dimension(Norb_tot,Norb_tot) :: fg,gdelta
    complex(8),allocatable,dimension(:,:,:) :: Gloc,Sigma,Zeta
    real(8)                                 :: wm(Lmats),wr(Lreal),dens(Norb_tot)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    Delta=zero
    !
    !MATSUBARA AXIS:
    allocate(Gloc(Norb_tot,Norb_tot,Lmats))
    allocate(Sigma(Norb_tot,Norb_tot,Lmats))
    allocate(Zeta(Norb_tot,Norb_tot,Lmats))
    Gloc  = zero
    Sigma = zero
    Zeta  = zero
    forall(i=1:LMats,iorb=1:Norb_d,jorb=1:Norb_d)Sigma(iorb,jorb,i)=impSmats(1,1,iorb,jorb,i)
    Zeta  = -Sigma
    forall(i=1:Lmats,iorb=1:Norb_tot)Zeta(iorb,iorb,i)=xi*wm(i)+xmu+zeta(iorb,iorb,i)
    do i=1,Lmats
       fg=zero
       do ik=1,Nkpts
          fg = fg + invert_gk(zeta(:,:,i),Hk(:,:,ik),Norb_tot)*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
       !
       if(Norb_tot>1)then
          call matrix_inverse(fg)
       else
          fg(1,1)=one/fg(1,1)
       endif
       if(cg_scheme=='weiss')then
          gdelta = fg + Sigma(:,:,i)
          call matrix_inverse(gdelta)
       else
          gdelta = zeta(:,:,i) - hloc(1,1,:,:) - fg
       endif
       delta(:,:,i)=gdelta(1:Norb_d,1:Norb_d)
    enddo
    !print
    do iorb=1,Norb_tot
       do jorb=1,Norb_tot
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
          call splot("Delta_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
       enddo
    enddo
    dens(1:Norb_d)=ed_dens(1:Norb_d)
    if(Norb_p > 0)then
       do iorb=1,Norb_p
          dens(Norb_d + iorb) =  get_density_fromFFT(Gloc(Norb_d + iorb,Norb_d + iorb,:),beta)
          write(*,"(A,F25.18)")"np      =",dens(Norb_d + iorb)
       enddo
    endif
    write(*,"(A,F25.18)")"n_total =",sum(dens)
    open(100,file="ndensity_all.ed",position='append')
    write(100,"(100F25.18)")(dens(iorb),iorb=1,Norb_tot),sum(dens)
    close(100)
    open(100,file="ndensity_last.ed")
    write(100,"(100F25.18)")(dens(iorb),iorb=1,Norb_tot),sum(dens)
    close(100)
    deallocate(Gloc,Sigma,Zeta)
    !
    !
    !REAL-AXIS:
    allocate(Gloc(Norb_tot,Norb_tot,Lreal))
    allocate(Sigma(Norb_tot,Norb_tot,Lreal))
    allocate(Zeta(Norb_tot,Norb_tot,Lreal))
    Gloc  = zero
    Sigma = zero
    Zeta  = zero
    forall(i=1:Lreal,iorb=1:Norb_d,jorb=1:Norb_d)Sigma(iorb,jorb,i)=impSreal(1,1,iorb,jorb,i)
    Zeta  = -Sigma
    forall(i=1:Lreal,iorb=1:Norb_tot)Zeta(iorb,iorb,i)=dcmplx(wr(i),eps)+xmu+zeta(iorb,iorb,i)
    do i=1,Lreal
       fg=zero
       do ik=1,Nkpts
          fg = fg + invert_gk(zeta(:,:,i),Hk(:,:,ik),Norb_tot)*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
    enddo
    do iorb=1,Norb_tot
       do jorb=1,Norb_tot
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi,dreal(gloc(iorb,jorb,:)))
       enddo
    enddo
    deallocate(Gloc,Sigma,Zeta)
    !
    !
    !SET GOAL DENSITY TO FIX
    if(Ntype==0)then
       nobj = sum(dens)
    else
       nobj = dens(Ntype)
    endif
  end subroutine get_delta







  !-------------------------------------------------------------------------------
  ! INVERT THE MATRIX (z-Hk) USED TO EVALUATE G_LOCAL
  !-------------------------------------------------------------------------------
  function invert_gk(zeta,Hk,Ndim) result(invg)
    integer                         :: Ndim
    complex(8),dimension(Ndim,Ndim) :: zeta,Hk
    complex(8),dimension(Ndim,Ndim) :: invg
    complex(8)                      :: h11,h22,h12
    invg=zeta-Hk
    if(Ndim==1)then
       invg(1,1)=one/invg(1,1)
    elseif(Ndim==2)then
       h11 = zeta(1,1) - hk(1,1)
       h22 = zeta(2,2) - hk(2,2)
       h12 = zeta(1,2) - hk(1,2)
       invg(1,1) = 1.d0/(h11 - abs(h12)**2/h22)
       invg(2,2) = 1.d0/(h22 - abs(h12)**2/h11)
       invg(1,2) = -h12/(h11*h22 - abs(h12)**2)
       invg(2,1) = conjg(invg(1,2))
    else
       call matrix_inverse(invg)
    endif
  end function invert_gk






  !-------------------------------------------------------------------------------
  ! GET THE LOCAL DENSITY FROM G(iw)
  !-------------------------------------------------------------------------------
  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT


end program ed_lda



