program ed_wsm_3d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                       :: iloop,Lk,Nso
  logical                       :: converged
  !Bath:
  integer                       :: Nb
  real(8),allocatable           :: Bath(:),Bath_Prev(:)
  !The local hybridization function:
  complex(8),allocatable        :: Weiss(:,:,:,:,:)
  complex(8),allocatable        :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable        :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable        :: Hk(:,:,:),wsmHloc(:,:),sigmaWSM(:,:)
  real(8),allocatable           :: Wtk(:)
  integer,allocatable           :: ik2ix(:),ik2iy(:),ik2iz(:)
  !variables for the model:
  integer                       :: Nk,Nkpath
  real(8)                       :: e0,mh,lambda,bx,by,bz,BIA
  real(8)                       :: wmixing
  character(len=16)             :: finput
  character(len=32)             :: hkfile
  logical                       :: spinsym,getpoles
  complex(8),dimension(4,4)     :: Gamma1,Gamma2,Gamma3,Gamma5
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED_WSM.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=30)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(e0,"E0",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.5d0)
  call parse_input_variable(bx,"BX",finput,default=0.3d0)
  call parse_input_variable(by,"BY",finput,default=0.d0)
  call parse_input_variable(bz,"BZ",finput,default=0.d0)
  call parse_input_variable(BIA,"BIA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma3=kron_pauli( pauli_tau_x,-pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  !
  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(SigmaWSM(Nso,Nso))
  call set_sigmaWSM()           !this set sigma_WSM(0) to zero
  !
  !
  !
  !
  !
  !Buil the Hamiltonian on a grid or on path
  call build_hk(trim(hkfile))
  !
  !TESTS
  !Find out if it is a semimetal
  call is_weyl_brutal( Nso, so2j(Smats(:,:,:,:,1),Nso) )
  !call chern_retriever([-3.1415926535897931d0,-0.16755166848237124d0,-0.60737458069855776d0])
  !call chern_retriever([pi*0.134188,pi*0.325236,pi*0.5436])
  STOP
  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_Prev(Nb))
  call ed_init_solver(bath,j2so(wsmHloc))
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     !
     !
     ! retrieve the self-energies
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call build_z2_indices( so2j(Smats(:,:,:,:,1),Nso) ) !so2j transforms a Nspin:Nspin:Norb:Norb into a Nso:Nso matrix
     !
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=3)
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc))
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc))
     endif
     !
     !
     !
     !Fit the new bath, starting from the old bath + the supplied Weiss/Delta
     if(ed_mode=="normal")then
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
        call ed_chi2_fitgf(Weiss,bath,ispin=2)
     else
        call ed_chi2_fitgf(Weiss,bath)
     endif
     !
     !
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
     Bath_Prev=Bath
     !
     !
     !
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call end_loop
     !
  enddo

  
  ! compute the local gf:
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=3)
  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Wtk,Smats)


  !!Get 3d Bands from Top. Hamiltonian
  call solve_hk_topological( so2j(Smats(:,:,:,:,1),Nso) )

  !!Find out if it is a semimetal
  call is_weyl_brutal( Nso, so2j(Smats(:,:,:,:,1),Nso) )
  !call chern_retriever([0.0d0,0.0d0,0.0d0])
  !call chern_retriever([pi*0.134188,pi*0.325236,pi*0.5436])
contains



  !---------------------------------------------------------------------
  !PURPOSE: GET WSM HAMILTONIAN IN THE FULL BZ
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,k,ik,iorb,jorb,io,ispin
    integer                             :: ix,iy,iz
    real(8)                             :: kx,ky,kz
    real(8)                             :: foo
    integer                             :: unit
    !
    !get H(k) and solve the non-interacting problem along a path in 3d:
    call build_hk_path()
    !
    !Get H(k) in the BZ:    
    Lk=Nk**3
    !
    write(LOGfile,*)"Build H(k) for WSM:"
    write(*,*)"# of k-points     :",Lk
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(Wtk(Lk))
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hk,hk_weyl,Nso,[Nk,Nk,Nk])
    Wtk = 1d0/Lk
    if(present(file))call TB_write_hk(Hk,trim(file),Nso,&
         Nd=Norb,Np=1,Nineq=1,&
         Nkvec=[Nk,Nk,Nk])
    !
    !
    !   
    !GET LOCAL PART OF THE HAMILTONIAN
    if(allocated(wsmHloc))deallocate(wsmHloc)
    allocate(wsmHloc(Nso,Nso))
    wsmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(wsmHloc))<1.d-9)wsmHloc=0d0
    call TB_write_Hloc(wsmHloc)
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE WSM HAMILTONIAN ALONG A 3D PATH IN THE BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path(kpath_)
    integer                                :: i,j
    integer                                :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable     :: kpath
    character(len=64)                      :: file
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    !
    sigmaWSM=zero
    if(present(kpath_))then
       write(LOGfile,*)"Build H(k) WSM along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))      
       kpath=kpath_
       file="Eig_path.nint"
    else
       write(LOGfile,*)"Build H(k) WSM along the high-symmetry path:"
       Npts = 8
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_m1
       kpath(2,:)=kpoint_x2
       kpath(3,:)=kpoint_gamma
       kpath(4,:)=kpoint_x1
       kpath(5,:)=kpoint_m2
       kpath(6,:)=kpoint_r
       kpath(7,:)=kpoint_x3
       kpath(8,:)=kpoint_gamma
       file="Eigenbands.nint"
    endif
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(Wtk(Lk))
    !
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hk,hk_weyl,Nso,kpath,Nkpath)
    Wtk = 1d0/Lk
    !
    call TB_solve_model(hk_weyl,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "M","X","G","X1","A","R","Z","G"],&
         file="Eigenband.nint")
  end subroutine build_hk_path


  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    write(LOGfile,*)"Build H_TOP(k) WSM along the X-G-M-R-Z-A-G-Z path:"
    !
    Npts = 8
    Lk=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_m1
    kpath(2,:)=kpoint_x2
    kpath(3,:)=kpoint_gamma
    kpath(4,:)=kpoint_x1
    kpath(5,:)=kpoint_m2
    kpath(6,:)=kpoint_r
    kpath(7,:)=kpoint_x3
    kpath(8,:)=kpoint_gamma
    call set_sigmaWSM(sigma)
    call TB_solve_model(hk_weyl,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "M","X","G","X1","A","R","Z","G"],&
         file="Eig_Htop.ed")
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !WSM HAMILTONIAN:
  !--------------------------------------------------------------------!
  subroutine set_SigmaWSM(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmaWSM = zero;if(present(sigma))sigmaWSM=sigma
  end subroutine set_SigmaWSM

  function hk_weyl(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_weyl: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - e0*(cos(kx) + cos(ky) + cos(kz)) )*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y +&
         by*pauli_tau_y + bz*pauli_tau_z
    Hk(3:4,3:4) = conjg( &
         (Mh - e0*(cos(-kx) + cos(-ky) + cos(-kz)) )*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y +&
         by*pauli_tau_y - bz*pauli_tau_z)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x -xi*BIA*pauli_tau_y + bx*pauli_tau_z
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x +xi*BIA*pauli_tau_y + bx*pauli_tau_z
    !
    !add the SigmaWSM term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaWSM)
    !
  end function hk_weyl

  !--------------------------------------------------------------------!
  !WSM REGION FINDER:
  !--------------------------------------------------------------------!

  subroutine is_weyl_brutal(N,sigma)
    integer                                 :: weyl_index, unit,test_brutal=0,i,j,k,mash_thickness
    real(8)                                    :: step
    complex(8),dimension(Nso,Nso)           :: sigma(Nso,Nso)
    integer                                 :: N
    real(8),dimension(:),allocatable        :: kpoint
    real(8),dimension(N)                    :: eval
    complex(8),dimension(N,N)               :: ham

    if(N/=4)stop "hk_weyl: error in N dimensions"
    !
    !
    weyl_index=0
    call set_sigmaWSM(sigma)
    !
    mash_thickness=10
    step=2*pi/(mash_thickness*Nk)
    !
    !write(*,*) "Starting Weyl point search"
    !
    kpoint=[-pi,-pi,-pi]
    xloop: do i=0,mash_thickness*Nk-1
      kpoint(1)=-pi+step*i
      yloop: do j=0,mash_thickness*Nk-1
        kpoint(2)=-pi+step*j
        zloop: do k=0,mash_thickness*Nk-1
          kpoint(3)=-pi+step*k
          test_brutal=test_brutal+1
          ham=hk_weyl(kpoint,N)
          call eigh(ham,Eval)
          if (abs(Eval(3))+abs(Eval(2)) .lt. 0.005) then
            weyl_index=1
            write(*,*) "----------------------------------------"
            write(*,*) "Brutal iterations: found Weyl point at ",kpoint, "value is ", abs(Eval(3))
            call chern_retriever(kpoint,step)
            write(*,*) "----------------------------------------"
            !exit xloop
          end if
        end do zloop
      end do yloop
    end do xloop
    unit=free_unit()
    open(unit,file="is_weyl.ed")
    write(unit,*)weyl_index
    close(unit)
  end subroutine is_weyl_brutal


  subroutine retrieve_third_band(n,x,fvec,iflag)
    integer                                 :: iflag,n, INFO
    real(8),dimension(N)                    :: x
    real(8),dimension(N)                    :: fvec
    complex(8),dimension(Nso,Nso)           :: ham
    complex(8),dimension(Nso,Nso)           :: sigma
    real(8),dimension(Nso)                  :: eval
    !
    sigma=so2j(Smats(:,:,:,:,1),Nso)
    call set_sigmaWSM(sigma)
    ham=hk_weyl(x,Nso)
    call eigh(ham,Eval)
    fvec=[Eval(3),0.0d0,0.0d0]
  end subroutine retrieve_third_band


  !--------------------------------------------------------------------!
  !CHERN NUMBER EVALUATION  (NEEDS FIXING AND OPTIMIZATION)
  !--------------------------------------------------------------------!

  subroutine chern_retriever(kpoint,cubesize)   !integrates Berry flux on a cubic surface around the point kpoint
    real(8)                         :: z2
    complex(8)                      :: phase
    real(8),dimension(3)            :: kpoint,kpoint_
    real(8)                         :: e,cubesize
    integer                         :: unit,TEST,perm,j,q,N,face_indx,face_sign,side_indx,side_sign,path_indx,border,varying_indx,run_direction
    integer,dimension(6,6)          :: permutations
    integer,dimension(4)            :: versor
    complex(8),dimension(4,4)       :: BlochOld,BlochNew
    complex(8),dimension(2,2)       :: OverlapMatrix
    real(8),dimension(4)            :: Eigval
    do TEST=2,2
    e=cubesize/(2*TEST)
    N=500
    !write(*,*) "Testing chirality"
    !write(*,*) "Integrating on a cube around ",kpoint
    z2=0.0d0
    !
    permutations(1,:)=[3,1,2,-1,-2,1] ![index fixing the face, 4*index fixing the path side, coordinate varying on side labelled by #5]
    permutations(2,:)=[-3,-2,-1,2,1,2] 
    permutations(3,:)=[1,2,3,-2,-3,2]
    permutations(4,:)=[-1,-3,-2,3,2,3]
    permutations(5,:)=[2,3,1,-3,-1,3]
    permutations(6,:)=[-2,-1,-3,1,3,1]
    versor=[1,-1,-1,1]
    !
    do perm=1,6
      face_indx=ABS(permutations(perm,1))   !which plane am I parallel to?
      face_sign=SIGN(1,permutations(perm,1)) !top-botton, left-right face, now I have selected one
      phase=1.0d0   !do path around this face
      !
      kpoint_=zero  
      BlochOld=zero
      BlochNew=zero
      OverlapMatrix=zero
      !
      kpoint_(face_indx)=kpoint(face_indx)+e*face_sign !fix one coordinate of the running point
      !
      do border=2,5 !run along border
        !
        side_indx=ABS(permutations(perm,border))      !which direction am I parallel to?
        side_sign=SIGN(1,permutations(perm,border))   !top-bottom, left-right side, now I've selected one
        varying_indx=ABS(permutations(perm,border+1)) !the point runs along the other coordinate
        run_direction=face_sign*versor(border-1)                !direction to run along to have a closed path
        !
        kpoint_(side_indx)=kpoint(side_indx)+e*side_sign !fix another coordinate of the running point
        !
        do path_indx=1,2*N !run on one side of the border
          !
          if (all(BlochOld(:,1) .eq. 0.0d0)) then !if I have yet to start, do one cycle
            kpoint_(varying_indx)=kpoint(varying_indx)-e*run_direction
            BlochOld=hk_weyl(kpoint_,4)
            call eigh(BlochOld,Eigval)  !store old eigenvectors
          endif
          kpoint_(varying_indx)=kpoint(varying_indx)-e*run_direction+(e*path_indx/N)*run_direction
          BlochNew=hk_weyl(kpoint_,4)
          call eigh(BlochNew,Eigval)  !store new eigenvectors
          do j=1,2
            do q=1,2
              OverlapMatrix(j,q)=dot_product(BlochOld(:,j),BlochNew(:,q))
            enddo
          enddo
          BlochOld=BlochNew !update eigenvectors
          phase=phase*det(OverlapMatrix)
        enddo !end run on one side of the border
      enddo !end run on the border
      z2=z2+IMAG(log(phase))!  !sum the phase to the face
    enddo !end run on face
    write(*,*) "Chirality is ",NINT(z2/(2*pi))
    unit=free_unit()
    open(unit,file="Chirality.ed",position="append")
    if (abs(z2)>0.5d0)then
      write(unit,'(4F16.9)')kpoint(1),kpoint(2),kpoint(3),z2/(2*pi)
    endif
    close(unit)
  enddo
  end subroutine chern_retriever






  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine build_z2_indices(sigma0)
    integer                                :: unit
    complex(8),dimension(Nso,Nso),optional :: sigma0(Nso,Nso)
    complex(8),dimension(Nso,Nso)          :: sigma0_
    !
    integer,dimension(4)                   :: z2
    !
    sigma0_=zero;if(present(sigma0))sigma0_=sigma0
    !Evaluate the Z2 index:
    !STRONG TI
    z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi,sigma0_)
    !WEAK TI
    !K=1: n_1=1, n_2,3=0,1
    z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    !K=2: n_2=1, n_1,2=0,1
    z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi,sigma0_)
    !k=3: n_3=1, n_1,2=0,1
    z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    unit=free_unit()
    open(unit,file="z2_invariant.ed")
    write(unit,*)z2
    close(unit)
  end subroutine build_z2_indices

  function z2_number(ktrims,sigma0) result(z2)
    real(8),dimension(:,:),intent(in)            :: ktrims
    complex(8),dimension(Nso,Nso)                :: sigma0
    complex(8),dimension(Nso,Nso,size(Ktrims,2)) :: Htrims
    complex(8),dimension(size(Ktrims,2))         :: Delta
    integer                                      :: z2
    integer                                      :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_weyl(Ktrims(:,itrim),Nso) + sigma0(:,:)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number



  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg,Nso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nso,Nso)               :: g
    integer                                     :: Nso,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nso,Nso)               :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so


end program ed_wsm_3d

