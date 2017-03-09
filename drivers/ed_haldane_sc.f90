program ed_haldane
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none

  integer                                         :: iloop,Lk,Nso,Nlso,Nlat
  logical                                         :: converged

  !Bath:
  integer                                         :: Nb
  real(8),allocatable,dimension(:,:)              :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)           :: halHloc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Hloc
  real(8),allocatable,dimension(:)                :: Wtk
  integer,allocatable,dimension(:)                :: ik2ix,ik2iy
  real(8)                                         :: a,a0,bklen

  real(8),dimension(2)                            :: d1,d2,d3
  real(8),dimension(2)                            :: a1,a2,a3
  real(8),dimension(2)                            :: bk1,bk2,pointK,pointKp

  !variables for the model:
  integer                                         :: Nk,Nkpath
  real(8)                                         :: ts,tsp,phi,delta,Mh,wmixing
  character(len=32)                               :: finput
  character(len=32)                               :: hkfile
  logical                                         :: spinsym
  !
  real(8),dimension(2)                            :: Eout
  real(8),allocatable,dimension(:)                :: dens
  !
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Sigma0,Self0



  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputHALDANE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputHALDANE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputHALDANE.conf",default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH","inputHALDANE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputHALDANE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  phi=phi*pi

  if(Nspin/=1.OR.Norb/=1)stop "Wrong setup from input file: Nspin=Norb=1"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso

  !see (ed_graphene for backupd version)
  !LATTICE BASIS:
  ! nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d2-d3                    !a*sqrt(3)[sqrt(3)/2,-1/2]
  a2 = d3-d1                    !a*sqrt(3)[-sqrt(3)/2,-1/2]
  a3 = d1-d2                    !a*sqrt(3)[0, 1]
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/sqrt(3d0)
  bk1=bklen*[ sqrt(3d0)/2d0 ,  1d0/2d0 ]
  bk2=bklen*[ sqrt(3d0)/2d0 , -1d0/2d0 ]

  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


  !Allocate Weiss Field:
  allocate(Weiss(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(Sigma0(Nlat,Nspin,Nspin,Norb,Norb));Sigma0=zero
  allocate(Self0(Nlat,Nspin,Nspin,Norb,Norb));Self0=zero


  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(halHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(Bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc,iprint=1)
     call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:,:),Nlat)
     call ed_get_sigma_real(Sreal(1,:,:,:,:,:,:),Nlat)
     call ed_get_self_matsubara(Smats(2,:,:,:,:,:,:),Nlat)
     call ed_get_self_real(Sreal(2,:,:,:,:,:,:),Nlat)
     Sigma0 = dreal(Smats(1,:,:,:,:,:,1))
     Self0  = dreal(Smats(2,:,:,:,:,:,1))

     ! compute the local gf:
     call dmft_gloc_matsubara_superc(Hk,Wtk,Gmats,Smats,iprint=4)


     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss_superc(Gmats,Smats,Weiss,Hloc,iprint=4)
     else
        call dmft_delta_superc(Gmats,Smats,Weiss,Hloc,iprint=4)
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     converged = check_convergence(Weiss(1,1,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  call dmft_gloc_realaxis_superc(Hk,Wtk,Greal,Sreal,iprint=4)
  
  !Eout = dmft_kinetic_energy(Hk,Wtk,Smats(1,:,:,:,:,:,:),Smats(2,:,:,:,:,:,:))
  !print*,Eout


  call build_EigenBands()
  call eval_Chern_number()


contains




  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                       :: h0,hx,hy,hz
    real(8)                       :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = -2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = -2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
  end function hk_haldane_model


  function hk_NambuHaldane_model(kpoint,N) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: N
    complex(8),dimension(N,N)       :: hk
    complex(8),dimension(Nlso,Nlso) :: hk11,hk22
    !
    if(2*Nlso/=N)stop "hk_NambuHaldane_model ERROR: 2*Nlso != N"
    !
    hk11  =     hk_haldane_model(kpoint,Nlso)  +  nnn2lso_reshape(Sigma0,Nlat,Nspin,Norb)
    hk22  =    -transpose(hk_haldane_model(-kpoint,Nlso)  +  nnn2lso_reshape(Sigma0,Nlat,Nspin,Norb))
    !
    Hk(1:Nlso,1:Nlso)               = hk11
    Hk(1:Nlso,Nlso+1:2*Nlso)        = nnn2lso_reshape(Self0,Nlat,Nspin,Norb)
    Hk(Nlso+1:2*Nlso,1:Nlso)        = nnn2lso_reshape(Self0,Nlat,Nspin,Norb)
    Hk(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = hk22
  end function hk_NambuHaldane_model




  


  !---------------------------------------------------------------------
  !PURPOSE: Get Haldane Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky    
    real(8),dimension(2)                  :: kvec
    write(LOGfile,*)"Build H(k) for Haldane model (now also Nobel Laureate):"
    Lk=Nk*Nk
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(2,Nlso,Nlso,Lk))
    allocate(Wtk(Lk))
    print*,"Build Hk(Nlso,Nlso) for the Graphene model"
    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          Hk(1,:,:,ik) = hk_haldane_model(kvec,Nlso)
          Hk(2,:,:,ik) = -transpose(hk_haldane_model(-kvec,Nlso))
       enddo
    enddo
    Wtk = 1d0/Lk
    !
    allocate(halHloc(Nlso,Nlso))
    halHloc = sum(Hk(1,:,:,:),dim=3)/Lk
    where(abs(dreal(halHloc))<1.d-9)halHloc=0d0
    call TB_write_hloc(halHloc)
    !
    !
  end subroutine build_hk





  subroutine build_EigenBands()
    real(8),dimension(4,2)                 :: kpath
    !
    KPath(1,:)=[0d0,0d0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_model(hk_NambuHaldane_model,2*Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,orange1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.dat")
  end subroutine build_EigenBands



  subroutine eval_Chern_number()
    complex(8),dimension(2*Nlso,2*Nlso,Lk) :: HkNambu ![2Nlso][2Nlso][Lk]
    integer                                :: ik,ix,iy
    real(8)                                :: kx,ky
    real(8),dimension(2)                   :: kvec
    real(8) :: chern
    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          Hknambu(:,:,ik) =  hk_NambuHaldane_model(kvec,2*Nlso)
       enddo
    enddo
    call get_Chern_number(HkNambu,[Nk,Nk],2,Nk/pi2*Nk/pi2,Chern)
  end subroutine eval_Chern_number

  
  ! calcola il numero di chern di un generico stato dipendente da k con il metodo di Resta
  subroutine Get_Chern_number(Hk,Nkvec,Noccupied,one_over_area,Chern)
    complex(8),intent(in),dimension(:,:,:)    :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)           :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,intent(in)                        :: Noccupied
    real(8),intent(in)                        :: one_over_area
    real(8),intent(out)                       :: Chern
    !
    integer                                   :: Nlso
    integer                                   :: Nktot
    integer                                   :: Nkx,Nky
    integer                                   :: ikx,iky
    integer                                   :: ikxP,ikyP
    integer                                   :: ik,iocc,i,j
    complex(8),dimension(:,:),allocatable     :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable          :: Eigval ![Nlso]
    complex(8),dimension(:,:),allocatable     :: Gmat
    complex(8),dimension(:,:,:,:),allocatable :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                   :: Ulink
    real(8),dimension(:,:),allocatable        :: BerryCurvature
    real(8)                                   :: berry_phase
    integer                                   :: unit
    !
    Nlso  = size(Hk,1)
    Nktot = size(Hk,3)
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    call assert_shape(Hk,[Nlso,Nlso,Nktot],"Get_Chern_NUmber_NEW","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number_NEW: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nlso,Nlso))
    allocate(Eigval(Nlso))
    allocate(BlochStates(Nkx,Nky,Noccupied,Nlso))
    allocate(BerryCurvature(Nkx,Nky))
    allocate(Gmat(Noccupied,Noccupied))
    ik=0
    do ikx=1,Nkx
       do iky=1,Nky
          ik=ik+1
          Eigvec = Hk(:,:,ik)
          call eigh(Eigvec,Eigval)
          do iocc=1,Noccupied
             BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
          enddo
       enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
       ikxP = modulo(ikx,Nkx) + 1
       !ikxM = modulo(ikx-2,Nkx) + 1
       do iky= 1, Nky
          ikyP = modulo(iky,Nky) + 1
          !ikyM = modulo(iky-2,Nky) + 1
          !
          if(Noccupied==1)then
             Ulink(1) = dot_product(BlochStates(ikx,iky,1,:)  , BlochStates(ikx,ikyP,1,:))
             Ulink(2) = dot_product(BlochStates(ikx,ikyP,1,:) , BlochStates(ikxP,ikyP,1,:))
             Ulink(3) = dot_product(BlochStates(ikxP,ikyP,1,:), BlochStates(ikxP,iky,1,:))
             Ulink(4) = dot_product(BlochStates(ikxP,iky,1,:) , BlochStates(ikx,iky,1,:))
             !
          else
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j)=dot_product(BlochStates(ikx,iky,i,:)  , BlochStates(ikx,ikyP,j,:))
             Ulink(1) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikx,ikyP,i,:) , BlochStates(ikxP,ikyP,j,:))
             Ulink(2) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikxP,ikyP,i,:), BlochStates(ikxP,iky,j,:))
             Ulink(3) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikxP,iky,i,:) , BlochStates(ikx,iky,j,:))
             Ulink(4) = det(gmat)
             !
          endif
          !
          berry_phase = -dimag(zlog( product(Ulink(:))  ))
          chern = chern + berry_phase
          BerryCurvature(ikx,iky) = berry_phase*one_over_area
          !
       enddo
    enddo
    !
    chern=chern/pi2/3
    !
    open(unit=free_unit(unit),file="Chern_Number.dat",position='append')
    write(unit,*)chern
    write(*,*)chern
    close(unit)
    !
    call splot3d("Berry_Curvature.dat",&
         linspace(0d0,pi2,Nkx,iend=.false.),&
         linspace(0d0,pi2,Nky,iend=.false.),&
         BerryCurvature(:Nkx,:Nky))
    !
  end subroutine Get_Chern_number



end program ed_haldane



