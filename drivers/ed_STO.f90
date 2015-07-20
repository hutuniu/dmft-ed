
program ed_STO
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb
  real(8),allocatable    :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:),Ti3dt2g_Hloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin
  real(8)                :: mh,lambda,wmixing,sumdens
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym
  !convergence function
  complex(8),allocatable :: delta_conv(:,:,:)
  !
  real(8),dimension(2)   :: Eout
#ifdef _MPI
  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call ed_read_input(trim(finput))

  !if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))

  allocate(delta_conv(Nso,Nso,Lmats))

  !Buil the Hamiltonian on a grid or on  path

  call build_hk(trim(hkfile))

  !Setup solver
  Nb=get_bath_size()

  write(*,*) "Nb dri",Nb

  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(bath)
  call set_hloc(my_reshape_A1_to_A2(Ti3dt2g_Hloc))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call ed_get_gloc(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=0)
     call ed_get_weiss(Gmats,Smats,Delta,Hloc=my_reshape_A1_to_A2(Ti3dt2g_Hloc),iprint=0)
     !Fit the new bath, starting from the old bath + the supplied delta 
     Bath_=bath
     if (ed_mode=="normal") then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        call spin_symmetrize_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(delta,bath)
     endif

     !MIXING:
     !if(iloop>1)
     Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     !Bath_=Bath

     do i=1,Lmats
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    delta_conv(io,jo,i)=delta(ispin,jspin,iorb,jorb,i)
                 enddo
              enddo
           enddo
        enddo
     enddo


     converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !if(ED_MPI_ID==0)converged = check_convergence_global(delta_conv(:,:,:),dmft_error,nsuccess,nloop)
     !#ifdef _MPI
     !call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     !#endif

     sumdens=sum(ed_get_dens())
     write(*,*) "sumdens",sumdens,"xmu",xmu,"converged",converged
     if(nread/=0.d0)call search_chemical_potential(xmu,sumdens,converged)
     write(*,*) "sumdens",sumdens,"xmu",xmu,"converged",converged

     if(ED_MPI_ID==0)call end_loop
  enddo
  !Get Kinetic Energy:
  Eout = ed_kinetic_energy(Hk,Wtk,Smats)
  !
#ifdef _MPI
  call MPI_FINALIZE(ED_MPI_ERR)
#endif
contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky,kz    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hk_reshaped
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    if(ED_MPI_ID==0)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**3
    if(ED_MPI_ID==0)write(*,*)"# of k-points     :",Lk
    if(ED_MPI_ID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nso,Nso,Lk));allocate(wtk(Lk));allocate(kxgrid(Nk),kygrid(Nk),kzgrid(Nk))
    kxgrid = kgrid(Nk)
    kygrid = kgrid(Nk)
    kzgrid = kgrid(Nk)

    Hk     = build_hk_model(hk_Ti3dt2g,Nso,kxgrid,kygrid,kzgrid)

    wtk = 1.0d0/Lk
    if(ED_MPI_ID==0.AND.present(file))then
       call write_hk_w90(file,Nso,Nd=Norb,Np=1,Nineq=1,hk=Hk,kxgrid=kxgrid,kygrid=kxgrid,kzgrid=kzgrid)
    endif

    allocate(Ti3dt2g_Hloc(Nso,Nso))
    Ti3dt2g_Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=0d0


    if(ED_MPI_ID==0) then
       call write_Hloc(Ti3dt2g_Hloc)
       write(*,*)
       write(*,*) "Sum over k of H(k) nella versione A1"
       write(*,*) "real"
       do i=1,Nso
          write(*,'(6F10.4)') (real(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(*,*) "complex"
       do i=1,Nso
          write(*,'(6F10.4)') (dimag(Ti3dt2g_Hloc(i,j)),j=1,Nso)
       enddo
       write(*,*)
    endif

    !Build the local GF:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu,Hk(:,:,ik))/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
       enddo
    enddo
    do iorb=1,Nso
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.ed",wm,Gmats(iorb,iorb,:))
       call splot("G0loc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.ed",wr,-dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
    enddo
    !
  end subroutine build_hk



  function hk_Ti3dt2g(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    complex(8),dimension(2,2) :: s_x,s_y,s_z
    complex(8),dimension(2,2) :: t_inter
    real(8)                   :: kx,ky,kz
    real(8)                   :: Eo,t1,t2,t3
    real(8)                   :: soc,ivb
    integer                   :: N,ndx
    complex(8)                :: oneI
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    oneI=cmplx(0.0d0,1.0d0)
    s_x=cmplx(0.0d0,0.0d0);s_y=cmplx(0.0d0,0.0d0);s_z=cmplx(0.0d0,0.0d0)
    s_x(1,2)=cmplx(1.0d0,0.0d0); s_x(2,1)=cmplx(1.0d0,0.0d0)
    s_y(1,2)=cmplx(0.0d0,-1.0d0);s_y(2,1)=cmplx(0.0d0,1.0d0)
    s_z(1,1)=cmplx(1.0d0,0.0d0); s_z(2,2)=cmplx(-1.0d0,0.0d0)
    kx=kvec(1);ky=kvec(2);kz=kvec(3)

    Eo = 3.31
    t1 = 0.277
    t2 = 0.031
    t3 = 0.076
    soc=0.0d0
    ivb=0.0d0

    Hk = zero

    do i=1,Norb
       ndx=2*i-1
       Hk(ndx:ndx+1,ndx:ndx+1) = band_cos_omo(kx,ky,kz,Eo,t1,t2,t3)
    enddo
    Hk(1,2)=0.5d0

    !Diagonale
   ! Hk(1:2,1:2) = band_cos_omo(kx,ky,kz,Eo,t1,t2,t3)
   ! Hk(3:4,3:4) = band_zx(kx,ky,kz,Eo,t1,t2,t3)
   ! Hk(5:6,5:6) = band_xy(kx,ky,kz,Eo,t1,t2,t3)

    !triangolo superiore
   ! Hk(1:2,3:4)=oneI*s_z*soc/2.
   ! Hk(1:2,5:6)=-1.*oneI*s_y*soc/2.+ivb*2*oneI*sin(kx)
   ! Hk(3:4,5:6)=oneI*s_x*soc/2.+ivb*2*oneI*sin(ky)

    !triangolo inferiore
    do i=1,N
       do j=1,N
          Hk(j,i)=conjg(Hk(i,j))
       enddo
    enddo

    Hk = my_reshape_Z_to_A1(Hk)
  
  end function hk_Ti3dt2g

  !2x2 band structure
  function band_cos_omo(kx,ky,kz,Eo,t1,t2,t3) result(hk)
    real(8)                   :: kx,ky,kz
    real(8)                   :: Eo,t1,t2,t3
    complex(8),dimension(2,2) :: hk
    hk = zero
    hk(1,1) = -2.*(cos(kx)+cos(ky)+cos(kz))
    hk(2,2) = hk(1,1)
  end function band_cos_omo
  function band_yz(kx,ky,kz,Eo,t1,t2,t3) result(hk)
    real(8)                   :: kx,ky,kz
    real(8)                   :: Eo,t1,t2,t3
    complex(8),dimension(2,2) :: hk
    hk = zero
    hk(1,1) = Eo-2.*t2*cos(kx)-2.*t1*cos(ky)-2.*t1*cos(kz)-4.*t3*cos(ky)*cos(kz)
    hk(2,2) = hk(1,1)
  end function band_yz
  function band_zx(kx,ky,kz,Eo,t1,t2,t3) result(hk)
    real(8)                   :: kx,ky,kz
    real(8)                   :: Eo,t1,t2,t3
    complex(8),dimension(2,2) :: hk
    hk = zero
    hk(1,1) = Eo-2.*t1*cos(kx)-2.*t2*cos(ky)-2.*t1*cos(kz)-4.*t3*cos(kz)*cos(kx)
    hk(2,2) = hk(1,1)
  end function band_zx
  function band_xy(kx,ky,kz,Eo,t1,t2,t3) result(hk)
    real(8)                   :: kx,ky,kz
    real(8)                   :: Eo,t1,t2,t3
    complex(8),dimension(2,2) :: hk
    hk = zero
    hk(1,1) = Eo-2.*t1*cos(kx)-2.*t1*cos(ky)-2.*t2*cos(kz)-4.*t3*cos(kx)*cos(ky)
    hk(2,2) = hk(1,1)
  end function band_xy




  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                                    :: iw
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: hk
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: g0k
    integer                                       :: i,ndx
    g0k=zero
    do i=1,Norb
      ndx=2*i-1
      g0k(ndx:ndx+1,ndx:ndx+1)=inverse_g0k2x2(iw,hk(ndx:ndx+1,ndx:ndx+1))
    enddo
  end function inverse_g0k
  !
  function inverse_g0k2x2(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k2x2



  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  function my_reshape_Z_to_A1(fg) result(g)
    complex(8),dimension(Nso,Nso)                   :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io1 = iorb + (ispin-1)*Norb
                   jo1 = jorb + (jspin-1)*Norb

                   io2 = ispin + (iorb-1)*Nspin
                   jo2 = jspin + (jorb-1)*Nspin

                   g(io1,jo1)  = fg(io2,jo2)
                enddo
             enddo
          enddo
       enddo
  end function my_reshape_Z_to_A1

  function my_reshape_A1_to_A2(fg) result(g)
    complex(8),dimension(Nso,Nso)                   :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb)     :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin,io,jo
       g = zero
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   g(ispin,jspin,iorb,jorb)  = fg(io,jo)
                enddo
             enddo
          enddo
       enddo
  end function my_reshape_A1_to_A2


end program ed_STO



