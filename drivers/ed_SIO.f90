
program ed_SIO
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb,unit,Lstart
  real(8),allocatable    :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:,:),Sreal(:,:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:,:),Greal(:,:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  complex(8),allocatable :: Ti3dt2g_Hloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:),kzgrid(:)
  !variables for the model:
  integer                :: Nk,Nkpath,i,j,iorb,jorb,io,jo,ispin,jspin,ilat
  real(8)                :: wmixing,dens_per_site,soc
  real(8),allocatable    :: orb_dens(:,:)
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  character(len=32)      :: HlocFILE
  logical                :: spinsym,IOfile
  !convergence function
  complex(8),allocatable :: delta_conv(:,:,:,:),delta_conv_avrg(:)
  !rotation on impHloc
  complex(8),allocatable                            :: impHloc_rot(:,:)
  real(8),allocatable                               :: impHloc_eig(:)

#ifdef _MPI_INEQ
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,     "FINPUT",           default='inputED_SIO.in')
  call parse_input_variable(hkfile,   "HKFILE",finput,    default="hkfile.in")
  call parse_input_variable(HlocFILE, "HlocFILE",finput,  default="impHloc.in")
  call parse_input_variable(nk,       "NK",finput,        default=10)
  call parse_input_variable(wmixing,  "WMIXING",finput,   default=0.5d0)
  call parse_input_variable(soc,      "SOC",finput,       default=1.d0)
  call parse_input_variable(Nlat,     "NLAT",finput,      default=8)
  call ed_read_input(trim(finput))

  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !Allocate convergence funct
  allocate(delta_conv(Nlat,Nso,Nso,Lmats),delta_conv_avrg(Lmats))
  allocate(orb_dens(Nlat,Norb))
  !
  !Read the Hamiltonian
  call read_hk(trim(hkfile))
  !
  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb),Bath_(Nlat,Nb))
  call ed_init_solver_lattice(bath)
  inquire(file="hamiltonian_site0008.restart",exist=IOfile)
  if(.not.IOfile) then
     if(mpiID==0)write(*,*)"Spin alignment"
     call break_symmetry_bath(bath(1,:),0.2d0,+1.d0,.false.)
     call break_symmetry_bath(bath(2,:),0.2d0,-1.d0,.false.)
     call break_symmetry_bath(bath(3,:),0.2d0,-1.d0,.false.)
     call break_symmetry_bath(bath(4,:),0.2d0,+1.d0,.false.)
     call break_symmetry_bath(bath(5,:),0.2d0,+1.d0,.false.)
     call break_symmetry_bath(bath(6,:),0.2d0,-1.d0,.false.)
     call break_symmetry_bath(bath(7,:),0.2d0,-1.d0,.false.)
     call break_symmetry_bath(bath(8,:),0.2d0,+1.d0,.false.)
  endif
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop")
     !
     call ed_solve_lattice(bath,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3)
     !
     if(mpiID==0)call rotate_Gimp()
     if(mpiID==0)call orbital_spin_mixture()
     !
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)
     call ed_get_sigma_real_lattice(Sreal,Nlat)
     !
     !if (iloop==1) then
        !if(mpiID==0) call spin_symmetrize_lattice(Smats,5)
        !if(mpiID==0) call spin_symmetrize_lattice(Sreal,5)
     !endif      
     !
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=3)
     call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3)
     Bath_=bath
     call ed_chi2_fitgf_lattice(bath,delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc))
     if(iloop>1) Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath
     delta_conv=zero
     delta_conv_avrg=zero
     do i=1,Lmats
        do ilat=1,Nlat
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       if((ispin.eq.jspin).and.(iorb.eq.jorb)) then
                          io = iorb + (ispin-1)*Norb
                          jo = jorb + (jspin-1)*Norb
                          delta_conv(ilat,io,jo,i)=delta(ilat,ispin,jspin,iorb,jorb,i)
                          delta_conv_avrg(i)=delta_conv_avrg(i)+delta_conv(ilat,io,jo,i)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     delta_conv_avrg=delta_conv_avrg/(Nso*Nlat)
     if(mpiID==0) converged = check_convergence(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(mpiID==0) converged = check_convergence_global(delta_conv_avrg,dmft_error,nsuccess,nloop)
     !if(mpiID==0) converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !if(mpiID==0)converged = check_convergence_global(delta_conv(:,:,:),dmft_error,nsuccess,nloop)
     !
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     !
     orb_dens=ed_get_dens_lattice(Nlat)
     dens_per_site=sum(orb_dens)/Nlat
     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged
     if(nread/=0.d0)call search_chemical_potential(xmu,dens_per_site,converged)
     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged
     !
     if(mpiID==0)call end_loop
  enddo
#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif
contains



  !_______________________________________________________________________
  !                            HAMILTONIAN
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: H(k) file for main program and write G0_loc
  !---------------------------------------------------------------------
  subroutine read_hk(file)
    character(len=*),optional                         :: file
    integer                                           :: i,j,ik=0
    integer                                           :: ix,iy
    real(8)                                           :: kx,ky,kz    
    integer                                           :: io,jo
    integer                                           :: iorb,jorb,ispin,jspin,ilat,jlat
    integer                                           :: unit
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lmats)     :: Gmats
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lreal)     :: Greal
    real(8)                                           :: wm(Lmats),wr(Lreal),dw,xmu_0
    real(8)                                           :: dumR(Nlat*Nso,Nlat*Nso),dumI(Nlat*Nso,Nlat*Nso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)  :: Hloc_dum
    logical                                           :: intersite
    complex(8),allocatable                            :: site_matrix_in(:,:),site_matrix_out(:,:)
    complex(8),allocatable                            :: impHloc_app(:,:)

    if(mpiID==0)write(LOGfile,*)"Read H(k) for SIO:"
    Lk=Nk
    if(mpiID==0)write(*,*)"# of k-points     :",Lk
    if(mpiID==0)write(*,*)"# of sites        :",Nlat
    if(mpiID==0)write(*,*)"# of SO-bands     :",Nso
    if(mpiID==0)write(*,*)"# of SOC factor   :",soc
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nlat*Nso,Nlat*Nso,Lk));allocate(wtk(Lk))
    wtk = 1.0d0/Lk
    !
    open(unit=123,file='hk_2_ED.dat',status='old',action='read')
    do ik=1,Lk
       do io=1,Nlat*Nspin*Norb
          read(123,'(50F10.5)') (dumR(io,jo),jo=1,Nlat*Nspin*Norb)
       enddo
       do io=1,Nlat*Nspin*Norb
          read(123,'(50F10.5)') (dumI(io,jo),jo=1,Nlat*Nspin*Norb)
       enddo
       Hk(:,:,ik)=cmplx(dumR,dumI)
    enddo
    if (soc .ne. 1.0d0) then
       write(*,*)"rescaling SOC"
       do io=1,Nlat*Nspin*Norb
          do jo=1,Nlat*Nspin*Norb
             if(io.ne.jo) Hk(io,jo,:)=Hk(io,jo,:)/soc
          enddo
       enddo
    endif
    close(123)
    !
    allocate(Ti3dt2g_Hloc(Nlat*Nso,Nlat*Nso))
    Ti3dt2g_Hloc = sum(Hk,dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=zero
    if(mpiID==0) then
       call write_Hloc(Ti3dt2g_Hloc,HlocFILE)
    endif
    !
    Hloc_dum=reshape_A1_to_A2_L(Ti3dt2g_Hloc)
    open(unit=100,file='impEDloc_R.dat',status='unknown',action='write',position='rewind')
    open(unit=101,file='impEDloc_I.dat',status='unknown',action='write',position='rewind')
    if(mpiID==0) then
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      write(100,'(5I3,10F15.10)')ilat,ispin,jspin,iorb,jorb,real(Hloc_dum(ilat,ispin,jspin,iorb,jorb))
                   enddo
                enddo
             enddo
          enddo
       enddo
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      write(101,'(5I3,10F15.10)')ilat,ispin,jspin,iorb,jorb,aimag(Hloc_dum(ilat,ispin,jspin,iorb,jorb))
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    close(100);close(101)
    !
    !rotation on impHloc
    if(allocated(impHloc_rot)) deallocate(impHloc_rot)
    allocate(impHloc_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb));impHloc_rot=zero
    if(allocated(impHloc_eig)) deallocate(impHloc_eig)
    allocate(impHloc_eig(Nlat*Nspin*Norb));impHloc_eig=0.d0
    !
    intersite=.false.
    ! con/senza termini inter-sito
    if (intersite) then
       !
       impHloc_rot=zero
       impHloc_rot=Ti3dt2g_Hloc
       call matrix_diagonalize(impHloc_rot,impHloc_eig,'V','U')
       !
       if(allocated(site_matrix_in)) deallocate(site_matrix_in)
       allocate(site_matrix_in(Nlat*Nspin*Norb,Nlat*Nspin*Norb));site_matrix_in=zero
       if(allocated(site_matrix_out)) deallocate(site_matrix_out)
       allocate(site_matrix_out(Nlat*Nspin*Norb,Nlat*Nspin*Norb));site_matrix_out=zero
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                site_matrix_in(io,io)=cmplx(float(ilat),0.0d0)
             enddo
          enddo
       enddo
       !
       site_matrix_out=matmul(transpose(conjg(impHloc_rot)),matmul(site_matrix_in,impHloc_rot))
       !
       open(unit=104,file='site_mixing.dat',status='unknown',action='write',position='rewind')
       write(104,*)"Sites before rotation"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100I12)')(nint(real(site_matrix_in(io,jo))),jo=1,Nlat*Nspin*Norb)
       enddo
       write(104,*)"Sites after rotation"
       write(104,*)"R:"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100F12.4)')(real(site_matrix_out(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(104,*)"I:"
       do io=1,Nlat*Nspin*Norb
          write(104,'(100F12.4)')(aimag(site_matrix_out(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       close(104)
       !
    else
       !
       if(allocated(impHloc_app)) deallocate(impHloc_app)
       allocate(impHloc_app(Nspin*Norb,Nspin*Norb));impHloc_app=zero
       !
       impHloc_rot=zero
       do ilat=1,Nlat
          io=1+(ilat-1)*Nso
          jo=Nso+(ilat-1)*Nso
          impHloc_app=Ti3dt2g_Hloc(io:jo,io:jo)
          call matrix_diagonalize(impHloc_app,impHloc_eig(io:jo),'V','U')
          impHloc_rot(io:jo,io:jo)=impHloc_app
       enddo
       !
    endif
    !
    open(unit=102,file='impHloc_eig.dat',status='unknown',action='write',position='rewind')
    do ilat=1,Nlat
       write(102,'(1I3,20F25.20)')ilat,(impHloc_eig(io),io=1+(ilat-1),Nspin*Norb+(ilat-1))
    enddo
    close(102)
    !
    open(unit=103,file='impHloc_rot.dat',status='unknown',action='write',position='rewind')
    write(103,*)"impHloc rotation, Real Part"
    do io=1,Nlat*Nspin*Norb
       write(103,'(100F12.4)')(real(impHloc_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(103,*)"impHloc rotation, Iaginary Part"
    do io=1,Nlat*Nspin*Norb
       write(103,'(100F12.4)')(aimag(impHloc_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    close(103)
    !
    go to 223
    !
    !Build the local GF in the spin-orbital Basis:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    xmu_0=Ti3dt2g_Hloc(1,1)
    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu_0,Hk(:,:,ik) )/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu_0,Hk(:,:,ik))/Lk
       enddo
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_iw.ed",wm,dimag(Gmats(io,jo,:)),real(Gmats(io,jo,:)) )
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,real(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    do i=1,Lmats
       Gmats(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(Gmats(:,:,i),impHloc_rot))
    enddo
    do i=1,Lreal
       Greal(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(Greal(:,:,i),impHloc_rot))
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_iw.ed",wm,dimag(Gmats(io,jo,:)),real(Gmats(io,jo,:)) )
                   call splot("G0loc_rot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,real(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    223 continue
    !
  end subroutine read_hk


  !_______________________________________________________________________
  !                                    Gfs
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: G0_loc functions DA RIFARE ATTENZIONE CHE H(k) è nella forma A1
  !---------------------------------------------------------------------
  function inverse_g0k(iw,hk) result(g0k)
    implicit none
    complex(8)                                              :: iw
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)   :: hk
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)   :: g0k
    !
    g0k=zero
    g0k=iw*eye(Nlat*Nspin*Norb)-hk
    !
    call inv(g0k)
    !
  end function inverse_g0k


  !---------------------------------------------------------------------
  !PURPOSE: rotaizione delle Gimp per portarle in una base simile a J,jz
  !---------------------------------------------------------------------
  subroutine rotate_Gimp()
    implicit none
    complex(8),allocatable             :: G_in(:,:,:),G_out(:,:,:),Gso(:,:,:,:,:,:)
    integer                            :: ilat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    !
    write(*,*) "A(w) rotation"
    !
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate(  Gso(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gso=zero
    allocate( G_in(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_in=zero
    allocate(G_out(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_out=zero
    !
    call ed_get_gimp_real_lattice(Gso,Nlat)
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   G_in(io,jo,:)=Gso(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(G_in(:,:,i),impHloc_rot))
    enddo
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("impGrot_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(G_out(io,jo,:))/pi,dreal(G_out(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine rotate_Gimp


  !---------------------------------------------------------------------
  !PURPOSE: 
  !---------------------------------------------------------------------
  subroutine spin_symmetrize_lattice(Self,lat)
    complex(8),allocatable,intent(inout)   :: Self(:,:,:,:,:,:)
    integer                                :: lat
    complex(8),allocatable                 :: Self_aux(:,:,:,:,:,:)

    if(mpiID==0)write(*,*)"Symmetrizing Sigma"
    allocate(Self_aux(Nlat,Nspin,Nspin,Norb,Norb,size(Self,dim=6)));Self_aux=zero
    Self_aux=Self
    do ilat=lat,Nlat
       if(mpiID==0)write(*,*)"site",ilat
       do iorb=1,Norb
          do jorb=1,Norb
             Self(ilat,1,1,iorb,jorb,:)=Self_aux(ilat,2,2,iorb,jorb,:)
             Self(ilat,1,2,iorb,jorb,:)=Self_aux(ilat,2,1,iorb,jorb,:)
             Self(ilat,2,1,iorb,jorb,:)=Self_aux(ilat,1,2,iorb,jorb,:)
             Self(ilat,2,2,iorb,jorb,:)=Self_aux(ilat,1,1,iorb,jorb,:)
          enddo
       enddo
    enddo
    deallocate(Self_aux)
  end subroutine spin_symmetrize_lattice


  !---------------------------------------------------------------------
  !PURPOSE: 
  !---------------------------------------------------------------------
  subroutine orbital_spin_mixture()
    implicit none
    complex(8),allocatable             :: Gso(:,:,:,:,:,:),Stot(:,:,:,:)
    integer                            :: ilat,io,jo
    integer                            :: ispin,jspin
    integer                            :: iorb,jorb
    real(8)                            :: wm(Lmats),wr(Lreal),dw
    real(8)                            :: site_mag(Nlat,Norb)
    !
    write(*,*) "Computing oprbital spin mixture per site"
    write(*,*) "Lmats used:",Lmats
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate( Gso(Nlat,Nspin,Nspin,Norb,Norb,Lmats)); Gso=zero
    allocate(Stot(Nlat,3,Norb,Norb));Stot=zero
    !
    call ed_get_gimp_matsubara_lattice(Gso,Nlat)
    !
    do ilat=1,Nlat
       do iorb=1,Norb
          do jorb=1,Norb
             !Sx
             Stot(ilat,1,iorb,jorb)=sum(    (Gso(ilat,1,2,iorb,jorb,:)+Gso(ilat,2,1,iorb,jorb,:) ))/beta
           !  Stot(ilat,1,iorb,jorb)=   (Gso(ilat,1,2,iorb,jorb,1)+Gso(ilat,2,1,iorb,jorb,1))/3.
             !Sy
             Stot(ilat,2,iorb,jorb)=sum( xi*(Gso(ilat,2,1,iorb,jorb,:)-Gso(ilat,1,2,iorb,jorb,:) ))/beta
           !  Stot(ilat,2,iorb,jorb)=xi*(Gso(ilat,2,1,iorb,jorb,1)-Gso(ilat,1,2,iorb,jorb,1))/3.
             !Sz
             Stot(ilat,3,iorb,jorb)=sum(    (Gso(ilat,1,1,iorb,jorb,:)-Gso(ilat,2,2,iorb,jorb,:) ))/beta
           !  Stot(ilat,3,iorb,jorb)=   (Gso(ilat,1,1,iorb,jorb,1)-Gso(ilat,2,2,iorb,jorb,1))/3.
          enddo
       enddo
    enddo
    !
    site_mag=0.d0
    site_mag=ed_get_mag_lattice(Nlat)
    !
    open(unit=105,file='Spin_mixture.dat',status='unknown',position='rewind',action='write',form='formatted')
    write(105,'(a100)') "#diagonal site, diagonal orbital"
    do ilat=1,Nlat
      ! write(105,'(2a8,30a20)') "#site","orbital","Re{Sx}","Re{Sy}","Re{Sz}","Im{Sx}","Im{Sy}","Im{Sz}","mag"
      ! do iorb=1,Norb
      !    write(105,'(2I8,30F20.12)') ilat, iorb, real(Stot(ilat,1,iorb,iorb)), real(Stot(ilat,2,iorb,iorb)), real(Stot(ilat,3,iorb,iorb)) &
      !                                          ,aimag(Stot(ilat,1,iorb,iorb)),aimag(Stot(ilat,2,iorb,iorb)),aimag(Stot(ilat,3,iorb,iorb)) &
      !                                          ,site_mag(ilat,iorb)/2.
      ! enddo
       write(105,'(I8,30F20.12)') ilat,  real(Stot(ilat,1,1,1)), real(Stot(ilat,1,2,2)), real(Stot(ilat,1,3,3)) &
                                      ,  real(Stot(ilat,2,1,1)), real(Stot(ilat,2,2,2)), real(Stot(ilat,2,3,3)) &
                                      ,  real(Stot(ilat,3,1,1)), real(Stot(ilat,3,2,2)), real(Stot(ilat,3,3,3)) &
                                      , aimag(Stot(ilat,1,1,1)),aimag(Stot(ilat,1,2,2)),aimag(Stot(ilat,1,3,3)) &       
                                      , aimag(Stot(ilat,2,1,1)),aimag(Stot(ilat,2,2,2)),aimag(Stot(ilat,2,3,3)) &
                                      , aimag(Stot(ilat,3,1,1)),aimag(Stot(ilat,3,2,2)),aimag(Stot(ilat,3,3,3)) &
                                      , site_mag(ilat,1)/2.,site_mag(ilat,2)/2.,site_mag(ilat,3)/2.
    enddo
    write(105,*)
    write(105,*)
    write(105,'(a100)') "#diagonal site, inter-orbital"
    do ilat=1,Nlat
       write(105,'(a8,I3)') "#site:",ilat
       write(105,'(30a20)') "#Sx(orb_1)","Sx(orb_2)","Sx(orb_3)","Sy(orb_1)","Sy(orb_2)","Sy(orb_3)","Sz(orb_1)","Sz(orb_2)","Sz(orb_3)"
       do iorb=1,Norb
          write(105,'(30F20.12)') (real(Stot(ilat,1,iorb,jorb)),jorb=1,Norb) &
                                 ,(real(Stot(ilat,2,iorb,jorb)),jorb=1,Norb) &
                                 ,(real(Stot(ilat,3,iorb,jorb)),jorb=1,Norb)
       enddo
       write(105,*)
       do iorb=1,Norb
          write(105,'(30F20.12)') (aimag(Stot(ilat,1,iorb,jorb)),jorb=1,Norb) &
                                 ,(aimag(Stot(ilat,2,iorb,jorb)),jorb=1,Norb) &
                                 ,(aimag(Stot(ilat,3,iorb,jorb)),jorb=1,Norb)
       enddo
    enddo
    close(105)
    deallocate(Gso,Stot)
    !
  end subroutine orbital_spin_mixture


  !_______________________________________________________________________
  !                            reshape functions
  !_______________________________________________________________________
  !---------------------------------------------------------------------
  !PURPOSE: reshape functions
  !  Z  = [Nspin,Nspin]*Norb
  !  A1 = [Norb*Norb]*Nspin
  !  A2 = [Nspin,Nspin,Norb,Norb]
  !---------------------------------------------------------------------
  function reshape_Z_to_A1(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = iorb + (ispin-1)*Norb
                jo1 = jorb + (jspin-1)*Norb
                !I-index
                io2 = ispin + (iorb-1)*Nspin
                jo2 = jspin + (jorb-1)*Nspin
                !switch
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function reshape_Z_to_A1

  function reshape_A1_to_Z(fg) result(g)
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: fg
    complex(8),dimension((Nspin*Norb),(Nspin*Norb)) :: g
    integer                                         :: i,j,iorb,jorb,ispin,jspin
    integer                                         :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = ispin + (iorb-1)*Nspin
                jo1 = jspin + (jorb-1)*Nspin
                !I-index
                io2 = iorb + (ispin-1)*Norb
                jo2 = jorb + (jspin-1)*Norb
                !switchHloc
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function reshape_A1_to_Z

  function reshape_A1_to_A2(fg) result(g)
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
  end function reshape_A1_to_A2

  function reshape_A1_to_A2_L(fg) result(g)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)          :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: g
    integer                                          :: i,j,iorb,jorb,ispin,jspin,io,jo,ilat
    g = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   g(ilat,ispin,jspin,iorb,jorb)  = fg(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function reshape_A1_to_A2_L

  !---------------------------------------------------------------------
  !PURPOSE: Inversion test
  !---------------------------------------------------------------------
  subroutine inversion_test(A,B,tol)
    implicit none
    complex (kind=8), intent(in)   ::   A(Nspin*Norb,Nspin*Norb)
    complex (kind=8), intent(in)   ::   B(Nspin*Norb,Nspin*Norb)
    real    (kind=4), intent(in)   ::   tol
    integer (kind=2)               ::   dime

    if (size(A).ne.size(B)) then
       write(*,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    if (abs(float(dime)-real(sum(matmul(A,B)))).gt.tol) write(*,'(A30)') "inversion test fail"
  end subroutine inversion_test


end program ed_SIO

