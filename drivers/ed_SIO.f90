
program ed_STO
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
  integer                :: Nb,unit
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
  real(8)                :: wmixing,dens_per_site
  real(8),allocatable    :: orb_dens(:,:)
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  character(len=32)      :: HlocFILE
  logical                :: spinsym
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
  call parse_input_variable(nkpath,   "NKPATH",finput,    default=500)
  call parse_input_variable(wmixing,  "WMIXING",finput,   default=0.5d0)
  call parse_input_variable(spinsym,  "SPINSYM",finput,   default=.true.)
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


  !Read the Hamiltonian on a grid
  call read_hk(trim(hkfile))

  !stop

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nlat,Nb),Bath_(Nlat,Nb))

  call ed_init_solver_lattice(bath)                !ok
  !call set_hloc(reshape_A1_to_A2_L(Ti3dt2g_Hloc))  !ok questo lo deve produrre read_hk

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath) |CONTROLLA INPUT
     call ed_solve_lattice(bath,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3)                  !ok
     if(mpiID==0)call rotate_Gimp()
     call ed_get_sigma_matsubara_lattice(Smats,Nlat)                                             !ok
     call ed_get_sigma_real_lattice(Sreal,Nlat)                                                  !ok
     call ed_get_gloc_lattice(Hk,Wtk,Gmats,Greal,Smats,Sreal,iprint=3)                           !ok
     call ed_get_weiss_lattice(Gmats,Smats,Delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc),iprint=3) !ok
     !Fit the new bath, starting from the old bath + the supplied delta 
     Bath_=bath
     call ed_chi2_fitgf_lattice(bath,delta,Hloc=reshape_A1_to_A2_L(Ti3dt2g_Hloc))
     !mixing:
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
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif

     orb_dens=ed_get_dens_lattice(Nlat)
     dens_per_site=sum(orb_dens)/Nlat

     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged
     if(nread/=0.d0)call search_chemical_potential(xmu,dens_per_site,converged)
     if (mpiID==0) write(*,*) "dens_per_site",dens_per_site,"xmu",xmu,"converged",converged

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
    real(8)                                           :: wm(Lmats),wr(Lreal),dw
    real(8)                                           :: dumR(Nlat*Nso,Nlat*Nso),dumI(Nlat*Nso,Nlat*Nso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)  :: Hloc_dum

    if(mpiID==0)write(LOGfile,*)"Read H(k) for STO:"
    Lk=Nk
    if(mpiID==0)write(*,*)"# of k-points     :",Lk
    if(mpiID==0)write(*,*)"# of sites        :",Nlat
    if(mpiID==0)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    allocate(Hk(Nlat*Nso,Nlat*Nso,Lk));allocate(wtk(Lk))
    wtk = 1.0d0/Lk

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
    close(123)

    allocate(Ti3dt2g_Hloc(Nlat*Nso,Nlat*Nso))
    Ti3dt2g_Hloc = sum(Hk,dim=3)/Lk
    where(abs((Ti3dt2g_Hloc))<1.d-9)Ti3dt2g_Hloc=zero
    if(mpiID==0) then
       call write_Hloc(Ti3dt2g_Hloc,HlocFILE)
    endif

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

    !Build the local GF in the spin-orbital Basis:
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)

    do ik=1,Lk
       do i=1,Lmats
          Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k( xi*wm(i)+xmu , Hk(:,:,ik) )/Lk
       enddo
       do i=1,Lreal
          Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))/Lk
       enddo
    enddo
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_iw.ed",wm,Gmats(io,jo,:))
                   call splot("G0loc_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_lat"//reg(txtfy(ilat))//"_realw.ed",wr,-dimag(Greal(io,jo,:))/pi,dreal(Greal(io,jo,:)))
                enddo
             enddo
          enddo
       enddo
    enddo

    !rotation on impHloc
    if(allocated(impHloc_rot)) deallocate(impHloc_rot)
    allocate(impHloc_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb));impHloc_rot=zero
    if(allocated(impHloc_eig)) deallocate(impHloc_eig)
    allocate(impHloc_eig(Nlat*Nspin*Norb));impHloc_eig=0.d0
    impHloc_rot=Ti3dt2g_Hloc

    !1) CON TERMINI INTER-SITO
    !call matrix_diagonalize(impHloc_rot,impHloc_eig,'V','U')

    !2) SENZA TERMINI INTER-SITO
    impHloc_rot=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   impHloc_rot(io,jo)=Ti3dt2g_Hloc(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    call matrix_diagonalize(impHloc_rot,impHloc_eig,'V','U')

    open(unit=102,file='impHloc_eig.dat',status='unknown',action='write',position='rewind')
    open(unit=103,file='impHloc_rot.dat',status='unknown',action='write',position='rewind')

    do ilat=1,Nlat
       write(102,'(1I3,20F25.20)')ilat,(impHloc_eig(io),io=1+(ilat-1),Nspin*Norb+(ilat-1))
       do io=1+(ilat-1),Nso+(ilat-1)
          write(103,'(1I3,20F25.20)')ilat,(impHloc_rot(io,jo),jo=1+(ilat-1),Nspin*Norb+(ilat-1))
       enddo
    enddo
    close(102);close(103)

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
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)   :: g0k,g0k_tmp
    integer                                                 :: i,ndx

    g0k=zero
    g0k_tmp=zero

    g0k=iw*eye(Nlat*Nspin*Norb)-hk
    g0k_tmp=g0k

    call inv(g0k)
    !call inversion_test(g0k,g0k_tmp,1.e-6)
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

    write(*,*) "A(w) rotation"

    wr = linspace(wini,wfin,Lreal,mesh=dw)
    allocate(  Gso(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gso=zero
    allocate( G_in(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_in=zero
    allocate(G_out(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));G_out=zero

    call ed_get_gimp_real_lattice(Gso,Nlat)

    !test
    write(*,*) sum(Gso(:,1,2,:,:,:))
    write(*,*) sum(Gso(:,2,1,:,:,:))
    write(*,*) sum(Gso(:,1,:,1,1,:))
    write(*,*) sum(Gso(:,1,2,:,:,:))

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

    do i=1,Lreal
       G_out(:,:,i)=matmul(transpose(conjg(impHloc_rot)),matmul(G_in(:,:,i),impHloc_rot))
    enddo

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

  end subroutine rotate_Gimp



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


end program ed_STO

