!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE SCIFOR
  implicit none
  integer                :: iloop,Lk
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  integer                :: Nporb,Ntorb
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=32)      :: hkfile
  integer                :: ntype
  real(8)                :: nobj
  real(8)                :: ts(3)
  character(len=16)      :: finput

  !parse additional variables && read input file && read Hk
  call parse_cmd_variable(finput,'FINPUT',default='inputED.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default='hkfile.in')
  call parse_input_variable(ntype,"NTYPE",finput,default=0)
  call parse_input_variable(Lk,"Lk",finput,default=1000)
  call parse_input_variable(Nporb,"NPORB",finput,default=0)
  call parse_input_variable(ts,"TS",finput,default=[0.5d0,0.5d0,0.5d0])
  !
  call ed_read_input(trim(finput))
  !

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))


  !setup solver
  !read Hloc
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)


  call read_hk(trim(hkfile))
  if(ntype>Ntorb)stop "Ntype > Ntorb"

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


contains


  subroutine read_hk(file)
    character(len=*) :: file
    integer          :: i,j,ik,iorb,jorb,Norb_d,Norb_p
    real(8)          :: de,e,kx,ky,kz,foo
    logical          :: ioexist
    
    inquire(file=trim(file),exist=ioexist)
    if(.not.ioexist)then
       Nporb=0
       Ntorb=Norb+Nporb
       allocate(Hk(Ntorb,Ntorb,Lk))
       allocate(dos_wt(Lk))
       de=2.d0/dble(Lk)
       do ik=1,Lk
          e = -1.d0 + dble(ik-1)*de
          Hk(:,:,ik) = Hloc(1,1,:,:) !set local part
          do iorb=1,Norb
             Hk(iorb,iorb,ik) = Hk(iorb,iorb,ik) - 2.d0*ts(iorb)*e
          enddo
          dos_wt(ik)=dens_bethe(e,1.d0)*de
       enddo
       !
    else
       !
       open(50,file=file,status='old')
       read(50,*)Lk,Norb_d,Nporb,foo,foo
       if(Norb_d/=Norb)stop "Can not read Hk.file: check Norb_d != Norb"
       Ntorb=Norb+Nporb
       allocate(Hk(Ntorb,Ntorb,Lk))
       allocate(dos_wt(Lk))
       do ik=1,Lk
          read(50,"(3(F10.7,1x))")kx,ky,kz
          do iorb=1,Ntorb
             read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Ntorb)
          enddo
       enddo
       close(50)
       dos_wt=1.d0/dble(Lk)
       Hloc(1,1,:Norb,:Norb) = sum(Hk(:Norb,:Norb,:),dim=3)/dble(Lk)
       !
    endif
    print*,"Hloc:"
    call print_Hloc(Hloc)
  end subroutine read_hk


  !+----------------------------------------+


  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8)                              :: iw
    complex(8),dimension(Ntorb,Ntorb)       :: fg,gdelta
    complex(8),allocatable,dimension(:,:,:) :: Gloc,Sigma,Zeta
    real(8)                                 :: wm(Lmats),wr(Lreal),dens(Ntorb)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta=zero
    !
    allocate(Gloc(Ntorb,Ntorb,Lmats),Sigma(Ntorb,Ntorb,Lmats),Zeta(Ntorb,Ntorb,Lmats))
    Sigma = zero
    Gloc  = zero
    Zeta  = zero
    forall(i=1:LMats,iorb=1:Norb,jorb=1:Norb)Sigma(iorb,jorb,i)=impSmats(1,1,iorb,jorb,i)
    Zeta  = -Sigma
    forall(i=1:Lmats,iorb=1:Ntorb)Zeta(iorb,iorb,i)=xi*wm(i)+xmu+zeta(iorb,iorb,i)
    do i=1,Lmats
       !
       fg=zero
       do ik=1,Lk
          fg = fg + invert_gk(zeta(:,:,i),Hk(:,:,ik),Ntorb)*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
       !
       if(Ntorb>1)then
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
       delta(:,:,i)=gdelta(1:Norb,1:Norb)
    enddo
    !print
    do iorb=1,Ntorb
       do jorb=1,Ntorb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
          call splot("Delta_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
       enddo
    enddo
    dens(1:Norb)=ed_dens(1:Norb)
    if(Nporb>0)then
       do iorb=1,Nporb
          dens(Norb+iorb) =  get_density_fromFFT(Gloc(Norb+iorb,Norb+iorb,:),beta)
          write(*,"(A,F25.18)")"np      =",dens(Norb+iorb)
       enddo
    endif
    write(*,"(A,F25.18)")"n_total =",sum(dens)
    open(100,file="nd_np_ntot_all.ed",position='append')
    write(100,"(100F25.18)")(dens(iorb),iorb=1,Ntorb),sum(dens)
    close(100)
    open(100,file="nd_np_ntot_last.ed")
    write(100,"(100F25.18)")(dens(iorb),iorb=1,Ntorb),sum(dens)
    close(100)
    deallocate(Gloc,Sigma,Zeta)


    allocate(Gloc(Ntorb,Ntorb,Lreal),Sigma(Ntorb,Ntorb,Lreal),Zeta(Ntorb,Ntorb,Lreal))
    Sigma = zero ; Gloc  = zero ; Zeta  = zero
    forall(i=1:Lreal,iorb=1:Norb,jorb=1:Norb)Sigma(iorb,jorb,i)=impSreal(1,1,iorb,jorb,i)
    Zeta  = -Sigma
    forall(i=1:Lreal,iorb=1:Ntorb)Zeta(iorb,iorb,i)=dcmplx(wr(i),eps)+xmu+zeta(iorb,iorb,i)
    do i=1,Lreal
       fg=zero
       do ik=1,Lk
          fg = fg + invert_gk(zeta(:,:,i),Hk(:,:,ik),Ntorb)*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
    enddo
    do iorb=1,Ntorb
       do jorb=1,Ntorb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi,dreal(gloc(iorb,jorb,:)))
       enddo
    enddo
    deallocate(Gloc,Sigma,Zeta)


    if(Ntype==0)then
       nobj = sum(dens)
    else
       nobj = dens(Ntype)
    endif

  end subroutine get_delta


  function invert_gk(zeta,Hk,Ndim) result(invg)
    integer                         :: Ndim
    complex(8),dimension(Ndim,Ndim) :: zeta,Hk
    complex(8),dimension(Ndim,Ndim) :: invg
    invg=zeta-Hk
    if(Ndim==1)then
       invg(1,1)=one/invg(1,1)
    else
       call matrix_inverse(invg)
    endif
  end function invert_gk

  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT


end program lancED



