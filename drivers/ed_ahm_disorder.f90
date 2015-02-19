!########################################################
!PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using  DMFT-ED
! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. 
!########################################################
program ed_ahm_disorder
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  complex(8),allocatable,dimension(:,:,:)  :: Smats,Sreal !self_energies
  complex(8),allocatable,dimension(:,:,:)  :: Gmats,Greal !local green's functions
  complex(8),allocatable,dimension(:,:,:)  :: Delta,errDelta
  real(8),allocatable,dimension(:)         :: erandom
  real(8),allocatable,dimension(:,:,:)     :: bath,bath_old,errBath
  logical                                  :: converged
  real(8)                                  :: wmixing,Wdis
  integer                                  :: i,is,iloop,nrandom,idum
  integer                                  :: Nb(2)
  real(8),dimension(:),allocatable         :: wm,wr
  real(8),dimension(:,:),allocatable       :: nii,dii,pii,eii ![Nlat][Norb]/[4]
  complex(8),dimension(:,:,:),allocatable  :: Hk

  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)


  ! READ INPUT FILES !
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(Wdis,"WDIS","inputRDMFT.in",default=0.d0)
  call parse_input_variable(idum,"IDUM","inputRDMFT.in",default=1234567)
  call ed_read_input("inputRDMFT.in")

  ! SET THE COMPRESSION THRESHOLD TO 1Mb (1024Kb)
  call set_store_size(1024)

  ! GET THE ACTUAL NUMBER OF LATTICE SITES !
  Nlat = Nside**2


  ! TO BE UPDATE WITH THE DMFT_TIGHT_BINDING ROUTINES !
  ! GET THE TIGHT BINDING HAMILTONIAN FOR THE SQUARE LATTICE !
  call get_lattice_hamiltonian(Nside,Nside)


  ! ALLOCATE MATSUBARA AND REAL FREQUENCIES !
  allocate(wm(Lmats),wr(Lreal))
  wini  =wini-Wdis
  wfin  =wfin+Wdis
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)


  ! RANDOM ENERGIES ! 
  allocate(erandom(Nlat))
  call random_seed(size=nrandom)  
  call random_seed(put=(/(idum,i=1,nrandom)/))
  call random_number(erandom)
  erandom=(2.d0*erandom-1.d0)*Wdis/2.d0
  !erandom(Nlat) = -sum(erandom(1:Nlat-1))
  call store_data("erandomVSisite.ed",erandom)
  call splot("erandom.ed",erandom,(/(1,i=1,size(erandom))/))


  ! ALLOCATE A BATH FOR EACH IMPURITY !
  Nb=get_bath_size()
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_old(Nlat,Nb(1),Nb(2)))

  ! ALLOCATE ALL THE REQUIRED VARIABLES !
  allocate(nii(Nlat,Norb))
  allocate(dii(Nlat,Norb))
  allocate(pii(Nlat,Norb))
  allocate(eii(Nlat,4))
  allocate(Smats(2,Nlat,Lmats))
  allocate(Sreal(2,Nlat,Lreal))
  allocate(Gmats(2,Nlat,Lmats))
  allocate(Greal(2,Nlat,Lreal))
  allocate(Delta(2,Nlat,Lmats))
  !
  allocate(Hk(Nlat,Nlat,1))
  Hk(:,:,1) = one*H0

  !+- initialize baths -+!
  call  ed_init_solver_lattice(bath)
  bath_old=bath

  !+- DMFT LOOP -+!
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop",unit=LOGfile)

     bath_old = bath
     call ed_solve_lattice(bath,Eloc=erandom)
     nii = ed_get_dens_lattice(Nlat)
     dii = ed_get_docc_lattice(Nlat)
     pii = ed_get_phisc_lattice(Nlat)
     eii = ed_get_eimp_lattice(Nlat)
     call ed_get_sigma_matsubara_lattice(Smats(1,:,:),Nlat)
     call ed_get_self_matsubara_lattice(Smats(2,:,:),Nlat)
     call ed_get_sigma_real_lattice(Sreal(1,:,:),Nlat)
     call ed_get_self_real_lattice(Sreal(2,:,:),Nlat)

     call ed_get_gloc_lattice(Hk,[1d0],Gmats=Gmats,Greal=Greal,Smats=Smats,Sreal=Sreal,iprint=0,Eloc=erandom)
     call ed_get_weiss_lattice(Nlat,Gmats,Smats,Delta,Eloc=erandom)
     call ed_chi2_fitgf_lattice(bath,Delta,Eloc=erandom)
     bath = wmixing*bath + (1.d0-wmixing)*bath_old
     if(rdmft_phsym)then
        do i=1,Nlat
           call ph_symmetrize_bath(bath(i,:,:))
        enddo
     endif
     if(mpiID==0)converged = check_convergence_local(bath(:,1,:),dmft_error,Nsuccess,nloop,index=2,total=3,id=0,file="BATHerror.err",reset=.false.)
     if(mpiID==0)converged = check_convergence(Delta(1,:,:),dmft_error,Nsuccess,nloop,index=3,total=3,id=0,file="DELTAerror.err",reset=.false.)
     if(mpiID==0)converged = check_convergence_local(pii,dmft_error,Nsuccess,nloop,index=1,total=3,id=0,file="error.err")
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     if(mpiID==0)call end_loop()
  enddo
  !+-------------------------------------+!

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains


  subroutine print_sc_out(converged)
    integer                        :: i,j,is,row,col,unit
    real(8)                        :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)        :: cdwii,rii,sii,zii
    real(8),dimension(Nside,Nside) :: dij,nij,cij,pij,eij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)      :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev,Eout
    logical                        :: converged
    complex(8),dimension(2,Lmats)  :: aGmats,aSmats
    complex(8),dimension(2,Lreal)  :: aGreal,aSreal
    character(len=50)              :: suffix,cloop,cfoo


    if(mpiID==0)then
       suffix=".ed"
       write(cfoo,'(I4.4)')iloop
       cloop = "_loop"//trim(cfoo)//".ed"

       !Get CDW "order parameter"
       do is=1,Nlat
          row=irow(is)
          col=icol(is)
          cdwii(is) = (-1.d0)**(row+col)*(nii(is,1)-1.d0)
       enddo

       nimp = sum(nii(:,1))/dble(Nlat)
       phi  = sum(pii(:,1))/dble(Nlat)
       docc = sum(dii(:,1))/dble(Nlat)
       ccdw = sum(cdwii)/dble(Nlat)
       print*,"<nimp>  =",nimp
       print*,"<phi>   =",phi
       print*,"<docc>  =",docc
       print*,"<ccdw>  =",ccdw

       call splot("nVSiloop.ed",iloop,nimp,append=.true.)
       call splot("phiVSiloop.ed",iloop,phi,append=.true.)
       call splot("doccVSiloop.ed",iloop,docc,append=.true.)
       call splot("ccdwVSiloop.ed",iloop,ccdw,append=.true.)
       call store_data("nVSisite"//trim(suffix),nii(:,1),(/(dble(i),i=1,Nlat)/))
       call store_data("phiVSisite"//trim(suffix),pii(:,1),(/(dble(i),i=1,Nlat)/))
       call store_data("doccVSisite"//trim(suffix),dii(:,1),(/(dble(i),i=1,Nlat)/))
       call store_data("epotVSisite"//trim(suffix),eii(:,1),(/(dble(i),i=1,Nlat)/))
       call store_data("cdwVSisite"//trim(suffix),cdwii,(/(dble(i),i=1,Nlat)/))


       !<DEBUG: to be removed or moved under converged section below
       ! call store_data("nVSisite"//trim(cloop),nii)
       ! call store_data("phiVSisite"//trim(cloop),pii)
       ! call store_data("LG_iw"//trim(suffix),wm(1:Lmats),Gmats(1,1:Nlat,1:Lmats))
       ! call store_data("LF_iw"//trim(suffix),wm(1:Lmats),Gmats(2,1:Nlat,1:Lmats))
       ! call store_data("LG_realw"//trim(suffix),wr(1:Lreal),Greal(1,1:Nlat,1:Lreal))
       ! call store_data("LF_realw"//trim(suffix),wr(1:Lreal),Greal(2,1:Nlat,1:Lreal))
       !>DEBUG


       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then

          Eout = ed_kinetic_energy_lattice(Hk,[1d0],Smats(1,:,:),Smats(2,:,:))
          unit=free_unit()
          open(unit,file='internal_energy.ed')
          write(unit,'(10(F18.10))')Eout(1),sum(eii(:,1))/Nlat,sum(eii(:,2))/Nlat,sum(eii(:,3))/Nlat,sum(eii(:,4))/Nlat,Eout(2)
          close(unit)
          open(unit,file='internal_energy_info.ed')
          write(unit,"(A,10A18)")"#","1<K>","2<Hi>","3<V=Hi-Ehf>","4<Eloc>","5<Ehf>","6<E0>"
          close(unit)


          call store_data("LG_iw"//trim(suffix),Gmats(1,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LF_iw"//trim(suffix),Gmats(2,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LG_realw"//trim(suffix),Greal(1,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LF_realw"//trim(suffix),Greal(2,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSigma_iw"//trim(suffix),Smats(1,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSelf_iw"//trim(suffix),Smats(2,1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSigma_realw"//trim(suffix),Sreal(1,1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSelf_realw"//trim(suffix),Sreal(2,1:Nlat,1:Lreal),wr(1:Lreal))

          ! call store_data("LDelta_iw"//trim(suffix),wm(1:Lmats),Delta(1,1:Nlat,1:Lmats))
          ! call store_data("LGamma_iw"//trim(suffix),wm(1:Lmats),Delta(2,1:Nlat,1:Lmats))

          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Nlat
             row=irow(is)
             col=icol(is)
             sii(is)   = dimag(Smats(1,is,1))-&
                  wm(1)*(dimag(Smats(1,is,2))-dimag(Smats(1,is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(Gmats(1,is,1))-&
                  wm(1)*(dimag(Gmats(1,is,2))-dimag(Gmats(1,is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(1,is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)
          do row=1,Nside
             grid_x(row)=row
             grid_y(row)=row
             do col=1,Nside
                i            = ij2site(row,col)
                nij(row,col) = nii(i,1)
                dij(row,col) = dii(i,1)
                pij(row,col) = pii(i,1)
                eij(row,col) = eii(i,1)
             enddo
          enddo

          call store_data("rhoVSisite"//trim(suffix),rii,(/(dble(i),i=1,Nlat)/))
          call store_data("sigmaVSisite"//trim(suffix),sii,(/(dble(i),i=1,Nlat)/))
          call store_data("zetaVSisite"//trim(suffix),zii,(/(dble(i),i=1,Nlat)/))
          call splot3d("3d_nVSij"//trim(suffix),grid_x,grid_y,nij)
          call splot3d("3d_doccVSij"//trim(suffix),grid_x,grid_y,dij)
          call splot3d("3d_phiVSij"//trim(suffix),grid_x,grid_y,pij)
          call splot3d("3d_epotVSij"//trim(suffix),grid_x,grid_y,eij)


          !Plot averaged local functions
          aGmats    = sum(Gmats,dim=2)/real(Nlat,8)
          aSmats = sum(Smats,dim=2)/real(Nlat,8)

          aGreal    = sum(Greal,dim=2)/real(Nlat,8)
          aSreal = sum(Sreal,dim=2)/real(Nlat,8)

          call splot("aSigma_iw"//trim(suffix),wm,aSmats(1,:))
          call splot("aSelf_iw"//trim(suffix),wm,aSmats(2,:))
          call splot("aSigma_realw"//trim(suffix),wr,aSreal(1,:))
          call splot("aSelf_realw"//trim(suffix),wr,aSreal(2,:))

          call splot("aG_iw"//trim(suffix),wm,aGmats(1,:))
          call splot("aF_iw"//trim(suffix),wm,aGmats(2,:))
          call splot("aG_realw"//trim(suffix),wr,aGreal(1,:))
          call splot("aF_realw"//trim(suffix),wr,aGreal(2,:))


          call get_moments(nii(:,1),mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot("statistics.n"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(dii(:,1),mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.docc"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(pii(:,1),mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.phi"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot("statistics.cdwn"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(zii,mean,sdev,var,skew,kurt)
          call splot("statistics.zeta"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(sii,mean,sdev,var,skew,kurt)
          call splot("statistics.sigma"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(rii,mean,sdev,var,skew,kurt)
          call splot("statistics.rho"//trim(suffix),mean,sdev,var,skew,kurt)

          data_covariance(1,:)=nii(:,1)
          data_covariance(2,:)=pii(:,1)
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file="covariance_n.phi"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          covariance_nd=0d0
          do i=1,2
             do j=1,2
                if(data_sdev(i)/=0d0.AND.data_sdev(j)/=0d0)then
                   covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
                endif
             enddo
          enddo
          open(10,file="correlation_n.phi"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
       end if

    end if
  end subroutine print_sc_out



  ! !******************************************************************
  ! !******************************************************************



  ! subroutine search_mu(convergence)
  !   integer, save         ::nindex
  !   integer               ::nindex1
  !   real(8)               :: naverage,ndelta1
  !   logical,intent(inout) :: convergence
  !   if(mpiID==0)then
  !      naverage=sum(nii)/dble(Nlat)
  !      nindex1=nindex
  !      ndelta1=rdmft_ndelta
  !      if((naverage >= rdmft_nread+rdmft_nerror))then
  !         nindex=-1
  !      elseif(naverage <= rdmft_nread-rdmft_nerror)then
  !         nindex=1
  !      else
  !         nindex=0
  !      endif
  !      if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !         rdmft_ndelta=ndelta1/2.d0
  !      else
  !         rdmft_ndelta=ndelta1
  !      endif
  !      xmu=xmu+real(nindex,8)*rdmft_ndelta
  !      write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",rdmft_nread,"| shift=",nindex*rdmft_ndelta,"| mu=",xmu
  !      write(*,"(A,f15.12,A,f15.12)")"Density Error:",abs(naverage-nread),'/',rdmft_nerror
  !      print*,""
  !      if(abs(naverage-rdmft_nread)>rdmft_nerror)convergence=.false.
  !      call splot("muVSiter.ed",xmu,abs(naverage-rdmft_nread),append=.true.)
  !   endif
  !   call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  ! end subroutine search_mu


end program ed_ahm_disorder
