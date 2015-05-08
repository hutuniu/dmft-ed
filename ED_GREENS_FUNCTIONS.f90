!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: free_unit,reg,free_units,txtfy
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv
  USE PLAIN_LANCZOS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  USE ED_MATVEC
  !
  implicit none
  private 

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer               :: state_vec
  complex(8),dimension(:),pointer            :: state_cvec
  real(8)                                    :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable           :: wm,tau,wr,vm


  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:)      :: auxGmats,auxGreal

  !Poles & Weights 
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:,:,:) :: GFpoles,GFweights


  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)         :: Chitau
  complex(8),allocatable,dimension(:,:)      :: Chiw,Chiiw


  public                                     :: buildgf_impurity
  public                                     :: buildchi_impurity

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Green's function calculation
  !+------------------------------------------------------------------+
  subroutine buildgf_impurity()
    call allocate_grids
    if(.not.allocated(impGmats))stop "build_gf_super: impGmats not allocated"
    if(.not.allocated(impGreal))stop "build_gf_super: impGreal not allocated"
    if(.not.allocated(impFmats))stop "build_gf_super: impFmats not allocated"
    if(.not.allocated(impFreal))stop "build_gf_super: impFreal not allocated"
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    if(.not.allocated(impSmats)) stop "build_gf_super: impSmats not allocated"
    if(.not.allocated(impSreal)) stop "build_gf_super: impSreal not allocated"
    if(.not.allocated(impSAmats))stop "build_gf_super: impSAmats not allocated"
    if(.not.allocated(impSAreal))stop "build_gf_super: impSAreal not allocated"    
    impSmats = zero
    impSreal = zero
    impSAmats = zero
    impSAreal = zero
    !
    if(.not.allocated(GFpoles))   allocate(GFpoles(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    if(.not.allocated(GFweights)) allocate(GFweights(Nspin,Nspin,Norb,Norb,2,lanc_nGFiter))
    GFpoles=zero
    GFweights=zero
    !
    select case(ed_mode)
    case default
       call build_gf_normal()
       call get_sigma_print_gf_normal()
    case ("superc")
       call build_gf_superc()
       call get_sigma_print_gf_superc()
    end select
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
    if(allocated(GFpoles))deallocate(GFpoles)
    if(allocated(GFweights))deallocate(GFweights)
  end subroutine buildgf_impurity
  !+------------------------------------------------------------------+
  !                    GREEN'S FUNCTIONS 
  !+------------------------------------------------------------------+
  include 'ed_greens_funcs_build_gf_normal.f90'
  include 'ed_greens_funcs_build_gf_superc.f90'



  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Susceptibility calculation
  !+------------------------------------------------------------------+
  subroutine buildchi_impurity()
    call allocate_grids()
    if(.not.allocated(Chitau))allocate(Chitau(Norb+1,0:Ltau))
    if(.not.allocated(Chiw))  allocate(Chiw(Norb+1,Lreal))
    if(.not.allocated(Chiiw)) allocate(Chiiw(Norb+1,0:Lmats))
    Chitau=zero
    Chiw=zero
    Chiiw=zero
    !
    call build_chi_spin()
    call print_chi_spin()
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
    if(allocated(Chitau))deallocate(Chitau)
    if(allocated(Chiw))deallocate(Chiw)
    if(allocated(Chiiw))deallocate(Chiiw)
  end subroutine buildchi_impurity
  !+------------------------------------------------------------------+
  !                    SPIN SUSCPTIBILITY
  !+------------------------------------------------------------------+
  include 'ed_greens_funcs_build_chi_spin.f90'





  !+------------------------------------------------------------------+
  !PURPOSE  : Print normal Green's functions
  !+------------------------------------------------------------------+
  subroutine get_sigma_print_gf_normal
    integer                                           :: i,j,ispin,isign,unit(7),iorb,jorb
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats,invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impG0real,invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp,impG0
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    select case(bath_type)
       !
       !
       !
    case default                !Diagonal in both spin and orbital
       !
       allocate(getIorb(Norb),getJorb(Norb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
       if(totNorb/=Norb)stop "get_sigma_print_gf_normal error counting the orbitals"
       !
       !Get G0^-1
       do i=1,Lmats
          invG0mats(:,:,:,:,i) = invg0_bath_mats(xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          invG0real(:,:,:,:,i) = invg0_bath_real(wr(i)+xi*eps,dmft_bath)
       enddo
       !Get Gimp^-1
       do ispin=1,Nspin
          do iorb=1,Norb
             invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
             invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
             impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get G0and:
       do i=1,Lmats
          impG0mats(:,:,:,:,i) = g0and_bath_mats(xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          impG0real(:,:,:,:,i) = g0and_bath_real(wr(i)+xi*eps,dmft_bath)
       enddo
       !
       !
       !
    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       !
       !
       !
       allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
       totNorb=l
       if(totNorb/=(Norb*(Norb+1)/2))stop "get_sigma_print_gf_normal error counting the orbitals"
       !
       !Get G0^-1
       do i=1,Lmats
          invG0mats(:,:,:,:,i) = invg0_bath_mats(xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          invG0real(:,:,:,:,i) = invg0_bath_real(wr(i)+xi*eps,dmft_bath)
       enddo
       !Get Gimp^-1
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGmats(ispin,ispin,:,:,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGreal(ispin,ispin,:,:,i)=invGimp
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
          !
          impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       enddo
       !Get G0and:
       do i=1,Lmats
          impG0mats(:,:,:,:,i) = g0and_bath_mats(xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          impG0real(:,:,:,:,i) = g0and_bath_real(wr(i)+xi*eps,dmft_bath)
       enddo
       !
       !
    end select
    !!
    !!
    !!
    !!
    !!
    !Print the impurity functions:
    if(ED_MPI_ID==0)then
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
          call open_units(reg(suffix))
          if(ed_verbose<4)then
             do i=1,Lmats
                write(unit(1),"(F26.15,6(F26.15))")wm(i),(dimag(impSmats(ispin,ispin,iorb,jorb,i)),dreal(impSmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(2),"(F26.15,6(F26.15))")wr(i),(dimag(impSreal(ispin,ispin,iorb,jorb,i)),dreal(impSreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do isign=1,2
                do i=1,lanc_nGFiter
                   write(unit(3),"(6(F26.15,1x))")(GFpoles(ispin,ispin,iorb,iorb,isign,i),GFweights(ispin,ispin,iorb,iorb,isign,i),ispin=1,Nspin)
                enddo
             enddo
          endif
          !
          if(ed_verbose<2)then
             do i=1,Lmats
                write(unit(4),"(F26.15,6(F26.15))")wm(i),(dimag(impGmats(ispin,ispin,iorb,jorb,i)),dreal(impGmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(5),"(F26.15,6(F26.15))")wr(i),(dimag(impGreal(ispin,ispin,iorb,jorb,i)),dreal(impGreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
          endif
          !
          if(ed_verbose<1)then
             do i=1,Lmats
                write(unit(6),"(F26.15,6(F26.15))")wm(i),(dimag(impG0mats(ispin,ispin,iorb,jorb,i)),dreal(impG0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(7),"(F26.15,6(F26.15))")wr(i),(dimag(impG0real(ispin,ispin,iorb,jorb,i)),dreal(impG0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
          endif
          call close_units()
       enddo
    endif
    !
  contains
    !
    subroutine open_units(string)
      character(len=*) :: string
      unit=free_units(size(unit))
      if(ed_verbose<4)then
         open(unit(1),file="impSigma"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(2),file="impSigma"//string//"_realw"//reg(ed_file_suffix)//".ed")
         open(unit(3),file="Gpoles_weights"//string//reg(ed_file_suffix)//".ed")
      endif
      if(ed_verbose<2)then
         open(unit(4),file="impG"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(5),file="impG"//string//"_realw"//reg(ed_file_suffix)//".ed")
      endif
      if(ed_verbose<1)then
         open(unit(6),file="impG0"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(7),file="impG0"//string//"_realw"//reg(ed_file_suffix)//".ed")
      endif
    end subroutine open_units
    !
    subroutine close_units()
      if(ed_verbose<4)then
         close(unit(1))
         close(unit(2))
         close(unit(3))
      endif
      if(ed_verbose<2)then
         close(unit(4))
         close(unit(5))
      endif
      if(ed_verbose<1)then
         close(unit(6))
         close(unit(7))
      endif
    end subroutine close_units
    !
  end subroutine get_sigma_print_gf_normal









  
  !+------------------------------------------------------------------+
  !PURPOSE  : Print Superconducting Green's functions
  !+------------------------------------------------------------------+
  subroutine get_sigma_print_gf_superc
    integer                                               :: i,j,ispin,unit(12),iorb,jorb
    complex(8)                                            :: iw
    complex(8)                                            :: det_mats(Lmats),det_real(Lreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)     :: impG0mats,impF0mats,invG0mats,invF0mats,invGmats,invFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)     :: impG0real,impF0real,invG0real,invF0real,invGreal,invFreal
    complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)       :: invGimp
    character(len=20)                                     :: suffix
    integer,dimension(:),allocatable                      :: getIorb,getJorb
    integer                                               :: totNorb,l
    !
    select case(bath_type)
       !
       !
    case default
       !
       allocate(getIorb(Norb),getJorb(Norb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
       if(totNorb/=Norb)stop "get_sigma_print_gf_superc error counting the orbitals"
       !
       !Get G0^-1,F0^-1
       ispin=1
       do i=1,Lmats
          invG0mats(ispin,ispin,:,:,i) = invg0_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
          invF0mats(ispin,ispin,:,:,i) = invf0_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          invG0real(ispin,ispin,:,:,i) = invg0_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
          invF0real(ispin,ispin,:,:,i) = invf0_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
       enddo
       !Get Gimp^-1
       do iorb=1,Norb
          det_mats  =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          invGmats(ispin,ispin,iorb,iorb,:) = conjg(impGmats(ispin,ispin,iorb,iorb,:))/det_mats
          invFmats(ispin,ispin,iorb,iorb,:) = impFmats(ispin,ispin,iorb,iorb,:)/det_mats
          !
          det_real  = impGreal(ispin,ispin,iorb,iorb,:)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1)) + impFreal(ispin,ispin,iorb,iorb,:)**2
          invGreal(ispin,ispin,iorb,iorb,:) =  conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1))/det_real(:)
          invFreal(ispin,ispin,iorb,iorb,:) =  impFreal(ispin,ispin,iorb,iorb,:)/det_real(:)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       do iorb=1,Norb
          impSmats(ispin,ispin,iorb,iorb,:)  = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
          impSAmats(ispin,ispin,iorb,iorb,:) = invF0mats(ispin,ispin,iorb,iorb,:) - invFmats(ispin,ispin,iorb,iorb,:)
          !
          impSreal(ispin,ispin,iorb,iorb,:)  = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          impSAreal(ispin,ispin,iorb,iorb,:) = invF0real(ispin,ispin,iorb,iorb,:) - invFreal(ispin,ispin,iorb,iorb,:)
       enddo
       !Get G0and:
       do i=1,Lmats
          impG0mats(ispin,ispin,:,:,i) = g0and_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
          impF0mats(ispin,ispin,:,:,i) = f0and_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          impG0real(ispin,ispin,:,:,i) = g0and_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
          impF0real(ispin,ispin,:,:,i) = f0and_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
       enddo
       !
       !
       !
    case ("hybrid")
       !
       !
       !
       allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
       totNorb=l
       if(totNorb/=(Norb*(Norb+1)/2))stop "get_sigma_print_gf_superc error counting the orbitals"
       !
       !
       !Get G0^-1,F0^-1
       ispin=1
       do i=1,Lmats
          invG0mats(ispin,ispin,:,:,i) = invg0_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
          invF0mats(ispin,ispin,:,:,i) = invf0_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          invG0real(ispin,ispin,:,:,i) = invg0_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
          invF0real(ispin,ispin,:,:,i) = invf0_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
       enddo
       !Get Gimp^-1
       do i=1,Lmats
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGmats(ispin,ispin,iorb,jorb,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFmats(ispin,ispin,iorb,jorb,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFmats(ispin,ispin,iorb,jorb,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGmats(ispin,ispin,iorb,jorb,i))
          call inv(invGimp)
          invGmats(ispin,ispin,:,:,i) = invGimp(1:Norb,1:Norb)
          invFmats(ispin,ispin,:,:,i) = invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       do i=1,Lreal
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGreal(ispin,ispin,iorb,jorb,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFreal(ispin,ispin,iorb,jorb,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFreal(ispin,ispin,iorb,jorb,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGreal(ispin,ispin,iorb,jorb,Lreal-i+1))
          call inv(invGimp)
          invGreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,1:Norb)
          invFreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       !
       impSmats(ispin,ispin,:,:,:)  = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
       impSAmats(ispin,ispin,:,:,:) = invF0mats(ispin,ispin,:,:,:) - invFmats(ispin,ispin,:,:,:)
       !
       impSreal(ispin,ispin,:,:,:)  = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       impSAreal(ispin,ispin,:,:,:) = invF0real(ispin,ispin,:,:,:) - invFreal(ispin,ispin,:,:,:)
       !
       !Get G0and:
       do i=1,Lmats
          impG0mats(ispin,ispin,:,:,i) = g0and_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
          impF0mats(ispin,ispin,:,:,i) = f0and_bath_mats(ispin,ispin,xi*wm(i),dmft_bath)
       enddo
       do i=1,Lreal
          impG0real(ispin,ispin,:,:,i) = g0and_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
          impF0real(ispin,ispin,:,:,i) = f0and_bath_real(ispin,ispin,wr(i)+xi*eps,dmft_bath)
       enddo
       !
       !
    end select
    !!
    !!
    !!
    !!
    !!PRINT OUT GF:
    if(ED_MPI_ID==0)then
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
          call open_units(reg(suffix))
          if(ed_verbose<4)then
             do i=1,Lmats
                write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSmats(ispin,ispin,iorb,jorb,i)),dreal(impSmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(2),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSAmats(ispin,ispin,iorb,jorb,i)),dreal(impSAmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(3),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSreal(ispin,ispin,iorb,jorb,i)),dreal(impSreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSAreal(ispin,ispin,iorb,jorb,i)),dreal(impSAreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
          endif
          if(ed_verbose<2)then
             do i=1,Lmats
                write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impGmats(ispin,ispin,iorb,jorb,i)),dreal(impGmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(6),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impFmats(ispin,ispin,iorb,jorb,i)),dreal(impFmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(7),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impGreal(ispin,ispin,iorb,jorb,i)),dreal(impGreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(8),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impFreal(ispin,ispin,iorb,jorb,i)),dreal(impFreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
          endif
          if(ed_verbose<1)then
             do i=1,Lmats
                write(unit(9),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impG0mats(ispin,ispin,iorb,jorb,i)),dreal(impG0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(10),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impF0mats(ispin,ispin,iorb,jorb,i)),dreal(impF0mats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(11),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impG0real(ispin,ispin,iorb,jorb,i)),dreal(impG0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(12),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impF0real(ispin,ispin,iorb,jorb,i)),dreal(impF0real(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
          endif
          call close_units
       enddo
    endif
    !
    !
  contains
    !
    subroutine open_units(string)
      character(len=*) :: string
      unit=free_units(size(unit))
      if(ed_verbose<4)then
         open(unit(1),file="impSigma"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(2),file="impSelf"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(3),file="impSigma"//string//"_realw"//reg(ed_file_suffix)//".ed")
         open(unit(4),file="impSelf"//string//"_realw"//reg(ed_file_suffix)//".ed")
      endif
      if(ed_verbose<2)then
         open(unit(5),file="impG"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(6),file="impF"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(7),file="impG"//string//"_realw"//reg(ed_file_suffix)//".ed")
         open(unit(8),file="impF"//string//"_realw"//reg(ed_file_suffix)//".ed")
      endif
      if(ed_verbose<1)then
         open(unit(9),file="impG0"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(10),file="impF0"//string//"_iw"//reg(ed_file_suffix)//".ed")
         open(unit(11),file="impG0"//string//"_realw"//reg(ed_file_suffix)//".ed")
         open(unit(12),file="impF0"//string//"_realw"//reg(ed_file_suffix)//".ed")
      endif
    end subroutine open_units
    !
    subroutine close_units()
      if(ed_verbose<4)then
         close(unit(1))
         close(unit(2))
         close(unit(3))
         close(unit(4))
      endif
      if(ed_verbose<2)then
         close(unit(5))
         close(unit(6))
         close(unit(7))
         close(unit(8))
      endif
      if(ed_verbose<1)then
         close(unit(9))
         close(unit(10))
         close(unit(11))
         close(unit(12))
      endif
    end subroutine close_units
    !
  end subroutine get_sigma_print_gf_superc










  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_chi_spin
    integer                               :: i,j,iorb
    integer                               :: unit(3)
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          unit(1)=free_unit()
          open(unit(1),file="Chi_orb"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
          unit(2)=free_unit()
          open(unit(2),file="Chi_orb"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
          unit(3)=free_unit()
          open(unit(3),file="Chi_orb"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
          do i=0,Ltau
             write(unit(1),*)tau(i),chitau(iorb,i)
          enddo
          do i=1,Lreal
             if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,i)),dreal(chiw(iorb,i))
          enddo
          do i=0,Lmats
             write(unit(3),*)vm(i),dimag(chiiw(iorb,i)),dreal(chiiw(iorb,i))
          enddo
          close(unit(1))
          close(unit(2))
          close(unit(3))
       enddo
       if(Norb>1)then
          iorb=Norb+1
          unit(1)=free_unit()
          open(unit(1),file="Chi_tot_tau"//reg(ed_file_suffix)//".ed")
          unit(2)=free_unit()
          open(unit(2),file="Chi_tot_realw"//reg(ed_file_suffix)//".ed")
          unit(3)=free_unit()
          open(unit(3),file="Chi_tot_iv"//reg(ed_file_suffix)//".ed")
          do i=0,Ltau
             write(unit(1),*)tau(i),chitau(iorb,i)
          enddo
          do i=1,Lreal
             if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,i)),dreal(chiw(iorb,i))
          enddo
          do i=0,Lmats
             write(unit(3),*)vm(i),dimag(chiiw(iorb,i)),dreal(chiiw(iorb,i))
          enddo
          close(unit(1))
          close(unit(2))
          close(unit(3))
       endif
    endif
  end subroutine print_chi_spin




  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids




end MODULE ED_GREENS_FUNCTIONS
