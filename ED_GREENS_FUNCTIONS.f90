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
  USE SF_LINALG,  only: matrix_inverse
  !<debug
  USE DMFT_FFTGF
  !>debug
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
  complex(8),allocatable,dimension(:,:)      :: Gaux_mats,Gaux_real


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
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
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
    if(.not.ed_supercond)then
       call build_gf_normal()
       call get_sigma_print_gf_normal()
    else
       call build_gf_superc()
       call get_sigma_print_gf_superc()
    end if
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(GFpoles))deallocate(GFpoles)
    if(allocated(GFweights))deallocate(GFweights)
  end subroutine buildgf_impurity
  !+------------------------------------------------------------------+
  !                    GREEN'S FUNCTIONS 
  !+------------------------------------------------------------------+
  include 'ed_build_gf_normal.f90'
  include 'ed_build_gf_superc.f90'



  !+------------------------------------------------------------------+
  !PURPOSE  : Interface routine for Susceptibility calculation
  !+------------------------------------------------------------------+
  subroutine buildchi_impurity()
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*(2*arange(1,Lmats)-1)
    do i=0,Lmats
       vm(i) = pi/beta*2*i
    enddo
    wr     = linspace(0d0,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)

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
  include 'ed_build_chi_spin.f90'





  !+------------------------------------------------------------------+
  !PURPOSE  : Print normal Green's functions
  !+------------------------------------------------------------------+
  subroutine get_sigma_print_gf_normal
    integer                                           :: i,j,ispin,isign,unit(7),iorb,jorb
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impG0real
    complex(8),dimension(Norb,Norb)                   :: invGimp,impG0
    character(len=20)                                 :: suffix
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !                        !this is ensured by the special *per impurity" bath structure
       !                        !no intra-orbital hoopings
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Lmats
                fg0 = xi*wm(i) + xmu - impHloc(ispin,ispin,iorb,iorb) - delta_bath_mats(ispin,iorb,xi*wm(i),dmft_bath)
                impSmats(ispin,ispin,iorb,iorb,i)= fg0 - one/impGmats(ispin,ispin,iorb,iorb,i)
                impG0mats(ispin,ispin,iorb,iorb,i) = one/fg0
             enddo
             do i=1,Lreal
                fg0 = wr(i) + xmu - impHloc(ispin,ispin,iorb,iorb) - delta_bath_real(ispin,iorb,wr(i)+xi*eps,dmft_bath)
                impSreal(ispin,ispin,iorb,iorb,i)= fg0 - one/impGreal(ispin,ispin,iorb,iorb,i)
                impG0real(ispin,ispin,iorb,iorb,i) = one/fg0
             enddo
          enddo
       enddo
       !
       if(ED_MPI_ID==0)then	
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))
             call open_units(reg(suffix))
             if(ed_verbose<4)then
                do i=1,Lmats
                   write(unit(1),"(F26.15,6(F26.15))")wm(i),(dimag(impSmats(ispin,ispin,iorb,iorb,i)),dreal(impSmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
                do i=1,Lreal
                   write(unit(2),"(F26.15,6(F26.15))")wr(i),(dimag(impSreal(ispin,ispin,iorb,iorb,i)),dreal(impSreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
                do isign=1,2
                   do i=1,lanc_nGFiter
                      write(unit(3),"(6(F26.15,1x))")(GFpoles(ispin,ispin,iorb,iorb,isign,i),GFweights(ispin,ispin,iorb,iorb,isign,i),ispin=1,Nspin)
                   enddo
                   write(unit(3),*)""
                enddo
             endif
             !
             if(ed_verbose<2)then
                do i=1,Lmats
                   write(unit(4),"(F26.15,6(F26.15))")wm(i),(dimag(impGmats(ispin,ispin,iorb,iorb,i)),dreal(impGmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
                do i=1,Lreal
                   write(unit(5),"(F26.15,6(F26.15))")wr(i),(dimag(impGreal(ispin,ispin,iorb,iorb,i)),dreal(impGreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
             endif
             !
             if(ed_verbose<1)then
                do i=1,Lmats
                   write(unit(6),"(F26.15,6(F26.15))")wm(i),(dimag(impG0mats(ispin,ispin,iorb,iorb,i)),dreal(impG0mats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
                do i=1,Lreal
                   write(unit(7),"(F26.15,6(F26.15))")wr(i),(dimag(impG0real(ispin,ispin,iorb,iorb,i)),dreal(impG0real(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
                enddo
             endif
             call close_units
          enddo
       endif

    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       !                        !intra-orbital hopping allow for mixed _ab GF
       do ispin=1,Nspin         !Spin diagona
          do iorb=1,Norb        !Orbital diagonal part GF_0=(iw+mu)_aa-hloc_aa-Delta_aa
             do i=1,Lmats
                impG0mats(ispin,ispin,iorb,iorb,i)= xi*wm(i)+xmu-impHloc(ispin,ispin,iorb,iorb)-delta_bath_mats(ispin,iorb,iorb,xi*wm(i),dmft_bath)
             enddo
             do i=1,Lreal
                impG0real(ispin,ispin,iorb,iorb,i)= wr(i)+xi*eps+xmu-impHloc(ispin,ispin,iorb,iorb)-delta_bath_real(ispin,iorb,iorb,wr(i)+xi*eps,dmft_bath)
             enddo
          enddo
          do iorb=1,Norb         !Orbital non-diagonal part
             do jorb=iorb+1,Norb !GF_0=-hloc_ab-Delta_ab
                do i=1,Lmats
                   impG0mats(ispin,ispin,iorb,jorb,i)= -impHloc(ispin,ispin,iorb,jorb)-delta_bath_mats(ispin,iorb,jorb,xi*wm(i),dmft_bath)
                   impG0mats(ispin,ispin,jorb,iorb,i)= -impHloc(ispin,ispin,jorb,iorb)-delta_bath_mats(ispin,jorb,iorb,xi*wm(i),dmft_bath)
                enddo
                do i=1,Lreal
                   impG0real(ispin,ispin,iorb,jorb,i)= -impHloc(ispin,ispin,iorb,jorb)-delta_bath_real(ispin,iorb,jorb,wr(i)+xi*eps,dmft_bath)
                   impG0real(ispin,ispin,jorb,iorb,i)= -impHloc(ispin,ispin,jorb,iorb)-delta_bath_real(ispin,jorb,iorb,wr(i)+xi*eps,dmft_bath)
                enddo
             enddo
          enddo
       enddo
       !
       !                         !Get Sigma and G_0 by matrix inversions:
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             impG0   = impG0mats(ispin,ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSmats(ispin,ispin,:,:,i) = impG0 - invGimp
             call matrix_inverse(impG0)
             impG0mats(ispin,ispin,:,:,i)=impG0
          enddo
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             impG0   = impG0real(ispin,ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSreal(ispin,ispin,:,:,i) = impG0 - invGimp
             call matrix_inverse(impG0)
             impG0real(ispin,ispin,:,:,i)=impG0
          enddo
       enddo
       !
       !Print the impurity functions:
       if(ED_MPI_ID==0)then
          do iorb=1,Norb
             do jorb=1,Norb
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
          enddo
       endif
    end select

  contains

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

  end subroutine get_sigma_print_gf_normal




  !+------------------------------------------------------------------+
  !PURPOSE  : Print Superconducting Green's functions
  !+------------------------------------------------------------------+
  subroutine get_sigma_print_gf_superc
    integer                                        :: i,j,ispin,unit(12),iorb,jorb
    complex(8)                                     :: iw
    complex(8),allocatable,dimension(:)            :: det
    complex(8),allocatable,dimension(:,:)          :: fg0,fg,sigma
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats,impF0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impG0real,impF0real
    character(len=20)                              :: suffix
    !
    !Diagonal in both spin and orbital
    !this is ensured by the special *per impurity" bath structure
    !no intra-orbital hoopings
    !THIS IS SUPERCONDUCTING CASE
    allocate(fg0(2,Lmats),fg(2,Lmats),det(Lmats))
    do ispin=1,Nspin
       do iorb=1,Norb
          det     =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          fg(1,:) =  conjg(impGmats(ispin,ispin,iorb,iorb,:))/det
          fg(2,:) =  impFmats(ispin,ispin,iorb,iorb,:)/det
          do i=1,LMats
             iw = xi*wm(i)
             fg0(1,i) = iw+xmu-impHloc(ispin,ispin,iorb,iorb)-delta_bath_mats(ispin,iorb,iw,dmft_bath)
             fg0(2,i) = -fdelta_bath_mats(ispin,iorb,iw,dmft_bath)
          enddo
          impSmats(ispin,ispin,iorb,iorb,:)= fg0(1,:) - fg(1,:)
          impSAmats(ispin,ispin,iorb,iorb,:)= fg0(2,:) - fg(2,:)
          det     =  abs(fg0(1,:))**2 + (fg0(2,:))**2
          impG0mats(ispin,ispin,iorb,iorb,:) = conjg(fg0(1,:))/det
          impF0mats(ispin,ispin,iorb,iorb,:) = fg0(2,:)/det
       enddo
    enddo
    deallocate(fg0,fg,det)


    allocate(fg0(2,Lreal),fg(2,Lreal),det(Lreal))
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Lreal
             iw=dcmplx(wr(i),eps)
             !TESTS SHOWS THAT THIS VERSION GIVES THE SAME RESULTS AS THE UNCOMMENTED LINES
             ! det(i)  = impGreal(ispin,ispin,iorb,iorb,i)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal+1-i)) + &
             !      impFreal(ispin,ispin,iorb,iorb,i)*conjg(impFreal(ispin,ispin,iorb,iorb,Lreal+1-i))
             ! fg(1,i) =  conjg(impGreal(ispin,ispin,iorb,iorb,Lreal+1-i))/det(i)
             ! fg(2,i) =  conjg(impFreal(ispin,ispin,iorb,iorb,Lreal+1-i))/det(i)
             det(i)  = -impGreal(ispin,ispin,iorb,iorb,i)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal+1-i)) - &
                  impFreal(ispin,ispin,iorb,iorb,i)*impFreal(ispin,ispin,iorb,iorb,i)
             fg(1,i) =  -conjg(impGreal(ispin,ispin,iorb,iorb,Lreal+1-i))/det(i)
             fg(2,i) =  -impFreal(ispin,ispin,iorb,iorb,i)/det(i)
             fg0(1,i) =  wr(i)+xmu-impHloc(ispin,ispin,iorb,iorb)-delta_bath_real(ispin,iorb,wr(i)+xi*eps,dmft_bath)
             fg0(2,i) = -fdelta_bath_real(ispin,iorb,wr(i)+xi*eps,dmft_bath)
          enddo
          impSreal(ispin,ispin,iorb,iorb,:) = fg0(1,:) - fg(1,:)
          impSAreal(ispin,ispin,iorb,iorb,:)= fg0(2,:) - fg(2,:)
          do i=1,Lreal
             det(i)     =  -fg0(1,i)*conjg(fg0(1,Lreal+1-i)) - fg0(2,i)*fg0(2,i)
             impG0real(ispin,ispin,iorb,iorb,i) = -conjg(fg0(1,Lreal+1-i))/det(i)
             impF0real(ispin,ispin,iorb,iorb,i) = -fg0(2,i)/det(i)
          enddo
       enddo
    enddo
    deallocate(fg0,fg,det)

    !
    if(ED_MPI_ID==0)then
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))
          call open_units(reg(suffix))
          if(ed_verbose<4)then
             do i=1,Lmats
                write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSmats(ispin,ispin,iorb,iorb,i)),dreal(impSmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(2),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSAmats(ispin,ispin,iorb,iorb,i)),dreal(impSAmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(3),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSreal(ispin,ispin,iorb,iorb,i)),dreal(impSreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSAreal(ispin,ispin,iorb,iorb,i)),dreal(impSAreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
          endif
          if(ed_verbose<2)then
             do i=1,Lmats
                write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impGmats(ispin,ispin,iorb,iorb,i)),dreal(impGmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(6),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impFmats(ispin,ispin,iorb,iorb,i)),dreal(impFmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(7),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impGreal(ispin,ispin,iorb,iorb,i)),dreal(impGreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(8),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impFreal(ispin,ispin,iorb,iorb,i)),dreal(impFreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
          endif
          if(ed_verbose<1)then
             do i=1,Lmats
                write(unit(9),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impG0mats(ispin,ispin,iorb,iorb,i)),dreal(impG0mats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lmats
                write(unit(10),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impF0mats(ispin,ispin,iorb,iorb,i)),dreal(impF0mats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(11),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impG0real(ispin,ispin,iorb,iorb,i)),dreal(impG0real(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Lreal
                write(unit(12),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impF0real(ispin,ispin,iorb,iorb,i)),dreal(impF0real(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             enddo
          endif
          call close_units
       enddo
    endif

  contains

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
