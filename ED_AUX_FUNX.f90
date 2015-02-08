!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none
  private

  interface print_state_vector
     module procedure print_state_vector_ivec,print_state_vector_int
  end interface print_state_vector

  interface set_Hloc
     module procedure set_Hloc_1,set_Hloc_2,set_Hloc_3d,set_Hloc_3c
  end interface set_Hloc

  interface get_Hloc
     module procedure get_Hloc_1,get_Hloc_2
  end interface get_Hloc

  interface get_sigma
     module procedure get_sigma_1,get_sigma_2,get_sigma_1_,get_sigma_2_
  end interface get_sigma

  interface set_sigma
     module procedure set_sigma_1,set_sigma_2,set_sigma_1_,set_sigma_2_
  end interface set_sigma


  public :: print_Hloc
  public :: set_Hloc
  public :: get_Hloc  
  !
  public :: get_sigma
  public :: set_sigma
  !
  public :: init_ed_structure
  public :: search_chemical_potential
  !
  public :: setup_pointers
  public :: setup_pointers_sc
  public :: build_sector
  public :: bdecomp
  public :: bjoin
  public :: flip_state
  public :: print_state_vector
  public :: c,cdg
  public :: binary_search
  public :: twin_sector_order
  public :: get_twin_sector


  !LATTICE:
  public :: blocks_to_matrix
  public :: matrix_to_blocks
  public :: select_block
  public :: stride_index
  public :: lat_orb2lo
  public :: get_independent_sites  

  !OBSOLETE (to be removed)
  public :: get_lattice_hamiltonian      !+- produce the 2D square/ slab-chain Hamiltonian -+! 


contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure(Hunit)
    character(len=64)                        :: Hunit
    logical                                  :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: imHloc         !local hamiltonian, imag part
    integer                                  :: i,NP,nup,ndw,iorb,jorb,ispin,jspin
    !
    !Norb=# of impurity orbitals
    !Nbath=# of bath sites (per orbital or not depending on bath_type)
    !Ns=total number of sites
    !Nbo=total number of bath sites (all sites - impurity sites)
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
    end select
    Nbo   = Ns-Norb
    Ntot  = 2*Ns
    NN    = 2**Ntot
    !
    nup=Ns/2
    ndw=Ns-nup
    if(.not.ed_supercond)then
       Nsect = (Ns+1)*(Ns+1)
       NP=get_sector_dimension(nup,ndw)
    else
       Nsect = Ntot+1
       NP=get_sc_sector_dimension(0)
    endif
    !
    if(ED_MPI_ID==0)then
       write(LOGfile,*)"Summary:"
       write(LOGfile,*)"--------------------------------------------"
       write(LOGfile,*)'Number of impurities         = ',Norb
       write(LOGfile,*)'Number of bath/impurity      = ',Nbath
       write(LOGfile,*)'Total # of Bath sites/spin   = ',Nbo
       write(LOGfile,*)'Total # of sites/spin        = ',Ns
       write(LOGfile,*)'Maximum dimension            = ',NP
       write(LOGfile,*)'Total size, Hilber space dim.= ',Ntot,NN
       write(LOGfile,*)'Number of sectors            = ',Nsect
       write(LOGfile,*)"--------------------------------------------"
    endif

    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0.d0
    imHloc = 0.d0

    inquire(file=Hunit,exist=control)
    if(control)then
       if(ED_MPI_ID==0)write(LOGfile,*)"Reading impHloc from file: "//Hunit
       open(50,file=Hunit,status='old')
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       close(50)
    else
       if(ED_MPI_ID==0)then
          write(LOGfile,*)"impHloc file not found."
          write(LOGfile,*)"impHloc should be defined elsewhere..."
       endif
    endif
    impHloc = dcmplx(reHloc,imHloc)
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")"H_local:"
       call print_Hloc(impHloc)
    endif



    allocate(impIndex(Norb,2))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect),getsz(Nsect),twin_mask(Nsect))
    if(.not.ed_supercond)then
       allocate(getsector(0:Ns,0:Ns))
    else
       allocate(getsector(-Ns:Ns,1))
    endif
    allocate(getCsector(2,Nsect))
    allocate(getCDGsector(2,Nsect))
    allocate(getBathStride(Norb,Nbath))
    allocate(neigen_sector(Nsect))


    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       lanc_nstates_sector=1            !set the required eigen per sector to 1 see later for neigen_sector
       finiteT=.false.          !set to do zero temperature calculations
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif


    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          if(ED_MPI_ID==0)write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          if(ED_MPI_ID==0)write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif

    endif

    if(finiteT)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
    else
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
    endif

    !Some check:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch..."
    if(Norb>3)stop "Norb > 3 ERROR. ask developer or develop your own on separate branch..." 
    if(nerr < dmft_error) nerr=dmft_error
    if(ed_supercond)then
       if(Nspin>1)stop "SC+AFM ERROR. ask developer or develop your own on separate branch..." 
       if(Norb>1)stop "SC Multi-Band not yet implemented. Wait for the developer to understand what to do..."
       if(ed_type=='c')stop "SC with Hermitian H not yet implemented. Wait for the developer to code it..."
       if(ed_twin)stop  "SC reduction with twim sectors not yet implemented. Ask developer or do it your own..."
    endif


    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
       !nloop=(i-1)*niter                !increase the max number of dmft loop allowed so to do threshold loop
       !write(LOGfile,"(A,I10)")"Increased Nloop to:",nloop
    endif
    if(Nspin>1.AND.ed_twin.eqv..true.)write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"

    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal))

    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impFmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impFreal(Nspin,Nspin,Norb,Norb,Lreal))


    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_phisc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))

  end subroutine init_ed_structure






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hloc(hloc,unit)
    integer,optional                            :: unit
    integer                                     :: iorb,jorb,ispin,jspin
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               (&
               (&
               '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
               jorb =1,Norb),&
               jspin=1,Nspin)
       enddo
    enddo
    if(present(unit))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine set_Hloc_1(hloc,ispin)
    complex(8),dimension(:,:) :: hloc
    integer                   :: ispin
    if(size(hloc,1)/=Norb.OR.size(hloc,2)/=Norb)stop "set_impHloc error: wrong dimensions of Hloc"
    impHloc(ispin,ispin,1:Norb,1:Norb) = Hloc
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_1
  !
  subroutine set_Hloc_2(hloc)
    complex(8),dimension(:,:,:,:) :: hloc
    if(size(hloc,1)/=Nspin.OR.size(hloc,2)/=Nspin)stop "set_impHloc error: wrong Nspin dimensions of Hloc"
    if(size(hloc,3)/=Norb.OR.size(hloc,4)/=Norb)stop "set_impHloc error: wrong Norb dimensions of Hloc"
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = Hloc
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_2
  !
  subroutine set_Hloc_3d(hloc)
    real(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_3d
  !
  subroutine set_Hloc_3c(hloc)
    complex(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_3c



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine get_Hloc_1(hloc,ispin)
    complex(8),dimension(:,:) :: hloc
    integer                   :: ispin
    if(size(hloc,1)/=Norb.OR.size(hloc,2)/=Norb)stop "set_impHloc error: wrong dimensions of Hloc"
    Hloc = impHloc(ispin,ispin,1:Norb,1:Norb)
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine get_Hloc_1
  !
  subroutine get_Hloc_2(hloc)
    complex(8),dimension(:,:,:,:) :: hloc
    if(size(hloc,1)/=Nspin.OR.size(hloc,2)/=Nspin)stop "set_impHloc error: wrong Nspin dimensions of Hloc"
    if(size(hloc,3)/=Norb.OR.size(hloc,4)/=Norb)stop "set_impHloc error: wrong Norb dimensions of Hloc"
    Hloc = impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb)
    write(LOGfile,*)"Updated impHloc:"
    call print_Hloc(impHloc)
  end subroutine get_Hloc_2




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine set_sigma_1(sigma,ispin)
    complex(8),dimension(:,:,:) :: sigma
    integer                     :: ispin
    if(size(sigma,1)/=Norb.OR.size(sigma,2)/=Norb)stop "set_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,3)/=Lmats.OR.size(sigma,3)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigma,3)==Lmats) then
       impSmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) = sigma(1:Norb,1:Norb,1:Lmats)
    else
       impSreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) = sigma(1:Norb,1:Norb,1:Lreal)
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine set_sigma_1
  !
  subroutine set_sigma_2(sigma)
    complex(8),dimension(:,:,:,:,:) :: sigma
    if(size(sigma,1)/=Nspin.OR.size(sigma,2)/=Nspin)stop "set_impSigma error: wrong Nspin dimensions of Sigma"
    if(size(sigma,3)/=Norb.OR.size(sigma,4)/=Norb)stop "set_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,5)/=Lmats.OR.size(sigma,5)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigma,5)==Lmats) then
       impSmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) = sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    else
       impSreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) = sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine set_sigma_2
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine set_sigma_1_(sigma,sigmaA,ispin)
    complex(8),dimension(:,:,:) :: sigma,sigmaA
    integer                     :: ispin
    if(size(sigma,1)/=Norb.OR.size(sigma,2)/=Norb)stop "set_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,3)/=Lmats.OR.size(sigma,3)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigmaA,1)/=Norb.OR.size(sigmaA,2)/=Norb)stop "set_impSigma error: wrong Norb dimensions of SigmaA"
    if(size(sigmaA,3)/=Lmats.OR.size(sigmaA,3)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of SigmaA"
    !
    if(size(sigma,3)==Lmats) then
       impSmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) = sigma(1:Norb,1:Norb,1:Lmats)
    else
       impSreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) = sigma(1:Norb,1:Norb,1:Lreal)
    end if
    !
    if(size(sigmaA,3)==Lmats) then
       impSAmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) = sigmaA(1:Norb,1:Norb,1:Lmats)
    else
       impSAreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) = sigmaA(1:Norb,1:Norb,1:Lreal)
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine set_sigma_1_
  !
  subroutine set_sigma_2_(sigma,sigmaA)
    complex(8),dimension(:,:,:,:,:) :: sigma,sigmaA
    if(size(sigma,1)/=Nspin.OR.size(sigma,2)/=Nspin) stop "set_impSigma error: wrong Nspin dimensions of Sigma"
    if(size(sigma,3)/=Norb.OR.size(sigma,4)/=Norb)   stop "set_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,5)/=Lmats.OR.size(sigma,5)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigmaA,1)/=Nspin.OR.size(sigmaA,2)/=Nspin) stop "set_impSigma error: wrong Nspin dimensions of SigmaA"
    if(size(sigmaA,3)/=Norb.OR.size(sigmaA,4)/=Norb)   stop "set_impSigma error: wrong Norb dimensions of SigmaA"
    if(size(sigmaA,5)/=Lmats.OR.size(sigmaA,5)/=Lreal) stop "set_impSigma error: wrong Lfreq dimensions of SigmaA"
    !
    if(size(sigma,5)==Lmats) then
       impSmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) = sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    else
       impSreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) = sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)
    end if
    !
    if(size(sigmaA,5)==Lmats) then
       impSAmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) = sigmaA(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    else
       impSAreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) = sigmaA(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine set_sigma_2_








  !+----------------------------------+!
  !+ PURPOSE
  !+----------------------------------+!
  subroutine get_sigma_1(sigma,ispin)
    complex(8),dimension(:,:,:) :: sigma
    integer                     :: ispin
    if(size(sigma,1)/=Norb.OR.size(sigma,2)/=Norb)stop "get_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,3)/=Lmats.OR.size(sigma,3)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigma,3)==Lmats) then
       sigma(1:Norb,1:Norb,1:Lmats)=impSmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) 
    else
       sigma(1:Norb,1:Norb,1:Lreal)=impSreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) 
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine get_sigma_1
  !
  subroutine get_sigma_2(sigma)
    complex(8),dimension(:,:,:,:,:) :: sigma
    if(size(sigma,1)/=Nspin.OR.size(sigma,2)/=Nspin)stop "get_impSigma error: wrong Nspin dimensions of Sigma"
    if(size(sigma,3)/=Norb.OR.size(sigma,4)/=Norb)stop "get_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,5)/=Lmats.OR.size(sigma,5)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigma,5)==Lmats) then
       sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)=impSmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) 
    else
       sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)=impSreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) 
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine get_sigma_2
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine get_sigma_1_(sigma,sigmaA,ispin)
    complex(8),dimension(:,:,:) :: sigma,sigmaA
    integer                     :: ispin
    if(size(sigma,1)/=Norb.OR.size(sigma,2)/=Norb)stop "get_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,3)/=Lmats.OR.size(sigma,3)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigmaA,1)/=Norb.OR.size(sigmaA,2)/=Norb)stop "get_impSigma error: wrong Norb dimensions of SigmaA"
    if(size(sigmaA,3)/=Lmats.OR.size(sigmaA,3)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of SigmaA"
    !
    if(size(sigma,3)==Lmats) then
       sigma(1:Norb,1:Norb,1:Lmats)=impSmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) 
    else
       sigma(1:Norb,1:Norb,1:Lreal)=impSreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) 
    end if
    !
    if(size(sigmaA,3)==Lmats) then
       sigmaA(1:Norb,1:Norb,1:Lmats)=impSAmats(ispin,ispin,1:Norb,1:Norb,1:Lmats) 
    else
       sigmaA(1:Norb,1:Norb,1:Lreal)=impSAreal(ispin,ispin,1:Norb,1:Norb,1:Lreal) 
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine get_sigma_1_
  !
  subroutine get_sigma_2_(sigma,sigmaA)
    complex(8),dimension(:,:,:,:,:) :: sigma,sigmaA
    if(size(sigma,1)/=Nspin.OR.size(sigma,2)/=Nspin) stop "get_impSigma error: wrong Nspin dimensions of Sigma"
    if(size(sigma,3)/=Norb.OR.size(sigma,4)/=Norb)   stop "get_impSigma error: wrong Norb dimensions of Sigma"
    if(size(sigma,5)/=Lmats.OR.size(sigma,5)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of Sigma"
    if(size(sigmaA,1)/=Nspin.OR.size(sigmaA,2)/=Nspin) stop "get_impSigma error: wrong Nspin dimensions of SigmaA"
    if(size(sigmaA,3)/=Norb.OR.size(sigmaA,4)/=Norb)   stop "get_impSigma error: wrong Norb dimensions of SigmaA"
    if(size(sigmaA,5)/=Lmats.OR.size(sigmaA,5)/=Lreal) stop "get_impSigma error: wrong Lfreq dimensions of SigmaA"
    !
    if(size(sigma,5)==Lmats) then
       sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)=impSmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) 
    else
       sigma(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)=impSreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) 
    end if
    !
    if(size(sigmaA,5)==Lmats) then
       sigmaA(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)=impSAmats(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats) 
    else
       sigmaA(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal)=impSAreal(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lreal) 
    end if
    write(LOGfile,*)"Updated impSigma"
  end subroutine get_sigma_2_








  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: i,in,dim,isector,jsector,dimup,dimdw
    integer                          :: nup,ndw,jup,jdw,iorb
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getsector(nup,ndw)=isector
          getnup(isector)=nup
          getndw(isector)=ndw
          dim = get_sector_dimension(nup,ndw)
          getdim(isector)=dim
          neigen_sector(isector) = min(dim,lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    enddo
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsect
          nup=getnup(isector)
          ndw=getndw(isector)
          if(nup<ndw)twin_mask(isector)=.false.
       enddo
       if(ED_MPI_ID==0)write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsect
    endif
    if(ED_MPI_ID==0)call stop_timer

    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo

    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)       = Norb + i
       enddo
    end select

    getCsector=0
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo

    getCDGsector=0
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers


  subroutine setup_pointers_sc
    integer                          :: i,isz,in,dim,isector,jsector
    integer                          :: sz,iorb,dim2,jsz
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do isz=-Ns,Ns
       sz=abs(isz)
       isector=isector+1
       getsector(isz,1)=isector
       getsz(isector)=isz
       dim = get_sc_sector_dimension(isz)
       getdim(isector)=dim
       neigen_sector(isector) = min(dim,lanc_nstates_sector)   !init every sector to required eigenstates
    enddo
    twin_mask=.true.
    !<TODO 
    !build the twin sector statements in the Superconducting channel.
    !>TODO
    if(ED_MPI_ID==0)call stop_timer

    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo

    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)      = Norb + i
       enddo
    end select

    getCsector=0
    !c_up
    do isector=1,Nsect
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCsector(1,isector)=jsector
    enddo
    !c_dw
    do isector=1,Nsect
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCsector(2,isector)=jsector
    enddo

    getCDGsector=0
    !cdg_up
    do isector=1,Nsect
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCDGsector(1,isector)=jsector
    enddo
    !cdg_dw
    do isector=1,Nsect
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_sc






  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP>|ImpDW,BathDW >
  subroutine build_sector(isector,map,dim2)
    integer              :: i,j,isector,iup,idw,mz,dim
    integer,optional     :: dim2
    integer              :: nup,ndw,sz
    integer              :: ivec(Ntot)
    integer,dimension(:) :: map
    !if(size(map)/=getdim(isector))stop "error in build_sector: wrong dimension of map"
    dim=0
    if(.not.ed_supercond)then
       nup = getnup(isector)
       ndw = getndw(isector)
       do i=1,NN
          call bdecomp(i,ivec)
          iup = sum(ivec(1:Ns))
          idw = sum(ivec(Ns+1:2*Ns))
          if(iup==nup.AND.idw==ndw)then
             dim           = dim+1 !count the states in the sector (n_up,n_dw)
             map(dim)      = i       !build the map to full space states
          endif
       enddo
    else
       sz = getsz(isector)
       do i=1,NN
          call bdecomp(i,ivec)
          mz = sum(ivec(1:Ns)) - sum(ivec(Ns+1:2*Ns))
          if(mz==sz)then
             dim             = dim+1 !count the states in the sector (n_up,n_dw)
             map(dim)        = i       !build the map to full space states
          endif
       enddo
    endif
    if(present(dim2))dim2=dim
  end subroutine build_sector







  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Ntot)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bdecomp(i,ivec)
    integer :: ivec(Ntot)         
    integer :: l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end subroutine bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Ntot) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bjoin(ivec,i)
    integer,dimension(ntot) :: ivec
    integer                 :: i,j
    i=1
    do j=1,Ntot
       i=i+ivec(j)*2**(j-1)
    enddo
  end subroutine bjoin












  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw. 
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                          :: isector
    integer,dimension(:)             :: order
    integer,dimension(:),allocatable :: Hmap
    integer                          :: i,dim
    dim = getdim(isector)
    if(size(Order)/=dim)stop "ED_AUX_FUNX/build_twin_sector: wrong dimensions."
    allocate(Hmap(dim))
    call build_sector(isector,Hmap) !build the map from the A-sector to \HHH
    do i=1,dim                      !
       Order(i)=flip_state(Hmap(i)) !get the list of states in \HHH corresponding to sector B twin of A
    enddo                           !
    call sort_array(Order)          !return the ordering of B-states in \HHH with respect to those of A
  end subroutine twin_sector_order
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array





  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==0)then
       j=0
    else
       if(m==1)then
          j=i-1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          !km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i-2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine c



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==1)then
       j=0
    else
       if(m==1)then
          j=i+1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          !km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i+2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine cdg





  !+------------------------------------------------------------------+
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ib,unit)
    integer,optional :: unit
    integer :: i,j,unit_
    integer :: ib(ntot)
    unit_=6;if(present(unit))unit_=unit
    call bjoin(ib,i)
    write(unit_,"(I3,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
    write(unit_,"(A1,A1)",advance="no")">","|"
    write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
    write(unit_,"(A1)",advance="yes")">"
  end subroutine print_state_vector_ivec
  !
  subroutine print_state_vector_int(i,unit)
    integer,optional :: unit
    integer :: i,j,unit_
    integer :: ib(ntot)
    unit_=6;if(present(unit))unit_=unit
    call bdecomp(i,ib)
    write(unit_,"(I3,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
    write(unit_,"(A2)",advance="no")">|"
    write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
    write(unit_,"(A1)",advance="yes")">"
  end subroutine print_state_vector_int




  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into 
  ! j=|{dw}>|{up}>
  !+------------------------------------------------------------------+
  function flip_state(m) result(j)
    integer :: m
    integer :: j
    integer :: ivec(Ntot),foo(Ntot)
    call bdecomp(m,ivec)
    foo(1:Ns)=Ivec(Ns+1:2*Ns)
    foo(Ns+1:2*Ns)=Ivec(1:Ns)
    call bjoin(foo,j)
  end function flip_state






  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into 
  ! j=|{dw}>|{up}>
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in) :: isector
    integer :: jsector
    integer :: nup,ndw
    nup=getnup(isector)
    ndw=getndw(isector)
    jsector=getsector(ndw,nup)
  end function get_twin_sector



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial
  !+------------------------------------------------------------------+
  function get_sector_dimension(nup,ndw) result(dim)
    integer :: nup,ndw,dim,dimup,dimdw
    dimup=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
    dimdw=(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
    dim=dimup*dimdw
  end function get_sector_dimension




  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial
  !+------------------------------------------------------------------+
  function get_sc_sector_dimension(mz) result(dim)
    integer :: mz
    integer :: i,dim,Nb
    dim=0
    Nb=Ns-mz
    do i=0,Nb/2 
       dim=dim + 2**(Nb-2*i)*nchoos(ns,Nb-2*i)*nchoos(ns-Nb+2*i,i)
    enddo
  end function get_sc_sector_dimension


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor
  !+------------------------------------------------------------------+
  function nchoos(n1,n2)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*real(n1+1-i,8)/real(i,8)
    enddo
    nchoos = int(xh + 0.5d0)
  end function nchoos




  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine search_chemical_potential(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    if(ED_MPI_ID==0)then
       ndiff=ntmp-nread
       nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
       !
       !check actual value of the density *ntmp* with respect to goal value *nread*
       count=count+1
       totcount=totcount+1
       if(count>2)then
          do i=1,2
             nindex_old(i+1)=nindex_old(i)
          enddo
       endif
       nindex_old(1)=nindex
       !
       if(ndiff >= nth)then
          nindex=-1
       elseif(ndiff <= -nth)then
          nindex=1
       else
          nindex=0
       endif
       !
       ndelta_old=ndelta
       bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
       !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       if(bool)then
          ndelta=ndelta_old*nratio !decreasing the step
       else
          ndelta=ndelta_old
       endif
       !
       if(ndelta_old<1.d-9)then
          ndelta_old=0.d0
          nindex=0
       endif
       !update chemical potential
       var=var+dble(nindex)*ndelta
       !xmu=xmu+dble(nindex)*ndelta
       !
       !Print information
       write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
       if(nindex>0)then
          write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
       elseif(nindex<0)then
          write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
       else
          write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
       endif
       write(LOGfile,"(A,f15.9)")"var  = ",var
       write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
       unit=free_unit()
       open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
       write(unit,*)var,ntmp,ndiff
       close(unit)
       !
       !check convergence within actual threshold
       !if reduce is activetd
       !if density is in the actual threshold
       !if DMFT is converged
       !if threshold is larger than nerror (i.e. this is not last loop)
       bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
       if(bool)then
          nth_magnitude_old=nth_magnitude        !save old threshold magnitude
          nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
          nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
          count=0                                !reset the counter
          converged=.false.                      !reset convergence
          ndelta=ndelta_old*nratio                  !reduce the delta step
          !
       endif
       !
       !if density is not converged set convergence to .false.
       if(abs(ntmp-nread)>nth)converged=.false.
       !
       !check convergence for this threshold
       !!---if smallest threshold-- NO MORE
       !if reduce is active (you reduced the treshold at least once)
       !if # iterations > max number
       !if not yet converged
       !set threshold back to the previous larger one.
       !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
       bool=ireduce.AND.(count>niter).AND.(.not.converged)
       if(bool)then
          ireduce=.false.
          nth=10.d0**(nth_magnitude_old)
       endif
       !
       write(LOGfile,"(A,I5)")"count= ",count
       write(LOGfile,"(A,L2)"),"Converged=",converged
       print*,""
       !
    endif
#ifdef _MPI
    call MPI_BCAST(xmu,1,MPI_Double_Precision,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
  end subroutine search_chemical_potential





  !-------------------------------------------------------------------------------------------
  !PURPOSE:  perform all the reduction and broadcast from Nlat blocks of size Nso=Nspin*Norb 
  ! (spin-orbital blocks) into a large matrix of order Nlat*Nso (just like Hk)
  !-------------------------------------------------------------------------------------------
  function blocks_to_matrix(Vblocks) result(Matrix)
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    integer                                               :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = (ip-1)*Nspin*Norb + 1
       j = ip*Nspin*Norb
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix) result(Vblocks)
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    integer                                               :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = (ip-1)*Nspin*Norb + 1
       j = ip*Nspin*Norb
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks

  function select_block(ip,Matrix) result(Vblock)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Vblock
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    integer                                               :: i,j,ip
    Vblock=zero
    i = (ip-1)*Nspin*Norb + 1
    j = ip*Nspin*Norb
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block





  !-------------------------------------------------------------------------------------------
  !PURPOSE: 
  ! - strid_index: given an array of indices ivec=[i_1,...,i_L] and the corresponding array of 
  !                ranges nvec=[N_1,...,N_L], evaluate the index of the stride corresponding the 
  !                actual ivec with respect to nvec as 
  !                I_stride = i1 + \sum_{k=2}^L (i_k-1)\prod_{l=1}^{k-1}N_k
  !                assuming that the  i_1 .inner. i_2 .... .inner i_L (i_1 fastest index, i_L
  !                slowest index)
  ! - lat_orb2lo: evaluate the stride for the simplest case of Norb and Nlat
  !-------------------------------------------------------------------------------------------
  function stride_index(Ivec,Nvec) result(i)
    integer,dimension(:)          :: Ivec
    integer,dimension(size(Ivec)) :: Nvec
    integer                       :: i,k,Ntmp
    I = Ivec(1)
    do k=2,size(Ivec)
       Ntmp=product(Nvec(1:k-1))
       I = I + (Ivec(k)-1)*Ntmp
    enddo
  end function stride_index

  function lat_orb2lo(ilat,iorb) result(io)
    integer :: ilat,iorb,io
    io=(ilat-1)*Norb + iorb
  end function lat_orb2lo





  !+-------------------------------------------------------+!
  ! PURPOSE: find all inequivalent sites with respect the   !
  ! user defined symmetry operations                        !
  ! Build and check maps from the full(independent) lattice !
  ! to the independent(full) lattice                        !
  !+-------------------------------------------------------+!
  subroutine get_independent_sites(symmetry_operations)
    integer                       :: i,row,col,unit,isymm
    integer,dimension(Nlat)       :: tmp_search
    integer,dimension(Nlat,Nsymm) :: tmp_map
    integer                       :: i_ind,check_maps
    character(len=5)              :: tmp_suffix
    interface
       function symmetry_operations(site_in) result(sites_out)
         implicit none
         integer                  :: site_in
         integer,allocatable      :: sites_out(:)
       end function symmetry_operations
    end interface
    !+- search and get number of independent sites -+!
    tmp_search=0
    i_ind=0
    do i=1,Nlat
       tmp_map(i,:)=symmetry_operations(i)
       if(tmp_search(i).ge.0) then
          i_ind=i_ind+1
          tmp_search(i)=i
          do isymm=1,Nsymm
             tmp_search(tmp_map(i,isymm))=-1
          end do
       end if
    end do
    !
    Nindep=i_ind
    ! (remember: each site is connected with Nsymm sites (+ 1 = identity)) !
    allocate(indep_list(Nindep),map_lat2ind(Nlat),map_ind2lat(Nindep,Nsymm+1))  
    !
    !+- get list of independent sites -+!
    i_ind=0
    unit=free_unit()    
    if(mpiID==0) open(unit,file='independent_sites.lattice')
    do i=1,Nlat
       if(tmp_search(i).ge.0) then
          i_ind=i_ind+1
          indep_list(i_ind) = tmp_search(i)
          if(mpiID==0) write(unit,*) dble(icol(indep_list(i_ind))),dble(irow(indep_list(i_ind))),dble(indep_list(i_ind)) 
       end if
    end do
    if(mpiID==0) close(unit)
    !+-  build maps -+!
    !
    do i_ind=1,Nindep
       map_lat2ind(indep_list(i_ind))=i_ind
       do isymm=1,Nsymm
          map_lat2ind(tmp_map(indep_list(i_ind),isymm))=i_ind
       end do
    end do
    ! 
    do i_ind=1,Nindep
       unit=free_unit()
       write(tmp_suffix,'(I4.4)') i_ind
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(mpiID==0) open(unit,file='equivalents'//trim(tmp_suffix)//'.lattice')
       map_ind2lat(i_ind,1) = indep_list(i_ind)
       if(mpiID==0) write(unit,*) icol(indep_list(i_ind)),irow(indep_list(i_ind))
       do isymm=1,Nsymm
          map_ind2lat(i_ind,isymm+1) = tmp_map(indep_list(i_ind),isymm)
          if(mpiID==0) write(unit,*) icol(tmp_map(indep_list(i_ind),isymm)),irow(tmp_map(indep_list(i_ind),isymm))
       end do
       if(mpiID==0) close(unit)
    end do
    !+- check maps +-!
    do i_ind=1,Nindep
       do isymm=1,Nsymm+1
          check_maps=map_ind2lat(i_ind,isymm)
          if(i_ind /= map_lat2ind(check_maps)) stop "WRONG MAPS"
       end do
    end do
  end subroutine get_independent_sites






  !+--------------------------------------------------------------------------------+!
  ! !!!!! OBSOLETE ROUTINE !!!!!!! Hlattice defined in the driver
  ! THIS IS LEFT ONLY FOR BACK COMPATIBILITY WITH OLD CODES (TO BE REMOVED)
  ! SUPERSEDED BY DMFT_TIGHT_BINDING PROCEDURES GENERATING SQUARE LATTICE TB HAMILTONIAN
  ! VIA DFT OF THE k-space FUNCTION.
  !+--------------------------------------------------------------------------------+!
  subroutine get_lattice_hamiltonian(Nrow,Ncol,pbc_row,pbc_col,ts)
    integer          :: Nrow
    integer,optional :: Ncol
    integer          :: Nsquare
    logical,optional :: pbc_row,pbc_col
    logical          :: pbc_row_,pbc_col_
    logical          :: symm
    integer          :: i,jj,j,k,row,col,link(4)
    integer          :: unit
    real(8),optional :: ts
    real(8)          :: ts_
    !
    ts_=1d0;if(present(ts))ts_=ts
    pbc_row_=.false.;pbc_col_=.false.
    if(pbcflag) then
       pbc_row_=.true.
       pbc_col_=.true.
    end if
    if(present(pbc_row)) pbc_row_=pbc_row
    if(present(pbc_col)) pbc_col_=pbc_col
    !
    allocate(H0(Nlat,Nlat))
    H0 = 0.d0
    unit=free_unit()
    if(mpiID==0) open(unit,file='rdmft_sites.lattice')
    if(present(Ncol)) then 
       !+- 2D LATTICE (NROW x NCOL) -+!
       if(Nlat /= Nrow*Ncol) stop "Nlat != Nrow*Ncol"
       allocate(icol(Nlat),irow(Nlat))
       allocate(ij2site(Nrow,Ncol))
       do row=0,Nrow-1
          do col=0,Ncol-1
             i=col+ 1 + row*Ncol
             !
             irow(i)=row+1
             icol(i)=col+1
             ij2site(row+1,col+1)=i
             !
             if(mpiID==0) write(unit,*) dble(col+1),dble(row+1)
             !right hop
             link(1)= i + 1     
             if((col+1)==Ncol) then
                if(pbc_col_) then
                   link(1)=1+row*Ncol  
                else
                   link(1)=0  
                end if
             end if
             !left  hop
             link(3)= i - 1    
             if((col-1)<0)     then
                if(pbc_col_) then
                   link(3)=Ncol+row*Ncol
                else
                   link(3)=0  
                end if
             end if
             !up    hop
             link(2)= i + Ncol 
             if((row+1)==Nrow) then
                if(pbc_row_) then
                   link(2)=col+1
                else
                   link(2)=0  
                end if
             end if
             !down  hop
             link(4)= i - Ncol 
             if((row-1)<0)     then
                if(pbc_row_) then
                   link(4)=col+1+(Nrow-1)*Ncol
                else
                   link(4)=0  
                end if
             end if
             !
             do jj=1,4
                if(link(jj)>0)H0(i,link(jj))=-ts_ !! ts_ must be negative.
             enddo
             !
          enddo
       enddo
    else
       !+- 1D LATTICE (CHAIN) -+!
       if(Nlat /= Nrow) stop "Nlat != Nrow"
       do i=1,Nrow-1
          H0(i,i+1)=-ts_
          H0(i+1,i)=-ts_
       end do
    end if
    if(mpiID==0) close(unit)
  end subroutine get_lattice_hamiltonian




  ! subroutine search_mu(ntmp,convergence)
  !   logical,intent(inout) :: convergence
  !   real(8)               :: ntmp
  !   logical               :: check
  !   integer,save          :: count=0
  !   integer,save          :: nindex=0
  !   real(8)               :: ndelta1,nindex1
  !   if(count==0)then
  !      inquire(file="searchmu_file.restart",exist=check)
  !      if(check)then
  !         open(10,file="searchmu_file.restart")
  !         read(10,*)ndelta,nindex
  !         close(10)
  !      endif
  !   endif
  !   count=count+1
  !   nindex1=nindex
  !   ndelta1=ndelta
  !   if((ntmp >= nread+nerr))then
  !      nindex=-1
  !   elseif(ntmp <= nread-nerr)then
  !      nindex=1
  !   else
  !      nindex=0
  !   endif
  !   if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !      ndelta=ndelta1/2.d0 !decreasing the step       
  !   else
  !      ndelta=ndelta1
  !   endif
  !   xmu=xmu+real(nindex,8)*ndelta
  !   if(abs(ntmp-nread)>nerr)convergence=.false.
  !   write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp," /",nread,&
  !        "| shift=",nindex*ndelta,"| xmu=",xmu
  !   write(*,"(A,f15.12)")"dn=",abs(ntmp-nread)
  !   print*,""
  !   print*,"Convergence:",convergence
  !   print*,""
  !   open(10,file="searchmu_file.restart.new")
  !   write(10,*)ndelta,nindex,xmu
  !   close(10)
  ! end subroutine search_mu

END MODULE ED_AUX_FUNX
