MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
  implicit none
  private

  interface print_state_vector
     module procedure print_state_vector_ivec
     module procedure print_state_vector_ivec_ud
     module procedure print_state_vector_int
  end interface print_state_vector


  public :: setup_ed_dimensions
  public :: init_ed_structure
  public :: setup_pointers_normal
  public :: setup_pointers_superc
  public :: setup_pointers_nonsu2
  public :: build_sector
  public :: bdecomp
  public :: bjoin
  public :: flip_state
  public :: print_state_vector
  public :: c,cdg
  public :: binary_search
  public :: twin_sector_order
  public :: get_twin_sector



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
  !  
  !+------------------------------------------------------------------+
  subroutine setup_ed_dimensions()
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
    case ('replica')
       Ns = Norb*(Nbath+1)
    end select
    Nlevels  = 2*Ns
    !
    select case(ed_mode)
    case default
       Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
       Nhel     = 1
    case ("superc")
       Nsectors = Nlevels+1     !sz=-Ns:Ns=2*Ns+1=Nlevels+1
       Nhel     = 1
    case("nonsu2")
       Nsectors = Nlevels+1     !n=0:2*Ns=2*Ns+1=Nlevels+1
       Nhel     = 2
    end select
  end subroutine setup_ed_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure(Hunit)
    character(len=64)                        :: Hunit
    logical                                  :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: imHloc         !local hamiltonian, imag part
    integer                                  :: i,dim_sector_max,iorb,jorb,ispin,jspin
    !
    call setup_ed_dimensions()
    select case(ed_mode)
    case default
       dim_sector_max=get_normal_sector_dimension(nup=Ns/2,ndw=Ns-Ns/2)
    case ("superc")
       dim_sector_max=get_superc_sector_dimension(0)
    case("nonsu2")
       dim_sector_max=get_nonsu2_sector_dimension(Ns)
    end select
    !
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I)")'# of levels/spin      = ',Ns
       write(LOGfile,"(A,I)")'Total size            = ',Nlevels
       write(LOGfile,"(A,I)")'# of impurities       = ',Norb
       write(LOGfile,"(A,I)")'# of bath/impurity    = ',Nbath
       write(LOGfile,"(A,I)")'# of Bath levels/spin = ',Ns-Norb
       write(LOGfile,"(A,I)")'Largest Sector        = ',dim_sector_max
       write(LOGfile,"(A,I)")'Number of sectors     = ',Nsectors
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    !
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0d0 ; imHloc = 0d0
    !
    !Search and read impHloc
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
          call sleep(2)
       endif
    endif
    impHloc = dcmplx(reHloc,imHloc)
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")"H_local:"
       call print_Hloc(impHloc)
    endif
    !
    !
    !Allocate indexing arrays
    allocate(impIndex(Norb,2))
    allocate(getdim(Nsectors),twin_mask(Nsectors))
    allocate(getnup(Nsectors),getndw(Nsectors),getsz(Nsectors),getn(Nsectors))
    select case(ed_mode)
    case default
       allocate(getsector(0:Ns,0:Ns))
    case ("superc")
       allocate(getsector(-Ns:Ns,1))
    case ("nonsu2")
       allocate(getsector(0:Nlevels,1))
    end select
    allocate(getCsector(2,Nsectors))
    allocate(getCDGsector(2,Nsectors))
    allocate(getBathStride(Norb,Nbath))
    allocate(neigen_sector(Nsectors))
    !
    !
    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.          !set to do zero temperature calculations
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !
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
    !
    !
    if(ED_MPI_ID==0)then
       if(finiteT)then
          write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
          write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
          write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
          call sleep(1)
       else
          write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
          call sleep(1)
       endif
    endif
    !
    !
    !#############################################################################################
    !CHECKS:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "Nspin > 2 ERROR"
    if(Norb>3)stop "Norb > 3 ERROR" 
    if(nerr < dmft_error) nerr=dmft_error
    if(ed_mode=="superc")then
       if(Nspin>1)stop "SC+AFM ERROR." 
       if(Norb>1)stop "SC Multi-Band IS NOT TESTED. remove this line in ED_AUX_FUNX to proceed."
       if(ed_type=='c')stop "SC with Hermitian H NOT IMPLEMENTED."
       if(ed_twin)stop  "SC + ED_TWIN NOT TESTED. remove this line in ED_AUX_FUNX to proceed."
    endif
    if(ed_mode=="nonsu2")then
       if(Nspin/=2)then
          write(LOGfile,"(A)")"ED msg: ed_mode=nonSU2 with Nspin!=1 is not allowed."
          write(LOGfile,"(A)")"        to enforce spin symmetry up-dw set ed_para=T."
          stop
       endif
       ! if(ed_twin)stop  "NONSU2 + ED_TWIN NOT TESTED. remove this line in ED_AUX_FUNX to proceed."
    endif
    !#############################################################################################
    !
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
       !nloop=(i-1)*niter                !increase the max number of dmft loop allowed so to do threshold loop
       !write(LOGfile,"(A,I10)")"Increased Nloop to:",nloop
    endif
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impFmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impFreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_phisc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))
  end subroutine init_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  !
  !NORMAL CASE
  !
  subroutine setup_pointers_normal
    integer                          :: i,in,dim,isector,jsector
    integer                          :: nup,ndw,jup,jdw,iorb
    integer                          :: unit,status,istate,counter
    logical                          :: IOfile
    integer                          :: anint
    real(8)                          :: adouble
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getsector(nup,ndw)=isector
          getnup(isector)=nup
          getndw(isector)=ndw
          dim = get_normal_sector_dimension(nup,ndw)
          getdim(isector)=dim
       enddo
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          read(unit,"(i3,f18.12,2x,ES19.12,1x,2i3,3x,i3,i10)",iostat=status)istate,adouble,adouble,nup,ndw,isector,anint
          list_sector(istate)=isector
          if(nup/=getnup(isector).OR.ndw/=getndw(isector))&
               stop "setup_pointers_normal error: nup!=getnup(isector).OR.ndw!=getndw(isector) "
       enddo
       close(unit)
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsectors
          nup=getnup(isector)
          ndw=getndw(isector)
          if(nup<ndw)twin_mask(isector)=.false.
       enddo
       if(ED_MPI_ID==0)write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    if(ED_MPI_ID==0)call stop_timer
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
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
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_normal

  ! 
  !SUPERCONDUCTING
  !
  subroutine setup_pointers_superc
    integer                          :: i,isz,in,dim,isector,jsector
    integer                          :: sz,iorb,jsz
    integer                          :: unit,status,istate
    logical                          :: IOfile
    integer                          :: anint
    real(8)                          :: adouble
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do isz=-Ns,Ns
       sz=abs(isz)
       isector=isector+1
       getsector(isz,1)=isector
       getsz(isector)=isz
       dim = get_superc_sector_dimension(isz)
       getdim(isector)=dim
       ! neigen_sector(isector) = min(dim,lanc_nstates_sector)   !init every sector to required eigenstates
    enddo
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          read(unit,"(i3,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)",iostat=status) istate,adouble,adouble,sz,isector,anint
          list_sector(istate)=isector
          if(sz/=getsz(isector))stop "setup_pointers_superc error: sz!=getsz(isector)."
       enddo
       close(unit)
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    twin_mask=.true.
    if(ed_twin)then
       write(LOGfile,*)"TWIN STATES IN SC CHANNEL: NOT TESTED!!"
       call sleep(3)
       do isector=1,Nsectors
          sz=getsz(isector)
          if(sz>0)twin_mask(isector)=.false.
       enddo
       if(ED_MPI_ID==0)write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    if(ED_MPI_ID==0)call stop_timer
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
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
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !    
    getCsector=0
    !c_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCsector(1,isector)=jsector
    enddo
    !c_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    !cdg_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCDGsector(1,isector)=jsector
    enddo
    !cdg_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_superc

  !
  !NON SU(2) SYMMETRIC
  !
  subroutine setup_pointers_nonsu2
    integer                          :: i,dim,isector,jsector
    integer                          :: in,jn,iorb
    integer                          :: unit,status,istate
    logical                          :: IOfile
    integer                          :: anint
    real(8)                          :: adouble
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do in=0,Nlevels
       isector=isector+1
       getsector(in,1)=isector
       getn(isector)=in
       dim = get_nonsu2_sector_dimension(in)
       getdim(isector)=dim
       neigen_sector(isector) = min(dim,lanc_nstates_sector)
    enddo
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          read(unit,"(i3,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)",iostat=status) istate,adouble,adouble,in,isector,anint
          list_sector(istate)=isector
          if(in/=getn(isector))stop "setup_pointers_superc error: n!=getn(isector)."
       enddo
       close(unit)
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    twin_mask=.true.
    if(ed_twin)then
       write(LOGfile,*)"TWIN STATES IN SC CHANNEL: NOT TESTED!!"
       call sleep(3)
       do isector=1,Nsectors
          in=getn(isector)
          if(in>Ns)twin_mask(isector)=.false.
          print*,twin_mask(isector),in
       enddo
       if(ED_MPI_ID==0)write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    if(ED_MPI_ID==0)call stop_timer
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
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
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector=0
    !c_{up,dw}
    do isector=1,Nsectors
       in=getn(isector);if(in==0)cycle
       jn=in-1
       jsector=getsector(jn,1)
       getCsector(1,isector)=jsector
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    !cdg_{up,dw}
    do isector=1,Nsectors
       in=getn(isector);if(in==Nlevels)cycle
       jn=in+1
       jsector=getsector(jn,1)
       getCDGsector(1,isector)=jsector
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_nonsu2





  !+------------------------------------------------------------------+
  !PURPOSE  : return the dimension of a sector
  !+------------------------------------------------------------------+
  !NORMAL
  function get_normal_sector_dimension(nup,ndw) result(dim)
    integer          :: nup
    integer,optional :: ndw
    integer          :: dim,dimup,dimdw
    if(present(ndw))then
       dimup = binomial(Ns,nup)    !this ensures better evaluation of the dimension
       dimdw = binomial(Ns,ndw)    !as it avoids large numbers
    else
       dimup = binomial(Ns,nup)
       dimdw = 1
    endif
    dim=dimup*dimdw
  end function get_normal_sector_dimension
  !SUPERC
  function get_superc_sector_dimension(mz) result(dim)
    integer :: mz
    integer :: i,dim,Nb
    dim=0
    Nb=Ns-mz
    do i=0,Nb/2 
       dim=dim + 2**(Nb-2*i)*binomial(ns,Nb-2*i)*binomial(ns-Nb+2*i,i)
    enddo
  end function get_superc_sector_dimension
  !NONSU2
  function get_nonsu2_sector_dimension(n) result(dim)
    integer :: n
    integer :: dim
    dim=binomial(2*Ns,n)
  end function get_nonsu2_sector_dimension










  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,Hup,Hdw)
    integer                   :: isector
    type(sector_map)          :: Hup
    type(sector_map),optional :: Hdw
    integer                   :: nup,ndw,sz,nt
    integer                   :: nup_,ndw_,sz_,nt_
    integer                   :: i
    integer                   :: iup,idw
    integer                   :: dim
    integer                   :: ivec(Ns),jvec(Ns)
    select case(ed_mode)
    case default
       nup = getNup(isector)
       ndw = getNdw(isector)
       if(present(Hdw))then
          !
          call map_allocate(Hup,get_normal_sector_dimension(nup))
          dim=0
          do i=0,2**Ns-1
             ivec  = bdecomp(i,Ns)
             nup_  = sum(ivec) 
             if(nup_ /= nup)cycle
             dim       = dim+1
             Hup%map(dim) = i
          enddo
          call map_allocate(Hdw,get_normal_sector_dimension(ndw))
          dim=0
          do i=0,2**Ns-1
             ivec  = bdecomp(i,Ns)
             ndw_  = sum(ivec) 
             if(ndw_ /= ndw)cycle
             dim       = dim+1
             Hdw%map(dim) = i
          enddo
          !
       else
          !
          call map_allocate(Hup,get_normal_sector_dimension(nup,ndw))
          dim=0
          do idw=0,2**Ns-1
             jvec  = bdecomp(idw,Ns)
             ndw_  = sum(jvec)
             if(ndw_ /= ndw)cycle
             do iup=0,2**Ns-1
                ivec  = bdecomp(iup,Ns)
                nup_  = sum(ivec)
                if(nup_ /= nup)cycle
                dim      = dim+1
                Hup%map(dim) = iup + idw*2**Ns
             enddo
          enddo
          !
       endif
       !
       !
    case ("superc")
       sz = getSz(isector)
       call map_allocate(Hup,get_superc_sector_dimension(sz))
       dim=0
       do idw=0,2**Ns-1
          jvec = bdecomp(idw,Ns)
          ndw_ = sum(jvec)
          do iup=0,2**Ns-1
             ivec = bdecomp(iup,Ns)
             nup_ = sum(ivec)
             sz_  = nup_ - ndw_
             if(sz_ == sz)then
                dim=dim+1
                Hup%map(dim)=iup + idw*2**Ns
             endif
          enddo
       enddo
       !
       !
    case ("nonsu2")
       nt = getN(isector)
       call map_allocate(Hup,get_nonsu2_sector_dimension(nt))
       dim=0
       do idw=0,2**Ns-1
          jvec = bdecomp(idw,Ns)
          do iup=0,2**Ns-1
             ivec = bdecomp(iup,Ns)
             nt_  = sum(ivec) + sum(jvec)
             if(nt_ == nt)then
                dim=dim+1
                Hup%map(dim)=iup + idw*2**Ns
             endif
          enddo
       enddo
    end select
  end subroutine build_sector






  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c
  ! subroutine c(m,i,j,sgn)
  !   integer :: ib(Nlevels)
  !   integer :: i,j,m,km,k
  !   integer :: isg
  !   real(8) :: sgn
  !   call bdecomp(i,ib)
  !   if (ib(m)==0)then
  !      j=0
  !   else
  !      if(m==1)then
  !         j=i-1
  !      else
  !         km=0
  !         do k=1,m-1
  !            km=km+ib(k)
  !         enddo
  !         !km=sum(ib(1:m-1))
  !         isg=(-1)**km
  !         j=(i-2**(m-1))*isg
  !      endif
  !   endif
  !   sgn=dfloat(j)/dfloat(abs(j))
  !   j=abs(j)
  ! end subroutine c



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg
  ! subroutine cdg(m,i,j,sgn)
  !   integer :: ib(Nlevels)
  !   integer :: i,j,m,km,k
  !   integer :: isg
  !   real(8) :: sgn
  !   call bdecomp(i,ib)
  !   if (ib(m)==1)then
  !      j=0
  !   else
  !      if(m==1)then
  !         j=i+1
  !      else
  !         km=0
  !         do k=1,m-1
  !            km=km+ib(k)
  !         enddo
  !         !km=sum(ib(1:m-1))
  !         isg=(-1)**km
  !         j=(i+2**(m-1))*isg
  !      endif
  !   endif
  !   sgn=dfloat(j)/dfloat(abs(j))
  !   j=abs(j)
  ! end subroutine cdg






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw. 
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                          :: isector
    integer,dimension(:)             :: order
    ! integer,dimension(:),allocatable :: Hmap
    type(sector_map)                 :: H
    integer                          :: i,dim
    dim = getdim(isector)
    if(size(Order)/=dim)stop "ED_AUX_FUNX/build_twin_sector: wrong dimensions."
    ! allocate(Hmap(dim))
    call map_allocate(H,dim)
    call build_sector(isector,H) !build the map from the A-sector to \HHH
    do i=1,dim                      !
       Order(i)=flip_state(H%map(i)) !get the list of states in \HHH corresponding to sector B twin of A
    enddo                           !
    call sort_array(Order)          !return the ordering of B-states in \HHH with respect to those of A
    deallocate(H%map)
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  ! superc: j=|{dw}>|{up}>  , sz  --> -sz
  ! nonsu2: j=|{!up}>|{!dw}>, n   --> 2*Ns-n
  !+------------------------------------------------------------------+
  function flip_state(m,n) result(j)
    integer          :: m
    integer,optional :: n
    integer          :: j
    integer          :: ivec(2*Ns),foo(2*Ns)

    select case(ed_mode)
    case default
       !Exchange UP-config |{up}> with DW-config |{dw}>
       foo(1:Ns)     =Ivec(Ns+1:2*Ns)
       foo(Ns+1:2*Ns)=Ivec(1:Ns)
    case("superc")
       !Invert the overall spin sign: |{up}> <---> |{dw}>
       Ivec = bdecomp(m,2*Ns)
       foo(1:Ns)     =Ivec(Ns+1:2*Ns)
       foo(Ns+1:2*Ns)=Ivec(1:Ns)
    case ("nonsu2")
       !Exchange Occupied sites (1) with Empty sites (0)
       Ivec = bdecomp(m,2*Ns)
       where(Ivec==1)foo=0
       where(Ivec==0)foo=1
    end select
    j = bjoin(foo,2*Ns)
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  ! sz      ==> -sz     (total spin flip)
  ! n       ==> 2*Ns-n  (particle hole)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in) :: isector
    integer :: jsector
    integer :: iup,idw,in,isz
    select case(ed_mode)
    case default
       iup=getnup(isector)
       idw=getndw(isector)
       jsector=getsector(idw,iup)
    case ("superc")
       isz=getsz(isector)
       jsector=getsector(-isz,1)
    case("nonsu2")
       in=getn(isector)
       jsector=getsector(Nlevels-in,1)
    end select
  end function get_twin_sector










  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp
  ! subroutine bdecomp(i,ivec)
  !   integer :: ivec(Nlevels)         
  !   integer :: l,i
  !   logical :: busy
  !   !this is the configuration vector |1,..,Ns,Ns+1,...,Nlevels>
  !   !obtained from binary decomposition of the state/number i\in 2^2*Ns
  !   do l=0,Nlevels-1
  !      busy=btest(i-1,l)
  !      ivec(l+1)=0
  !      if(busy)ivec(l+1)=1
  !   enddo
  ! end subroutine bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin
  ! subroutine bjoin(ivec,i)
  !   integer,dimension(Nlevels) :: ivec
  !   integer                 :: i,j
  !   i=1
  !   do j=1,Nlevels
  !      i=i+ivec(j)*2**(j-1)
  !   enddo
  ! end subroutine bjoin

  function Imap_ud(iup,idw,Ns) result(i)
    integer :: iup,idw,Ns
    integer :: i
    i   = iup + idw*2**Ns
  end function Imap_ud



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
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
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
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



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
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
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




  !+------------------------------------------------------------------+
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ivec,unit)
    integer,intent(in) :: ivec(:)
    integer,optional   :: unit
    integer            :: unit_
    integer            :: dim,i,j,Ntot
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    i= bjoin(ivec,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_ivec
  !
  subroutine  print_state_vector_ivec_ud(ivec,jvec,unit)
    integer,intent(in) :: ivec(:),jvec(size(ivec))
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,iup,idw,Ntot
    character(len=2)   :: fbt
    character(len=20)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//",1x,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    iup = bjoin(ivec,Ntot)
    idw = bjoin(jvec,Ntot)
    i = bjoin([ivec,jvec],2*Ntot)
    write(unit_,"(I9,1x,I4,1x,A1)",advance="no")i,iup,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A1,I4,A2)",advance="no")">",idw," |"
    write(unit_,"(10I1)",advance="no")(jvec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")ibits(i,0,Ntot),ibits(i,Ntot,2*Ntot)
  end subroutine print_state_vector_ivec_ud
  !
  subroutine print_state_vector_int(i,Ntot,unit)
    integer,intent(in) :: i
    integer,intent(in) :: Ntot
    integer,optional   :: unit
    integer            :: unit_
    integer            :: dim,j
    integer            :: ivec(Ntot)
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    ivec = bdecomp(i,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_int

  ! subroutine print_state_vector_ivec(ib,unit)
  !   integer,optional :: unit
  !   integer :: i,j,unit_
  !   integer :: ib(Nlevels)
  !   unit_=6;if(present(unit))unit_=unit
  !   call bjoin(ib,i)
  !   write(unit_,"(I3,1x,A1)",advance="no")i,"|"
  !   write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
  !   write(unit_,"(A1,A1)",advance="no")">","|"
  !   write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
  !   write(unit_,"(A1)",advance="yes")">"
  ! end subroutine print_state_vector_ivec
  !
  ! subroutine print_state_vector_int(i,unit)
  !   integer,optional :: unit
  !   integer :: i,j,unit_
  !   integer :: ib(Nlevels)
  !   unit_=6;if(present(unit))unit_=unit
  !   call bdecomp(i,ib)
  !   write(unit_,"(I3,1x,A1)",advance="no")i,"|"
  !   write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
  !   write(unit_,"(A2)",advance="no")">|"
  !   write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
  !   write(unit_,"(A1)",advance="yes")">"
  ! end subroutine print_state_vector_int




end MODULE ED_SETUP
