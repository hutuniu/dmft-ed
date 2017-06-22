!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit
  USE SF_STAT
  USE SF_SP_LINALG
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN_MATVEC
  !
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private


  public :: diagonalize_impurity

  public :: ed_diag_set_MPI
  public :: ed_diag_del_MPI


#ifdef _MPI
  integer :: MpiComm=MPI_UNDEFINED
#else
  integer :: MpiComm=0
#endif
  logical :: MpiStatus=.false.
  integer :: MPI_SIZE=1
  logical :: MPI_MASTER=.true.  !
  integer :: unit


contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the MPI-Parallel environment for ED_DIAG
  !+------------------------------------------------------------------+
  subroutine ed_diag_set_MPI(comm)
#ifdef _MPI
    integer :: comm
    MpiComm  = comm
    MpiStatus = .true.
    MPI_SIZE  = get_Size_MPI(MpiComm)
    MPI_MASTER= get_Master_MPI(MpiComm)
#else
    integer,optional :: comm
#endif
  end subroutine ed_diag_set_MPI

  subroutine ed_diag_del_MPI()
#ifdef _MPI
    MpiComm  = MPI_UNDEFINED
    MpiStatus = .false.
#endif
  end subroutine ed_diag_del_MPI





  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine diagonalize_impurity()
    call ed_diag_c
    call ed_analysis
  end subroutine diagonalize_impurity







  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine ed_diag_c
    integer                :: nup,ndw,isector,dim
    integer                :: isect,izero,sz,nt
    integer                :: i,j,iter,unit
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Ei,Jz
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:)
    logical                :: lanc_solve,Tflag,lanc_verbose
    !
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    ! if(MPI_MASTER)then
    write(LOGfile,"(A)")"Diagonalize impurity H:"
    call start_timer()
    ! endif
    !
    lanc_verbose=.false.
    if(ed_verbose>2)lanc_verbose=.true.
    !
    iter=0
    sector: do isector=1,Nsectors
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       !DEBUG>>
       if(Jz_basis.and.Jz_max.and.abs(gettwoJz(isector))>int(2.*Jz_max_value))cycle
       !>>DEBUG
       iter=iter+1
       Tflag    = twin_mask(isector).AND.ed_twin
       select case(ed_mode)
       case default
          Tflag = Tflag.AND.(getnup(isector)/=getndw(isector))
       case ("superc")
          Tflag = Tflag.AND.(getsz(isector)/=0)
       case("nonsu2")
          Tflag = Tflag.AND.(getn(isector)/=Ns)
       end select
       Dim      = getdim(isector)
       Neigen   = min(dim,neigen_sector(isector))
       Nitermax = min(dim,lanc_niter)
       Nblock   = min(dim,lanc_ncv_factor*Neigen + lanc_ncv_add)!min(dim,5*Neigen+10)
       !
       lanc_solve  = .true.
       if(Neigen==dim)lanc_solve=.false.
       if(dim<=max(lanc_dim_threshold,MPI_SIZE))lanc_solve=.false.
       !
       ! if(MPI_MASTER)then
       if(ed_verbose==3)then
          select case(ed_mode)
          case default
             nup  = getnup(isector)
             ndw  = getndw(isector)
             write(LOGfile,"(1X,I4,A,I4,A6,I2,A6,I2,A6,I15,A12,3I6)")&
                  iter,"-Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock
          case ("superc")
             sz   = getsz(isector)
             write(LOGfile,"(1X,I4,A,I4,A5,I4,A6,I15,A12,3I6)")&
                  iter,"-Solving sector:",isector," sz:",sz," dim=",getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock
          case ("nonsu2")
             if(Jz_basis)then
                nt   = getn(isector)
                Jz   = gettwoJz(isector)/2.
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,F5.1,A6,I15,A12,3I6)")&
                     iter,"-Solving sector:",isector," n:",nt," Jz:",Jz," dim=",getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock

             else
                nt   = getn(isector)
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,I15,A12,3I6)")&
                     iter,"-Solving sector:",isector," n:",nt," dim=",getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock
             endif
          end select
       elseif(ed_verbose==1.OR.ed_verbose==2)then
          call eta(iter,count(twin_mask),LOGfile)
       endif
       ! endif
       !
       if(lanc_solve)then
          if(allocated(eig_values))deallocate(eig_values)
          if(allocated(eig_basis))deallocate(eig_basis)
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0d0 ; eig_basis=zero
          !
          call setup_Hv_sector(isector)
          if(ed_sparse_H)call ed_buildH_c()
          !
          if(MpiStatus)then
             call sp_eigh(MpiComm,spHtimesV_cc,Dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,tol=lanc_tolerance)
          else
             call sp_eigh(spHtimesV_cc,Dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,tol=lanc_tolerance)
          endif
          call delete_Hv_sector()
       else
          if(allocated(eig_values))deallocate(eig_values)
          if(allocated(eig_basis))deallocate(eig_basis)
          allocate(eig_values(Dim),eig_basis(Dim,dim))
          eig_values=0d0 ; eig_basis=zero
          call setup_Hv_sector(isector)
          call ed_buildH_c(eig_basis)
          call delete_Hv_sector()
          call eigh(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=one
       endif
       !
       if(spH0%status)call sp_delete_matrix(spH0)
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,twin=Tflag,size=lanc_nstates_total)
          enddo
       else
          do i=1,Neigen
             enemin = eig_values(i)
             if (enemin < oldzero-10.d0*gs_threshold)then
                oldzero=enemin
                call es_free_espace(state_list)
                call es_add_state(state_list,enemin,eig_basis(1:dim,i),isector,twin=Tflag)
             elseif(abs(enemin-oldzero) <= gs_threshold)then
                oldzero=min(oldzero,enemin)
                call es_add_state(state_list,enemin,eig_basis(1:dim,i),isector,twin=Tflag)
             endif
          enddo
       endif
       if(MPI_MASTER)then
          unit=free_unit()
          open(unit,file="eigenvalues_list"//reg(ed_file_suffix)//".ed",position='append',action='write')
          call print_eigenvalues_list(isector,eig_values(1:Neigen),unit)
          close(unit)
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    !if(MPI_MASTER)
    call stop_timer(LOGfile)
  end subroutine ed_diag_c










  !###################################################################################################
  !
  !    POST-PROCESSING ROUTINES
  !
  !###################################################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : analyse the spectrum and print some information after 
  !lanczos diagonalization. 
  !+------------------------------------------------------------------+
  subroutine ed_analysis()
    integer             :: nup,ndw,sz,n,isector,dim
    integer             :: istate
    integer             :: i,unit
    integer             :: Nsize,NtoBremoved,nstates_below_cutoff
    integer             :: numgs
    real(8)             :: Egs,Ei,Ec,Etmp
    type(histogram)     :: hist
    real(8)             :: hist_a,hist_b,hist_w
    integer             :: hist_n
    integer,allocatable :: list_sector(:),count_sector(:)    
    !POST PROCESSING:
    if(MPI_MASTER)then
       unit=free_unit()
       open(unit,file="state_list"//reg(ed_file_suffix)//".ed")
       call print_state_list(unit)
       close(unit)
    endif
    if(ed_verbose>=2)call print_state_list(LOGfile)
    !
    zeta_function=0d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei   = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       zeta_function=real(state_list%size,8)
    end if
    !
    !
    numgs=es_return_gs_degeneracy(state_list,gs_threshold)
    if(numgs>Nsectors)stop "ed_diag: too many gs"
    ! if(MPI_MASTER.AND.ed_verbose>=2)then
    if(ed_verbose>=2)then
       do istate=1,numgs
          isector = es_return_sector(state_list,istate)
          Egs     = es_return_energy(state_list,istate)
          select case(ed_mode)
          case default
             nup  = getnup(isector)
             ndw  = getndw(isector)
             write(LOGfile,"(A,F20.12,2I4)")'Egs =',Egs,nup,ndw
          case("superc")
             sz  = getsz(isector)
             write(LOGfile,"(A,F20.12,I4)")'Egs =',Egs,sz
          case("nonsu2")
             n  = getn(isector)
             write(LOGfile,"(A,F20.12,I4)")'Egs =',Egs,n
          end select
       enddo
       write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    endif
    !
    !
    !
    !get histogram distribution of the sector contributing to the evaluated spectrum:
    !go through states list and update the neigen_sector(isector) sector-by-sector
    if(finiteT)then
       if(MPI_MASTER)then
          unit=free_unit()
          open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",position='append')
          hist_n = Nsectors
          hist_a = 1d0
          hist_b = dble(Nsectors)
          hist_w = 1d0
          hist = histogram_allocate(hist_n)
          call histogram_set_range_uniform(hist,hist_a,hist_b)
          do i=1,state_list%size
             isector = es_return_sector(state_list,i)
             call histogram_accumulate(hist,dble(isector),hist_w)
          enddo
          call histogram_print(hist,unit)
          write(unit,*)""
          close(unit)
       endif
       !
       !
       !
       allocate(list_sector(state_list%size),count_sector(Nsectors))
       !get the list of actual sectors contributing to the list
       do i=1,state_list%size
          list_sector(i) = es_return_sector(state_list,i)
       enddo
       !count how many times a sector appears in the list
       do i=1,Nsectors
          count_sector(i) = count(list_sector==i)
       enddo
       !adapt the number of required Neig for each sector based on how many
       !appeared in the list.
       do i=1,Nsectors
          if(any(list_sector==i))then !if(count_sector(i)>1)then
             neigen_sector(i)=neigen_sector(i)+1
          else
             neigen_sector(i)=neigen_sector(i)-1
          endif
          !prevent Neig(i) from growing unbounded but 
          !try to put another state in the list from sector i
          if(neigen_sector(i) > count_sector(i))neigen_sector(i)=count_sector(i)+1
          if(neigen_sector(i) <= 0)neigen_sector(i)=1
       enddo
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-beta(Ec-Egs)) < \epsilon_c
       ! if this condition is violated then required number of states is increased
       ! if number of states is larger than those required to fullfill the cutoff: 
       ! trim the list and number of states.
       Egs  = state_list%emin
       Ec   = state_list%emax
       Nsize= state_list%size
       if(exp(-beta*(Ec-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + lanc_nstates_step
          !if(MPI_MASTER)
          write(LOGfile,"(A,I4)")"Increasing lanc_nstates_total:",lanc_nstates_total
       else
          ! !Find the energy level beyond which cutoff condition is verified & cut the list to that size
          write(LOGfile,*)
          isector = es_return_sector(state_list,state_list%size)
          Ei      = es_return_energy(state_list,state_list%size)
          do while ( exp(-beta*(Ei-Egs)) <= cutoff )
             ! if(ed_verbose>=1.AND.MPI_MASTER)then
             if(ed_verbose>=1)then
                select case(ed_mode)
                case default
                   write(LOGfile,"(A,I4,2x,2I4,2x,I5)")&
                        "Trimming state:",isector,getnup(isector),getndw(isector),state_list%size
                case("superc")
                   write(LOGfile,"(A,I4,2x,I4,2x,I5)")&
                        "Trimming state:",isector,getsz(isector),state_list%size
                case("nonsu2")
                   write(LOGfile,"(A,I4,2x,I4,2x,I5)")&
                        "Trimming state:",isector,getn(isector),state_list%size
                end select
             endif
             call es_pop_state(state_list)
             isector = es_return_sector(state_list,state_list%size)
             Ei      = es_return_energy(state_list,state_list%size)
          enddo
          ! if(ed_verbose>=1.AND.MPI_MASTER)then
          if(ed_verbose>=1)then
             write(LOGfile,*)"Trimmed state list:"          
             call print_state_list(LOGfile)
          endif
          !
          lanc_nstates_total=max(state_list%size,lanc_nstates_step)+lanc_nstates_step
          ! if(ed_verbose>=1.AND.MPI_MASTER)write(*,"(A,I4)")"Adjusting lanc_nstates_total to:",lanc_nstates_total
          write(*,"(A,I4)")"Adjusting lanc_nstates_total to:",lanc_nstates_total
          !
       endif
    endif
  end subroutine ed_analysis


  subroutine print_state_list(unit)
    integer :: nup,ndw,sz,n,isector
    integer :: istate
    integer :: unit
    real(8) :: Estate,Jz
    ! if(MPI_MASTER)then
    select case(ed_mode)
    case default
       write(unit,"(A)")"# i       E_i           exp(-(E-E0)/T)       nup ndw  Sect     Dim"
    case ("superc")
       write(unit,"(A)")"# i       E_i           exp(-(E-E0)/T)       Sz     Sect     Dim"
    case ("nonsu2")
       if(Jz_basis)then
          write(unit,"(A3,A18,2x,A19,1x,A3,3x,A4,3x,A3,A10)")"# i","E_i","exp(-(E-E0)/T)","n","Jz","Sect","Dim"
       else
          write(unit,"(A3,A18,2x,A19,1x,A3,3x,A3,A10)")"# i","E_i","exp(-(E-E0)/T)","n","Sect","Dim"
       endif
    end select
    do istate=1,state_list%size
       Estate  = es_return_energy(state_list,istate)
       isector = es_return_sector(state_list,istate)
       select case(ed_mode)
       case default
          nup   = getnup(isector)
          ndw   = getndw(isector)
          write(unit,"(i3,f18.12,2x,ES19.12,1x,2i3,3x,i3,i10)")&
               istate,Estate,exp(-beta*(Estate-state_list%emin)),nup,ndw,isector,getdim(isector)
       case("superc")
          sz   = getsz(isector)
          write(unit,"(i3,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)")&
               istate,Estate,exp(-beta*(Estate-state_list%emin)),sz,isector,getdim(isector)
       case("nonsu2")
          n    = getn(isector)
          if(Jz_basis)then
             Jz   = gettwoJz(isector)/2.
             write(unit,"(i3,f18.12,2x,ES19.12,1x,i3,3x,F4.1,3x,i3,i10)")&
                  istate,Estate,exp(-beta*(Estate-state_list%emin)),n,Jz,isector,getdim(isector)
          else
             write(unit,"(i3,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)")&
                  istate,Estate,exp(-beta*(Estate-state_list%emin)),n,isector,getdim(isector)
          endif
       end select
    enddo
    ! endif
  end subroutine print_state_list

  subroutine print_eigenvalues_list(isector,eig_values,unit)
    integer              :: isector
    real(8),dimension(:) :: eig_values
    integer              :: unit,i
    ! if(MPI_MASTER)then
    select case(ed_mode)
    case default
       write(unit,"(A7,A3,A3)")" # Sector","Nup","Ndw"
       write(unit,"(I4,2x,I3,I3)")isector,getnup(isector),getndw(isector)
    case ("superc")
       write(unit,"(A7,A4)")" # Sector","Sz"
       write(unit,"(I4,2x,I4)")isector,getsz(isector)
    case ("nonsu2")
       write(unit,"(A7,A3)")" # Sector","N"
       write(unit,"(I4,2x,I4)")isector,getn(isector)
    end select
    do i=1,size(eig_values)
       write(unit,*)eig_values(i)
    enddo
    write(unit,*)""
    ! endif
  end subroutine print_eigenvalues_list





end MODULE ED_DIAG









