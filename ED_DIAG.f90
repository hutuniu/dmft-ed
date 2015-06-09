!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE SF_CONSTANTS
  USE SF_LINALG, only: matrix_diagonalize
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit
  USE SF_STAT
  !
  USE ARPACK_LANCZOS
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  USE ED_MATVEC
  implicit none
  private

  public :: diagonalize_impurity

  real(8),allocatable,dimension(:) :: egs_values


contains

  !                    LANCZOS DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine diagonalize_impurity
    select case(ed_type)
    case default
       call ed_diag_d
    case('c')
       call ed_diag_c
    end select
    call ed_analysis
  end subroutine diagonalize_impurity



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_diag_d
    integer             :: nup,ndw,isector,dim
    integer             :: sz,nt
    integer             :: i,iter
    integer             :: numgs
    integer             :: Nitermax,Neigen,Nblock
    real(8)             :: oldzero,enemin
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    logical             :: lanc_solve,Tflag
    if(allocated(egs_values))deallocate(egs_values)
    allocate(egs_values(Nsectors))
    egs_values=0d0
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    numgs=0
    if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer()
    if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Diagonalize impurity H:"
    iter=0
    sector: do isector=1,Nsectors
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       iter=iter+1
       if(ED_MPI_ID==0)then
          if(ed_verbose<0)then
             dim      = getdim(isector)
             select case(ed_mode)
             case default
                nup  = getnup(isector)
                ndw  = getndw(isector)
                write(LOGfile,"(1X,I4,A,I4,A6,I2,A6,I2,A6,I15)")iter,"-Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
             case ("superc")
                sz   = getsz(isector)
                write(LOGfile,"(1X,I4,A,I4,A5,I4,A6,I15)")iter,"-Solving sector:",isector," sz:",sz," dim=",getdim(isector)
             case ("nonsu2")
                nt   = getn(isector)
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,I15)")iter,"-Solving sector:",isector," n:",nt," dim=",getdim(isector)
             end select
          elseif(ed_verbose<2)then
             call eta(iter,count(twin_mask),LOGfile)
          endif
       endif
       ! Tflag    = twin_mask(isector).AND.ed_twin.AND.(getnup(isector)/=getndw(isector))
       Tflag    = twin_mask(isector).AND.ed_twin
       select case(ed_mode)
       case default
          Tflag=Tflag.AND.(getnup(isector)/=getndw(isector))
       case ("superc")
          Tflag=Tflag.AND.(getsz(isector)/=0)
       case("nonsu2")
          Tflag=Tflag.AND.(getn(isector)/=Ns)
       end select
       dim      = getdim(isector)
       Neigen   = min(dim,neigen_sector(isector))
       Nitermax = min(dim,lanc_niter)
       Nblock   = min(dim,5*Neigen+10)
       !
       lanc_solve  = .true.
       if(Neigen==dim)lanc_solve=.false.
       if(dim<=max(512,ED_MPI_SIZE))lanc_solve=.false.
       !
       if(lanc_solve)then
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0.d0 ; eig_basis=0.d0
          call ed_buildH_d(isector)
#ifdef _MPI
          call lanczos_parpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_dd)
#else
          call lanczos_arpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_dd)
#endif
          if(spH0%status)call sp_delete_matrix(spH0)
       else
          allocate(eig_values(dim),eig_basis(dim,dim))
          eig_values=0.d0 ; eig_basis=0.d0 
          call ed_buildH_d(isector,eig_basis)
          if(spH0%status)call sp_delete_matrix(spH0)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=1.d0
       endif
       !
       egs_values(isector)=eig_values(1)
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,twin=Tflag,size=lanc_nstates_total)
          enddo
       else
          enemin = eig_values(1)
          if (enemin < oldzero-10.d0*gs_threshold)then
             numgs=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          elseif(abs(enemin-oldzero) <= gs_threshold)then
             numgs=numgs+1
             if (numgs > Nsectors)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          endif
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  end subroutine ed_diag_d



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine ed_diag_c
    integer                :: nup,ndw,isector,dim
    integer                :: sz,nt
    integer                :: i,iter
    integer                :: numgs
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:)
    logical                :: lanc_solve,Tflag
    if(allocated(egs_values))deallocate(egs_values)
    allocate(egs_values(Nsectors))
    egs_values=0.d0
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    numgs=0
    if(ed_verbose<3.AND.ED_MPI_ID==0)call start_timer()
    if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A)")"Diagonalize impurity H:"
    iter=0
    sector: do isector=1,Nsectors
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       iter=iter+1
       if(ED_MPI_ID==0)then
          if(ed_verbose<0)then
             dim      = getdim(isector)
             select case(ed_mode)
             case default
                nup  = getnup(isector)
                ndw  = getndw(isector)
                write(LOGfile,"(1X,I4,A,I4,A6,I2,A6,I2,A6,I15)")iter,"-Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
             case ("superc")
                sz   = getsz(isector)
                write(LOGfile,"(1X,I4,A,I4,A5,I4,A6,I15)")iter,"-Solving sector:",isector," sz:",sz," dim=",getdim(isector)
             case ("nonsu2")
                nt   = getn(isector)
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,I15)")iter,"-Solving sector:",isector," n:",nt," dim=",getdim(isector)
             end select
          elseif(ed_verbose<2)then
             call eta(iter,count(twin_mask),LOGfile)
          endif
       endif
       !Tflag    = twin_mask(isector).AND.ed_twin.AND.(getnup(isector)/=getndw(isector))
       Tflag    = twin_mask(isector).AND.ed_twin
       select case(ed_mode)
       case default
          Tflag=Tflag.AND.(getnup(isector)/=getndw(isector))
       case ("superc")
          Tflag=Tflag.AND.(getsz(isector)/=0)
       case("nonsu2")
          Tflag=Tflag.AND.(getn(isector)/=Ns)
       end select
       dim      = getdim(isector)
       Neigen   = min(dim,neigen_sector(isector))
       Nitermax = min(dim,lanc_niter)
       Nblock   = min(dim,5*Neigen+10)
       !
       lanc_solve  = .true.
       if(Neigen==dim)lanc_solve=.false.
       if(dim<=max(512,ED_MPI_SIZE))lanc_solve=.false.
       !
       if(lanc_solve)then
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0.d0 ; eig_basis=zero
          call ed_buildH_c(isector)
#ifdef _MPI
          call lanczos_parpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_cc)
#else
          call lanczos_arpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_cc)
#endif
          if(spH0%status)call sp_delete_matrix(spH0)
       else
          allocate(eig_values(Dim),eig_basis(Dim,dim))
          eig_values=0.d0 ; eig_basis=zero
          call ed_buildH_c(isector,eig_basis)
          if(spH0%status)call sp_delete_matrix(spH0)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=one
       endif
       !
       egs_values(isector)=eig_values(1)
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,twin=Tflag,size=lanc_nstates_total)
          enddo
       else
          enemin = eig_values(1)
          if (enemin < oldzero-10.d0*gs_threshold)then
             numgs=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          elseif(abs(enemin-oldzero) <= gs_threshold)then
             numgs=numgs+1
             if (numgs > Nsectors)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          endif
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    if(ed_verbose<3.AND.ED_MPI_ID==0)call stop_timer
  end subroutine ed_diag_c





  subroutine print_state_list(unit)
    integer :: nup,ndw,sz,n,isector
    integer :: istate
    integer :: unit
    real(8) :: Estate
    if(ed_twin)then
       do isector=1,Nsectors
          if(.not.twin_mask(isector))egs_values(get_twin_sector(isector))=egs_values(isector)
       enddo
    endif
    if(ED_MPI_ID==0)then
       select case(ed_mode)
       case default
          write(unit,"(A)")"# i       E_i            exp(-(E-E0)/T)   nup ndw  Sect     Dim"
       case ("superc")
          write(unit,"(A)")"# i       E_i            exp(-(E-E0)/T)     Sz     Sect     Dim"
       case ("nonsu2")
          write(unit,"(A)")"# i       E_i            exp(-(E-E0)/T)      n     Sect     Dim"
       end select
       do istate=1,state_list%size
          Estate  = es_return_energy(state_list,istate)
          isector = es_return_sector(state_list,istate)
          select case(ed_mode)
          case default
             nup   = getnup(isector)
             ndw   = getndw(isector)
             write(unit,"(i3,f18.12,E18.9,1x,2i3,3x,i3,i10)")istate,Estate,exp(-beta*(Estate-state_list%emin)),nup,ndw,isector,getdim(isector)
          case("superc")
             sz   = getsz(isector)
             write(unit,"(i3,f18.12,E18.9,1x,i3,3x,i3,i10)")istate,Estate,exp(-beta*(Estate-state_list%emin)),sz,isector,getdim(isector)
          case("nonsu2")
             n    = getn(isector)
             write(unit,"(i3,f18.12,E18.9,1x,i3,3x,i3,i10)")istate,Estate,exp(-beta*(Estate-state_list%emin)),n,isector,getdim(isector)
          end select
       enddo
    endif
  end subroutine print_state_list



  !+-------------------------------------------------------------------+
  !PURPOSE  : analyse the spectrum and print some information after 
  !lanczos  diagonalization. 
  !+------------------------------------------------------------------+
  subroutine ed_analysis()
    integer             :: nup,ndw,sz,n,isector,dim
    integer             :: istate
    integer             :: i,unit
    integer             :: Nsize,NtoBremoved,nstates_below_cutoff
    integer             :: numgs
    real(8)             :: Egs,Ei,Ec
    type(histogram)     :: hist
    real(8)             :: hist_a,hist_b,hist_w
    integer             :: hist_n
    integer,allocatable :: list_sector(:),count_sector(:)    
    !POST PROCESSING:
    unit=free_unit()
    open(unit,file="state_list"//reg(ed_file_suffix)//".ed")
    call print_state_list(unit)
    close(unit)
    if(ed_verbose<=3)call print_state_list(LOGfile)
    !
    zeta_function=0d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei            = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       zeta_function=real(state_list%size,8)
    end if
    !
    numgs=es_return_gs_degeneracy(state_list,gs_threshold)
    if(numgs>Nsectors)stop "ed_diag: too many gs"
    do istate=1,numgs
       isector = es_return_sector(state_list,istate)
       Egs     = es_return_energy(state_list,istate)
       dim     = getdim(isector)
       select case(ed_mode)
       case default
          nup  = getnup(isector)
          ndw  = getndw(isector)
          if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A,F20.12,2I4)")'Egs =',Egs,nup,ndw
       case("superc")
          sz  = getsz(isector)
          if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A,F20.12,I4)")'Egs =',Egs,sz
       case("nonsu2")
          n  = getn(isector)
          if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A,F20.12,I4)")'Egs =',Egs,n
       end select
    enddo
    if(ed_verbose<3.AND.ED_MPI_ID==0)write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    !
    !
    !Get histogram distribution of the sector contributing to the evaluated spectrum:
    !Go thru states list and update the neigen_sector(isector) sector-by-sector
    if(finiteT)then
       if(ED_MPI_ID==0)then
          unit=free_unit()
          open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",access='append')
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
          if(ED_MPI_ID==0)write(LOGfile,"(A,I4)")"Increasing lanc_nstates_total:",lanc_nstates_total
       else
          !Find the energy level beyond which cutoff condition is verified & cut the list to that size
          do i=1,Nsize
             Ei = es_return_energy(state_list,i)
             if(exp(-beta*(Ei-Egs)) <= cutoff)exit
          enddo
          nstates_below_cutoff=i
          if(trim_state_list)lanc_nstates_total=max(nstates_below_cutoff,lanc_nstates_step)
          NtoBremoved=state_list%size - nstates_below_cutoff + 1
          do i=1,NtoBremoved
             call es_pop_state(state_list)
          enddo
          if(ed_verbose<4.AND.ED_MPI_ID==0.AND.trim_state_list)write(*,"(A,I4)")"Adjusting lanc_nstates_total to:",lanc_nstates_total
          if(ed_verbose<4.AND.ED_MPI_ID==0)write(*,"(A,I4)")"Trim list_size to         :",state_list%size
       endif
    endif
  end subroutine ed_analysis




end MODULE ED_DIAG
