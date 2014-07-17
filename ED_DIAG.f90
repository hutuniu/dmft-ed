!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE CONSTANTS
  USE MATRIX, only: matrix_diagonalize
  USE TIMER
  USE IOTOOLS, only:reg,free_unit
  USE STATISTICS
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

  public :: lanc_ed_diag

  real(8),allocatable,dimension(:) :: egs_values

contains

  !                    LANCZOS DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag
    logical :: iverbose_
    select case(ed_type)
    case default
       call lanc_ed_diag_d
    case('c')
       call lanc_ed_diag_c
    end select
    call lanc_ed_analysis
  end subroutine lanc_ed_diag



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag_d
    integer             :: nup,ndw,sz,n,isector,dim
    integer             :: i,j,iter
    integer             :: gs_count
    integer             :: Nitermax,Neigen,Nblock
    real(8)             :: oldzero,enemin,Egs,Ei,Ec
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    logical             :: lanc_solve,Tflag
    if(allocated(egs_values))deallocate(egs_values)
    allocate(egs_values(Nsectors))
    egs_values=0.d0
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    gs_count=0
    if(ed_verbose<2.AND.ED_MPI_ID==0)call start_progress(LOGfile)
    iter=0
    sector: do isector=1,Nsectors
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       iter=iter+1
       if(ED_MPI_ID==0)then
          if(ed_verbose==0)then
             call progress(iter,count(twin_mask))
          elseif(ed_verbose==-1)then
             select case(ed_mode)
             case default
                nup  = getnup(isector)
                ndw  = getndw(isector)
                write(LOGfile,"(1X,I4,A,I4,A6,I2,A6,I2,A6,I15)")iter,"-Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
             case ("superc")
                sz   = getsz(isector)
                write(LOGfile,"(1X,I4,A,I4,A5,I4,A6,I15)")iter,"-Solving sector:",isector," sz:",sz," dim=",getdim(isector)
             case ("nonsu2")
                n    = getn(isector)
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,I15)")iter,"-Solving sector:",isector," n:",n," dim=",getdim(isector)
             end select
          endif
       endif
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
       if(Neigen==dim.OR.dim<=max(16,ED_MPI_SIZE))lanc_solve=.false.
       !
       if(lanc_solve)then
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0.d0 ; eig_basis=0.d0
          call ed_buildH_d(isector)
#ifdef _MPI
          call lanczos_parpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_dd)
#else
          call lanczos_arpack( dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_dd)
#endif
          if(spH0%status)call sp_delete_matrix(spH0)
       else
          allocate(eig_values(dim),eig_basis(dim,dim))
          eig_values=0.d0 ; eig_basis=0.d0 
          call ed_buildH_d(isector,eig_basis)
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
             gs_count=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          elseif(abs(enemin-oldzero) <= gs_threshold)then
             gs_count=gs_count+1
             if (gs_count > Nsectors)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          endif
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    if(ed_verbose<2.AND.ED_MPI_ID==0)call stop_progress
  end subroutine lanc_ed_diag_d





  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag_c
    integer                :: nup,ndw,sz,n,isector,dim
    integer                :: i,j,iter
    integer                :: gs_count
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Egs,Ei,Ec
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:)
    logical                :: lanc_solve,Tflag
    if(allocated(egs_values))deallocate(egs_values)
    allocate(egs_values(Nsectors))
    egs_values=0.d0
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    gs_count=0
    if(ed_verbose<2.AND.ED_MPI_ID==0)call start_progress(LOGfile)
    iter=0
    sector: do isector=1,Nsectors
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       iter=iter+1
       if(ED_MPI_ID==0)then
          if(ed_verbose==0)then
             call progress(iter,count(twin_mask))
          elseif(ed_verbose==-1)then
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
                write(LOGfile,"(1X,I4,A,I4,A4,I4,A6,I15)")iter,"-Solving sector:",isector," n:",n," dim=",getdim(isector)
             end select
          endif
       endif
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
       if(Neigen==dim.OR.dim<=max(16,ED_MPI_SIZE))lanc_solve=.false.
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
             gs_count=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          elseif(abs(enemin-oldzero) <= gs_threshold)then
             gs_count=gs_count+1
             if (gs_count > Nsectors)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector,twin=Tflag)
          endif
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    if(ed_verbose<2.AND.ED_MPI_ID==0)call stop_progress
  end subroutine lanc_ed_diag_c









  !+-------------------------------------------------------------------+
  !PURPOSE  : analyse the spectrum and print some information after 
  !lanczos  diagonalization. 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_analysis()
    integer             :: nup,ndw,sz,n,isector,jsector,dim
    integer             :: istate
    integer             :: i,j,unit
    integer             :: gs_degeneracy
    real(8)             :: Egs,Estate,Emax
    logical             :: lanc_solve
    type(histogram)     :: hist
    real(8)             :: hist_a,hist_b,hist_w
    integer             :: hist_n
    integer,allocatable :: list_sector(:),count_sector(:)    
    !POST PROCESSING:
    if(ed_twin)then
       do isector=1,Nsectors
          if(.not.twin_mask(isector))egs_values(get_twin_sector(isector))=egs_values(isector)
       enddo
    endif
    if(ed_verbose<2.AND.ED_MPI_ID==0)then
       unit=free_unit()
       open(unit,file="state_list"//reg(ed_file_suffix)//".ed")
       select case(ed_mode)
       case default
          write(unit,"(A)")"#i       E_i           exp(-(E-E0)/T)    nup ndw Sect  Dim"
       case ("superc")
          write(unit,"(A)")"#i       E_i           exp(-(E-E0)/T)     Sz    Sect    Dim"
       case ("nonsu2")
          write(unit,"(A)")"#i       E_i           exp(-(E-E0)/T)      n    Sect    Dim"
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
       close(unit)
       !
       unit=free_unit()
       open(unit,file="egs_values"//reg(ed_file_suffix)//".ed")
       do isector=1,Nsectors
          write(unit,*)isector,egs_values(isector)
       end do
       close(unit)
       !
    endif
    zeta_function = 0.d0
    Egs           = state_list%emin
    if(finiteT)then
       do istate=1,state_list%size
          Estate        = es_return_energy(state_list,istate)
          zeta_function = zeta_function + exp(-beta*(Estate-Egs))
       enddo
    else
       zeta_function=real(state_list%size,8)
    end if
    !
    gs_degeneracy=es_return_gs_degeneracy(state_list,gs_threshold)
    if(gs_degeneracy>Nsectors)stop "ed_diag: too many gs"
    do istate=1,gs_degeneracy
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
    !Get histogram distribution of the sector contributing to the evaluated spectrum:
    !Go thru states list and update the neigen_sector(isector) sector-by-sector
    if(ed_verbose<2.AND.finiteT.AND.ED_MPI_ID==0)then
       unit=free_unit()
       open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",access='append')
       hist_n = Nsectors
       hist_a = 1.d0
       hist_b = real(Nsectors,8)
       hist_w = 1.d0
       hist = histogram_allocate(hist_n)
       call histogram_set_range_uniform(hist,hist_a,hist_b)
       do i=1,state_list%size
          isect0 = es_return_sector(state_list,i)
          call histogram_accumulate(hist,dble(isect0),hist_w)
       enddo
       call histogram_print(hist,unit)
       write(unit,*)""
       close(unit)
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
       Egs   = state_list%emin
       Emax  = state_list%emax
       if(exp(-beta*(Emax-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + 2!*lanc_nincrement
          if(ED_MPI_ID==0)write(*,"(A,I4)")"Increasing lanc_nstates_total+2:",lanc_nstates_total
       endif
    endif
  end subroutine lanc_ed_analysis




end MODULE ED_DIAG
