!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE CONSTANTS, only:zero,pi,xi
  USE IOTOOLS, only:free_unit,reg,txtfy
  USE ARRAYS, only: arange
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_MATVEC
  implicit none
  private
  !
  public                             :: observables_impurity
  public                             :: energy_impurity
  !
  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:),allocatable   :: phisc
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimensioN(:,:),allocatable :: zimp,simp
  real(8)                            :: s2tot
  real(8)                            :: Egs
  real(8)                            :: Ehartree
  real(8)                            :: Eknot
  real(8)                            :: Epot
  real(8)                            :: Eint
  real(8)                            :: Dust,Dund,Dse,Dph
  real(8)                            :: Ekin
  real(8)                            :: Etot


contains 


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer,dimension(Ntot)          :: ib
    integer                          :: i,j
    integer                          :: izero
    integer                          :: isector,jsector
    integer                          :: idim,jdim
    integer                          :: isz,jsz
    integer                          :: iorb,jorb,ispin
    integer                          :: numstates
    integer                          :: r,m
    real(8)                          :: sgn
    real(8)                          :: gs_weight
    real(8)                          :: Ei
    real(8)                          :: peso
    real(8)                          :: norm
    real(8),dimension(Norb)          :: nup,ndw,Sz,nt
    real(8),dimension(:),pointer     :: gsvec
    complex(8),dimension(:),pointer  :: gscvec
    integer,allocatable,dimension(:) :: Hmap,HJmap
    real(8),allocatable              :: vvinit(:)
    !
    !LOCAL OBSERVABLES:
    ! density, 
    ! double occupancy, 
    ! magnetization, 
    ! orbital//spin correlations  
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    if(ed_supercond) allocate(phisc(Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc   = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    numstates=state_list%size
    do izero=1,numstates
       isector = es_return_sector(state_list,izero)
       Ei      = es_return_energy(state_list,izero)
       idim     = getdim(isector)
       !
       if(ed_type=='d')then
          gsvec  => es_return_vector(state_list,izero)
          norm=sqrt(dot_product(gsvec,gsvec))
       elseif(ed_type=='c')then
          gscvec  => es_return_cvector(state_list,izero)
          norm=sqrt(dot_product(gscvec,gscvec))
       endif
       if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       allocate(Hmap(idim))
       call build_sector(isector,Hmap)
       !
       do i=1,idim
          m=Hmap(i)
          call bdecomp(m,ib)
          !
          if(ed_type=='d')then
             gs_weight=peso*gsvec(i)**2
          elseif(ed_type=='c')then
             gs_weight=peso*abs(gscvec(i))**2
          endif
          !
          !Get operators:
          do iorb=1,Norb
             nup(iorb)= dble(ib(iorb))
             ndw(iorb)= dble(ib(iorb+Ns))
             sz(iorb) = (nup(iorb) - ndw(iorb))/2.d0
             nt(iorb) =  nup(iorb) + ndw(iorb)
          enddo
          !
          !Evaluate averages of observables:
          do iorb=1,Norb
             dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
             dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
             dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
             docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
             magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
             sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
             n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
             do jorb=iorb+1,Norb
                sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
             enddo
          enddo
          s2tot = s2tot  + (sum(sz))**2*gs_weight
       enddo
       if(associated(gsvec))nullify(gsvec)
       if(associated(gscvec))nullify(gscvec)
       deallocate(Hmap)
    enddo
    !
    !SUPERCONDUCTING ORDER PARAMETER
    if(ed_supercond) then
       phisc = 0.d0
       do ispin=1,Nspin
          do iorb=1,Norb
             numstates=state_list%size
             do izero=1,numstates
                !
                isector = es_return_sector(state_list,izero)
                Ei      = es_return_energy(state_list,izero)
                idim    = getdim(isector)
                if(ed_type=='d')then
                   gsvec  => es_return_vector(state_list,izero)
                   norm=sqrt(dot_product(gsvec,gsvec))
                elseif(ed_type=='c')then
                   gscvec  => es_return_cvector(state_list,izero)
                   norm=sqrt(dot_product(gscvec,gscvec))
                endif
                if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
                !
                peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
                peso = peso/zeta_function
                !
                allocate(Hmap(idim))
                call build_sector(isector,Hmap)
                !GET <(C_UP + CDG_DW)(CDG_UP + C_DW)> = 
                !<C_UP*CDG_UP> + <CDG_DW*C_DW> + <C_UP*C_DW> + <CDG_DW*CDG_UP> = 
                !<N_UP> + < 1 - N_DW> + 2*<PHI>
                isz = getsz(isector)
                if(isz<Ns)then
                   jsz     = isz+1
                   jsector = getsector(jsz,1)
                   jdim    = getdim(jsector)
                   allocate(HJmap(jdim),vvinit(jdim))
                   call build_sector(jsector,HJmap)
                   vvinit=0.d0
                   do i=1,idim
                      m=Hmap(i)
                      call bdecomp(m,ib)
                      if(ib(iorb)==0)then
                         call cdg(iorb,m,r,sgn)
                         j=binary_search(HJmap,r)
                         vvinit(j) = sgn*gsvec(m)
                      endif
                   enddo
                   do i=1,idim
                      m=Hmap(i)
                      call bdecomp(m,ib)
                      if(ib(iorb+Ns)==1)then
                         call c(iorb+Ns,m,r,sgn)
                         j=binary_search(HJmap,r)
                         vvinit(j) = vvinit(j) + sgn*gsvec(m)
                      endif
                   enddo
                   deallocate(HJmap)
                   phisc(iorb) = phisc(iorb) + dot_product(vvinit,vvinit)*peso
                   deallocate(vvinit)
                endif
                if(associated(gsvec)) nullify(gsvec)
                deallocate(Hmap)
                !
             enddo
             phisc(iorb) = 0.5d0*(phisc(iorb) - dens_up(iorb) - (1.d0-dens_dw(iorb)))
          enddo
       enddo
    end if
    !
    !
    call get_szr
    if(ED_MPI_ID==0)then
       if(iolegend)call write_legend
       call write_observables()
       write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
       if(ed_supercond)then
          write(LOGfile,"(A,20f18.12,A)")    "phi "//reg(ed_file_suffix)//"=",(phisc(iorb),iorb=1,Norb),(uloc(iorb)*phisc(iorb),iorb=1,Norb)
       else
          write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)       
       endif
       if(ed_verbose<0)then
          write(LOGfile,"(A,20f18.12,A)")    "sz2 "//reg(ed_file_suffix)//"=",((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
          if(Nspin==2)then
             write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
          endif
       endif
    endif
    !
    do iorb=1,Norb
       ed_dens(iorb)=dens(iorb)
       ed_docc(iorb)=docc(iorb)
       if(ed_supercond)ed_phisc(iorb)=phisc(iorb)
    enddo
    !
#ifdef _MPI
    call MPI_BCAST(ed_dens,Norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
    call MPI_BCAST(ed_docc,Norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
    if(ed_supercond)call MPI_BCAST(ed_phisc,Norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ED_MPI_ERR)
#endif
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
    deallocate(simp,zimp)
    if(ed_supercond)deallocate(phisc)
  end subroutine observables_impurity






  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine energy_impurity()
    integer,dimension(Ntot)                           :: ib
    integer                                           :: i,j
    integer                                           :: izero
    integer                                           :: isector
    integer                                           :: idim
    integer                                           :: iorb,jorb,ispin
    integer                                           :: numstates
    integer                                           :: m,k1,k2,k3,k4
    real(8)                                           :: sg1,sg2,sg3,sg4
    real(8)                                           :: gs_weight,spin_degenereacy
    real(8)                                           :: Ei
    real(8)                                           :: peso
    real(8)                                           :: norm
    real(8),dimension(Norb)                           :: nup,ndw
    real(8),dimension(Nspin,Norb)                     :: eloc
    real(8),dimension(:),pointer                      :: gsvec
    complex(8),dimension(:),pointer                   :: gscvec
    integer,allocatable,dimension(:)                  :: Hmap
    logical                                           :: Jcondition
    real(8),dimension(Lmats)                          :: wm
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impDelta,impKappa
    !
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    !
    Egs     = state_list%emin
    Ehartree= 0.d0
    Eknot   = 0.d0
    Epot    = 0.d0
    Eint    = 0.d0
    Dust    = 0.d0
    Dund    = 0.d0
    Dse     = 0.d0
    Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=Hloc(ispin,ispin,iorb,iorb)
       enddo
    enddo

    numstates=state_list%size
    do izero=1,numstates
       isector = es_return_sector(state_list,izero)
       Ei      = es_return_energy(state_list,izero)
       idim    = getdim(isector)
       !
       if(ed_type=='d')then
          gsvec  => es_return_vector(state_list,izero)
          norm=sqrt(dot_product(gsvec,gsvec))
       elseif(ed_type=='c')then
          gscvec  => es_return_cvector(state_list,izero)
          norm=sqrt(dot_product(gscvec,gscvec))
       endif
       if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       allocate(Hmap(idim))
       call build_sector(isector,Hmap)
       !
       do i=1,idim
          m=Hmap(i)
          call bdecomp(m,ib)
          !
          if(ed_type=='d')then
             gs_weight=peso*gsvec(i)**2
          elseif(ed_type=='c')then
             gs_weight=peso*abs(gscvec(i))**2
          endif
          !
          !Get operators:
          do iorb=1,Norb
             nup(iorb)= dble(ib(iorb))
             ndw(iorb)= dble(ib(iorb+Ns))
          enddo
          !
          !start evaluating the Tr(H_loc) to estimate potential energy
          !
          !LOCAL ENERGY
          Eknot = Eknot + dot_product(eloc(1,:),nup)*gs_weight + dot_product(eloc(Nspin,:),ndw)*gs_weight
          !
          !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
          !Euloc=\sum=i U_i*(n_u*n_d)_i
          if(.not.ed_supercond) then
             Epot = Epot + dot_product(uloc,nup*ndw)*gs_weight
          else
             Epot = Epot - uloc(1)*(nup(1)-0.5d0)*(ndw(1)-0.5d0)*gs_weight
          end if
          !
          !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
          !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
          !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   Epot = Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   Dust = Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                enddo
             enddo
          endif
          !
          !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
          !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   Epot = Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   Dund = Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                enddo
             enddo
          endif
          !
          !SPIN-EXCHANGE (S-E) TERMS
          !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
          if(Norb>1.AND.Jhflag)then
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=((iorb/=jorb).AND.&
                        (ib(jorb)==1)      .AND.&
                        (ib(iorb+Ns)==1)   .AND.&
                        (ib(jorb+Ns)==0)   .AND.&
                        (ib(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,m,k1,sg1)
                      call c(iorb+Ns,k1,k2,sg2)
                      call cdg(jorb+Ns,k2,k3,sg3)
                      call cdg(iorb,k3,k4,sg4)
                      j=binary_search(Hmap,k4)
                      Epot = Epot + Jh*sg1*sg2*sg3*sg4*gs_weight
                      Dse  = Dse  + sg1*sg2*sg3*sg4*gs_weight
                   endif
                enddo
             enddo
          endif
          !
          !PAIR-HOPPING (P-H) TERMS
          !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
          !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
          if(Norb>1.AND.Jhflag)then
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=((iorb/=jorb).AND.&
                        (ib(jorb)==1)      .AND.&
                        (ib(jorb+Ns)==1)   .AND.&
                        (ib(iorb+Ns)==0)   .AND.&
                        (ib(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,m,k1,sg1)
                      call c(jorb+Ns,k1,k2,sg2)
                      call cdg(iorb+Ns,k2,k3,sg3)
                      call cdg(iorb,k3,k4,sg4)
                      j=binary_search(Hmap,k4)
                      Epot = Epot + Jh*sg1*sg2*sg3*sg4*gs_weight
                      Dph  = Dph  + sg1*sg2*sg3*sg4*gs_weight
                   endif
                enddo
             enddo
          endif
          !
          !HARTREE-TERMS CONTRIBUTION:
          if(hfmode.and..not.ed_supercond)then
             Ehartree=Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      Ehartree=Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*Ust*gs_weight
                      Ehartree=Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*(Ust-Jh)*gs_weight
                   enddo
                enddo
             endif
          endif
       enddo
       if(associated(gsvec))nullify(gsvec)
       if(associated(gscvec))nullify(gscvec)
       deallocate(Hmap)
    enddo
    Epot = Epot + Ehartree

    Ekin = 0.d0
    !Get Delta function: we only need the orbital diagonal part as Ekin = Tr(Delta*G)
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Lmats
                impDelta(ispin,ispin,iorb,iorb,i) = delta_bath_mats(ispin,iorb,xi*wm(i),dmft_bath)
             enddo
          enddo
       enddo
    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Lmats
                impDelta(ispin,ispin,iorb,iorb,i) = delta_bath_mats(ispin,iorb,iorb,xi*wm(i),dmft_bath)
             enddo
          enddo
       enddo
    end select

    spin_degenereacy=dble(3-Nspin)            !spin_degeneracy=2 if Nspin=1, spin_degeneracy=1 if Nspin=2
    do ispin=1,Nspin
       do iorb=1,Norb
          Ekin = Ekin + spin_degenereacy*2.0*sum(dreal(impDelta(ispin,ispin,iorb,iorb,:)*impGmats(ispin,ispin,iorb,iorb,:) ))/beta
       enddo
    enddo

    Etot = Ekin +  Epot + Eknot

    if(ed_verbose<0)then
       write(LOGfile,"(A,10f18.12)")"<K>     =",Ekin
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",Epot!,Uloc(1)*(ed_docc(1)-ed_dens(1)/2.d0+0.25d0)
       write(LOGfile,"(A,10f18.12)")"<H>     =",Etot
       write(LOGfile,"(A,10f18.12)")"<V>     =",Epot-Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",Dund
       write(LOGfile,"(A,10f18.12)")"Dse     =",Dse
       write(LOGfile,"(A,10f18.12)")"Dph     =",Dph
    endif

    call write_energy()
  end subroutine energy_impurity









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3.d0*pi/beta
    do ispin=1,Nspin
       do iorb=1,Norb
          simp(iorb,ispin) = dimag(impSmats(ispin,ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,ispin,iorb,iorb,2))-dimag(impSmats(ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="columns_info.ed")
    if(.not.ed_supercond)then
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
            reg(txtfy(5*Norb+1))//"s2",&
            reg(txtfy(5*Norb+2))//"egs",&
            ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((5+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            ((reg(txtfy((6+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    else
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(Norb+iorb))//"phi_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(2*Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(3*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(4*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(5*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
            reg(txtfy(6*Norb+1))//"s2",&
            reg(txtfy(6*Norb+2))//"egs",&
            ((reg(txtfy(6*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((6+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((6+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            ((reg(txtfy((7+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    endif
    close(unit)
    !
    unit = free_unit()
    open(unit,file="control_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         (reg(txtfy(iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(Norb+1))//"U'",&
         reg(txtfy(Norb+2))//"Jh",&
         reg(txtfy(Norb+3))//"<K>",&
         reg(txtfy(Norb+4))//"<Hi>",&
         reg(txtfy(Norb+5))//"<H>",&
         reg(txtfy(Norb+6))//"<V>=<Hi-Ehf>",&
         reg(txtfy(Norb+7))//"<E0>",&
         reg(txtfy(Norb+8))//"<Ehf>",&
         reg(txtfy(Norb+9))//"<Dst>",&
         reg(txtfy(Norb+10))//"<Dnd>",&
         reg(txtfy(Norb+11))//"<Dse>",&
         reg(txtfy(Norb+12))//"<Dph>"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    if(ed_verbose<2)then
       unit = free_unit()
       open(unit,file="control_vars"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh
       close(unit)
       !
       unit = free_unit()
       open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
       if(.not.ed_supercond)then
          write(unit,"(90(F15.9,1X))")&
               (dens(iorb),iorb=1,Norb),&
               (docc(iorb),iorb=1,Norb),&
               (dens_up(iorb),iorb=1,Norb),&
               (dens_dw(iorb),iorb=1,Norb),&
               (magz(iorb),iorb=1,Norb),&
               s2tot,egs,&
               ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
               ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
               ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
               ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       else
          write(unit,"(90(F15.9,1X))")&
               (dens(iorb),iorb=1,Norb),&
               (phisc(iorb),iorb=1,Norb),&
               (docc(iorb),iorb=1,Norb),&
               (dens_up(iorb),iorb=1,Norb),&
               (dens_dw(iorb),iorb=1,Norb),&
               (magz(iorb),iorb=1,Norb),&
               s2tot,egs,&
               ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
               ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
               ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
               ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       endif
       close(unit)    
    endif
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    if(.not.ed_supercond)then
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    else
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (phisc(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    endif
    close(unit)         
  end subroutine write_observables




  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_energy()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="energy_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90F15.9)")(uloc(iorb),iorb=1,Norb),Ust,Jh,Ekin,Epot,Etot,Epot-Ehartree,Eknot,Ehartree,Dust,Dund,Dse,Dph
    close(unit)
    !
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")(uloc(iorb),iorb=1,Norb),Ust,Jh,Ekin,Epot,Etot,Epot-Ehartree,Eknot,Ehartree,Dust,Dund,Dse,Dph
    close(unit)
  end subroutine write_energy

end MODULE ED_OBSERVABLES
