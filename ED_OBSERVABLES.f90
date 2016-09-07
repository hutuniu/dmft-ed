!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_MATVEC
  USE ED_BATH
  implicit none
  private
  !
  public                             :: observables_impurity
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


contains 


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer,dimension(Nlevels)       :: ib
    integer                          :: i,j
    integer                          :: izero
    integer                          :: isector,jsector
    integer                          :: idim,jdim
    integer                          :: isz,jsz
    integer                          :: iorb,jorb,ispin,jspin,isite,jsite
    integer                          :: numstates
    integer                          :: r,m,k
    real(8)                          :: sgn,sgn1,sgn2
    real(8)                          :: gs_weight
    real(8)                          :: Ei
    real(8)                          :: peso
    real(8)                          :: norm
    real(8),dimension(Norb)          :: nup,ndw,Sz,nt
    real(8),dimension(:),pointer     :: gsvec
    complex(8),dimension(:),pointer  :: gscvec
    integer,allocatable,dimension(:) :: Hmap,HJmap
    real(8),allocatable              :: vvinit(:)
    !!<DEBUG
    !logical :: converged
    !real(8) :: pdens
    !>DEBUG

    !
    !LOCAL OBSERVABLES:
    ! density, 
    ! double occupancy, 
    ! magnetization, 
    ! orbital//spin correlations  
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(phisc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    phisc   = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
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
       !pdens=0d0
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
          !pdens     = pdens      +  nt(1)*gs_weight*zeta_function
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
       !!<DEBUG  comment
       !print*,"sectors contribution to dens:"
       !select case(ed_mode)
       !case default
       !   print*,isector,getnup(isector),getndw(isector),pdens       
       !case ("superc")
       !   print*,isector,getsz(isector),pdens       
       !case("nonsu2")
       !   print*,isector,getn(isector),pdens       
       !end select
       !!>DEBUG
       if(associated(gsvec))nullify(gsvec)
       if(associated(gscvec))nullify(gscvec)
       deallocate(Hmap)
    enddo
    !
    !SUPERCONDUCTING ORDER PARAMETER
    if(ed_mode=="superc")then
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
                         vvinit(j) = sgn*gsvec(i)
                      endif
                   enddo
                   do i=1,idim
                      m=Hmap(i)
                      call bdecomp(m,ib)
                      if(ib(iorb+Ns)==1)then
                         call c(iorb+Ns,m,r,sgn)
                         j=binary_search(HJmap,r)
                         vvinit(j) = vvinit(j) + sgn*gsvec(i)
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
    !<<DEBUG
    !IMPURITY DENSITY MATRIX
    if(allocated(imp_density_matrix)) deallocate(imp_density_matrix);allocate(imp_density_matrix(Nspin,Nspin,Norb,Norb))
    imp_density_matrix=zero
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
       !Diagonal densities
       do ispin=1,Nspin
          do iorb=1,Norb
             isite=impIndex(iorb,ispin)
             do m=1,idim
                i=Hmap(m)
                call bdecomp(i,ib)
                if(ed_type=='d')imp_density_matrix(ispin,ispin,iorb,iorb) = imp_density_matrix(ispin,ispin,iorb,iorb) + peso*ib(isite)*gsvec(m)*gsvec(m)
                if(ed_type=='c')imp_density_matrix(ispin,ispin,iorb,iorb) = imp_density_matrix(ispin,ispin,iorb,iorb) + peso*ib(isite)*conjg(gscvec(m))*gscvec(m)
             enddo
          enddo
       enddo
       !off-diagonal
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ed_mode=="normal").and.(ispin/=jspin))cycle !  ed_mode=="normal" ==>    spin off-dig term not calculated
                   if((bath_type=="normal").and.(iorb/=jorb))cycle !bath_type=="normal" ==> orbital off-dig term not calculated
                   isite=impIndex(iorb,ispin)
                   jsite=impIndex(jorb,jspin)
                   do m=1,idim
                      i=Hmap(m)
                      call bdecomp(i,ib)
                      if((ib(isite)==1).and.(ib(jsite)==0))then
                         call c(isite,i,r,sgn1)
                         call cdg(jsite,r,k,sgn2)
                         j=binary_search(Hmap,k)
                         if(ed_type=='d')imp_density_matrix(ispin,jspin,iorb,jorb) = imp_density_matrix(ispin,jspin,iorb,jorb) + peso*sgn1*gsvec(m)*sgn2*gsvec(j)
                         if(ed_type=='c')imp_density_matrix(ispin,jspin,iorb,jorb) = imp_density_matrix(ispin,jspin,iorb,jorb) + peso*sgn1*gscvec(m)*sgn2*conjg(gscvec(j))
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       deallocate(Hmap)
    enddo
    imp_density_matrix = imp_density_matrix/float(numstates)
    !IMPURITY DENSITY OPERATORS
    if((Nspin/=1).and.(Norb==3))then
    if(allocated(impStot))    deallocate(impStot);   allocate(impStot(3,Norb,Norb));  impStot=zero
    if(allocated(impLtot))    deallocate(impLtot);   allocate(impLtot(3,Nspin,Nspin));impLtot=zero
    if(allocated(impj_aplha)) deallocate(impj_aplha);allocate(impj_aplha(3));         impj_aplha=zero
    impLdotS=zero
    !
    !#####################################################
    !#                    S(iorb,jorb)                   #
    !#####################################################
    !
    !Sx =    [ <c+_up,c_dw> + <c+_dw,c_up> ]_(iorb,jorb)
    !Sy = xi*[ <c+_dw,c_up> - <c+_up,c_dw> ]_(iorb,jorb)
    !Sz =    [ <c+_up,c_up> - <c+_dw,c_dw> ]_(iorb,jorb)
    !
    do iorb=1,Norb
       do jorb=1,Norb
          if(ed_mode=="normal")cycle                     !  ed_mode=="normal" ==>    spin off-dig term not calculated
          if((bath_type=="normal").and.(iorb/=jorb))cycle!bath_type=="normal" ==> orbital off-dig term not calculated
          impStot(1,iorb,jorb) = 0.5d0*( imp_density_matrix(1,2,iorb,jorb) + imp_density_matrix(2,1,iorb,jorb) )
          impStot(2,iorb,jorb) = 0.5d0*( imp_density_matrix(2,1,iorb,jorb) - imp_density_matrix(1,2,iorb,jorb) )*xi
          impStot(3,iorb,jorb) = 0.5d0*( imp_density_matrix(1,1,iorb,jorb) - imp_density_matrix(2,2,iorb,jorb) )
       enddo
    enddo
    !
    !#####################################################
    !#                   L(ispin,jspin)                  #
    !#####################################################
    !1=yz 2=zx 3=xy
    !Lx = xi*[ <c+_3,c_2> - <c+_2,c_3> ]_(ispin,jspin)
    !Ly = xi*[ <c+_1,c_3> - <c+_3,c_1> ]_(ispin,jspin)
    !Lz = xi*[ <c+_2,c_1> - <c+_1,c_2> ]_(ispin,jspin)
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          if((ed_mode=="normal").and.(ispin/=jspin))cycle!  ed_mode=="normal" ==>    spin off-dig term not calculated
          if(bath_type=="normal")cycle                   !bath_type=="normal" ==> orbital off-dig term not calculated
          impLtot(1,ispin,jspin) = ( imp_density_matrix(ispin,jspin,3,2) - imp_density_matrix(ispin,jspin,2,3) )*xi
          impLtot(2,ispin,jspin) = ( imp_density_matrix(ispin,jspin,1,3) - imp_density_matrix(ispin,jspin,3,1) )*xi
          impLtot(3,ispin,jspin) = ( imp_density_matrix(ispin,jspin,2,1) - imp_density_matrix(ispin,jspin,1,2) )*xi
       enddo
    enddo
    !
    !#####################################################
    !#                        LdotS                      #
    !#####################################################
    ! 1=yz 2=zx 3=xy
    ! + xi*[ <c+_3up,c_1up> - <c+_1up,c_3up> ]  ==  spin - diagonal 1
    ! + xi*[ <c+_2dw,c_1dw> - <c+_1dw,c_2dw> ]  ==  spin - diagonal 2
    ! -    [ <c+_2up,c_2dw> + <c+_2dw,c_2up> ]  ==  orb  - diagonal 1
    ! + xi*[ <c+_3dw,c_3up> - <c+_3up,c_3dw> ]  ==  orb  - diagonal 2
    ! -    [ <c+_3dw,c_1up> + <c+_1up,c_3dw> ]  ==  full off - diagonal 1
    ! + xi*[ <c+_2up,c_1dw> - <c+_1dw,c_2up> ]  ==  full off - diagonal 2
    !
    impLdotS =  (imp_density_matrix(1,1,3,1)-imp_density_matrix(1,1,1,3))*xi &
               +(imp_density_matrix(2,2,2,1)-imp_density_matrix(2,2,1,2))*xi &
               -(imp_density_matrix(1,2,2,2)+imp_density_matrix(2,1,2,2))    &
               +(imp_density_matrix(2,1,3,3)-imp_density_matrix(1,2,3,3))*xi &
               -(imp_density_matrix(2,1,3,1)+imp_density_matrix(1,2,1,3))    &
               +(imp_density_matrix(1,2,2,1)-imp_density_matrix(2,1,1,2))*xi
    impLdotS = impLdotS/2.d0
    !
    !#####################################################
    !#              ja=trace{La}+trace{Sa}               #
    !#####################################################
    !
    impj_aplha(1)=trace(impStot(1,:,:))+trace(impLtot(1,:,:))
    impj_aplha(2)=trace(impStot(2,:,:))+trace(impLtot(2,:,:))
    impj_aplha(3)=trace(impStot(3,:,:))+trace(impLtot(3,:,:))
    !
    endif
    !<<DEBUG
    !
    call get_szr
    if(ED_MPI_ID==0)then
       if(iolegend)call write_legend
       call write_observables()
       write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
       select case(ed_mode)
       case default
          write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
       case("superc")
          write(LOGfile,"(A,20f18.12,A)")    "phi "//reg(ed_file_suffix)//"=",(phisc(iorb),iorb=1,Norb),(abs(uloc(iorb))*phisc(iorb),iorb=1,Norb)
       end select
       if(ed_verbose<3)then
          if(Nspin==2)then
             write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
             write(LOGfile,"(A,10f18.12,A)") " Ji "//reg(ed_file_suffix)//"=",(real(impj_aplha(i)),i=1,3)
          endif
       endif
    endif
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_phisc(iorb)  =phisc(iorb)
    enddo
    !
    deallocate(dens,docc,phisc,dens_up,dens_dw,magz,sz2,n2)
    deallocate(simp,zimp)    
  end subroutine observables_impurity

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
    open(unit,file="observables_info.ed")
    select case(ed_mode)
    case default
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
    case ("superc")
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
    end select
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
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
       open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh
       close(unit)
       !
       unit = free_unit()
       open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
       select case(ed_mode)
       case default
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
       case ("superc")
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
       end select
       close(unit)    
    endif
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    select case(ed_mode)
    case default
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
    case ("superc")
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
    end select
    close(unit)         
  end subroutine write_observables




end MODULE ED_OBSERVABLES
