!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine build_gf_normal()
  integer :: iorb,jorb,ispin,i
  !
  !
  !NORMAL: (default)
  do ispin=1,Nspin
     do iorb=1,Norb
        if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
        select case(ed_type)
        case default
           call lanc_build_gf_normal_d(iorb,ispin)
        case ('c')
           call lanc_build_gf_normal_c(iorb,ispin)
        end select
     enddo
  enddo
  !
  !
  !HYBRID:
  if(bath_type=="hybrid")then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
              select case(ed_type)
              case default
                 call lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
              case ('c')
                 call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)                    
              end select
           enddo
        enddo
     enddo
     !Put here off-diagonal manipulation by symmetry:
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           impGmats(:,:,iorb,jorb,:) = 0.5d0*(impGmats(:,:,iorb,jorb,:) &
                - (one-xi)*impGmats(:,:,iorb,iorb,:) - (one-xi)*impGmats(:,:,jorb,jorb,:))
           impGreal(:,:,iorb,jorb,:) = 0.5d0*(impGreal(:,:,iorb,jorb,:) &
                - (one-xi)*impGreal(:,:,iorb,iorb,:) - (one-xi)*impGreal(:,:,jorb,jorb,:))
           impGmats(:,:,jorb,iorb,:) = impGmats(:,:,iorb,jorb,:)
           impGreal(:,:,jorb,iorb,:) = impGreal(:,:,iorb,jorb,:)
        enddo
     enddo
  endif
  !
  !
  !REPLICA:
  if(bath_type=="replica")then
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if((dmft_bath%mask(ispin,ispin,iorb,jorb,1).eqv. .true.).or.(dmft_bath%mask(ispin,ispin,iorb,jorb,2).eqv. .true.))then
                 if(ed_verbose<3.AND.MPI_MASTER)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                 select case(ed_type)
                 case default
                    call lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
                 case ('c')
                    call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)                    
                 end select
              endif
           enddo
        enddo
     enddo
     !Put here off-diagonal manipulation by symmetry:
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              if((dmft_bath%mask(ispin,ispin,iorb,jorb,1).eqv. .true.).or.(dmft_bath%mask(ispin,ispin,iorb,jorb,2).eqv. .true.))then
                 impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                      - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                 impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                      - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                 impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                 impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
              endif
           enddo
        enddo
     enddo
  endif
  !
end subroutine build_gf_normal






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_d(iorb,ispin)
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)  
  integer                          :: iorb,ispin,isite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map)                  :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_d



!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine nu_lanc_build_gf_normal_d(iorb,ispin)
  real(8),allocatable :: vvinit(:)
  integer             :: iorb,ispin,isite,istate
  integer             :: idim,isector
  integer             :: jdim,jsector
  integer             :: ib(Nlevels)
  integer             :: m,i,j,r
  real(8)             :: sgn,norm2
  complex(8)          :: cnorm2
  integer             :: Nlanc,Nitermax,Neigen,Nblock
  real(8),allocatable :: eig_val(:)
  real(8),allocatable :: eig_vec(:,:)
  type(sector_map)    :: HI,HJ
  !
  isite=impIndex(iorb,ispin)
  !
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  !
  do istate=1,state_list%size
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     !
     if(abs(sqrt(dot_product(state_vec,state_vec))-1d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !
     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then
        !
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        !
        call build_sector(jsector,HJ)
        !
        jdim  = getdim(jsector)
        allocate(vvinit(jdim));vvinit=0d0
        !
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        !Setup the Eigvalues and Eigvector to store info:
        Nlanc    = min(jdim,lanc_nGFiter)
        Nitermax = Nlanc
        Neigen   = 1
        Nblock   = min(jdim,2*Neigen + 1)
        if(allocated(eig_val))deallocate(eig_val)
        if(allocated(eig_vec))deallocate(eig_vec)
        allocate(eig_val(Neigen))      ; eig_val=0d0
        allocate(eig_vec(Nlanc,Neigen)) ; eig_vec=0d0
        if(MpiStatus)then
           call sp_eigh(MpiComm,spHtimesV_dd,Nlanc,Neigen,Nblock,Nitermax,eig_val,eig_vec,tol=lanc_tolerance,v0=vvinit)
        else
           call sp_eigh(spHtimesV_dd,Nlanc,Neigen,Nblock,Nitermax,eig_val,eig_vec,tol=lanc_tolerance,v0=vvinit)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal_bis(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)




        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        !Build H in the Isector
        call ed_buildH_d(isector)

        !Setup the Eigvalues and Eigvector to store info:
        Neigen   = 1
        Nitermax = min(dim,lanc_nGFiter)
        Nblock   = min(dim,2*Neigen + 1)
        if(allocated(eig_values))deallocate(eig_values)
        if(allocated(eig_vector))deallocate(eig_vector)
        allocate(eig_values(Neigen),eig_vector(iDim,Neigen))
        eig_values=0d0 ; eig_vector=0d0
        if(MpiStatus)then
           call sp_eigh(MpiComm,spHtimesV_dd,iDim,Neigen,Nblock,Nitermax,eig_values,eig_vector,tol=lanc_tolerance,v0=vvinit)
        else
           call sp_eigh(spHtimesV_dd,iDim,Neigen,Nblock,Nitermax,eig_values,eig_vector,tol=lanc_tolerance,v0=vvinit)
        endif
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal_bis(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)

        ! alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        ! call ed_buildH_d(jsector)
        ! call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        ! cnorm2=one*norm2
        ! call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)



        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine nu_lanc_build_gf_normal_d






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_c(iorb,ispin)
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: iorb,ispin,isite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  isite=impIndex(iorb,ispin)
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ) !note that here you are doing twice the map building...
        vvinit=0.d0
        do m=1,idim                     !loop over |gs> components m
           i=HI%map(m)                    !map m to Hilbert space state i
           ib = bdecomp(i,2*Ns)            !i into binary representation
           if(ib(isite)==0)then          !if impurity is empty: proceed
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)      !map r back to  jsector
              vvinit(j) = sgn*state_cvec(m)  !build the cdg_up|gs> state
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(LOGfile,"(A,2I3,I15)")' del particle:',&
             getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_c






!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE PRECISION
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_mix_d(iorb,jorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  real(8),allocatable              :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_vec  => es_return_vector(state_list,istate)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dd,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_vec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_d(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_dc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_vec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_mix_d



!+------------------------------------------------------------------+
!PURPOSE  : DOUBLE COMPLEX
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  isite=impIndex(iorb,ispin)  !orbital 1
  jsite=impIndex(jorb,ispin)  !orbital 2
  !
  numstates=state_list%size
  !   
  if(ed_verbose<3.AND.MPI_MASTER)call start_timer
  do istate=1,numstates
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !

     !EVALUATE (c^+_iorb + c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb + c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=0.d0
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,vvinit,alfa_,beta_,nlanc)
        cnorm2=one*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(vvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==0)then
              call cdg(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==0)then
              call cdg(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif

     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        if(ed_verbose<1.AND.MPI_MASTER)write(*,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(isite)==1)then
              call c(isite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(jsite)==1)then
              call c(jsite,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        alfa_=0.d0 ; beta_=0.d0 ; nlanc=min(jdim,nitermax)
        call ed_buildH_c(jsector)
        call sp_lanc_tridiag(lanc_spHtimesV_cc,cvinit,alfa_,beta_,nlanc)
        cnorm2=-xi*norm2
        call add_to_lanczos_gf_normal(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
        deallocate(cvinit)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  if(ed_verbose<3.AND.MPI_MASTER)call stop_timer
  deallocate(alfa_,beta_)
end subroutine lanc_build_gf_normal_mix_c








