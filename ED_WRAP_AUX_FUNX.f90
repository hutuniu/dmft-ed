module ED_WRAP_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_IOTOOLS, only: free_unit
  implicit none
  private

  interface extract_Hloc
     module procedure extract_Hloc_1
     module procedure extract_Hloc_2
  end interface extract_Hloc

  interface reshape_Hloc
     module procedure reshape_Hloc_1
     module procedure reshape_Hloc_2
  end interface reshape_Hloc

  !LATTICE:
  public :: blocks_to_matrix
  public :: matrix_to_blocks
  public :: select_block
  public :: stride_index
  public :: extract_Hloc
  public :: reshape_Hloc
  public :: get_independent_sites  

  !OBSOLETE (to be removed)
  public :: get_lattice_hamiltonian      !+- produce the 2D square/ slab-chain Hamiltonian -+! 

contains


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! perform all the reduction and broadcast from Nlat blocks of size Nso=Nspin*Norb 
  ! (spin-orbital blocks) into a large matrix of order Nlat*Nso (just like Hk)
  !+-----------------------------------------------------------------------------+!
  function blocks_to_matrix(Vblocks) result(Matrix)
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    integer                                               :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nspin*Norb
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
       i = 1 + (ip-1)*Nspin*Norb
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




  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! - strid_index: given an array of indices ivec=[i_1,...,i_L] and the corresponding 
  !                array of ranges nvec=[N_1,...,N_L], evaluate the index of the 
  !                stride corresponding the actual ivec with respect to nvec as 
  !                I_stride = i1 + \sum_{k=2}^L (i_k-1)\prod_{l=1}^{k-1}N_k
  !                assuming that the  i_1 .inner. i_2 .... .inner i_L (i_1 fastest index, 
  !                i_L slowest index)
  ! - lat_orb2lo: evaluate the stride for the simplest case of Norb and Nlat
  !+-----------------------------------------------------------------------------+!
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



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  
  ! get the block diagonal local part of a given Hamiltonian H(k_\perp;Ri,Rj)
  !+-----------------------------------------------------------------------------+!
  function extract_Hloc_1(Hk,Nlat,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nlat,Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: i,iorb,ispin,ilat,is
    integer                                     :: j,jorb,jspin,jlat,js
    Hloc = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc_1

  function extract_Hloc_2(Hk,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: i,iorb,ispin,ilat,is
    integer                                     :: j,jorb,jspin,jlat,js
    Hloc = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i = stride_index([iorb,ispin],[Norb,Nspin])
                j = stride_index([jorb,jspin],[Norb,Nspin])
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                print*,i,is,j,js
                Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc_2


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape the block diagonal local part of a Hamiltonian from the [Nlso][Nlso]
  ! shape to the [Nlat][Nspin][Nspin][Norb][Norb] shape used in the ED code.
  !+-----------------------------------------------------------------------------+!
  function reshape_Hloc_1(Hloc,Nlat,Nspin,Norb) result(edHloc)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: edHloc
    !
    integer                                               :: i,iorb,ispin,ilat,is
    integer                                               :: j,jorb,jspin,jlat,js
    edHloc=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   edHloc(ilat,ispin,jspin,iorb,jorb) = Hloc(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function reshape_Hloc_1

  function reshape_Hloc_2(Hloc,Nspin,Norb) result(edHloc)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hloc
    integer                                     :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: edHloc
    !
    integer                                     :: i,iorb,ispin,ilat,is
    integer                                     :: j,jorb,jspin,jlat,js
    edHloc=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                edHloc(ispin,jspin,iorb,jorb) = Hloc(is,js)
             enddo
          enddo
       enddo
    enddo
  end function reshape_Hloc_2


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  
  ! find all inequivalent sites with respect the user defined symmetry operations 
  ! Build and check maps from the full(independent) lattice to the independent 
  ! (full) lattice                        !
  !+-----------------------------------------------------------------------------+!
  subroutine get_independent_sites(symmetry_operations,Nsymm,Nindep)
    integer,intent(in)                 :: Nsymm
    integer,intent(inout)              :: Nindep
    integer                            :: i,row,col,unit,isymm
    integer,dimension(Nlat)            :: tmp_search
    integer,dimension(Nlat,Nsymm)      :: tmp_map
    integer                            :: i_ind,check_maps
    character(len=5)                   :: tmp_suffix
    integer,dimension(:),allocatable   :: map_lat2ind
    integer,dimension(:,:),allocatable :: map_ind2lat
    integer,allocatable,dimension(:)   :: indep_list
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
    allocate(Ineq_sites_list(Nindep))
    Ineq_sites_list = indep_list
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
    real(8),optional :: ts
    real(8)          :: ts_
    logical          :: symm
    integer          :: i,jj,j,k,row,col,link(4)
    integer          :: unit
    !
    pbc_row_=.false.;pbc_col_=.false.
    if(pbcflag) then
       pbc_row_=.true.
       pbc_col_=.true.
    end if
    if(present(pbc_row)) pbc_row_=pbc_row
    if(present(pbc_col)) pbc_col_=pbc_col
    ts_=0.5d0;if(present(ts))ts_=ts
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
                if(link(jj)>0)H0(i,link(jj))=-ts_ !! ts must be negative.
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



end module ED_WRAP_AUX_FUNX
