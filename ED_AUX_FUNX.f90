MODULE ED_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg
  implicit none
  private


  interface set_Hloc
     module procedure set_Hloc_1
     module procedure set_Hloc_2
     module procedure set_Hloc_3d
     module procedure set_Hloc_3c
  end interface set_Hloc

  interface get_Hloc
     module procedure get_Hloc_1
     module procedure get_Hloc_2
  end interface get_Hloc

  interface print_Hloc
     module procedure print_Hloc_2
     module procedure print_Hloc_4
  end interface print_Hloc

  interface select_block
     module procedure select_block_Nlso
     module procedure select_block_NNN
  end interface select_block

  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape

  interface extract_Hloc
     module procedure extract_Hloc_1
     module procedure extract_Hloc_2
  end interface extract_Hloc



  public :: set_Hloc
  public :: get_Hloc  
  public :: print_Hloc
  public :: blocks_to_matrix
  public :: matrix_to_blocks
  public :: select_block
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  public :: extract_Hloc
  public :: stride_index
  public :: get_independent_sites  
  !OBSOLETE (to be removed)
  public :: search_chemical_potential

contains




  !##################################################################
  !##################################################################
  !                   HLOC ROUTINES
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hloc_2(hloc,file)
    character(len=*),optional :: file
    integer                   :: Ni,Nj,iorb,jorb,unit
    complex(8),dimension(:,:) :: hloc
    unit=LOGfile;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    Ni = size(hloc,1)
    Nj = size(hloc,2)
    if(present(file))then
       do iorb=1,Ni
          write(unit,"(90F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(90F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc_2
  !
  subroutine print_Hloc_4(hloc,file)
    character(len=*),optional                   :: file
    integer                                     :: iorb,jorb,ispin,jspin,unit
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    unit=LOGfile;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    if(present(file))then
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
    else
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
    endif
  end subroutine print_Hloc_4






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine set_Hloc_1(hloc,ispin)
    complex(8),dimension(:,:) :: hloc
    integer                   :: ispin
    if(size(hloc,1)/=Norb.OR.size(hloc,2)/=Norb)stop "set_impHloc error: wrong dimensions of Hloc"
    impHloc(ispin,ispin,1:Norb,1:Norb) = Hloc
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")""
       write(LOGfile,"(A)")"Updated impHloc:"
       if(ed_verbose<4)call print_Hloc(impHloc)
    endif
  end subroutine set_Hloc_1
  !
  subroutine set_Hloc_2(hloc)
    complex(8),dimension(:,:,:,:) :: hloc
    if(size(hloc,1)/=Nspin.OR.size(hloc,2)/=Nspin)stop "set_impHloc error: wrong Nspin dimensions of Hloc"
    if(size(hloc,3)/=Norb.OR.size(hloc,4)/=Norb)stop "set_impHloc error: wrong Norb dimensions of Hloc"
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = Hloc
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")""
       write(LOGfile,"(A)")"Updated impHloc:"
       if(ed_verbose<4)call print_Hloc(impHloc)
    endif
  end subroutine set_Hloc_2
  !
  subroutine set_Hloc_3d(hloc)
    real(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")""
       write(LOGfile,"(A)")"Updated impHloc:"
       if(ed_verbose<4)call print_Hloc(impHloc)
    endif
  end subroutine set_Hloc_3d
  !
  subroutine set_Hloc_3c(hloc)
    complex(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")""
       write(LOGfile,"(A)")"Updated impHloc:"
       if(ed_verbose<4)call print_Hloc(impHloc)
    endif
  end subroutine set_Hloc_3c



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine get_Hloc_1(hloc,ispin)
    complex(8),dimension(:,:) :: hloc
    integer                   :: ispin
    if(size(hloc,1)/=Norb.OR.size(hloc,2)/=Norb)stop "set_impHloc error: wrong dimensions of Hloc"
    Hloc = impHloc(ispin,ispin,1:Norb,1:Norb)
  end subroutine get_Hloc_1
  !
  subroutine get_Hloc_2(hloc)
    complex(8),dimension(:,:,:,:) :: hloc
    if(size(hloc,1)/=Nspin.OR.size(hloc,2)/=Nspin)stop "set_impHloc error: wrong Nspin dimensions of Hloc"
    if(size(hloc,3)/=Norb.OR.size(hloc,4)/=Norb)stop "set_impHloc error: wrong Norb dimensions of Hloc"
    Hloc = impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb)
  end subroutine get_Hloc_2



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
    integer                                     :: j,jorb,jspin,js
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
    integer                                     :: i,iorb,ispin,is
    integer                                     :: j,jorb,jspin,js
    Hloc = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i = stride_index([iorb,ispin],[Norb,Nspin])
                j = stride_index([jorb,jspin],[Norb,Nspin])
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc_2










  !##################################################################
  !##################################################################
  !                   TRASFORMATION ROUTINES
  !##################################################################
  !##################################################################

  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! bcast the Blocks vector [Nlat][Nspin*Norb][Nspin*Norb]
  ! into a large matrix [Nlat*Nspin*Norb][Nlat*Nspin*Norb]
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


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! bcast the diagonal part of a large matrix [Nlat*Nspin*Norb][Nlat*Nspin*Norb]
  ! into Blocks vector [Nlat][Nspin*Norb][Nspin*Norb]
  !+-----------------------------------------------------------------------------+!
  function matrix_to_blocks(Matrix) result(Vblocks)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb)      :: Vblocks
    integer                                               :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nspin*Norb
       j = ip*Nspin*Norb
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! select a single block of the diagonal from a large matrix.
  !   + _Nlso = matrix has dimensions Nlso*Nlso
  !   + _NNN  = matrix has dimensions Nlat,Nspin,Nspin,Norb,Norb
  !+-----------------------------------------------------------------------------+!
  function select_block_Nlso(ip,Matrix) result(Vblock)
    integer                                               :: ip
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Vblock
    integer                                               :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nspin*Norb
    j =       ip*Nspin*Norb
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso




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
  ! find all inequivalent sites with respect the user defined symmetry operations 
  ! Build and check maps from the full(independent) lattice to the independent 
  ! (full) lattice                        !
  !+-----------------------------------------------------------------------------+!
  subroutine get_independent_sites(symmetry_operations,Nsymm,Nindep)
    integer,intent(in)                 :: Nsymm
    integer,intent(inout)              :: Nindep
    integer                            :: i,unit,isymm
    integer,dimension(Nlat)            :: tmp_search
    integer,dimension(Nlat,Nsymm)      :: tmp_map
    integer                            :: i_ind,check_maps
    character(len=5)                   :: tmp_suffix
    ! integer,dimension(:),allocatable   :: map_lat2ind
    ! integer,dimension(:,:),allocatable :: map_ind2lat
    !integer,allocatable,dimension(:)   :: indep_list
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
    write(*,*) Nlat
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
    write(*,*) "NINDEP",Nindep
    write(*,*) indep_list
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
    ! allocate(Ineq_sites_list(Nindep))
    ! Ineq_sites_list = indep_list
  end subroutine get_independent_sites









  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################

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
