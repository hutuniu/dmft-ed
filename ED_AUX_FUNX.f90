MODULE ED_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_LINALG
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
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

!  interface orbital_Lz_rotation
!     module procedure orbital_Lz_rotation_Norb
!     module procedure orbital_Lz_rotation_NorbNspin
!  end interface orbital_Lz_rotation


  public :: set_Hloc
  public :: get_Hloc  
  public :: print_Hloc
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  public :: so2os_reshape
  public :: os2so_reshape
  public :: stride_index
  !OBSOLETE (to be removed)
  public :: search_chemical_potential
  public :: search_chempot
  public :: SOC_jz_symmetrize
  public :: atomic_SOC
  public :: atomic_SOC_rotation
  public :: orbital_Lz_rotation_NorbNspin
  public :: orbital_Lz_rotation_Norb
  public :: atomic_j
  public :: tql2
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
          write(unit,"(2000F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(2000F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(2000(A1,F7.3,A1,F7.3,A1,2x))")&
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
             write(unit,"(2000F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(2000F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    else
       do ispin=1,Nspin
          do iorb=1,Norb
             write(LOGfile,"(2000(A1,F7.3,A1,F7.3,A1,2x))")&
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
    write(LOGfile,"(A)")""
    write(LOGfile,"(A)")"Updated impHloc:"
    if(ed_verbose<4)call print_Hloc(impHloc)
  end subroutine set_Hloc_1
  !
  subroutine set_Hloc_2(hloc)
    complex(8),dimension(:,:,:,:) :: hloc
    if(size(hloc,1)/=Nspin.OR.size(hloc,2)/=Nspin)stop "set_impHloc error: wrong Nspin dimensions of Hloc"
    if(size(hloc,3)/=Norb.OR.size(hloc,4)/=Norb)stop "set_impHloc error: wrong Norb dimensions of Hloc"
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = Hloc
    write(LOGfile,"(A)")""
    write(LOGfile,"(A)")"Updated impHloc:"
    if(ed_verbose<4)call print_Hloc(impHloc)
  end subroutine set_Hloc_2
  !
  subroutine set_Hloc_3d(hloc)
    real(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    write(LOGfile,"(A)")""
    write(LOGfile,"(A)")"Updated impHloc:"
    if(ed_verbose<4)call print_Hloc(impHloc)
  end subroutine set_Hloc_3d
  !
  subroutine set_Hloc_3c(hloc)
    complex(8) :: hloc
    impHloc(1:Nspin,1:Nspin,1:Norb,1:Norb) = hloc
    write(LOGfile,"(A)")""
    write(LOGfile,"(A)")"Updated impHloc:"
    if(ed_verbose<4)call print_Hloc(impHloc)
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
    integer                                            :: jorb,jspin,js
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
    integer                                               :: jorb,jspin,js
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
    integer                                            :: jorb,jspin,js
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
    integer                                               :: jorb,jspin,js
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


  function so2os_reshape(fg,Nspin,Norb) result(g)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: iorb,jorb,ispin,jspin
    integer                                     :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = iorb + (ispin-1)*Norb
                jo1 = jorb + (jspin-1)*Norb
                !I-index
                io2 = ispin + (iorb-1)*Nspin
                jo2 = jspin + (jorb-1)*Nspin
                !switch
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function so2os_reshape

  function os2so_reshape(fg,Nspin,Norb) result(g)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: iorb,jorb,ispin,jspin
    integer                                     :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = ispin + (iorb-1)*Nspin
                jo1 = jspin + (jorb-1)*Nspin
                !I-index
                io2 = iorb + (ispin-1)*Norb
                jo2 = jorb + (jspin-1)*Norb
                !switch
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function os2so_reshape




  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! - stride_index: given an array of indices ivec=[i_1,...,i_L] and the corresponding 
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
    write(LOGfile,"(A,L2)")"Converged=",converged
    print*,""
    !
  end subroutine search_chemical_potential







  subroutine search_chempot(xmu_tmp,dens_tmp,converged_)!,write_6)
    real(8),intent(in)          ::   dens_tmp
    real(8),intent(inout)       ::   xmu_tmp
    logical,intent(inout)       ::   converged_
    !logical,intent(in)          ::   write_6
    !internal
    real(8)                     ::   diffdens,delta_xmu,xmu_shift
    real(8)                     ::   denslarge,denssmall,xmu_old
    real(8),save                ::   diffdens_old
    real(8),save                ::   xmularge
    real(8),save                ::   xmusmall
    integer,save                ::   ilarge=0
    integer,save                ::   ismall=0
    integer,save                ::   inotbound
    integer,save                ::   ibound=0
    integer,save                ::   iattempt=1
    integer                     ::   unit
    !
    diffdens=dens_tmp-nread
    delta_xmu=0.2d0
    !
    if ((dabs(diffdens)).le.nerr) then
       converged_=.TRUE.
       inotbound=0
       !if(write_6)then
          write(LOGfile,*)
          write(LOGfile,*) "   ------------------- search chempot -----------------"
          write(LOGfile,'(A30,I3)')    "   Density ok in attempt: ",iattempt
          write(LOGfile,'(A30,F10.6)') "   tolerance: ",nerr
          write(LOGfile,'(A30,F10.6)') "   density: ",dens_tmp
          write(LOGfile,'(A30,F10.6)') "   error: ",diffdens
          write(LOGfile,'(A30,F10.6)') "   target desity: ",nread
          write(LOGfile,'(A30,F10.6)') "   xmu: ",xmu_tmp
          write(LOGfile,"(A30,L3)")    "   Converged(n): ",converged_
          write(LOGfile,*) "   ----------------------------------------------------"
          write(LOGfile,*)
       !endif
       unit=free_unit()
       open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
       write(unit,'(3F25.12)')xmu_tmp,dens_tmp,diffdens
       close(unit)
    else
       converged_=.FALSE.
       !if(write_6)then
          write(LOGfile,*)
          write(LOGfile,*) "   ------------------- search chempot -----------------"
          write(LOGfile,'(A30,2I5)')    "   Adjusting xmu #",iattempt,inotbound
          write(LOGfile,'(A10,F10.6,A8,F10.6,A8,F10.6)') "    n:",dens_tmp,"!= n:",nread,"error:",abs(dens_tmp-nread)
       !endif
       !vedo se la densità è troppa o troppo poca
       if (diffdens.gt.0.d0) then  
  !        ilarge=1
          xmularge=xmu_tmp
          denslarge=dens_tmp
       elseif (diffdens.lt.0.d0) then
  !        ismall=1
          xmusmall=xmu_tmp
          denssmall=dens_tmp
       endif
       if (ilarge*ismall.ne.0) ibound=ibound+1
       if (ilarge*ismall.eq.0) then
          !non ho ancora trovato un xmu per cui diffdens cambia segno
          inotbound=inotbound+1
          if (inotbound>=8)  delta_xmu = 0.5d0
          if (inotbound>=12) delta_xmu = 1.0d0
          if (inotbound>=16) delta_xmu = 1.5d0
          if (inotbound>=20) delta_xmu = 2.0d0
          xmu_shift = delta_xmu * diffdens
          xmu_old = xmu_tmp
          xmu_tmp = xmu_tmp - xmu_shift
          !
          write(LOGfile,*) "   Delta xmu: ",delta_xmu
          write(LOGfile,*) "   Old xmu: ",xmu_old
          write(LOGfile,*) "   Try xmu: ",xmu_tmp
          write(LOGfile,*) "   ----------------------------------------------------"
          write(LOGfile,*)
          !
       else
       !elseif((ilarge*ismall.ne.0).and.(ibound.ne.1))then
          !ho trovato un xmu per cui diffdens cambia segno
          write(LOGfile,*)"   xmu is bound",xmularge,"-",xmusmall
          xmu_shift =  sign(1.0d0,diffdens)*abs((xmusmall-xmularge)/2.)
          xmu_old = xmu_tmp
          xmu_tmp = xmu_tmp - xmu_shift
          !
          write(LOGfile,*) "   Old xmu: ",xmu_old
          write(LOGfile,*) "   Try xmu =",xmu_tmp
          write(LOGfile,*) "   ----------------------------------------------------"
          write(LOGfile,*)
          !
       endif
       !
       unit=free_unit()
       open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
       write(unit,'(3F25.12,1I5)')xmu_tmp,dens_tmp,diffdens,iattempt
       close(unit)
       !
    endif
    iattempt=iattempt+1
    diffdens_old=diffdens
    !
  end subroutine search_chempot



  subroutine SOC_jz_symmetrize_old(funct)
    !passed
    complex(8),allocatable,intent(inout)         ::  funct(:,:,:,:,:)
    complex(8),allocatable                       ::  symmetrized_funct(:,:,:,:,:)
    complex(8),allocatable                       ::  a_funct(:),b_funct(:),c_funct(:),d_funct(:),e_funct(:),f_funct(:)
    integer                                      ::  ispin,iorb
    integer                                      ::  ifreq,Lfreq
    if(size(funct,dim=1)/=Nspin)stop "wrong size 1 in SOC symmetrize input f"
    if(size(funct,dim=2)/=Nspin)stop "wrong size 2 in SOC symmetrize input f"
    if(size(funct,dim=3)/=Norb) stop "wrong size 3 in SOC symmetrize input f"
    if(size(funct,dim=4)/=Norb) stop "wrong size 4 in SOC symmetrize input f"
    Lfreq=size(funct,dim=5)
    allocate(symmetrized_funct(Nspin,Nspin,Norb,Norb,Lfreq));symmetrized_funct=zero
    allocate(a_funct(Lfreq));a_funct=zero
    allocate(b_funct(Lfreq));b_funct=zero
    allocate(c_funct(Lfreq));c_funct=zero
    allocate(d_funct(Lfreq));d_funct=zero
    allocate(e_funct(Lfreq));e_funct=zero
    allocate(f_funct(Lfreq));f_funct=zero
    !
    !diagonal
    do ispin=1,Nspin
       do iorb=1,Norb
          a_funct = a_funct + funct(ispin,ispin,iorb,iorb,:)
       enddo
    enddo
    a_funct = a_funct / ( Nspin * Norb )
    !
    !cerchi bianchi
    b_funct = ( funct(1,2,1,3,:) + funct(2,1,3,1,:) )/2.d0
    !cerchi neri
    c_funct = ( funct(1,2,3,1,:) + funct(2,1,1,3,:) )/2.d0
    !quadri bianchi
    e_funct = ( funct(1,1,2,1,:) + funct(1,2,3,2,:) + funct(2,2,1,2,:) + funct(2,1,3,2,:) )/4.d0
    !quadri neri
    d_funct = ( funct(1,1,1,2,:) + funct(1,2,2,3,:) + funct(2,2,2,1,:) + funct(2,1,2,3,:) )/4.d0
    !media
    f_funct = ( b_funct - c_funct - xi*e_funct + xi*d_funct )/4.d0
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          symmetrized_funct(ispin,ispin,iorb,iorb,:) = a_funct
       enddo
    enddo
    !
    do ifreq=1,Lfreq
       symmetrized_funct(:,:,:,:,ifreq)=symmetrized_funct(:,:,:,:,ifreq)-2.d0*f_funct(ifreq)*so2nn_reshape(atomic_SOC(),Nspin,Norb)
    enddo
    !
    funct = zero
    funct = symmetrized_funct
    !
    deallocate(symmetrized_funct)
  end subroutine SOC_jz_symmetrize_old




  subroutine SOC_jz_symmetrize(funct)
    !passed
    complex(8),allocatable,intent(inout)         ::  funct(:,:,:,:,:)
    complex(8),allocatable                       ::  funct_in(:,:,:),funct_out(:,:,:)
    complex(8),allocatable                       ::  a_funct(:),b_funct(:)
    integer                                      ::  ispin,io
    integer                                      ::  ifreq,Lfreq
    complex(8),allocatable                       ::  U(:,:),Udag(:,:)
    if(size(funct,dim=1)/=Nspin)stop "wrong size 1 in SOC symmetrize input f"
    if(size(funct,dim=2)/=Nspin)stop "wrong size 2 in SOC symmetrize input f"
    if(size(funct,dim=3)/=Norb) stop "wrong size 3 in SOC symmetrize input f"
    if(size(funct,dim=4)/=Norb) stop "wrong size 4 in SOC symmetrize input f"
    Lfreq=size(funct,dim=5)
    allocate(funct_in(Nspin*Norb,Nspin*Norb,Lfreq)); funct_in=zero
    allocate(funct_out(Nspin*Norb,Nspin*Norb,Lfreq));funct_out=zero
    allocate(U(Nspin*Norb,Nspin*Norb));U=zero
    allocate(Udag(Nspin*Norb,Nspin*Norb));Udag=zero
    allocate(a_funct(Lfreq));a_funct=zero
    allocate(b_funct(Lfreq));b_funct=zero
    !
    !function intake
    do ifreq=1,Lfreq
       funct_in(:,:,ifreq)=nn2so_reshape(funct(:,:,:,:,ifreq),Nspin,Norb)
    enddo
    !
    !function diagonalization
    if(Jz_basis)then
       U=matmul(transpose(conjg(orbital_Lz_rotation_NorbNspin())),atomic_SOC_rotation())
       Udag=transpose(conjg(U))
    else
       U=atomic_SOC_rotation()
       Udag=transpose(conjg(U))
    endif
    !
    do ifreq=1,Lfreq
       funct_out(:,:,ifreq)=matmul(Udag,matmul(funct_in(:,:,ifreq),U))
    enddo
    !
    !function symmetrization in the rotated basis
    do io=1,2
       a_funct(:)=a_funct(:)+funct_out(io,io,:)
    enddo
    a_funct = a_funct/2.d0
    do io=3,6
       b_funct(:)=b_funct(:)+funct_out(io,io,:)
    enddo
    b_funct = b_funct/4.d0
    funct_out=zero
    do io=1,2
       funct_out(io,io,:)=a_funct(:)
    enddo
    do io=3,6
       funct_out(io,io,:)=b_funct(:)
    enddo
    !
    !function rotation in the non-diagonal basis
    funct_in=zero
    do ifreq=1,Lfreq
       funct_in(:,:,ifreq)=matmul(U,matmul(funct_out(:,:,ifreq),Udag))
    enddo
    !
    !founction out
    funct=zero
    do ifreq=1,Lfreq
       funct(:,:,:,:,ifreq)=so2nn_reshape(funct_in(:,:,ifreq),Nspin,Norb)
    enddo
  end subroutine SOC_jz_symmetrize




  !+-------------------------------------------------------------------+
  !PURPOSE  : Atomic SOC and j vector components
  !+-------------------------------------------------------------------+
  function atomic_SOC() result (LS)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: LS,LS_
    integer                                      :: i,j
    LS_=zero;LS=zero
    LS_(1:2,3:4) = +Xi * pauli_z / 2.
    LS_(1:2,5:6) = -Xi * pauli_y / 2.
    LS_(3:4,5:6) = +Xi * pauli_x / 2.
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          LS_(j,i)=conjg(LS_(i,j))
       enddo
    enddo
    LS=so2os_reshape(LS_,Nspin,Norb)
  end function atomic_SOC

  function atomic_SOC_rotation() result (LS_rot)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: LS_rot,LS_rot_
    integer                                      :: i,j
    LS_rot_=zero;LS_rot=zero
    !
    ! {a,Sz}-->{J}
    !
    ![Norb*Norb]*Nspin notation
    !J=1/2 jz=-1/2
    LS_rot_(1,1)=+1.d0
    LS_rot_(2,1)=-Xi
    LS_rot_(6,1)=-1.d0
    LS_rot_(:,1)=LS_rot_(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    LS_rot_(4,2)=+1.d0
    LS_rot_(5,2)=+Xi
    LS_rot_(3,2)=+1.d0
    LS_rot_(:,2)=LS_rot_(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    LS_rot_(4,3)=+1.d0
    LS_rot_(5,3)=-Xi
    LS_rot_(:,3)=LS_rot_(:,3)/sqrt(2.)
    !J=3/2 jz=+3/2
    LS_rot_(1,4)=-1.d0
    LS_rot_(2,4)=-Xi
    LS_rot_(:,4)=LS_rot_(:,4)/sqrt(2.)
    !J=3/2 jz=-1/2
    LS_rot_(1,5)=+1.d0
    LS_rot_(2,5)=-Xi
    LS_rot_(6,5)=+2.d0
    LS_rot_(:,5)=LS_rot_(:,5)/sqrt(6.)
    !J=3/2 jz=+1/2
    LS_rot_(4,6)=-1.d0
    LS_rot_(5,6)=-Xi
    LS_rot_(3,6)=+2.d0
    LS_rot_(:,6)=LS_rot_(:,6)/sqrt(6.)
    !
    LS_rot=LS_rot_
    !
  end function atomic_SOC_rotation

  function orbital_Lz_rotation_Norb() result (U_rot)
    complex(8),dimension(Norb,Norb)              :: U_rot,U_rot_
    integer                                      :: i,j
    U_rot=zero;U_rot_=zero
    !
    ! {a}-->{Lz}
    !
    ![Norb*Norb] notation
    U_rot_(1,1)=-Xi/sqrt(2.)
    U_rot_(2,2)=+1.d0/sqrt(2.)
    U_rot_(3,3)=+Xi
    U_rot_(1,2)=-Xi/sqrt(2.)
    U_rot_(2,1)=-1.d0/sqrt(2.)
    !
    U_rot=U_rot_
    !
  end function orbital_Lz_rotation_Norb

  function orbital_Lz_rotation_NorbNspin() result (U_rot)
    complex(8),dimension(Norb,Norb)              :: U_rot_
    complex(8),dimension(Norb*Nspin,Norb*Nspin)  :: U_rot
    integer                                      :: i,j
    U_rot=zero;U_rot_=zero
    !
    ! {a,Sz}-->{Lz,Sz}
    !
    ![Norb*Norb]*Nspin notation
    U_rot_(1,1)=-Xi/sqrt(2.)
    U_rot_(2,2)=+1.d0/sqrt(2.)
    U_rot_(3,3)=+Xi
    U_rot_(1,2)=-Xi/sqrt(2.)
    U_rot_(2,1)=-1.d0/sqrt(2.)
    !
    U_rot(1:Norb,1:Norb)=U_rot_
    U_rot(1+Norb:2*Norb,1+Norb:2*Norb)=U_rot_
    !
  end function orbital_Lz_rotation_NorbNspin

  function atomic_j(component) result (ja)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: ja,ja_
    character(len=1)                             :: component
    integer                                      :: i,j
    ja_=zero;ja=zero
    if    (component=="x")then
       ja_(1:2,1:2) = pauli_x / 2.
       ja_(3:4,3:4) = pauli_x / 2.
       ja_(5:6,5:6) = pauli_x / 2.
       ja_(3:4,5:6) = -Xi * eye(2)
    elseif(component=="y")then
       ja_(1:2,1:2) = pauli_y / 2.
       ja_(3:4,3:4) = pauli_y / 2.
       ja_(5:6,5:6) = pauli_y / 2.
       ja_(1:2,5:6) = +Xi * eye(2)
    elseif(component=="z")then
       ja_(1:2,1:2) = pauli_z / 2.
       ja_(3:4,3:4) = pauli_z / 2.
       ja_(5:6,5:6) = pauli_z / 2.
       ja_(1:2,3:4) = -Xi * eye(2)
    endif
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          ja_(j,i)=conjg(ja_(i,j))
       enddo
    enddo
    ja=so2os_reshape(ja_,Nspin,Norb)
  end function atomic_j




  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
  !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
  !    tridiagonal matrix by the QL method.  The eigenvectors of a full
  !    symmetric matrix can also be found if TRED2 has been used to reduce this
  !    full matrix to tridiagonal form.
  !  Parameters:
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
  !    the matrix.  On output, the eigenvalues in ascending order.  If an error
  !    exit is made, the eigenvalues are correct but unordered for indices
  !    1,2,...,IERR-1.
  !
  !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
  !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
  !    On output, E has been destroyed.
  !
  !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
  !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
  !    the tridiagonal matrix are desired, Z must contain the identity matrix.
  !    On output, Z contains the orthonormal eigenvectors of the symmetric
  !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
  !    the eigenvectors associated with the stored eigenvalues.
  !
  !    Output, integer ( kind = 4 ) IERR, error flag.
  !    0, normal return,
  !    J, if the J-th eigenvalue has not been determined after
  !    30 iterations.
  !
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)
       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap


END MODULE ED_AUX_FUNX
