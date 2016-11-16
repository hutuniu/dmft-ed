!ANOMALous, REAL GREEN'S FUNCTION
subroutine ed_get_fimp_real_1(Freal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  Freal(:,:,:,:,:) = impFreal(:,:,:,:,:)
end subroutine ed_get_fimp_real_1

subroutine ed_get_fimp_real_2(Freal)
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Freal(io,jo,:) = impFreal(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_real_2

subroutine ed_get_fimp_real_3(Freal,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lreal),intent(inout) :: Freal
  integer                                   :: iorb,jorb,ispin,jspin
  Freal(:) = impFreal(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_real_3







subroutine ed_get_fimp_real_lattice_1(Freal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  Freal(1:Nsites,:,:,:,:,:) = Frealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_fimp_real_lattice_1

subroutine ed_get_fimp_real_lattice_2(Freal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Freal(ilat,io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_real_lattice_2

subroutine ed_get_fimp_real_lattice_3(Freal,Nsites,ispin,jspin,iorb,jorb)
  integer                                          :: Nsites
  complex(8),dimension(Nsites,Lreal),intent(inout) :: Freal
  integer                                          :: iorb,jorb,ispin,jspin
  Freal(1:Nsites,:) = Frealii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_real_lattice_3





subroutine ed_get_fimp_real_lattice_11(Freal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  Freal(:,:,:,:,:) = Frealii(ilat,:,:,:,:,:)
end subroutine ed_get_fimp_real_lattice_11

subroutine ed_get_fimp_real_lattice_21(Freal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Freal
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Freal(io,jo,:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_real_lattice_21

subroutine ed_get_fimp_real_lattice_31(Freal,ilat,ispin,jspin,iorb,jorb)
  integer                                   :: ilat
  complex(8),dimension(Lreal),intent(inout) :: Freal
  integer                                   :: iorb,jorb,ispin,jspin
  Freal(:) = Frealii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_real_lattice_31
