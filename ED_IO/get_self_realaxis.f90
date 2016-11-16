!ANOMALous, REAL SELF-ENERGY
subroutine ed_get_self_real_1(SAreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
  SAreal(:,:,:,:,:) = impSAreal(:,:,:,:,:)
end subroutine ed_get_self_real_1

subroutine ed_get_self_real_2(SAreal)
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              SAreal(io,jo,:) = impSAreal(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_real_2

subroutine ed_get_self_real_3(SAreal,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lreal),intent(inout) :: SAreal
  integer                                   :: iorb,jorb,ispin,jspin
  SAreal(:) = impSAreal(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_real_3











subroutine ed_get_self_real_lattice_1(SAreal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
  SAreal(1:Nsites,:,:,:,:,:) = SArealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_self_real_lattice_1

subroutine ed_get_self_real_lattice_2(SAreal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 SAreal(ilat,io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_real_lattice_2

subroutine ed_get_self_real_lattice_3(SAreal,Nsites,ispin,jspin,iorb,jorb)
  integer                                          :: Nsites
  complex(8),dimension(Nsites,Lreal),intent(inout) :: SAreal
  integer                                          :: iorb,jorb,ispin,jspin
  SAreal(1:Nsites,:) = SArealii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_real_lattice_3







subroutine ed_get_self_real_lattice_11(SAreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
  SAreal(:,:,:,:,:) = SArealii(ilat,:,:,:,:,:)
end subroutine ed_get_self_real_lattice_11

subroutine ed_get_self_real_lattice_21(SAreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: SAreal
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              SAreal(io,jo,:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_real_lattice_21

subroutine ed_get_self_real_lattice_31(SAreal,ilat,ispin,jspin,iorb,jorb)
  integer                                   :: ilat
  complex(8),dimension(Lreal),intent(inout) :: SAreal
  integer                                   :: iorb,jorb,ispin,jspin
  SAreal(:) = SArealii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_real_lattice_31
