!ANOMALous, MATSUBARA SELF-ENERGY
subroutine ed_get_self_matsubara_1(SAmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
  SAmats(:,:,:,:,:) = impSAmats(:,:,:,:,:)
end subroutine ed_get_self_matsubara_1

subroutine ed_get_self_matsubara_2(SAmats)
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              SAmats(io,jo,:) = impSAmats(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_matsubara_2

subroutine ed_get_self_matsubara_3(SAmats,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lmats),intent(inout) :: SAmats
  integer                                   :: iorb,jorb,ispin,jspin
  SAmats(:) = impSAmats(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_matsubara_3





subroutine ed_get_self_matsubara_lattice_1(SAmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
  SAmats(1:Nsites,:,:,:,:,:) = SAmatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_self_matsubara_lattice_1

subroutine ed_get_self_matsubara_lattice_2(SAmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 SAmats(ilat,io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_matsubara_lattice_2

subroutine ed_get_self_matsubara_lattice_3(SAmats,Nsites,ispin,jspin,iorb,jorb)
  integer                                          :: Nsites
  complex(8),dimension(Nsites,Lmats),intent(inout) :: SAmats
  integer                                          :: iorb,jorb,ispin,jspin
  SAmats(1:Nsites,:) = SAmatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_matsubara_lattice_3






subroutine ed_get_self_matsubara_lattice_11(SAmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
  SAmats(:,:,:,:,:) = SAmatsii(ilat,:,:,:,:,:)
end subroutine ed_get_self_matsubara_lattice_11

subroutine ed_get_self_matsubara_lattice_21(SAmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: SAmats
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              SAmats(io,jo,:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_self_matsubara_lattice_21

subroutine ed_get_self_matsubara_lattice_31(SAmats,ilat,ispin,jspin,iorb,jorb)
  integer                                   :: ilat
  complex(8),dimension(Lmats),intent(inout) :: SAmats
  integer                                   :: iorb,jorb,ispin,jspin
  SAmats(:) = SAmatsii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_self_matsubara_lattice_31
