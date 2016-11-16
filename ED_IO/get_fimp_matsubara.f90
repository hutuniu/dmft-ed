!ANOMALous, MATSUBARA GREEN'S FUNCTION
subroutine ed_get_fimp_matsubara_1(Fmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  Fmats(:,:,:,:,:) = impFmats(:,:,:,:,:)
end subroutine ed_get_fimp_matsubara_1

subroutine ed_get_fimp_matsubara_2(Fmats)
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
  integer  :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Fmats(io,jo,:) = impFmats(ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_matsubara_2

subroutine ed_get_fimp_matsubara_3(Fmats,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lmats),intent(inout) :: Fmats
  integer                                   :: iorb,jorb,ispin,jspin
  Fmats(:) = impFmats(ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_matsubara_3







subroutine ed_get_fimp_matsubara_lattice_1(Fmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  Fmats(1:Nsites,:,:,:,:,:) = Fmatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_fimp_matsubara_lattice_1

subroutine ed_get_fimp_matsubara_lattice_2(Fmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  do ilat=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Fmats(ilat,io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_matsubara_lattice_2

subroutine ed_get_fimp_matsubara_lattice_3(Fmats,Nsites,ispin,jspin,iorb,jorb)
  integer                                          :: Nsites
  complex(8),dimension(Nsites,Lmats),intent(inout) :: Fmats
  integer                                          :: iorb,jorb,ispin,jspin
  Fmats(1:Nsites,:) = Fmatsii(1:Nsites,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_matsubara_lattice_3






subroutine ed_get_fimp_matsubara_lattice_11(Fmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  Fmats(:,:,:,:,:) = Fmatsii(ilat,:,:,:,:,:)
end subroutine ed_get_fimp_matsubara_lattice_11

subroutine ed_get_fimp_matsubara_lattice_21(Fmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Fmats
  integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Fmats(io,jo,:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_fimp_matsubara_lattice_21

subroutine ed_get_fimp_matsubara_lattice_31(Fmats,ilat,ispin,jspin,iorb,jorb)
  integer                                   :: ilat
  complex(8),dimension(Lmats),intent(inout) :: Fmats
  integer                                   :: iorb,jorb,ispin,jspin
  Fmats(:) = Fmatsii(ilat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_fimp_matsubara_lattice_31
