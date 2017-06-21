!+------------------------------------------------------------------+
!                         SPIN-SPIN
!+------------------------------------------------------------------+
subroutine print_chi_spin
  integer                               :: i,iorb
  integer                               :: unit(3)
  do iorb=1,Norb
     call splot("spinChi_l"//str(iorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(iorb,0:))
     call splot("spinChi_l"//str(iorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,spinChi_w(iorb,:))
     call splot("spinChi_l"//str(iorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(iorb,:))
  enddo
  if(Norb>1)then
     iorb=Norb+1
     call splot("spinChi_tot"//str(iorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(iorb,0:))
     call splot("spinChi_tot"//str(iorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,spinChi_w(iorb,:))
     call splot("spinChi_tot"//str(iorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(iorb,:))
  endif
end subroutine print_chi_spin







!+------------------------------------------------------------------+
!                     DENSITY-DENSITY
!+------------------------------------------------------------------+
subroutine print_chi_dens
  integer                               :: i,j,iorb,jorb
  integer                               :: unit(3),unit_mix
  do iorb=1,Norb
     do jorb=iorb,Norb
        call splot("densChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tau(iorb,jorb,0:))
        call splot("densChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,densChi_w(iorb,jorb,:))
        call splot("densChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_iv(iorb,jorb,:))
     enddo
  enddo
end subroutine print_chi_dens


subroutine print_chi_dens_mix
  integer                               :: i,j,iorb,jorb
  integer                               :: unit(3),unit_mix
  do iorb=1,Norb
     do jorb=1,Norb
        call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_mix_tau(iorb,jorb,0:))
        call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,densChi_mix_w(iorb,jorb,:))
        call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_mix_iv(iorb,jorb,:))
     enddo
  enddo
end subroutine print_chi_dens_mix


subroutine print_chi_dens_tot
  integer                               :: i,j,iorb,jorb
  integer                               :: unit(3),unit_mix
  call splot("densChi_tot_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tot_tau(0:))
  call splot("densChi_tot_realw"//reg(ed_file_suffix)//".ed",wr,densChi_tot_w(:))
  call splot("densChi_tot_iw"//reg(ed_file_suffix)//".ed",vm,densChi_tot_iv(:))
end subroutine print_chi_dens_tot










!+------------------------------------------------------------------+
!                             PAIR
!+------------------------------------------------------------------+
subroutine print_chi_pair
  integer                               :: i,iorb
  integer                               :: unit(3)
  do iorb=1,Norb
     call splot("pairChi_orb"//str(iorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,pairChi_tau(iorb,0:))
     call splot("pairChi_orb"//str(iorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,pairChi_w(iorb,:))
     call splot("pairChi_orb"//str(iorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,pairChi_iv(iorb,:))
  enddo
end subroutine print_chi_pair
