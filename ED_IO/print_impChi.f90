!+------------------------------------------------------------------+
!                         SPIN-SPIN
!+------------------------------------------------------------------+
subroutine print_chi_spin
  integer                               :: i,iorb
  integer                               :: unit(3)
  do iorb=1,Norb
     unit=free_units(3)
     open(unit(1),file="spinChi_orb"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
     open(unit(2),file="spinChi_orb"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
     open(unit(3),file="spinChi_orb"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
     do i=0,Ltau
        write(unit(1),*)tau(i),spinChi_tau(iorb,i)
     enddo
     do i=1,Lreal
        if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(spinChi_w(iorb,i)),dreal(spinChi_w(iorb,i))
     enddo
     do i=0,Lmats
        write(unit(3),*)vm(i),dimag(spinChi_iv(iorb,i)),dreal(spinChi_iv(iorb,i))
     enddo
     do i=1,3
        close(unit(i))
     enddo
  enddo
  if(Norb>1)then
     iorb=Norb+1
     unit=free_units(3)
     open(unit(1),file="spinChi_tot"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
     open(unit(2),file="spinChi_tot"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
     open(unit(3),file="spinChi_tot"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
     do i=0,Ltau
        write(unit(1),*)tau(i),spinChi_tau(iorb,i)
     enddo
     do i=1,Lreal
        if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(spinChi_w(iorb,i)),dreal(spinChi_w(iorb,i))
     enddo
     do i=0,Lmats
        write(unit(3),*)vm(i),dimag(spinChi_iv(iorb,i)),dreal(spinChi_iv(iorb,i))
     enddo
     do i=1,3
        close(unit(i))
     enddo
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
        unit=free_units(3)
        open(unit(1),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_tau"//reg(ed_file_suffix)//".ed")
        open(unit(2),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw"//reg(ed_file_suffix)//".ed")
        open(unit(3),file="densChi_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw"//reg(ed_file_suffix)//".ed")
        do i=0,Ltau
           write(unit(1),*)tau(i),densChi_tau(iorb,jorb,i)
        enddo
        do i=1,Lreal
           if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_w(iorb,jorb,i)),dreal(densChi_w(iorb,jorb,i))
        enddo
        do i=0,Lmats
           write(unit(3),*)vm(i),dimag(densChi_iv(iorb,jorb,i)),dreal(densChi_iv(iorb,jorb,i))
        enddo
        do i=1,3
           close(unit(i))
        enddo
     enddo
  enddo
end subroutine print_chi_dens


subroutine print_chi_dens_mix
  integer                               :: i,j,iorb,jorb
  integer                               :: unit(3),unit_mix
  do iorb=1,Norb
     do jorb=1,Norb
        unit=free_units(3)
        open(unit(1),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_tau"//reg(ed_file_suffix)//".ed")
        open(unit(2),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw"//reg(ed_file_suffix)//".ed")
        open(unit(3),file="densChi_mix_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw"//reg(ed_file_suffix)//".ed")
        do i=0,Ltau
           write(unit(1),*)tau(i),densChi_mix_tau(iorb,jorb,i)
        enddo
        do i=1,Lreal
           if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_mix_w(iorb,jorb,i)),dreal(densChi_mix_w(iorb,jorb,i))
        enddo
        do i=0,Lmats
           write(unit(3),*)vm(i),dimag(densChi_mix_iv(iorb,jorb,i)),dreal(densChi_mix_iv(iorb,jorb,i))
        enddo
        do i=1,3
           close(unit(i))
        enddo
     enddo
  enddo
end subroutine print_chi_dens_mix


subroutine print_chi_dens_tot
  integer                               :: i,j,iorb,jorb
  integer                               :: unit(3),unit_mix
  unit=free_units(3)
  open(unit(1),file="densChi_tot_tau"//reg(ed_file_suffix)//".ed")
  open(unit(2),file="densChi_tot_realw"//reg(ed_file_suffix)//".ed")
  open(unit(3),file="densChi_tot_iw"//reg(ed_file_suffix)//".ed")
  do i=0,Ltau
     write(unit(1),*)tau(i),densChi_tot_tau(i)
  enddo
  do i=1,Lreal
     if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(densChi_tot_w(i)),dreal(densChi_tot_w(i))
  enddo
  do i=0,Lmats
     write(unit(3),*)vm(i),dimag(densChi_tot_iv(i)),dreal(densChi_tot_iv(i))
  enddo
  do i=1,3
     close(unit(i))
  enddo
end subroutine print_chi_dens_tot










!+------------------------------------------------------------------+
!                             PAIR
!+------------------------------------------------------------------+
subroutine print_chi_pair
  integer                               :: i,iorb
  integer                               :: unit(3)
  do iorb=1,Norb
     unit=free_units(3)
     open(unit(1),file="pairChi_orb"//reg(txtfy(iorb))//"_tau"//reg(ed_file_suffix)//".ed")
     open(unit(2),file="pairChi_orb"//reg(txtfy(iorb))//"_realw"//reg(ed_file_suffix)//".ed")
     open(unit(3),file="pairChi_orb"//reg(txtfy(iorb))//"_iw"//reg(ed_file_suffix)//".ed")
     do i=0,Ltau
        write(unit(1),*)tau(i),pairChi_tau(iorb,i)
     enddo
     do i=1,Lreal
        if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(pairChi_w(iorb,i)),dreal(pairChi_w(iorb,i))
     enddo
     do i=0,Lmats
        write(unit(3),*)vm(i),dimag(pairChi_iv(iorb,i)),dreal(pairChi_iv(iorb,i))
     enddo
     do i=1,3
        close(unit(i))
     enddo
  enddo
end subroutine print_chi_pair
