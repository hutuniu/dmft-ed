
!+------------------------------------------------------------------+
!PURPOSE  : Build the Self-energy functions, NONSU2 case
!+------------------------------------------------------------------+
subroutine get_sigma_nonsu2
  integer                                           :: i,j,isign,unit(7),iorb,jorb,ispin,jspin,io,jo
  complex(8)                                        :: fg0
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real
  complex(8),dimension(Nspin*Norb,Nspin*Norb)       :: invGimp,Foo
  character(len=20)                                 :: suffix
  !
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  impG0mats=zero
  impG0real=zero
  invG0mats = zero
  invG0real = zero
  !
  !Get G0^-1
  invG0mats(:,:,:,:,:)=invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
  invG0real(:,:,:,:,:)=invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
  !Get impDelta_anderson
  impDeltamats(:,:,:,:,:)=delta_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
  impDeltareal(:,:,:,:,:)=delta_bath_real(dcmplx(wr(:),eps),dmft_bath)
  !Get inverse functions
  invimpG0mats=invG0mats
  invimpG0real=invG0real
  !
  select case(bath_type)
     !
  case ("normal")
     !
     !Get Gimp^-1 - Matsubara freq.
     do i=1,Lmats
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    if (iorb.eq.jorb) then
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                    endif
                 enddo
              enddo
           enddo
        enddo
        call inv(invGimp)!<--- get [G_{imp}]^-1
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    if (iorb.eq.jorb) then
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !Get Gimp^-1 - Real freq.
     do i=1,Lreal
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    if (iorb.eq.jorb) then
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                    endif
                 enddo
              enddo
           enddo
        enddo
        call inv(invGimp)!<--- get [G_{imp}]^-1
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    if (iorb.eq.jorb) then
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (jspin-1)*Norb
                       impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
  case ("hybrid","replica")
     !
     !Get Gimp^-1 - Matsubara freq.
     do i=1,Lmats
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                 enddo
              enddo
           enddo
        enddo
        call inv(invGimp)!<--- get [G_{imp}]^-1
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                 enddo
              enddo
           enddo
        enddo
     enddo
     !Get Gimp^-1 - Real freq.
     do i=1,Lreal
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                 enddo
              enddo
           enddo
        enddo
        call inv(invGimp)!<--- get [G_{imp}]^-1
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
  end select
  !
  !Get G0and:
  impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
  impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
  !
  !
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine get_sigma_nonsu2
