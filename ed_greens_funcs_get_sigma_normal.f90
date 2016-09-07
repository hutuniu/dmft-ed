subroutine get_sigma_normal
  integer                                           :: i,ispin,iorb
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
  complex(8),dimension(Norb,Norb)                   :: invGimp
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  invG0mats = zero
  invGmats  = zero
  invG0real = zero
  invGreal  = zero
  !
  !Get G0^-1
  invG0mats(:,:,:,:,:) = invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
  invG0real(:,:,:,:,:) = invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
  !
  select case(bath_type)
  case default                !Diagonal in both spin and orbital
     !
     !Get Gimp^-1
     do ispin=1,Nspin
        do iorb=1,Norb
           invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
           invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
        enddo
     enddo
     !Get Sigma functions: Sigma= G0^-1 - G^-1
     impSmats=zero
     impSreal=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
           impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
        enddo
     enddo
     !
  case ("hybrid","replica")   !Diagonal in spin only. Full Orbital structure
     !
     !Get Gimp^-1
     do ispin=1,Nspin
        do i=1,Lmats
           invGimp = impGmats(ispin,ispin,:,:,i)
           call inv(invGimp)
           invGmats(ispin,ispin,:,:,i)=invGimp
        enddo
        !
        do i=1,Lreal
           invGimp = impGreal(ispin,ispin,:,:,i)
           call inv(invGimp)
           invGreal(ispin,ispin,:,:,i)=invGimp
        enddo
     enddo
     !Get Sigma functions: Sigma= G0^-1 - G^-1
     impSmats=zero
     impSreal=zero
     do ispin=1,Nspin
        impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
        !
        impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
     enddo
     !
  end select
  !
  !Get G0and:
  impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
  impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
  !!
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine get_sigma_normal
