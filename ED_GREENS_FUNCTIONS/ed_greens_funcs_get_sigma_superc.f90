!+------------------------------------------------------------------+
!PURPOSE  : Build the Self-energy functions, SUPERC case
!+------------------------------------------------------------------+
subroutine get_sigma_superc
  integer                                               :: i,ispin,iorb
  real(8)                                               :: det_mats(Lmats)
  complex(8)                                            :: det_real(Lreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)     :: invG0mats,invF0mats,invGmats,invFmats
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)     :: invG0real,invF0real,invGreal,invFreal
  complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)       :: invGimp
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  invG0mats = zero
  invF0mats = zero
  invGmats  = zero
  invFmats  = zero
  invG0real = zero
  invF0real = zero
  invGreal  = zero
  invFreal  = zero
  !
  !Get G0^-1,F0^-1
  ispin=1
  invG0mats(ispin,ispin,:,:,:) = invg0_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
  invF0mats(ispin,ispin,:,:,:) = invf0_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
  !
  invG0real(ispin,ispin,:,:,:) = invg0_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
  invF0real(ispin,ispin,:,:,:) = invf0_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
  !
  select case(bath_type)
  case default
     !      
     !Get Gimp^-1
     do iorb=1,Norb
        det_mats  =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
        invGmats(ispin,ispin,iorb,iorb,:) = conjg(impGmats(ispin,ispin,iorb,iorb,:))/det_mats
        invFmats(ispin,ispin,iorb,iorb,:) = impFmats(ispin,ispin,iorb,iorb,:)/det_mats
        !
        det_real  = -impGreal(ispin,ispin,iorb,iorb,:)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1)) - impFreal(ispin,ispin,iorb,iorb,:)**2
        invGreal(ispin,ispin,iorb,iorb,:) =  -conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1))/det_real(:)
        invFreal(ispin,ispin,iorb,iorb,:) =  -impFreal(ispin,ispin,iorb,iorb,:)/det_real(:)
     enddo
     !Get Sigma functions: Sigma= G0^-1 - G^-1
     impSmats=zero
     impSAmats=zero
     impSreal=zero
     impSAreal=zero
     do iorb=1,Norb
        impSmats(ispin,ispin,iorb,iorb,:)  = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
        impSAmats(ispin,ispin,iorb,iorb,:) = invF0mats(ispin,ispin,iorb,iorb,:) - invFmats(ispin,ispin,iorb,iorb,:)
        !
        impSreal(ispin,ispin,iorb,iorb,:)  = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
        impSAreal(ispin,ispin,iorb,iorb,:) = invF0real(ispin,ispin,iorb,iorb,:) - invFreal(ispin,ispin,iorb,iorb,:)
     enddo
     !
  case ("hybrid")
     !
     !Get Gimp^-1
     do i=1,Lmats
        invGimp=zero
        invGimp(1:Norb,1:Norb)               = impGmats(ispin,ispin,:,:,i)
        invGimp(1:Norb,Norb+1:2*Norb)        = impFmats(ispin,ispin,:,:,i)
        invGimp(Norb+1:2*Norb,1:Norb)        = impFmats(ispin,ispin,:,:,i)
        invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGmats(ispin,ispin,:,:,i))
        call inv(invGimp)
        invGmats(ispin,ispin,:,:,i) = invGimp(1:Norb,1:Norb)
        invFmats(ispin,ispin,:,:,i) = invGimp(1:Norb,Norb+1:2*Norb)
     enddo
     do i=1,Lreal
        invGimp=zero
        invGimp(1:Norb,1:Norb)               = impGreal(ispin,ispin,:,:,i)
        invGimp(1:Norb,Norb+1:2*Norb)        = impFreal(ispin,ispin,:,:,i)
        invGimp(Norb+1:2*Norb,1:Norb)        = impFreal(ispin,ispin,:,:,i)
        invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGreal(ispin,ispin,:,:,Lreal-i+1))
        call inv(invGimp)
        invGreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,1:Norb)
        invFreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,Norb+1:2*Norb)
     enddo
     !Get Sigma functions: Sigma= G0^-1 - G^-1
     impSmats=zero
     impSAmats=zero
     impSreal=zero
     impSAreal=zero
     !
     impSmats(ispin,ispin,:,:,:)  = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
     impSAmats(ispin,ispin,:,:,:) = invF0mats(ispin,ispin,:,:,:) - invFmats(ispin,ispin,:,:,:)
     !
     impSreal(ispin,ispin,:,:,:)  = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
     impSAreal(ispin,ispin,:,:,:) = invF0real(ispin,ispin,:,:,:) - invFreal(ispin,ispin,:,:,:)
     !
  end select
  !
  !Get G0and:
  impG0mats(ispin,ispin,:,:,:) = g0and_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
  impF0mats(ispin,ispin,:,:,:) = f0and_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
  !
  impG0real(ispin,ispin,:,:,:) = g0and_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
  impF0real(ispin,ispin,:,:,:) = f0and_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
  !!
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
end subroutine get_sigma_superc
