  states: do i=first_state,last_state
     m = H%map(i)
     impi = i-ishift
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     !
     !
     !IMPURITY  HAMILTONIAN
     include "ED_HAMILTONIAN_MATVEC/HxVimp.f90"


     !LOCAL INTERACTION
     include "ED_HAMILTONIAN_MATVEC/HxVint.f90"


     !BATH HAMILTONIAN
     include "ED_HAMILTONIAN_MATVEC/HxVbath.f90"


     !IMPURITY- BATH HYBRIDIZATION
     include "ED_HAMILTONIAN_MATVEC/HxVimp_bath.f90"

  enddo states

