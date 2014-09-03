MODULE DMFT_ED
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL, only: Hloc,&
       impSmats,impSAmats,impSreal,impSAreal,&
       ed_dens,ed_docc,ed_phisc,ed_Ekin,ed_Epot,ed_Ehartree,&
       ED_MPI_ID,ED_MPI_SIZE,ED_MPI_ERR
  USE ED_AUX_FUNX, only: print_Hloc,search_chemical_potential
  USE ED_BATH, only: get_bath_size,spin_symmetrize_bath,ph_symmetrize_bath,break_symmetry_bath
  USE ED_MAIN, only: init_ed_solver,ed_solver,ed_kinetic_energy
  !
  USE ED_CHI2FIT
END MODULE DMFT_ED

