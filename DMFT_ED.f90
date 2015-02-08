MODULE DMFT_ED
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL, only: &
       impSmats,impSAmats,impSreal,impSAreal,&
       ed_dens,ed_dens_up,ed_dens_dw,ed_docc,ed_phisc,&
       ed_Ekin,ed_Epot,ed_Ehartree,ed_Eknot,&
       ed_Dust,ed_Dund,ed_Dse,ed_Dph, &
       ED_MPI_ID,ED_MPI_SIZE,ED_MPI_ERR
  USE ED_AUX_FUNX, only: &
       print_Hloc,set_Hloc,get_Hloc,&
       get_Sigma,set_sigma,&
       search_chemical_potential
  USE ED_BATH, only: &
       allocate_bath,deallocate_bath,set_bath,copy_bath,get_bath_size,spin_symmetrize_bath,ph_symmetrize_bath,break_symmetry_bath,enforce_normal_bath,check_bath_dimension
  USE ED_BATH_TYPE
  USE ED_MAIN, only: init_ed_solver,ed_solver
  USE ED_ENERGY, only: ed_kinetic_energy,ed_kinetic_energy_sc
  USE ED_CHI2FIT
END MODULE DMFT_ED

