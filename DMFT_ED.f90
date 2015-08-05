MODULE DMFT_ED
  USE ED_INPUT_VARS

  USE ED_VARS_GLOBAL, only: &
       Nlat        ,&
       icol        ,&
       irow        ,&
       ij2site     ,&
       H0          ,&
       impSmats    ,&
       impSAmats   ,&
       impSreal    ,&
       impSAreal   ,&
       ed_dens     ,&
       indep_list  ,&
       map_lat2ind ,&
       map_ind2lat ,&
       ED_MPI_ID,ED_MPI_SIZE,ED_MPI_ERR

  USE ED_AUX_FUNX, only:&
       set_Hloc,&
       get_Hloc,&
       extract_Hloc          ,&
       blocks_to_matrix      ,&
       matrix_to_blocks      ,&
       select_block          ,&
       stride_index          ,&       
       lso2nnn_reshape       ,&
       nnn2lso_reshape       ,&
       get_independent_sites ,&  
       get_lattice_hamiltonian,& !OBSOLETE (to be removed)
       search_chemical_potential


  USE ED_BATH, only: &
       get_bath_size        ,&
       check_bath_dimension ,&
       spin_symmetrize_bath ,&
       ph_symmetrize_bath   ,&
       ph_trans_bath        ,&
       break_symmetry_bath  ,&
       enforce_normal_bath


  USE ED_MAIN,      only: &
       ed_init_solver         ,&
       ed_solve               ,&
       ed_get_sigma_matsubara ,&
       ed_get_self_matsubara  ,&
       ed_get_sigma_real      ,&
       ed_get_self_real       ,&
       ed_get_gimp_matsubara  ,&
       ed_get_fimp_matsubara  ,&
       ed_get_gimp_real       ,&
       ed_get_fimp_real       ,&
       ed_get_dens            ,&
       ed_get_mag             ,&
       ed_get_docc            ,&
       ed_get_phisc           ,&
       ed_get_eimp            ,&
       ed_get_epot            ,&
       ed_get_eint            ,&
       ed_get_ehartree        ,&
       ed_get_eknot           ,&
       ed_get_doubles         ,&
       ed_get_dust            ,&
       ed_get_dund            ,&
       ed_get_dse             ,&
       ed_get_dph             ,&
       ed_init_solver_lattice         ,&
       ed_solve_lattice               ,&
       ed_get_sigma_matsubara_lattice ,&
       ed_get_self_matsubara_lattice  ,&
       ed_get_sigma_real_lattice      ,&
       ed_get_self_real_lattice       ,&
       ed_get_gimp_matsubara_lattice  ,&
       ed_get_fimp_matsubara_lattice  ,&
       ed_get_gimp_real_lattice       ,&
       ed_get_fimp_real_lattice       ,&
       ed_get_dens_lattice            ,&
       ed_get_mag_lattice             ,&
       ed_get_docc_lattice            ,&
       ed_get_phisc_lattice           ,&
       ed_get_eimp_lattice            ,&
       ed_get_epot_lattice            ,&
       ed_get_eint_lattice            ,&
       ed_get_ehartree_lattice        ,&
       ed_get_eknot_lattice           ,&
       ed_get_doubles_lattice         ,&
       ed_get_dust_lattice            ,&
       ed_get_dund_lattice            ,&
       ed_get_dse_lattice             ,&
       ed_get_dph_lattice


  USE ED_GLOC,     only:&
       ed_get_gloc,&
       ed_get_gloc_lattice,&
       ed_get_gij_lattice


  USE ED_WEISS,    only: ed_get_weiss,ed_get_weiss_lattice


  USE ED_ENERGY,   only: ed_kinetic_energy,ed_kinetic_energy_lattice


  USE ED_CHI2FIT,  only: ed_chi2_fitgf,ed_chi2_fitgf_lattice


END MODULE DMFT_ED

