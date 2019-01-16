require 'mesa_script'

# 1.586e-16 = 1e10 g/s
log_mdots = (0...17).map{|i| -22 + i*0.5 }
puts log_mdots

# Model to load for this Teff
loadsave = '../../../sections2-3/make_more_models/DA_11000K.mod'
# Diffusion timescale for this Teff (from table 2 based on diffusion only calcs)
logtaudiff = 1.2236
# Surface conv mass for this Teff (also from table 2)
log_surf_cvz_dq = -11.87
# Parameters for resolving diffusion near the surface
diffMinT = 3e4 # Needs to be low for Teffs where there's no/thin surface convection.
diffVmax = 1e0 # Needs to be high (1e2) for Teffs where there's no/thin surface convection.


# Makes sure our reported surface abundances really reflect surface value.
avg_dq = (10**log_surf_cvz_dq)/2 

log_mdots.each do |log_mdot|
  logmdot_string = sprintf("%.2f",log_mdot)
  logdir = "LOGS"+logmdot_string
  puts logdir
  logfile = "logfile" + logmdot_string + ".txt"
  mdot = 10**log_mdot
  if(log_mdot > -15) then
    max_timestep = (10**logtaudiff)/(50.0)
  else
    max_timestep = (10**logtaudiff)/10.0
  end

  # 100 diffusion timescales
  maxage = (10**logtaudiff)*100
  
  Inlist.make_inlist('inlist_accrete') do
  
  # &star_job
    pgstar_flag true
  
    load_saved_model true
    saved_model_name loadsave
  
    set_initial_age true
    initial_age 0
    set_initial_model_number true
    initial_model_number 0
    set_initial_dt true
    years_for_initial_dt 1e-05
  
    set_to_this_tau_factor 0.1
    set_tau_factor true
  
    change_net true
    new_net_name 'pollution3.net'
  
  # / ! end of star_job namelist
  
    # &controls

        mass_change mdot
        surface_avg_abundance_dq avg_dq
        max_age maxage # Very many diffusion timescales
        max_years_for_timestep max_timestep

        diffusion_min_t_at_surface diffMinT
        diffusion_v_max diffVmax
  
        t_mix_limit 50000
  
          use_ledoux_criterion true
          thermohaline_coeff 1
          thermohaline_option 'Brown_Garaud_Stellmach_13'
  
          log_directory logdir
  
    use_gold_tolerances false
    use_eosDT2 true
    tiny_fuzz_for_eos 1e-12          

    d_mix_ignore_diffusion 1e15
    max_abar_for_burning -1

    photo_interval 100
  
    history_interval 10
  
    profile_interval 1000
  
  
    mixing_length_alpha 0.8
  
    mlt_option 'ML2'
  
    which_atm_option 'grey_and_kap'
  
    accrete_same_as_surface false
    accretion_dump_missing_metals_into_heaviest false
    accrete_given_mass_fractions true  
    num_accretion_species 10
    # Bulk Earth (McDonough 2000)
    accretion_species_id[1] = 'fe56'
    accretion_species_xa[1] = 0.319
    accretion_species_id[2] = 'o16'
    accretion_species_xa[2] = 0.297
    accretion_species_id[3] = 'si28'
    accretion_species_xa[3] = 0.161
    accretion_species_id[4] = 'mg24'
    accretion_species_xa[4] = 0.154
    accretion_species_id[5] = 'ni58'
    accretion_species_xa[5] = 0.01822
    accretion_species_id[6] = 'ca40'
    accretion_species_xa[6] = 0.0171
    accretion_species_id[7] = 'al27'
    accretion_species_xa[7] = 0.0159
    accretion_species_id[8] = 's32'
    accretion_species_xa[8] = 0.00635
    accretion_species_id[9] = 'cr52'
    accretion_species_xa[9] = 0.0047
    accretion_species_id[10] = 'na23'
    accretion_species_xa[10] = 0.0018
  
  
    mesh_delta_coeff 0.75
  
    min_dq 1e-21
  
    max_surface_cell_dq 1e-21
  
    show_diffusion_info true
  
    do_element_diffusion true
    diffusion_dt_limit 1e-15
  
    diffusion_min_dq_at_surface avg_dq
    diffusion_min_dq_ratio_at_surface 3
  
    diffusion_use_full_net true
  
    use_type2_opacities false
  
    varcontrol_target 0.0001
  
    diffusion_steps_limit 5
    diffusion_steps_hard_limit 50
  
    x_ctrl[10] = 9300
  
  
  # / ! end of controls namelist
  
  # &pgstar
  
  # / ! end of pgstar namelist
  end

  execute_command = './rn | tee ' + logfile
  puts execute_command
  system(execute_command)

end
