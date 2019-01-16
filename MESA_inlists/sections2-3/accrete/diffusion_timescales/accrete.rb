require 'mesa_script'

teffs = (1...31).map{|i| 21000 - i*500}
puts teffs

koester_dts_ca = [4.2924,3.7125,3.3303,2.8725,1.9997,0.4845,
                  -1.6941,-2.3590,-2.6305,-2.6277,-2.6220,-2.4688,
                  -2.4213,-2.3804,-2.3077,-2.2713,-2.2821,-2.2527,
                  -2.2271,-2.1953]
koester_teffs = [6000,7000,8000,9000,10000,11000,12000,
                 13000,14000,15000,16000,17000,18000,
                 19000,20000,21000,22000,23000,24000,25000]

koester_dts = Hash[koester_teffs.zip(koester_dts_ca)]

teffs.each do |teff|

  # Math to setup the things that vary for each run.
  teff_string = sprintf("%.0f",teff)
  filename = '../../make_more_models/DA_' + teff_string + 'K.mod'
  puts "loading file " + filename
  logdir = "LOGS_"+teff_string
  puts "saving to " + logdir
  logfile = teff_string + "_logfile.txt"
  puts "logfile " + logfile
  teff_rounded = (teff/1000)*1000 # abusing integer division to round
  koester_dt = koester_dts[teff_rounded] # rounded to nearest thousand to get a ballpark number from Koester's table
  puts "Temperature: " + teff_string
  puts "Koester timescale: " + sprintf("%.2f",koester_dt)
  maxdt = (10**koester_dt)/20.0
  maxage = maxdt*20.0*30.0 # Make sure going for many diffusion timescales.
  if(teff < 10500) then
    diffminT = 3e4
    maxv = 1e-2
  elsif(teff < 12500) then
    diffminT = 1.5e4
    maxv = 1e0
  else
    # Effectively turning off this setting.
    diffminT = 5e3
    maxv = 1e2
  end

  Inlist.make_inlist('inlist_accrete') do
    # &star_job
  
  
      load_saved_model true
      saved_model_name filename
  
      set_tau_factor true  
      set_to_this_tau_factor 0.1
  
      set_initial_dt true
      years_for_initial_dt 1e-5
  
      set_initial_age true
      initial_age 0
      set_initial_model_number true
      initial_model_number 0
  
      change_net true
      new_net_name 'pollution.net'
  
    # display on-screen plots
      pgstar_flag true
  
  # / !end of star_job namelist
  
  
  # &controls
  
        log_directory logdir
        max_years_for_timestep maxdt
        max_age maxage
        # If the convection is deep and we want to assume fully mixed,
        # convergence is better when you merge a reasonable number of surface cells.
        diffusion_min_t_at_surface diffminT
        t_mix_limit 2e4
        diffusion_v_max maxv 
  
        # MUST be smaller than the surface convection region at least.
        # Note that this doesn't actually matter for diffusion timescales or evolution,
        # it just applies to reported surface abundances.
        surface_avg_abundance_dq 1e-15 # A little deeper than the photosphere.

        max_abar_for_burning -1
  
        history_interval 1
        profile_interval 1000
        photo_interval 10000
        mesh_delta_coeff 0.75
        varcontrol_target 1e-4
        # varcontrol_target = 3d-4
  
        # High surface resolution
        min_dq 1e-21
        max_surface_cell_dq 1e-21
  
        use_type2_opacities false
  
        # Chosen to match Koester
        mlt_option 'ML2'
        mixing_length_alpha 0.8 # For DAs.
        # use_Ledoux_criterion = .true.
  
        which_atm_option 'grey_and_kap'
        # Switch atm to WD_tau_25_tables for Teff < x_ctrl(10)
        x_ctrl[10] = 9e3
  
        accrete_same_as_surface false
        accrete_given_mass_fractions true
        num_accretion_species 1
        accretion_species_id[1] = 'ISO_HERE'
        accretion_species_xa[1] = 1
        mass_change 1.586e-19 # Accretion rate in msun/yr
        # 1 msun/yr = 6.305e25 g/s
        # 1.586d-19 msun/yr = 1d7 g/s
  
        # element diffusion
        do_element_diffusion true
        diffusion_use_full_net true
  
        # Keep this smaller than the surface convection zone for pollution
        diffusion_min_dq_at_surface 2e-18
        diffusion_min_dq_ratio_at_surface 3
  
        show_diffusion_info true # terminal output for diffusion
        # show_diffusion_substep_info = .true. ! terminal output for diffusion
        # show_diffusion_timing = .true.
  
        diffusion_dt_limit 1e-15 # in seconds
        diffusion_steps_limit 5
        diffusion_steps_hard_limit 50
  
        # diffusion_use_isolve = .true.
        # diffusion_rtol_for_isolve = 1d-4
        # diffusion_atol_for_isolve = 1d-5
        # diffusion_maxsteps_for_isolve = 1000
        # diffusion_isolve_solver = 'ros2_solver'
  
  
     # / ! end of controls namelist
  end

  execute_command = './rn | tee ' + logfile
  system(execute_command)
  
end
