require 'mesa_script'


old_filename = 'DA_0.6.mod'

teffs = (1...31).map{|i| 21000 - i*500}
puts teffs

teffs.each do |teff|

  teff_string = sprintf("%.0f",teff)
  filename = 'DA_' + teff_string + 'K.mod'

  puts "loading file " + old_filename
  puts "saving file  " + filename

  Inlist.make_inlist('inlist_coolDA') do
    # &star_job  
  
      load_saved_model true
      saved_model_name old_filename

      save_model_when_terminate true
      save_model_filename filename
  
      set_tau_factor true  
      set_to_this_tau_factor 1.0
  
      set_initial_dt true
      # years_for_initial_dt = 1d-4
      years_for_initial_dt 1e2
  
      change_net true
      new_net_name 'pollution.net'
  
      # display on-screen plots
      pgstar_flag true
  
  
      # &controls
  
          history_interval 1

          # Make sure Teffs are accurate to at least 3 digits. (e.g. 20500 +- 50 at most)
          delta_lgTeff_limit 0.0005

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
  
        teff_lower_limit teff
  
        # turn off burning so no diffusion induced flashes
        # keep hydrogen layer thick
        max_abar_for_burning -1
        
        # element diffusion
        do_element_diffusion true
        diffusion_use_full_net true
  
        # Change this to be smaller than the surface convection zone for pollution
        diffusion_min_dq_at_surface 1e-16
        diffusion_min_dq_ratio_at_surface 3
        diffusion_min_t_at_surface 5e3
  
        show_diffusion_info true # terminal output for diffusion
        # show_diffusion_substep_info = .true. ! terminal output for diffusion
        # show_diffusion_timing = .true.
  
        diffusion_dt_limit 1e-10 # in seconds
        diffusion_v_max 1e-2
        diffusion_steps_limit 50
        diffusion_steps_hard_limit 100
  
        # diffusion_use_isolve = .true.
        # diffusion_rtol_for_isolve = 1d-4
        # diffusion_atol_for_isolve = 1d-5
        # diffusion_maxsteps_for_isolve = 1000
        # diffusion_isolve_solver = 'ros2_solver'
  
  
  
        # / ! end of controls namelist
  end

  system('./rn')
  old_filename = filename
  
end
