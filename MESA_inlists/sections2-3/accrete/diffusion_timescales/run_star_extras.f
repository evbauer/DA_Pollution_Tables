! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! Switch atmosphere option below a certain temperature.
         if(s% Teff < s% x_ctrl(10)) then
            s% which_atm_option = 'WD_tau_25_tables'
            s% tau_factor = 1
         end if

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 14
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         use ionization_lib, only: eval_typical_charge
         use chem_lib, only: chem_get_iso_id
         use chem_def
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: k, k_cz_bot, kphot, k_mass_lim, itrace
         real(dp) :: cz_xm, edv_iso
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Id in the NET for the iso we wish to trace.
         itrace = s% net_iso(ISO_HERE) ! e.g. s% net_iso(ica40) for Calcium
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         ! Mass exterior to tau_eff = 5
         ! Turn this one into a dummy by default, since compiling requires star_utils module.
         names(1) = 'photosphere5_xm_eff'
         vals(1) = -1
         !do k = 2,s% nz ! first k is not a tau_eff
         !   if(tau_eff(s,k) > 5d0) then
         !      vals(1) = s% xmstar*sum(s% dq(1:k-1))/Msun
         !      exit
         !   end if
         !end do

         ! Mass exterior to tau = 5
         names(2) = 'photosphere5_xm'
         do k = 2,s% nz ! first k is not a tau_eff
            if(s% tau(k) > 5d0) then
               vals(2) = s% xmstar*sum(s% dq(1:k-1))/Msun
               exit
            end if
         end do
         
         ! Mass exterior to tau = 2/3
         names(3) = 'photosphere23_xm'
         do k = 2,s% nz
            if(s% tau(k) > 2d0/3d0) then
               vals(3) = s% xmstar*sum(s% dq(1:k-1))/Msun
               kphot = k
               exit
            end if
         end do

         ! If Teff > 14000, then set a limit over which to look for k_cz_bot
         ! Otherwise, can look deeper. (Prevents finding the small zone below the surface at hotter temps)
         k_mass_lim = s% nz
         if(s% Teff > 14d3) then
            do k = 2, s% nz
               if(s% xmstar*sum(s% dq(1:k-1))/Msun > 1d-10) then
                  k_mass_lim = k
                  exit
               end if
            end do
         end if
         
         ! Mass exterior to tau = 1
         names(4) = 'photosphere1_xm'
         do k = 2,s% nz
            if(s% tau(k) > 1d0) then
               vals(4) = s% xmstar*sum(s% dq(1:k-1))/Msun
               exit
            end if
         end do

         names(5) = 'evan_cz_xm'
         ! Find cell where switches from convective to not
         k_cz_bot = 1
         do k=2, k_mass_lim
            if(s% mixing_type(k-1) == convective_mixing .and. &
                 s% mixing_type(k) /= convective_mixing) then
               k_cz_bot = k - 1
               exit
            end if
         end do
         ! print *, k_cz_bot
         ! print *, s% cz_bdy_dq(k_cz_bot)
         ! vals(5) = s% xmstar*(sum(s% dq(1:k_cz_bot-1)) - s% cz_bdy_dq(k_cz_bot))/Msun

         if(k_cz_bot < kphot) then
            ! Surface convection doesn't exist or is even shallower
            ! than the photosphere, so just use photosphere.
            k_cz_bot = kphot
         end if

         ! k_cz_bot is the last convective cell,
         ! the face of k_cz_bot + 1 is the face at the boundary of the convective zone
         ! (roughly, can also try to interpolate within cells to find exact boundary)
         cz_xm = s% xmstar*sum(s% dq(1:k_cz_bot))/Msun
         vals(5) = cz_xm

         ! Also use location of base of the convection zone to evaluate tyipcal charges there.
         names(6) = 'base_charge_Na'
         names(7) = 'base_charge_Mg'
         names(8) = 'base_charge_Si'
         names(9) = 'base_charge_Ca'
         names(10) = 'base_charge_Fe'
         vals(6) = s% typical_charge(s% net_iso(ina23),k_cz_bot)
         vals(7) = s% typical_charge(s% net_iso(img24),k_cz_bot)
         vals(8) = s% typical_charge(s% net_iso(isi28),k_cz_bot)
         vals(9) = s% typical_charge(s% net_iso(ica40),k_cz_bot)
         vals(10) = s% typical_charge(s% net_iso(ife56),k_cz_bot)

         ! Most quantities here are evaluated at the face of k_cz_bot+1,
         ! the outer edge of the first non-convective cell.
         ! See MESA1 fig 9.
         edv_iso = abs(s% edv(itrace,k_cz_bot+1))
         names(11) = 'inst_diffusion_timescale'
         vals(11) = cz_xm*Msun / ( pi4*pow2(s% r(k_cz_bot+1))*s% rho_face(k_cz_bot+1)*edv_iso)

         names(12) = 'steady_diffusion_timescale'
         ! Backs up a bit into convective zone to evaluate mass fraction.
         vals(12) = &
              s% xa(itrace,k_cz_bot-1) * &
              cz_xm*Msun / (s% star_mdot * Msun/secyer)

         names(13) = 'cz_bot_zone'
         vals(13) = k_cz_bot

         names(14) = 'edv_base'
         vals(14) = edv_iso

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info

      end module run_star_extras
      
