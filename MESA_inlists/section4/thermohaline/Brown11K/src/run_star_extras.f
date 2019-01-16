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

      ! use star_utils, only: tau_eff
      
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
         how_many_extra_history_columns = 18
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
         itrace = s% net_iso(ica40) ! e.g. s% net_iso(ica40) for Calcium
         
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
         do k = 2,s% nz ! first k is not a tau_eff
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
         do k = 2,s% nz ! first k is not a tau_eff
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

         ! If the convection is down in the helium or deeper,
         ! this is not the comnvection we want.
         if(s% xmstar*sum(s% dq(1:k_cz_bot))/Msun > 1d-5) then
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
         vals(12) = -1
         ! vals(12) = &
         !      s% xa(itrace,k_cz_bot-2) * &
         !      cz_xm*Msun / (s% star_mdot * Msun/secyer)

         names(13) = 'cz_bot_zone'
         vals(13) = k_cz_bot

         names(14) = 'edv_base'
         vals(14) = edv_iso

         itrace = s% net_iso(img24) 
         edv_iso = abs(s% edv(itrace,k_cz_bot+1))
         names(15) = 'inst_Mg'
         vals(15) = cz_xm*Msun / ( pi4*pow2(s% r(k_cz_bot+1))*s% rho_face(k_cz_bot+1)*edv_iso)

         itrace = s% net_iso(isi28) 
         edv_iso = abs(s% edv(itrace,k_cz_bot+1))
         names(16) = 'inst_Si'
         vals(16) = cz_xm*Msun / ( pi4*pow2(s% r(k_cz_bot+1))*s% rho_face(k_cz_bot+1)*edv_iso)
         
         itrace = s% net_iso(io16) 
         edv_iso = abs(s% edv(itrace,k_cz_bot+1))
         names(17) = 'inst_O'
         vals(17) = cz_xm*Msun / ( pi4*pow2(s% r(k_cz_bot+1))*s% rho_face(k_cz_bot+1)*edv_iso)

         itrace = s% net_iso(ife56) 
         edv_iso = abs(s% edv(itrace,k_cz_bot+1))
         names(18) = 'inst_Fe'
         vals(18) = cz_xm*Msun / ( pi4*pow2(s% r(k_cz_bot+1))*s% rho_face(k_cz_bot+1)*edv_iso)


      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         use chem_def
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, k, nc, m, hi, oi, fei
         ! These must have dimension matching the net (pollution3.net) + 1 for electrons.
         real(dp), dimension(18) :: A, charge, nd
         real(dp), dimension(18,18) :: Kdiff,zdiff,zdiff1,zdiff2
         real(dp) :: kappa, DHFe, DHO
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


         ! print *, "saving extra profile_columns" 
         ! Number of species in pollution3.net
         nc = 17
         m = nc+1
         if(nc /= s% species) then
            stop "net doesn't match assumed size in run_star_extras"
         end if
         
         names(1) = 'kappa_SM'
         names(2) = 'D12_H_Fe'
         names(3) = 'D12_H_O'

         A(1:nc) = chem_isos% Z_plus_N(s% chem_id(1:nc))
         A(m) = me/amu

         ! indices for various species
         hi = s% net_iso(ih1)
         oi = s% net_iso(io16)
         fei = s% net_iso(ife56)
         ! print *, "Iron index:", fei

         do k=1,s%nz

            do j=1,nc
               charge(j) = max( 1d0, s% typical_charge(j,k) )
               nd(j) = s% rho(k) * s% xa(j,k) / (A(j)*amu)
            end do
            ! print *, "Charge for 5th net iso: ", charge(5)
            ! Electron entries
            charge(m) = -1d0
            nd(m) = 0d0
            do j=1,nc
               nd(m) = nd(m) + charge(j)*nd(j)
            end do
            
            call get_SM_coeffs( &
                 nc,m, &
                 s% rho(k), s% T(k), &
                 A,charge,nd, &
                 Kdiff,zdiff,zdiff1,zdiff2,kappa)

            ! D12 = ( n1 n2/(n1 + n2) ) k T / K12 
            DHFe = nd(hi)*nd(fei)*boltzm*s% T(k) / ( (nd(hi) + nd(fei)) * Kdiff(hi,fei))
            DHO = nd(hi)*nd(oi)*boltzm*s% T(k) / ( (nd(hi) + nd(oi)) * Kdiff(hi,oi))

            vals(k,1) = kappa
            vals(k,2) = DHFe
            vals(k,3) = DHO
         end do
         
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

       ! Calculate coefficients given in Appendix C.3 of Stanton & Murillo, PR E 93, 043203 (2016)
      subroutine get_SM_coeffs(nc,m,rho,T,A,Z,nd,Kdiff,zdiff,zdiff1,zdiff2,kappa)
        integer, intent(in) :: nc, m
        real(dp), intent(in) :: rho, T
        real(dp), dimension(:), intent(in) :: A, Z, nd ! m
        real(dp), intent(out) :: kappa ! ion separation relative to electron screening length.
        real(dp), dimension(:,:), intent(inout) :: Kdiff,zdiff,zdiff1,zdiff2 ! (m,m) but electron entries shouldn't be used.
        ! Only the ion terms are modified by this routine right now.

        real(dp), dimension(nc,nc) :: g_plasma, mu
        real(dp), dimension(nc,nc,2,3) :: Keff, Omega ! Indexed by all species and different orders
        ! Fitting coefficients from Stanton & Murillo
        real(dp), dimension(2,3) :: a1,a2,a3,a4,a5,b0,b1,b2,b3,b4
        real(dp) :: lambda ! The screening length
        integer :: i,j,k,facmo,no,mo ! Last two are used for different orders of collisions.

        ! Initialize all to 0, since some entries never get set or used.
        a1(:,:) = 0d0
        a2(:,:) = 0d0
        a3(:,:) = 0d0
        a4(:,:) = 0d0
        a5(:,:) = 0d0
        b0(:,:) = 0d0
        b1(:,:) = 0d0
        b2(:,:) = 0d0
        b3(:,:) = 0d0
        b4(:,:) = 0d0

        ! Table IV, coefficients for fits (C23)-(C24)
        a1(1,1) = 1.4660d0
        a1(1,2) = 0.52094d0
        a1(1,3) = 0.30346d0
        a1(2,2) = 0.85401d0

        a2(1,1) = -1.7836d0
        a2(1,2) = 0.25153d0
        a2(1,3) = 0.23739d0
        a2(2,2) = -0.22898d0

        a3(1,1) = 1.4313d0
        a3(1,2) = -1.1337d0
        a3(1,3) = -0.62167d0
        a3(2,2) = -0.60059d0

        a4(1,1) = -0.55833d0
        a4(1,2) = 1.2155d0
        a4(1,3) = 0.56110d0
        a4(2,2) = 0.80591d0

        a5(1,1) = 0.061162d0
        a5(1,2) = -0.43784d0
        a5(1,3) = -0.18046d0
        a5(2,2) = -0.30555d0

        b0(1,1) = 0.081033d0
        b0(1,2) = 0.20572d0
        b0(1,3) = 0.68375d0
        b0(2,2) = 0.43475d0

        b1(1,1) = -0.091336d0
        b1(1,2) = -0.16536d0
        b1(1,3) = -0.38459d0
        b1(2,2) = -0.21147d0

        b2(1,1) = 0.051760d0
        b2(1,2) = 0.061572d0
        b2(1,3) = 0.10711d0
        b2(2,2) = 0.11116d0

        b3(1,1) = -0.50026d0
        b3(1,2) = -0.12770d0
        b3(1,3) = 0.10649d0
        b3(2,2) = 0.19665d0

        b4(1,1) = 0.17044d0
        b4(1,2) = 0.066993d0
        b4(1,3) = 0.028760d0
        b4(2,2) = 0.15195d0

        ! Get the screening length
        call lam_SM(nc,m,rho,T,Z,nd,lambda,kappa)

        ! Calculate g_plasma
        do i=1,nc
           do j=1,nc
              g_plasma(i,j) = Z(i)*Z(j)*qe*qe/(lambda*boltzm*T)
           end do
        end do

        ! Calculate Keff using Eqns (C22-C24)
        do i=1,nc
           do j=1,nc
              if( g_plasma(i,j) < 0d0) then ! Don't calculate for attractive potentials, set to 0
                 Keff(i,j,:,:) = 0d0
              else if( g_plasma(i,j) < 1d0) then ! Use eqn C23 for weakly coupled
                 do no=1,2
                    do mo=1,3 ! Implementing the (m-1)! term with a simple if statement.
                       if(mo .eq. 3) then
                          facmo = 2 ! (3-1)!
                       else
                          facmo = 1 ! (1-1)! and (2-1)!
                       end if
                       Keff(i,j,no,mo) = (-1d0*no/4d0)*facmo*safe_log_cr( &
                            a1(no,mo)*g_plasma(i,j) &
                            + a2(no,mo)*pow2(g_plasma(i,j)) &
                            + a3(no,mo)*pow3(g_plasma(i,j)) &
                            + a4(no,mo)*pow4(g_plasma(i,j)) &
                            + a5(no,mo)*pow5(g_plasma(i,j)) )
                    end do
                 end do
              else ! Use eqn C24 for strongly coupled
                 do no=1,2
                    do mo=1,3
                       Keff(i,j,no,mo) = &
                            (b0(no,mo) + b1(no,mo)*safe_log_cr(g_plasma(i,j)) &
                            + b2(no,mo)*pow2(safe_log_cr(g_plasma(i,j))) ) / &
                            (1d0 + b3(no,mo)*g_plasma(i,j) + b4(no,mo)*pow2(g_plasma(i,j)) )
                    end do
                 end do
              end if
           end do
        end do

        ! Calculate the collision integrals using (C19)
        do i=1,nc
           do j=1,nc
              mu(i,j) = amu*A(i)*A(j)/(A(i) + A(j)) ! Reduced mass for collision
              do no=1,2
                 do mo=1,3
                    Omega(i,j,no,mo) = &
                         sqrt(2d0*pi/(mu(i,j))) * &
                         (pow2(Z(i)*Z(j)*qe*qe)/pow_cr(boltzm*T,1.5d0)) * &
                         Keff(i,j,no,mo)
                 end do
              end do
           end do
        end do

        ! Resistance coefficient is (16/3) ni nj mu Omega11
        ! (compare Paquette eqn 22 with SM eq B8)
        ! The zdiffs are given in terms of collision integrals in Paquette eqns (14-16), (23-25)
        do i=1,nc
           do j=1,nc
              Kdiff(i,j) = (16d0/3d0)*nd(i)*nd(j)*mu(i,j)*Omega(i,j,1,1)
              zdiff(i,j) = 1d0 - 2d0*Omega(i,j,1,2)/(5d0*Omega(i,j,1,1))
              zdiff1(i,j) = 2.5d0 + (2d0*Omega(i,j,1,3) - 10d0*Omega(i,j,1,2))/(5d0*Omega(i,j,1,1))
              zdiff2(i,j) = Omega(i,j,2,2)/Omega(i,j,1,1)
           end do
        end do

        ! Note that SM don't give fits for attractive potentials, so we haven't
        ! touched the electron entries. They exit this routine unchanged, so they
        ! either need to be initialized before this routine is called or somehow
        ! calculated later.
        
      end subroutine get_SM_coeffs

      ! Screening Length according to Stanton & Murillo
      subroutine lam_SM(nc,m,rho,T,Z,nd,lam_eff,kappa)
        integer, intent(in) :: nc, m
        real(dp), intent(in) :: rho, T
        real(dp), dimension(:), intent(in) :: Z, nd ! m. charges, number densities of all species
        real(dp), intent(out) :: lam_eff
        real(dp), intent(out) :: kappa ! ion separation relative to electron screening length.

        real(dp) :: tiny_n, ne, EF, lam_e, lam_sum, rhotot, ni_tot
        real(dp), dimension(nc) :: ai, gam_is, lam_i ! ion sphere radius, coupling, screening for each type of ion
        integer :: i

        tiny_n = 1d-20 ! g/cc
        ne = nd(m)
        ! Electron Fermi energy
        EF = (hbar*hbar*pow_cr(3d0*pi*pi*ne,2d0/3d0))/(2d0*me)
        ! Electron screening length accounting for degeneracy correction
        lam_e = pow_cr(pi4*qe*qe*ne/sqrt(pow2(boltzm*T) + pow2((2d0/3d0)*EF)),-0.5d0)

        ! Compute kappa
        ni_tot = 0d0
        do i = 1,nc
           ni_tot = ni_tot + nd(i)
        end do
        kappa = pow_cr(3d0/(pi4*ni_tot),1d0/3d0)/lam_e

        rhotot = 0d0
        do i = 1,nc
           ! rhotot is a CHARGE density. Just trying to follow the notation of SM eq (34)
           rhotot = rhotot + Z(i)*qe*nd(i) ! = qe*ne?
        end do
        do i = 1,nc
           ai(i) = pow_cr(3d0*Z(i)*qe/(pi4*rhotot),1d0/3d0)
           gam_is(i) = Z(i)*Z(i)*qe*qe/(ai(i)*boltzm*T)
           ! Number densities that are 0 or tiny cause screening length to diverge.
           ! This is physical; nothing is there to screen anything.
           ! But numerically, I don't want to rely on fortran to handle division by zero and infinity,
           ! so just set these screening lengths to be huge and they won't
           ! contribute anything to overall screening length.
           if(nd(i) < tiny_n) then
              lam_i(i) = 1d99 ! cm. This won't contribute to any screening.
           else
              lam_i(i) = sqrt(boltzm*T/( pi4*pow2(Z(i)*qe)*nd(i) ))
           end if
        end do

        lam_sum = 1d0/pow2(lam_e) ! The electron part of the screening length.
        do i = 1,nc ! Sum over all the ions.
           lam_sum = lam_sum + pow_cr( pow2(lam_i(i))*(1d0+3d0*gam_is(i)), -1d0)
        end do
        lam_eff = pow_cr(lam_sum,-0.5d0)
      end subroutine lam_SM


      end module run_star_extras
      
