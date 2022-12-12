! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
      use chem_def
      use binary_def
      use math_lib

      implicit none
      
    contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns 
         
         !By SG 
         s% other_alpha_mlt => alpha_mlt_routine  
         !New routines to add

	      s% other_am_mixing => TSF

         ! s% other_wind => brott_wind
      end subroutine extras_controls
      
      !By SG
      subroutine alpha_mlt_routine(id, ierr)
         use chem_def, only: ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, h1
         real(dp) :: alpha_H, alpha_other, H_limit
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         alpha_H = s% x_ctrl(21)
         alpha_other = s% x_ctrl(22)
         H_limit = s% x_ctrl(23)
         h1 = s% net_iso(ih1)
         !write(*,1) 'alpha_H', alpha_H
         !write(*,1) 'alpha_other', alpha_other
         !write(*,1) 'H_limit', H_limit
         !write(*,2) 'h1', h1
         !write(*,2) 's% nz', s% nz
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
            !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k), 
         end do
         !stop
      end subroutine alpha_mlt_routine

      

      subroutine TSF(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k,j,op_err,nsmooth,nsmootham
	      real(dp) :: alpha,shearsmooth,nu_tsf,nu_tsf_t,omegac,omegag,omegaa,omegat,difft,diffm,brunts,bruntsn2

         call star_ptr(id,s,ierr)
         if (ierr /= 0) return

         nsmooth=1!
         nsmootham=1
         shearsmooth=1d-30
               op_err = 0
         alpha = 2.5d-1

         !Calculate shear at each zone, then calculate TSF torque
         do k=nsmooth+1,s% nz-(nsmooth+1)

            nu_tsf=1d-30
            nu_tsf_t=1d-30
            !Calculate smoothed shear, q= dlnOmega/dlnr
            shearsmooth = s% omega_shear(k)/(2.*nsmooth+1.)
            do j=1,nsmooth
               shearsmooth = shearsmooth + (1./(2.*nsmooth+1.))*( s% omega_shear(k-j) + s% omega_shear(k+j) )
            end do

                     diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) !Magnetic diffusivity
            difft = 16d0*5.67d-5*(s% T(k))**3/(3d0*s% opacity(k)*(s% rho(k))**2*s% Cp(k)) !Thermal diffusivity
            omegaa = s% omega(k)*(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))))**(1./3.) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(omegaa*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)+(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+(s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*abs(shearsmooth*s% omega(k)/brunts)**(1./3.) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0.) then
               if (pow2(brunts) > 2.*pow2(shearsmooth)*pow2(s% omega(k))) then
                  omegac = 1d0*s% omega(k)*((brunts/s% omega(k))**0.5)*(diffm/(pow2(s% r(k))*s% omega(k)))**0.25  !Critical field strength
                  !nu_tsf = 5d-1+5d-1*tanh(10d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                  !nu_tsf = nu_tsf*alpha**3*s% omega(k)*pow2(s% r(k))*(s% omega(k)/brunts)**2 !nu_omega for revised Tayler instability
                  nu_tsf = alpha**3*pow2(s% r(k))*s% omega(k)*(s% omega(k)/brunts)**2 !nu_omega for revised Tayler instability

               end if
               ! Add TSF enabled by thermal diffusion
               if (pow2(brunts) < 2.*pow2(shearsmooth)*pow2(s% omega(k))) then
                  nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
               end if
               s% am_nu_omega(k) = s% am_nu_omega(k) + max(nu_tsf,nu_tsf_t) + 1d-1
            end if
         end do

      end subroutine TSF

      


      ! integer function how_many_extra_profile_columns(id, id_extra)
      !    use star_def, only: star_info
      !    integer, intent(in) :: id, id_extra
      !    integer :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      !    how_many_extra_profile_columns = 8
      ! end function how_many_extra_profile_columns
      
      
      ! subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
      !    use star_def, only: star_info, maxlen_profile_column_name
      !    use const_def, only: dp
      !    integer, intent(in) :: id, id_extra, n, nz
      !    character (len=maxlen_profile_column_name) :: names(n)

      !    real(dp) :: vals(nz,n)
      !    real(dp) :: brunts,alpha,nu_tsf,nu_tsf_t,diffm,difft,omegac,omegaa,omegat,lognuomega,bruntsn2,shearsmooth
	   !    real(dp) :: xdiff,bruntt,bruntsold
      !    real(dp) :: xmagfmu, xmagft, xmagfdif, xmagfnu, &
      !       xkap, xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig, &
      !       xeta

      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    integer :: k,j,op_err,nsmooth,nsmootham
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
         
      !    !note: do NOT add the extra names to profile_columns.list
      !    ! the profile_columns.list is only for the built-in profile column options.
      !    ! it must not include the new column names you are adding here.

      !    !if (n /= 1) stop 'data_for_extra_profile_columns'
      !    names(1) = 'log_N' !log Brunt
      !    names(2) = 'q_smooth' !dimensionless shear, dlnOmega/dlnr
      !    names(3) = 'nu_TSF' !new AM diffusivity
      !    names(4) = 'log_nu_TSF'
      !    names(5) = 'diff_mag' !Magnetic diffusivity
      !    names(6) = 'log_diff_mag'
      !    names(7) = 'J_in' !Angular momentum contained within
      !    names(8) = 'log_Nef' !Effective Brunt

      !    vals(:,1) = 1d-50
      !    vals(:,2) = 1d-30
      !    vals(:,3) = 1d-50
      !    vals(:,4) = -50
      !    vals(:,5) = 1d-50
      !    vals(:,6) = -50
      !    vals(:,7) = 1d-50
      !    vals(:,8) = -50

      !          op_err = 0
      !    alpha = 2.5d-1
      !    nsmooth=1
      !    nsmootham=1	 
      !    do k=nsmooth+1,s% nz-(nsmooth+1)

      !       vals(k,1) = safe_log10(sqrt(abs(s% brunt_N2(k))))
      !       vals(k,2) = s% omega_shear(k)/(2.*nsmooth+1.)
      !       nu_tsf=1d-30
      !       nu_tsf_t=1d-30

      !       do j=1,nsmooth
      !          vals(k,2) = vals(k,2) + (1./(2.*nsmooth+1.))*( s% omega_shear(k-j) + s% omega_shear(k+j) )
      !       end do
      !       shearsmooth = vals(k,2)

      !                diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
      !       difft = 16d0*5.67d-5*(s% T(k))**3/(3d0*s% opacity(k)*(s% rho(k))**2*s% Cp(k))
      !       vals(k,5) = diffm
      !       vals(k,6) =  safe_log10(vals(k,5))

      !       omegaa = s% omega(k)*(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))))**(1./3.) !Alfven frequency at saturation, assuming adiabatic instability
      !       omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
      !       brunts = sqrt(abs( s% brunt_N2_composition_term(k)+(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
      !       bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+(s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
      !       brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
      !       brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
      !       omegaa = s% omega(k)*abs(shearsmooth*s% omega(k)/brunts)**(1./3.) !Recalculate omegaa

      !       ! Calculate nu_TSF
      !       if (s% brunt_N2(k) > 0.) then
      !          if (pow2(brunts) > 2.*pow2(vals(k,2))*pow2(s% omega(k))) then
      !             omegac = 1d0*s% omega(k)*((brunts/s% omega(k))**0.5)*(diffm/(pow2(s% r(k))*s% omega(k)))**0.25  !Critical field strength
      !             !nu_tsf = 5d-1+5d-1*tanh(10d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
      !             !nu_tsf = nu_tsf*alpha**3*s% omega(k)*pow2(s% r(k))*(s% omega(k)/brunts)**2 !nu_omega for revised Tayler instability
      !             nu_tsf = alpha**3*s% omega(k)*pow2(s% r(k))*(s% omega(k)/brunts)**2 !nu_omega for revised Tayler instability
      !          end if
      !          ! Add TSF enabled by thermal diffusion
      !          if (pow2(brunts) < 2.*pow2(vals(k,2))*pow2(s% omega(k))) then
      !             nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
      !          end if
      !          vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
      !          vals(k,4) = safe_log10(vals(k,3))
      !          vals(k,8) = safe_log10(brunts)
      !       end if
      !    end do
      !    do k=s% nz-1,1,-1
      !       vals(k,7) = vals(k+1,7) + s% j_rot(k)*s% dm_bar(k) 
      !    end do

      ! end subroutine data_for_extra_profile_columns


      ! ! integer function how_many_extra_history_columns(id, id_extra)
      ! !    integer, intent(in) :: id, id_extra
      ! !    integer :: ierr
      ! !    type (star_info), pointer :: s
      ! !    ierr = 0
      ! !    call star_ptr(id, s, ierr)
      ! !    if (ierr /= 0) return
      ! !    how_many_extra_history_columns = 9
      ! ! end function how_many_extra_history_columns
      
      
      ! subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
      !    integer, intent(in) :: id, id_extra, n
      !    character (len=maxlen_history_column_name) :: names(n)
      !    real(dp) :: vals(n)
      !    real(dp) :: ominteg,ninteg,brunts,dr,J1,J2,J3,JHe,JSi,JFe,nu_pulse
      !    integer :: j,k
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return


      !    names(1) = 'nu_max' !Aseteroseismic nu_max
      !    names(2) = 'om_g' !Core rotation as sensed by g-modes	
      !    names(3) = 'AM_1p4' !Angular momentum in inner 0.1 Msun
      !    names(4) = 'AM_1p6' !Angular momentum in inner 0.2 Msun
      !    names(5) = 'AM_1p8' !Angular momentum in inner 0.3 Msun
      !    names(6) = 'AM_He' !Angular momentum in helium core
      !    names(7) = 'AM_Si' !Angular momentum in silicon core
      !    names(8) = 'AM_Fe' !Angular momentum in iron core
      !    names(9) = 'mom_inertia' !Angular momentum in iron core

      !    vals(1) = 1d-50
      !    vals(2) = 1d-50
      !    vals(3) = 1d-50
      !    vals(4) = 1d-50
      !    vals(5) = 1d-50
      !    vals(6) = 1d-50
      !    vals(7) = 1d-50
      !    vals(8) = 1d-50
      !    vals(9) = 1d-50

      !    vals(1) = s% nu_max
      !    nu_pulse = s% nu_max !Set pulsation period via scaling relations
      !    if (s% log_surface_radius < -1d0) nu_pulse = 1000. !If white dwarf, set pulsation period = 1000 microhertz

      !    ominteg=1d-99
      !    ninteg=1d-50
      !    J1=1d-50
      !    J2=1d-50
      !    J3=1d-50
      !    JHe=1d-50
      !    JSi=1d-50
      !    JFe=1d-50
      !    do k = 2, s% nz-1

      !       if (s% brunt_N2(k) < 1d-14) brunts = 1d-14 
      !       if (s% brunt_N2(k) > 1d-14) brunts = s% brunt_N2(k)

      !       if (2.*pi*nu_pulse/1d6<sqrt(brunts)) then
      !       if (2.*pi*nu_pulse/1d6<sqrt(2.)*s% csound(k)/s% r(k)) then
      !          dr = (s% r(k+1)-s% r(k-1))/2.
      !          ominteg = ominteg + sqrt(brunts)*(dr/s% r(k))*s% omega(k)
      !             ninteg = ninteg + sqrt(brunts)*(dr/s% r(k))
      !       end if
      !       end if

      !       if (s% m(k)/1.99d33 < 1.4d0) then
      !       J1 = J1 + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if 

      !       if (s% m(k)/1.99d33 < 1.6d0) then
      !       J2 = J2 + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if

      !       if (s% m(k)/1.99d33 < 1.8d0) then
      !       J3 = J3 + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if  

      !       if (s% m(k)/1.99d33 < s% he_core_mass) then
      !       JHe = JHe + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if 
      !       if (s% m(k)/1.99d33 < s% si_core_mass) then
      !       JSi = JSi + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if 
      !       if (s% m(k)/1.99d33 < s% fe_core_mass) then
      !       JFe = JFe + (2./3.)*pow2(s% r(k))*s% dm(k)*s% omega(k)
      !       end if 

      !    end do

      !    vals(2) = ominteg/ninteg
      !    vals(3) = J1
      !    vals(4) = J2
      !    vals(5) = J3
      !    vals(6) = JHe
      !    vals(7) = JSi
      !    vals(8) = JFe
      !    vals(9) = dot_product(s% i_rot(:s% nz), s% dm_bar(:s%nz))

      !    !print*, s% nu_max, 2.*pi/(ominteg/ninteg)/86400.
      !    !print*, s% log_surface_radius

      ! end subroutine data_for_extra_history_columns




      real(dp) function diffmag(rho,T,abar,zbar,ierr)

         ! Written by S.-C. Yoon, Oct. 10, 2003
         ! Electrical conductivity according to Spitzer 1962
         ! See also Wendell et al. 1987, ApJ 313:284
         real(dp), intent(in) :: rho, T, abar, zbar
	      integer, intent(out) :: ierr
         real(dp) :: xmagfmu, xmagft, xmagfdif, xmagfnu, &
            xkap, xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig, &
            xeta

            xgamma = 0.2275d0*zbar*zbar*pow(rho*1.d-6/abar,1d0/3d0)*1.d8/T
            xlg = log10(xgamma)
            if (xlg < -1.5d0) then
               xsig1 = sige1(zbar,T,xgamma)
               xsig = xsig1
            else if (xlg >= -1.5d0 .and. xlg <= 0d0) then
               xxx = (xlg + 0.75d0)*4d0/3d0
               ffff = 0.25d0*(2d0-3d0*xxx + xxx*xxx*xxx)
               xsig1 = sige1(zbar,T,xgamma)
               xsig2 = sige2(T,rho,zbar,ierr)
               if (ierr /= 0) return
               xsig = (1d0-ffff)*xsig2 + ffff*xsig1
            else if (xlg > 0d0 .and. xlg < 0.5d0) then
               xsig2 = sige2(T,rho,zbar,ierr)
               if (ierr /= 0) return
               xsig = xsig2
            else if (xlg >= 0.5d0 .and. xlg < 1d0) then
               xxx = (xlg-0.75d0)*4d0
               ffff = 0.25d0*(2d0-3d0*xxx + xxx*xxx*xxx)
               xsig2 = sige2(T,rho,zbar,ierr)
               if (ierr /= 0) return
               xsig3 = sige3(zbar,T,xgamma)
               xsig = (1d0-ffff)*xsig3 + ffff*xsig2
            else
               xsig3 = sige3(zbar,T,xgamma)
               xsig = xsig3
            endif

            diffmag = 7.1520663d19/xsig ! magnetic diffusivity

      end function diffmag



      real(dp) function sige1(z,t,xgamma)
         ! Written by S.-C. Yoon, Oct. 10, 2003
         ! Electrical conductivity according to Spitzer 1962
         ! See also Wendell et al. 1987, ApJ 313:284
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: etan, xlambda,f
         if (t >= 4.2d5) then
            f = sqrt(4.2d5/t)
         else
            f = 1.
         end if
         xlambda = sqrt(3.*z*z*z)*pow(xgamma,-1.5d0)*f + 1.
         etan = 3.d11*z*log(xlambda)*pow(t,-1.5d0)             ! magnetic diffusivity
         etan = etan/(1.-1.20487*exp(-1.0576*pow(z,0.347044d0))) ! correction: gammae
         sige1 = clight*clight/(4.*pi*etan)                    ! sigma = c^2/(4pi*eta)
      end function sige1


      real(dp) function sige2(T,rho,zbar,ierr)
         ! writen by S.-C. YOON Oct. 10, 2003
         ! electrical conductivity using conductive opacity
         ! see Wendell et al. 1987 ApJ 313:284
         use kap_lib, only: kap_get_elect_cond_opacity
         real(dp), intent(in) :: t,rho,zbar
         integer, intent(out) :: ierr
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         call kap_get_elect_cond_opacity( &
            zbar, log10(rho), log10(T),  &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         sige2 = 1.11d9*T*T/(rho*kap)
      end function sige2


      real(dp) function sige3(z,t,xgamma)
         ! writen by S.-C. YOON Oct. 10, 2003
         ! electrical conductivity in degenerate matter,
         ! according to Nandkumar & Pethick (1984)
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: rme, rm23, ctmp, xi
         rme = 8.5646d-23*t*t*t*xgamma*xgamma*xgamma/pow5(z)  ! rme = rho6/mue
         rm23 = pow(rme,2d0/3d0)
         ctmp = 1.+ 1.018*rm23
         xi= sqrt(3.14159d0/3.)*log(z)/3. + 2.*log(1.32+2.33/sqrt(xgamma))/3.d0-0.484*rm23/ctmp
         sige3 = 8.630d21*rme/(z*ctmp*xi)
      end function sige3

      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      ! integer function extras_startup(id, restart, ierr)
      !    integer, intent(in) :: id
      !    logical, intent(in) :: restart
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      !    extras_startup = 0
      !    if (.not. restart) then
      !       call alloc_extra_info(s)
      !    else ! it is a restart
      !       call unpack_extra_info(s)
      !    end if
      ! end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      ! integer function extras_check_model(id, id_extra)
      !    integer, intent(in) :: id, id_extra
      !    integer :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      !    extras_check_model = keep_going         
      !    if (.false. .and. s% star_mass_h1 < 0.35d0) then
      !       ! stop when star hydrogen mass drops to specified level
      !       extras_check_model = terminate
      !       write(*, *) 'have reached desired hydrogen mass'
      !       return
      !    end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
      !    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      ! end function extras_check_model


      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      ! integer function extras_finish_step(id, id_extra)
      !    integer, intent(in) :: id, id_extra
      !    integer :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      !    extras_finish_step = keep_going
      !    call store_extra_info(s)

      !    ! to save a profile, 
      !       ! s% need_to_save_profiles_now = .true.
      !    ! to update the star log,
      !       ! s% need_to_update_history_now = .true.

      !    ! see extras_check_model for information about custom termination codes
      !    ! by default, indicate where (in the code) MESA terminated
      !    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      ! end function extras_finish_step
      
      
      ! subroutine extras_after_evolve(id, id_extra, ierr)
      !    integer, intent(in) :: id, id_extra
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      ! end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      ! subroutine alloc_extra_info(s)
      !    integer, parameter :: extra_info_alloc = 1
      !    type (star_info), pointer :: s
      !    call move_extra_info(s,extra_info_alloc)
      ! end subroutine alloc_extra_info
      
      
      ! subroutine unpack_extra_info(s)
      !    integer, parameter :: extra_info_get = 2
      !    type (star_info), pointer :: s
      !    call move_extra_info(s,extra_info_get)
      ! end subroutine unpack_extra_info
      
      
      ! subroutine store_extra_info(s)
      !    integer, parameter :: extra_info_put = 3
      !    type (star_info), pointer :: s
      !    call move_extra_info(s,extra_info_put)
      ! end subroutine store_extra_info
      
      
      ! subroutine move_extra_info(s,op)
      !    integer, parameter :: extra_info_alloc = 1
      !    integer, parameter :: extra_info_get = 2
      !    integer, parameter :: extra_info_put = 3
      !    type (star_info), pointer :: s
      !    integer, intent(in) :: op
         
      !    integer :: i, j, num_ints, num_dbls, ierr
         
      !    i = 0
      !    ! call move_int or move_flg    
      !    num_ints = i
         
      !    i = 0
      !    ! call move_dbl       
         
      !    num_dbls = i
         
      !    if (op /= extra_info_alloc) return
      !    if (num_ints == 0 .and. num_dbls == 0) return
         
      !    ierr = 0
      !    call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
      !    if (ierr /= 0) then
      !       write(*,*) 'failed in star_alloc_extras'
      !       write(*,*) 'alloc_extras num_ints', num_ints
      !       write(*,*) 'alloc_extras num_dbls', num_dbls
      !       stop 1
      !    end if
         
      !    contains
         
      !    subroutine move_dbl(dbl)
      !       real(dp) :: dbl
      !       i = i+1
      !       select case (op)
      !       case (extra_info_get)
      !          dbl = s% extra_work(i)
      !          !dbl = s% xtra1_array(i)
      !       case (extra_info_put)
      !          s% extra_work(i) = dbl
	   !     !s% xtra1_array(i) = dbl
      !       end select
      !    end subroutine move_dbl
         
      !    subroutine move_int(int)
      !       integer :: int
      !       i = i+1
      !       select case (op)
      !       case (extra_info_get)
      !          int = s% extra_iwork(i)
      !       case (extra_info_put)
      !          s% extra_iwork(i) = int
      !       end select
      !    end subroutine move_int
         
      !    subroutine move_flg(flg)
      !       logical :: flg
      !       i = i+1
      !       select case (op)
      !       case (extra_info_get)
      !          flg = (s% extra_iwork(i) /= 0)
      !       case (extra_info_put)
      !          if (flg) then
      !             s% extra_iwork(i) = 1
      !          else
      !             s% extra_iwork(i) = 0
      !          end if
      !       end select
      !    end subroutine move_flg
      
      ! end subroutine move_extra_info

   

      ! end module run_star_extras
      

      
      
      
      ! subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
      !    use star_def
         
      !    integer, intent(in) :: id
      !    real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
      !    ! NOTE: surface is outermost cell. not necessarily at photosphere.
      !    ! NOTE: don't assume that vars are set at this point.
      !    ! so if you want values other than those given as args,
      !    ! you should use values from s% xh(:,:) and s% xa(:,:) only.
      !    ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
      !    real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
      !    integer, intent(out) :: ierr

      !    integer :: h1, he4
      !    real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
      !       vink_wind, nieu_wind, hamann_wind
      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return

      !    L1 = Lsurf
      !    M1 = Msurf
      !    R1 = Rsurf
      !    T1 = Tsurf

      !    h1 = s% net_iso(ih1)
      !    he4 = s% net_iso(ihe4)
      !    Xs = s% xa(h1,1)
      !    Ys = s% xa(he4,1)
      !    ! Z=0.0142 is Z from Asplund et al. 2009
      !    !Z = 0.02 changed by SG in the below line
      !    Z_div_Z_solar = s% kap_rq% Zbase/0.02d0 !0.0142d0
      !    ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
      !    Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

      !    vink_wind = 0d0
      !    nieu_wind = 0d0
      !    hamann_wind = 0d0
      !    w = 0

      !    call eval_Vink_wind(vink_wind)
      !    call eval_Nieuwenhuijzen_wind(nieu_wind)
      !    call eval_Hamann_wind(hamann_wind)

      !    ! use 1/10 hamann
      !    hamann_wind = hamann_wind/10d0

      !    if (T1 < Teff_jump) then
      !       ! low T wind
      !       w = max(vink_wind, nieu_wind)
      !    else
      !       ! high T wind
      !       alfa = 0d0
      !       if (Xs > 0.7d0) then
      !          alfa = 1d0
      !       else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
      !          alfa = (Xs - 0.4d0)/0.3d0
      !       end if
      !       w = alfa * vink_wind + (1d0-alfa) * hamann_wind
      !    end if

      !    ierr = 0

      !    contains

      !    subroutine eval_Vink_wind(w)
      !       real(dp), intent(inout) :: w
      !       real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

      !       ! alfa = 1 for hot side, = 0 for cool side
      !       if (T1 > 27500d0) then
      !          alfa = 1
      !       else if (T1 < 22500d0) then
      !          alfa = 0
      !       else
      !          dT = 100d0
      !          if (T1 > Teff_jump + dT) then
      !             alfa = 1
      !          else if (T1 < Teff_jump - dT) then
      !             alfa = 0
      !          else
      !             alfa = (T1 - (Teff_jump - dT)) / (2*dT)
      !          end if
      !       end if
            
      !       if (alfa > 0) then ! eval hot side wind (eqn 24)
      !          vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
      !          vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
      !          logMdot = &
      !             - 6.697d0 &
      !             + 2.194d0*log10(L1/Lsun/1d5) &
      !             - 1.313d0*log10(M1/Msun/30) &
      !             - 1.226d0*log10(vinf_div_vesc/2d0) &
      !             + 0.933d0*log10(T1/4d4) &
      !             - 10.92d0*pow2(log10(T1/4d4)) &
      !             + 0.85d0*log10(Z_div_Z_solar)
      !          w1 = exp10(logMdot)
      !       else
      !          w1 = 0
      !       end if
            
      !       if (alfa < 1) then ! eval cool side wind (eqn 25)
      !          vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
      !          vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
      !          logMdot = &
      !             - 6.688d0 &
      !             + 2.210d0*log10(L1/Lsun/1d5) &
      !             - 1.339d0*log10(M1/Msun/30) &
      !             - 1.601d0*log10(vinf_div_vesc/2d0) &
      !             + 1.07d0*log10(T1/2d4) &
      !             + 0.85d0*log10(Z_div_Z_solar)
      !          w2 = exp10(logMdot)
      !       else
      !          w2 = 0
      !       end if
            
      !       w = alfa*w1 + (1 - alfa)*w2
            
      !    end subroutine eval_Vink_wind

      !    subroutine eval_Nieuwenhuijzen_wind(w)
      !       ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
      !       real(dp), intent(out) :: w
      !       real(dp) :: log10w
      !       include 'formats'
      !       log10w = -14.02d0 &
      !                +1.24d0*log10(L1/Lsun) &
      !                +0.16d0*log10(M1/Msun) &
      !                +0.81d0*log10(R1/Rsun) &
      !                +0.85d0*log10(Z_div_Z_solar)
      !       w = exp10(log10w)
      !    end subroutine eval_Nieuwenhuijzen_wind

      !    subroutine eval_Hamann_wind(w)
      !       ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
      !       real(dp), intent(out) :: w
      !       real(dp) :: log10w
      !       include 'formats'
      !       log10w = -11.95d0 &
      !                +1.5d0*log10(L1/Lsun) &
      !                -2.85d0*Xs &
      !                + 0.85d0*log10(Z_div_Z_solar)
      !       w = exp10(log10w)
      !    end subroutine eval_Hamann_wind

      ! end subroutine brott_wind
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      integer function extras_check_model(id)
         integer, intent(in) :: id
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      

      end module run_star_extras
      
