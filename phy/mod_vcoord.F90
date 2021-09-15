! ------------------------------------------------------------------------------
! Copyright (C) 2021 Mats Bentsen
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_vcoord
! ------------------------------------------------------------------------------
! This module contains parameter, variables and procedures related to the
! vertical coordinate.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_config, only: inst_suffix
   use mod_constants, only: g, epsil, spval, onem
   use mod_xc
   use mod_eos, only: sig
   use mod_state, only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu, pv
   use mod_hor3map, only: reconstruction_struct, remap_struct, &
                          hor3map_plm, hor3map_ppm, hor3map_monotonic, &
                          hor3map_non_oscillatory, &
                          hor3map_non_oscillatory_posdef, &
                          prepare_reconstruction, reconstruct, regrid, &
                          prepare_remapping, remap, &
                          hor3map_noerr, hor3map_errstr
   use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
   use mod_tracers, only: ntr, trc
#endif

   implicit none

   private

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      vcoord_type = 'isopyc_bulkml', &
      reconstruction_method = 'ppm', &
      density_limiting = 'monotonic', &
      tracer_limiting = 'monotonic', &
      velocity_limiting = 'monotonic'
   logical :: &
      density_pc_upper_bndr = .false., &
      density_pc_lower_bndr = .false., &
      tracer_pc_upper_bndr = .true., &
      tracer_pc_lower_bndr = .false., &
      velocity_pc_upper_bndr = .true., &
      velocity_pc_lower_bndr = .false.
   real(r8) :: &
      dpmin_surface  = 1.5_r8, &
      dpmin_interior = .1_r8

   ! Options derived from string options.
   integer :: &
      vcoord_type_tag, &
      reconstruction_method_tag, &
      density_limiting_tag, &
      tracer_limiting_tag, &
      velocity_limiting_tag

   ! Parameters:
   integer, parameter :: &
      isopyc_bulkml = 1, &        ! Vertical coordinate type: bulk surface mixed
                                  ! layer with isopycnic layers below.
      cntiso_hybrid = 2           ! Vertical coordinate type: Hybrid coordinate
                                  ! with pressure coordinates towards the
                                  ! surface and continuous isopycnal below.
   real(r8), parameter :: &
      bfsq_min      = 1.e-7_r8, & ! Minimum buoyancy frequency squared in
                                  ! monotonized potential density to be used in
                                  ! regridding [s-2].
      regrid_mval   = - 1.e33_r8  ! Missing value for regridding.


   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      sigmar     ! Reference potential density [g cm-3].

   public :: vcoord_type_tag, isopyc_bulkml, cntiso_hybrid, sigmar, &
             readnml_vcoord, inivar_vcoord, cntiso_hybrid_regrid_remap, &
             remap_velocity

contains

   subroutine readnml_vcoord
   ! ---------------------------------------------------------------------------
   ! Read variables in the namelist group 'vcoord' and resolve options.
   ! ---------------------------------------------------------------------------

      character(len = 80) :: nml_fname
      integer :: ios
      logical :: fexist

      namelist /vcoord/ &
         vcoord_type, reconstruction_method, &
         density_limiting, tracer_limiting, velocity_limiting, &
         density_pc_upper_bndr, density_pc_lower_bndr, &
         tracer_pc_upper_bndr, tracer_pc_lower_bndr, &
         velocity_pc_upper_bndr, velocity_pc_lower_bndr, &
         dpmin_surface, dpmin_interior

      ! Read variables in the namelist group 'vcoord'.
      if (mnproc == 1) then
         nml_fname = 'ocn_in'//trim(inst_suffix)
         inquire(file = nml_fname, exist = fexist)
         if (fexist) then
            open (unit = nfu, file = nml_fname, status = 'old', action = 'read')
         else
            nml_fname = 'limits'//trim(inst_suffix)
            inquire(file = nml_fname, exist = fexist)
            if (fexist) then
               open (unit = nfu, file = nml_fname, status = 'old', &
                     action = 'read')
            else
               write (lp,*) 'readnml_vcoord: could not find namelist file!'
               call xchalt('(readnml_vcoord)')
                      stop '(readnml_vcoord)'
            endif
         endif
         read (unit = nfu, nml = vcoord, iostat = ios)
         close (unit = nfu)
      endif
      call xcbcst(ios)
      if (ios /= 0) then
         if (mnproc == 1) &
            write (lp,*) 'readnml_vcoord: No vertical coordinate variable '//  &
                         'group found in namelist. Using defaults.'
      else
        call xcbcst(vcoord_type)
        call xcbcst(reconstruction_method)
        call xcbcst(density_limiting)
        call xcbcst(tracer_limiting)
        call xcbcst(velocity_limiting)
        call xcbcst(density_pc_upper_bndr)
        call xcbcst(density_pc_lower_bndr)
        call xcbcst(tracer_pc_upper_bndr)
        call xcbcst(tracer_pc_lower_bndr)
        call xcbcst(velocity_pc_upper_bndr)
        call xcbcst(velocity_pc_lower_bndr)
        call xcbcst(dpmin_surface)
        call xcbcst(dpmin_interior)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_vcoord: vertical coordinate variables:'
         write (lp,*) '  vcoord_type =            ', &
                      trim(vcoord_type)
         write (lp,*) '  reconstruction_method =  ', &
                      trim(reconstruction_method)
         write (lp,*) '  density_limiting =       ', &
                      trim(density_limiting)
         write (lp,*) '  tracer_limiting =        ', &
                      trim(tracer_limiting)
         write (lp,*) '  velocity_limiting =      ', &
                      trim(velocity_limiting)
         write (lp,*) '  density_pc_upper_bndr =  ', density_pc_upper_bndr
         write (lp,*) '  density_pc_lower_bndr =  ', density_pc_lower_bndr
         write (lp,*) '  tracer_pc_upper_bndr =   ', tracer_pc_upper_bndr
         write (lp,*) '  tracer_pc_lower_bndr =   ', tracer_pc_lower_bndr
         write (lp,*) '  velocity_pc_upper_bndr = ', velocity_pc_upper_bndr
         write (lp,*) '  velocity_pc_lower_bndr = ', velocity_pc_lower_bndr
         write (lp,*) '  dpmin_surface =          ', dpmin_surface
         write (lp,*) '  dpmin_interior =         ', dpmin_interior
      endif

      ! Resolve options.
      select case (trim(vcoord_type))
         case ('isopyc_bulkml')
            vcoord_type_tag = isopyc_bulkml
         case ('cntiso_hybrid')
            vcoord_type_tag = cntiso_hybrid
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: vcoord_type = ', &
                  trim(vcoord_type), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(reconstruction_method))
         case ('plm')
            reconstruction_method_tag = hor3map_plm
         case ('ppm')
            reconstruction_method_tag = hor3map_ppm
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: reconstruction_method = ', &
                  trim(reconstruction_method), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(density_limiting))
         case ('monotonic')
            density_limiting_tag = hor3map_monotonic
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: density_limiting = ', &
                  trim(density_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(tracer_limiting))
         case ('monotonic')
            tracer_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            tracer_limiting_tag = hor3map_non_oscillatory
         case ('non_oscillatory_posdef')
            tracer_limiting_tag = hor3map_non_oscillatory_posdef
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: tracer_limiting = ', &
                  trim(tracer_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(velocity_limiting))
         case ('monotonic')
            velocity_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            velocity_limiting_tag = hor3map_non_oscillatory
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: velocity_limiting = ', &
                  trim(velocity_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      
      ! Change units from [m] to [g cm-1 s-2] of depth interval variables.
      dpmin_surface = dpmin_surface*onem
      dpmin_interior = dpmin_interior*onem

   end subroutine readnml_vcoord

   subroutine inivar_vcoord
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k

   !$omp parallel do private(i, k)
      do j = 1 - nbdy, jj + nbdy
         do k = 1, kk
            do i = 1 - nbdy, ii + nbdy
               sigmar(i, j, k) = spval
            enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_vcoord

   subroutine cntiso_hybrid_regrid_remap(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      type(reconstruction_struct) :: rcs
      type(remap_struct) :: rms

      real(r8), dimension(kdm + 1) :: p_1d, prgrd_1d, sigmar_1d
      real(r8), dimension(kdm) :: temp_1d, saln_1d, sigma_1d
      real(r8) :: beta, sdpsum, smean, dpmin_max, dpmin, pku, pku_test, &
                  pmin, dpt, pt, ptu1, ptl1, ptu2, ptl2, w1, x
      integer :: i, j, k, l, kn, nt, ks, ke, kl, ku, errstat
      logical :: thin_layers, layer_added
#ifdef TRC
      real(r8), dimension(kdm, ntr) :: trc_1d
#endif

      ! Minimum potential density difference with respect to pressure for
      ! potential density to be used in regridding.
      beta = bfsq_min/(g*g)

      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
         
            ! Copy variables into 1D arrays.
            p_1d(1) = p(i, j, 1)
            do k = 1, kk
               kn = k + nn
               temp_1d(k) = temp(i, j, kn)
               saln_1d(k) = saln(i, j, kn)
               p_1d(k + 1) = p_1d(k) + dp(i, j, kn)
               sigma_1d(k) = sigma(i, j, kn)
               sigmar_1d(k) = sigmar(i, j, k)
#ifdef TRC
               do nt = 1, ntr
                 trc_1d(k, nt) = trc(i, j, kn, nt)
               enddo
#endif
            enddo
            sigmar_1d(kk + 1) = sigmar_1d(kk)

            ! Make sure potential density to be used in regridding is
            ! monotonically increasing with depth.
            kl = kk
            ku = kl - 1
            do while (ku > 0)
               thin_layers = p_1d(kl + 1) - p_1d(ku) < epsil
               if (thin_layers .or. &
                   sigma_1d(kl) - sigma_1d(ku) &
                   < .5_r8*beta*(p_1d(kl + 1) - p_1d(ku))) then
                  sdpsum = sigma_1d(ku)*(p_1d(ku + 1) - p_1d(ku)) &
                         + sigma_1d(kl)*(p_1d(kl + 1) - p_1d(kl))
                  if (.not. thin_layers) &
                     smean = sdpsum/(p_1d(kl + 1) - p_1d(ku))
                  do
                     layer_added = .false.
                     if (ku > 1) then
                        if (thin_layers) then
                           ku = ku - 1
                           sdpsum = sdpsum &
                                  + sigma_1d(ku)*(p_1d(ku + 1) - p_1d(ku))
                           thin_layers = p_1d(kl + 1) - p_1d(ku) < epsil
                           if (.not. thin_layers) &
                              smean = sdpsum/(p_1d(kl + 1) - p_1d(ku))
                           layer_added = .true.
                        else
                           if (smean - sigma_1d(ku - 1) &
                               < .5_r8*beta*(p_1d(kl + 1) - p_1d(ku - 1))) then
                              ku = ku - 1
                              sdpsum = sdpsum &
                                     + sigma_1d(ku)*(p_1d(ku + 1) - p_1d(ku))
                              smean = sdpsum/(p_1d(kl + 1) - p_1d(ku))
                              layer_added = .true.
                           endif
                        endif
                     endif
                     if (kl < kk) then
                        if (thin_layers) then
                           kl = kl + 1
                           sdpsum = sdpsum &
                                  + sigma_1d(kl)*(p_1d(kl + 1) - p_1d(kl))
                           thin_layers = p_1d(kl + 1) - p_1d(ku) < epsil
                           if (.not. thin_layers) &
                              smean = sdpsum/(p_1d(kl + 1) - p_1d(ku))
                           layer_added = .true.
                        else
                           if (sigma_1d(kl + 1) - smean &
                               < .5_r8*beta*(p_1d(kl + 2) - p_1d(ku))) then
                              kl = kl + 1
                              sdpsum = sdpsum &
                                     + sigma_1d(kl)*(p_1d(kl + 1) - p_1d(kl))
                              smean = sdpsum/(p_1d(kl + 1) - p_1d(ku))
                              layer_added = .true.
                           endif
                        endif
                     endif
                     if (.not. layer_added) exit
                  enddo
                  do k = ku, kl
                     sigma_1d(k) = smean &
                                 + .5_r8*beta*( p_1d(k ) + p_1d(k  + 1) &
                                              - p_1d(ku) - p_1d(kl + 1))
                  enddo
               endif
               kl = ku
               ku = kl - 1
            enddo

            ! Prepare reconstruction with current interface pressures.
            errstat = prepare_reconstruction(p_1d, reconstruction_method_tag, &
                                             rcs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

            ! Monotonically reconstruct potential density.
            errstat = reconstruct(rcs, sigma_1d, density_limiting_tag, &
                                  density_pc_upper_bndr, density_pc_lower_bndr)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

            ! On the basis of the reconstructed potential density, regrid
            ! interface pressures so interface potential densities match target
            ! values.
            errstat = regrid(rcs, sigmar_1d, prgrd_1d, regrid_mval)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

            ! Modify regridded interface pressures to ensure the water column is
            ! properly bounded.
            k = 1
            do
               ks = k
               if (prgrd_1d(k) /= regrid_mval) exit
               prgrd_1d(k) = p_1d(1)
               if (k > kk) exit
               k = k + 1
            enddo
            k = kk + 1
            do
               ke = k
               if (prgrd_1d(k) /= regrid_mval) exit
               prgrd_1d(k) = p_1d(kk + 1)
               if (k == 1) exit
               k = k - 1
            enddo
            prgrd_1d(1) = p_1d(1)
            prgrd_1d(kk + 1) = p_1d(kk + 1)

            ! If no regrid interface is found in the water column, try to place
            ! all water in the layer with potential density bounds that include
            ! the column mean potential density.
            if (ks == ke) then
               sdpsum = 0._r8
               do k = 1, kk
                  sdpsum = sdpsum + sigma_1d(k)*(p_1d(k + 1) - p_1d(k))
               enddo
               smean = sdpsum/(p_1d(kk + 1) - p_1d(1))
               ks = 2
               do while (ks <= kk)
                  if (smean < sigmar_1d(ks)) exit
                  ks = ks + 1
               enddo
               do k = ks, kk
                  prgrd_1d(k) = p_1d(kk + 1)
               enddo
               ke = ks - 1
            endif

            ! Modify interface pressures so that layer thicknesses are
            ! above a specified threshold.
            dpmin_max = (p_1d(kk + 1) - p_1d(1))/kk
            dpmin_max = dpmin_surface
            dpmin = min(dpmin_max, dpmin_surface, dpmin_interior)
            ks = max(2, ks)
            ke = min(kk, ke)
            k = ks
            do while (k <= ke)
               if (prgrd_1d(k + 1) - prgrd_1d(k) < dpmin) then
                  if (k == ke) then
                     prgrd_1d(k) = prgrd_1d(ke + 1)
                  else
                     ku = k
                     kl = k + 1
                     pku = .5_r8*(prgrd_1d(kl) + prgrd_1d(ku) - dpmin)
                     do
                        layer_added = .false.
                        kl = kl + 1
                        pku_test = ((pku - dpmin)*(kl - ku) + prgrd_1d(kl)) &
                                   /(kl - ku + 1)
                        if (pku_test + (kl - ku)*dpmin > prgrd_1d(kl)) then
                           if (kl == ke + 1) exit
                           pku = pku_test
                           layer_added = .true.
                        else
                           kl = kl - 1
                        endif
                        ku = ku - 1
                        pku_test = ((pku - dpmin)*(kl - ku) + prgrd_1d(ku)) &
                                   /(kl - ku + 1)
                        if (pku_test < prgrd_1d(ku)) then
                           if (ku == 1) exit
                           pku = pku_test
                           layer_added = .true.
                        else
                           ku = ku + 1
                        endif
                        if (.not. layer_added) exit
                     enddo
                     if     (ku == 1) then
                        do k = 2, kl
                           prgrd_1d(k) = min(prgrd_1d(ke + 1), &
                                             prgrd_1d(k - 1) + dpmin)
                        enddo
                        do k = kl + 1, ke
                           prgrd_1d(k) = &
                              min(prgrd_1d(ke + 1), &
                                  max(prgrd_1d(k), prgrd_1d(1) + dpmin*(k - 1)))
                        enddo
                     elseif (kl == ke + 1) then
                        do k = ku, kl
                           prgrd_1d(k) = prgrd_1d(ke + 1)
                        enddo
                     else
                        prgrd_1d(ku) = pku
                        do k = ku + 1, kl
                           prgrd_1d(k) = prgrd_1d(k - 1) + dpmin
                        enddo
                     endif
                     k = kl
                  endif
               endif
               k = k + 1
            enddo

            ! Modify regridded interface pressures to ensure that a minimum
            ! layer thickness towards the surface is maintained. A smooth
            ! transition between modified and unmodified interfaces is sought.
            dpmin = min(dpmin_max, dpmin_surface)
            dpt = dpmin
            do k = 2, ke
               pmin = p_1d(1) + dpmin*(k - 1)
               dpt = max(prgrd_1d(k + 1) - prgrd_1d(k), dpt)
               pt = max(prgrd_1d(k), pmin)
               ptu1 = pmin - dpt
               ptl1 = pmin + dpt
               ptu2 = pmin
               ptl2 = pmin + 2._r8*dpt
               w1 = min(1._r8,(prgrd_1d(k) - p_1d(1))/(pmin - p_1d(1)))
               if (prgrd_1d(k) > ptu1 .and. prgrd_1d(k) < ptl1) then
                  x = .5_r8*(prgrd_1d(k) - ptu1)/dpt
                  pt = pmin + dpt*x*x
               endif
               if (prgrd_1d(k + 1) > ptu2 .and. prgrd_1d(k + 1) < ptl2) then
                  x = .5_r8*(prgrd_1d(k + 1) - ptu2)/dpt
                  pt = w1*pt + (1._r8 - w1)*(pmin + dpt*x*x)
               endif
               prgrd_1d(k) = min(p_1d(ke + 1), pt)
            enddo
            
            ! Prepare remapping to layer structure with regridded interface
            ! pressures.
            errstat = prepare_remapping(rcs, prgrd_1d, rms)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

            ! Reconstruct and remap potential temperature.
            errstat = reconstruct(rcs, temp_1d, tracer_limiting_tag, &
                                  tracer_pc_upper_bndr, tracer_pc_lower_bndr)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif
            errstat = remap(rcs, rms, temp_1d)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

            ! Reconstruct and remap salinity.
            errstat = reconstruct(rcs, saln_1d, tracer_limiting_tag, &
                                  tracer_pc_upper_bndr, tracer_pc_lower_bndr)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif
            errstat = remap(rcs, rms, saln_1d)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_hybrid_regrid_remap)')
                      stop '(cntiso_hybrid_regrid_remap)'
            endif

#ifdef TRC
            ! Reconstruct and remap tracers.
            do nt = 1, ntr
               errstat = reconstruct(rcs, trc_1d(:, nt), tracer_limiting_tag, &
                                     tracer_pc_upper_bndr, tracer_pc_lower_bndr)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(cntiso_hybrid_regrid_remap)')
                         stop '(cntiso_hybrid_regrid_remap)'
               endif
               errstat = remap(rcs, rms, trc_1d(:, nt))
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(cntiso_hybrid_regrid_remap)')
                         stop '(cntiso_hybrid_regrid_remap)'
               endif
            enddo
#endif

            ! Update 3D arrays
            do k = 1, kk
               kn = k + nn
               temp(i, j, kn) = temp_1d(k)
               saln(i, j, kn) = saln_1d(k)
               dp(i, j, kn) = prgrd_1d(k + 1) - prgrd_1d(k)
               sigma(i, j, kn) = sig(temp_1d(k), saln_1d(k))
#ifdef TRC
               do nt = 1,ntr
                 trc(i, j, kn, nt) = trc_1d(k, nt)
               enddo
#endif
            enddo

         enddo
         enddo
      enddo

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'cntiso_hybrid_regrid_remap:'
         endif
         call chksummsk(dp(1 - nbdy, 1 - nbdy, k1n), ip, kk, 'dp')
         call chksummsk(temp(1 - nbdy, 1 - nbdy, k1n), ip, kk, 'temp')
         call chksummsk(saln(1 - nbdy, 1 - nbdy, k1n), ip, kk, 'saln')
         call chksummsk(sigma(1 - nbdy, 1 - nbdy, k1n), ip,  kk, 'sigma')
         call chksummsk(sigmar, ip,  kk, 'sigmar')
#ifdef TRC
         do nt = 1, ntr
            call chksummsk(trc(1-nbdy, 1-nbdy, k1n, nt), ip, kk, 'trc')
         enddo
#endif
      endif

   end subroutine cntiso_hybrid_regrid_remap

   subroutine remap_velocity(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      type(reconstruction_struct) :: rcs
      type(remap_struct) :: rms

      real(r8), dimension(kdm + 1) :: p_1d, prgrd_1d
      real(r8), dimension(kdm) :: u_1d, v_1d
      real(r8) :: q
      integer :: i, j, k, l, kn, errstat
#ifdef TRC
      real(r8), dimension(kdm, ntr) :: trc_1d
#endif

      call xctilr(dp(1 - nbdy, 1 - nbdy, k1n), 1, kk, 3, 3, halo_ps)

   !$omp parallel do private(k, kn, l, i)
      do j = -2, jj + 2
         do k = 1, kk
            kn = k + nn
            do l = 1, isp(j)
            do i = max(- 2, ifp(j, l)), min(ii + 2, ilp(j, l))
               p(i, j, k + 1) = p(i, j, k) + dp(i, j, kn)
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(k,kn,l,i,q)
      do j = - 1, jj + 2
         do k = 1, kk
            kn = k + nn
            do l = 1, isu(j)
            do i = max(- 1, ifu(j, l)), min(ii + 2, ilu(j, l))
               q = min(p(i, j, kk + 1), p(i - 1, j, kk + 1))
               dpu(i, j, kn) = &
                 .5_r8*( (min(q, p(i - 1, j, k + 1)) - min(q, p(i - 1, j, k))) &
                       + (min(q, p(i    , j, k + 1)) - min(q, p(i    , j, k))))
            enddo
            enddo
            do l = 1, isv(j)
            do i = max(- 1, ifv(j, l)), min(ii + 2, ilv(j, l))
               q = min(p(i, j, kk + 1), p(i, j - 1, kk + 1))
               dpv(i, j, kn) = &
                 .5_r8*( (min(q, p(i, j - 1, k + 1)) - min(q, p(i, j - 1, k))) &
                       + (min(q, p(i, j    , k + 1)) - min(q, p(i, j    , k))))
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      do j = 1, jj

         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
         
            ! Copy variables into 1D arrays.
            prgrd_1d(1) = pu(i, j, 1)
            do k = 1, kk
               kn = k + nn
               u_1d(k) = u(i, j, kn)
               p_1d(k) = pu(i, j, k)
               prgrd_1d(k + 1) = prgrd_1d(k) + dpu(i, j, kn)
            enddo
            p_1d(kk + 1) = pu(i, j, kk + 1)

            ! Prepare reconstruction with current interface pressures.
            errstat = prepare_reconstruction(p_1d, reconstruction_method_tag, &
                                             rcs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Prepare remapping to layer structure with regridded interface
            ! pressures.
            errstat = prepare_remapping(rcs, prgrd_1d, rms)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Reconstruct and remap u-component of velocity.
            errstat = reconstruct(rcs, u_1d, velocity_limiting_tag, &
                                  velocity_pc_upper_bndr, &
                                  velocity_pc_lower_bndr)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif
            errstat = remap(rcs, rms, u_1d)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Update 3D arrays
            do k = 1, kk
               kn = k + nn
               u(i, j, kn) = u_1d(k)
            enddo

         enddo
         enddo

         do l = 1, isv(j)
         do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
         
            ! Copy variables into 1D arrays.
            prgrd_1d(1) = pv(i, j, 1)
            do k = 1, kk
               kn = k + nn
               v_1d(k) = v(i, j, kn)
               p_1d(k) = pv(i, j, k)
               prgrd_1d(k + 1) = prgrd_1d(k) + dpv(i, j, kn)
            enddo
            p_1d(kk + 1) = pv(i, j, kk + 1)

            ! Prepare reconstruction with current interface pressures.
            errstat = prepare_reconstruction(p_1d, reconstruction_method_tag, &
                                             rcs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Prepare remapping to layer structure with regridded interface
            ! pressures.
            errstat = prepare_remapping(rcs, prgrd_1d, rms)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Reconstruct and remap v-component of velocity.
            errstat = reconstruct(rcs, v_1d, velocity_limiting_tag, &
                                  velocity_pc_upper_bndr, &
                                  velocity_pc_lower_bndr)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif
            errstat = remap(rcs, rms, v_1d)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Update 3D arrays
            do k = 1, kk
               kn = k + nn
               v(i, j, kn) = v_1d(k)
            enddo

         enddo
         enddo

      enddo

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'remap_velocity:'
         endif
         call chksummsk(dpu(1 - nbdy, 1 - nbdy, k1n), iu, kk, 'dpv')
         call chksummsk(dpv(1 - nbdy, 1 - nbdy, k1n), iv, kk, 'dpv')
         call chksummsk(u(1 - nbdy, 1 - nbdy, k1n), iu, kk, 'u')
         call chksummsk(v(1 - nbdy, 1 - nbdy, k1n), iv, kk, 'v')
      endif

   end subroutine remap_velocity

end module mod_vcoord
