! ------------------------------------------------------------------------------
! Copyright (C) 2021-2022 Mats Bentsen, Mehmet Ilicak
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

module mod_vdiff
! ------------------------------------------------------------------------------
! This module contains procedures for solving vertical diffusion equations.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
  use mod_constants, only: g, spcifh, alpha0, onem
   use mod_time, only: delt1
   use mod_xc
   use mod_eos, only: sig
   use mod_state, only: u, v, dp, dpu, dpv, temp, saln, sigma
   use mod_checksum, only: csdiag, chksummsk
   use mod_diffusion, only: Kvisc_m, Kdiff_t, Kdiff_s, t_ns_nonloc, s_nonloc
   use mod_forcing, only: surflx, sswflx, surrlx, salflx, salrlx, t_sw_nonloc
#ifdef TRC
   use mod_tracers, only: ntr, trc, trflx
#endif

   implicit none

   private

   real(r8), parameter :: &
      dpmin_vdiff  = 0.1_r8*onem

   public :: cntiso_hybrid_vdifft, cntiso_hybrid_vdiffm

contains

   subroutine cntiso_hybrid_vdifft(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm) :: dp_1d, temp_1d, saln_1d, &
                                  nut_1d, nus_1d, nutrc_1d
      real(r8), dimension(2:kdm) :: fpbase, fp, gam
      real(r8) :: cpi, dtg, c, bei, rhs
      integer :: i, j, k, l, kn, nt
#ifdef TRC
      real(r8), dimension(kdm, ntr) :: trc_1d
#endif

      cpi = 1._r8/spcifh    ! Multiplicative inverse of specific heat capacity.
      dtg = delt1*g
      c = g*g*delt1/(alpha0*alpha0)

      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
         
            ! Copy variables into 1D arrays.
            do k = 1, kk
               kn = k + nn
               dp_1d(k) = dp(i, j, kn)
               temp_1d(k) = temp(i, j, kn)
               saln_1d(k) = saln(i, j, kn)
               nut_1d(k) = Kdiff_t(i, j, k)
               nus_1d(k) = Kdiff_s(i, j, k)
#ifdef TRC
               do nt = 1, ntr
                 trc_1d(k, nt) = trc(i, j, kn, nt)
               enddo
               nutrc_1d(k) = Kdiff_t(i, j, k)
#endif
            enddo

            ! Vertical diffusion equations are solved by backward integration
            ! forming a tridiagonal set of equations:
            !
            !   - fp(k)*U(k-1) + (dp(k) + fp(k) + fp(k+1))*U(k) - fp(k+1)*U(k+1)
            !   = dp(k)*(u(k) + Q_nonloc(k))
            !
            ! Here u and U is the variable to be diffused at old and new
            ! time-level, respectively, and Q_nonloc is the divergence of
            ! non-local transport of surface flux.

            ! Diffusive interface fluxes, before multiplying with diffusivity.
            do k = 2, kk
               fpbase(k) = c/max(dpmin_vdiff, .5_r8*(dp_1d(k - 1) + dp_1d(k)))
            enddo

            ! Diffusion of potential temperature.
            do k = 2, kk
               fp(k) = nut_1d(k)*fpbase(k)
            enddo
            bei = 1._r8/(dp_1d(1) + fp(2))
            rhs = dp_1d(1)*temp_1d(1) &
                - ( (1._r8 - t_ns_nonloc(i,j,2))*(surflx(i,j) - sswflx(i,j)) &
                  + (1._r8 - t_sw_nonloc(i,j,2))*sswflx(i,j) &
                  + surrlx(i,j))*dtg*cpi
            temp_1d(1) = rhs*bei
            do k = 2, kk - 1
               gam(k) = - fp(k)*bei
               bei = 1._r8/(dp_1d(k) + fp(k)*(1._r8 + gam(k)) + fp(k + 1))
               rhs = dp_1d(k)*temp_1d(k) &
                   - ( (t_ns_nonloc(i,j,k) - t_ns_nonloc(i,j,k+1)) &
                       *(surflx(i,j) - sswflx(i,j)) &
                     + (t_sw_nonloc(i,j,k) - t_sw_nonloc(i,j,k+1)) &
                       *sswflx(i,j))*dtg*cpi
               temp_1d(k) = (rhs + fp(k)*temp_1d(k - 1))*bei
            enddo
            gam(kk) = - fp(kk)*bei
            bei = 1._r8/(dp_1d(kk) + fp(kk)*(1._r8 + gam(kk)))
            rhs = dp_1d(kk)*temp_1d(kk) &
                - ( (t_ns_nonloc(i,j,kk) - t_ns_nonloc(i,j,kk+1)) &
                    *(surflx(i,j) - sswflx(i,j)) &
                  + (t_sw_nonloc(i,j,kk) - t_sw_nonloc(i,j,kk+1)) &
                    *sswflx(i,j))*dtg*cpi
            temp_1d(kk) = (rhs + fp(kk)*temp_1d(kk - 1))*bei
            do k = kk - 1, 1, - 1
               temp_1d(k) = temp_1d(k) - gam(k + 1)*temp_1d(k + 1)
            enddo

            ! Diffusion of salinity.
            do k = 2, kk
               fp(k) = nus_1d(k)*fpbase(k)
            enddo
            bei = 1._r8/(dp_1d(1) + fp(2))
            rhs = dp_1d(1)*saln_1d(1) &
                - ((1._r8 - s_nonloc(i,j,2))*salflx(i,j) &
                  + salrlx(i,j))*dtg
            saln_1d(1) = rhs*bei
            do k = 2, kk - 1
               gam(k) = - fp(k)*bei
               bei = 1._r8/(dp_1d(k) + fp(k)*(1._r8 + gam(k)) + fp(k + 1))
               rhs = dp_1d(k)*saln_1d(k) &
                   - (s_nonloc(i,j,k) - s_nonloc(i,j,k+1))*salflx(i,j)*dtg
               saln_1d(k) = (rhs + fp(k)*saln_1d(k - 1))*bei
            enddo
            gam(kk) = - fp(kk)*bei
            bei = 1._r8/(dp_1d(kk) + fp(kk)*(1._r8 + gam(kk)))
            rhs = dp_1d(kk)*saln_1d(kk) &
                - (s_nonloc(i,j,kk) - s_nonloc(i,j,kk+1))*salflx(i,j)*dtg
            saln_1d(kk) = (rhs + fp(kk)*saln_1d(kk - 1))*bei
            do k = kk - 1, 1, - 1
               saln_1d(k) = saln_1d(k) - gam(k + 1)*saln_1d(k + 1)
            enddo

#ifdef TRC
            ! Diffusion of tracers.
            do k = 2, kk
               fp(k) = nutrc_1d(k)*fpbase(k)
            enddo
            bei = 1._r8/(dp_1d(1) + fp(2))
            do nt = 1, ntr
               rhs = dp_1d(1)*trc_1d(1,nt) &
                   - (1._r8 - s_nonloc(i,j,2))*trflx(nt,i,j)*dtg
               trc_1d(1, nt) = rhs*bei
            enddo
            do k = 2, kk - 1
               gam(k) = - fp(k)*bei
               bei = 1._r8/(dp_1d(k) + fp(k)*(1._r8 + gam(k)) + fp(k + 1))
               do nt = 1, ntr
                  rhs = dp_1d(k)*trc_1d(k,nt) &
                      - (s_nonloc(i,j,k) - s_nonloc(i,j,k+1))*trflx(nt,i,j)*dtg
                  trc_1d(k, nt) = (rhs + fp(k)*trc_1d(k - 1, nt))*bei
               enddo
            enddo
            gam(kk) = - fp(kk)*bei
            bei = 1._r8/(dp_1d(kk) + fp(kk)*(1._r8 + gam(kk)))
            do nt = 1, ntr
               rhs = dp_1d(kk)*trc_1d(kk,nt) &
                   - (s_nonloc(i,j,kk) - s_nonloc(i,j,kk+1))*trflx(nt,i,j)*dtg
               trc_1d(kk, nt) = (rhs + fp(kk)*trc_1d(kk - 1, nt))*bei
            enddo
            do k = kk - 1, 1, - 1
               do nt = 1, ntr
                  trc_1d(k, nt) = trc_1d(k, nt) - gam(k + 1)*trc_1d(k + 1, nt)
               enddo
            enddo
#endif

            ! Update 3D arrays
            do k = 1, kk
               kn = k + nn
               temp(i, j, kn) = temp_1d(k)
               saln(i, j, kn) = saln_1d(k)
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
            write (lp,*) 'cntiso_hybrid_vdifft:'
         endif
         call chksummsk(temp, ip, 2*kk, 'temp')
         call chksummsk(saln, ip, 2*kk, 'saln')
         call chksummsk(sigma, ip, 2*kk, 'sigma')
#ifdef TRC
         do nt = 1, ntr
            call chksummsk(trc(1-nbdy, 1-nbdy, 1, nt), ip, 2*kk, 'trc')
         enddo
#endif
      endif

   end subroutine cntiso_hybrid_vdifft

   subroutine cntiso_hybrid_vdiffm(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm) :: dp_1d, u_1d, v_1d, nuv_1d
      real(r8), dimension(2:kdm) :: fpbase, fp, gam
      real(r8) :: c, bei
      integer :: i, j, k, l, kn

      c = g*g*delt1/(alpha0*alpha0)

      call xctilr(Kvisc_m, 1, kk, 1, 1, halo_ps)

      do j = 1, jj

         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
         
            ! Copy variables into 1D arrays.
            do k = 1, kk
               kn = k + nn
               dp_1d(k) = dpu(i, j, kn)
               u_1d(k) = u(i, j, kn)
               nuv_1d(k) = .5_r8*(Kvisc_m(i-1, j, k) + Kvisc_m(i, j, k))
            enddo

            ! Vertical diffusion equations are solved by backward integration
            ! forming a tridiagonal set of equations:
            !
            !   - fp(k)*U(k-1) + (dp(k) + fp(k) + fp(k+1))*U(k) - fp(k+1)*U(k+1)
            !   = dp(k)*u(k)
            !
            ! Here u and U is the variable to be diffused at old and new
            ! time-level, respectively.

            ! Diffusive interface fluxes, before multiplying with diffusivity.
            do k = 2, kk
               fpbase(k) = c/max(dpmin_vdiff, .5_r8*(dp_1d(k - 1) + dp_1d(k)))
            enddo

            ! Diffusion of u-component of baroclinic velocity.
            do k = 2, kk
               fp(k) = nuv_1d(k)*fpbase(k)
            enddo
            bei = 1._r8/(dp_1d(1) + fp(2))
            u_1d(1) = dp_1d(1)*u_1d(1)*bei
            do k = 2, kk - 1
               gam(k) = - fp(k)*bei
               bei = 1._r8/(dp_1d(k) + fp(k)*(1._r8 + gam(k)) + fp(k + 1))
               u_1d(k) = (dp_1d(k)*u_1d(k) + fp(k)*u_1d(k - 1))*bei
            enddo
            gam(kk) = - fp(kk)*bei
            bei = 1._r8/(dp_1d(kk) + fp(kk)*(1._r8 + gam(kk)))
            u_1d(kk) = (dp_1d(kk)*u_1d(kk) + fp(kk)*u_1d(kk - 1))*bei
            do k = kk - 1, 1, - 1
               u_1d(k) = u_1d(k) - gam(k + 1)*u_1d(k + 1)
            enddo

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
            do k = 1, kk
               kn = k + nn
               dp_1d(k) = dpv(i, j, kn)
               v_1d(k) = v(i, j, kn)
               nuv_1d(k) = .5_r8*(Kvisc_m(i, j-1, k) + Kvisc_m(i, j, k))
            enddo

            ! Vertical diffusion equations are solved by backward integration
            ! forming a tridiagonal set of equations:
            !
            !   - fp(k)*U(k-1) + (dp(k) + fp(k) + fp(k+1))*U(k) - fp(k+1)*U(k+1)
            !   = dp(k)*u(k)
            !
            ! Here u and U is the variable to be diffused at old and new
            ! time-level, respectively.

            ! Diffusive interface fluxes, before multiplying with diffusivity.
            do k = 2, kk
               fpbase(k) = c/max(dpmin_vdiff, .5_r8*(dp_1d(k - 1) + dp_1d(k)))
            enddo

            ! Diffusion of v-component of baroclinic velocity.
            do k = 2, kk
               fp(k) = nuv_1d(k)*fpbase(k)
            enddo
            bei = 1._r8/(dp_1d(1) + fp(2))
            v_1d(1) = dp_1d(1)*v_1d(1)*bei
            do k = 2, kk - 1
               gam(k) = - fp(k)*bei
               bei = 1._r8/(dp_1d(k) + fp(k)*(1._r8 + gam(k)) + fp(k + 1))
               v_1d(k) = (dp_1d(k)*v_1d(k) + fp(k)*v_1d(k - 1))*bei
            enddo
            gam(kk) = - fp(kk)*bei
            bei = 1._r8/(dp_1d(kk) + fp(kk)*(1._r8 + gam(kk)))
            v_1d(kk) = (dp_1d(kk)*v_1d(kk) + fp(kk)*v_1d(kk - 1))*bei
            do k = kk - 1, 1, - 1
               v_1d(k) = v_1d(k) - gam(k + 1)*v_1d(k + 1)
            enddo

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
            write (lp,*) 'cntiso_hybrid_vdiffm:'
         endif
         call chksummsk(u, iu, 2*kk, 'u')
         call chksummsk(v, iv, 2*kk, 'v')
      endif

   end subroutine cntiso_hybrid_vdiffm

end module mod_vdiff
