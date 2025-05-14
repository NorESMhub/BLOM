! ------------------------------------------------------------------------------
! Copyright (C) 2025 Mats Bentsen
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

module mod_cppm
! ------------------------------------------------------------------------------
! This module contains routines for advection of layer pressure thickness and
! tracers using Compatible Piecewise Parabolic Method (CPPM). If monotonic
! limiting is selected, the method ensures monotonicity of layer thickness for
! non-divergent flow and monotonicity of tracers for all flows.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_time, only: nstep
   use mod_xc, only: idm, jdm, kdm, nbdy, ii, jj, ip, &
                     halo_ps, halo_us, halo_vs, halo_uv, halo_vv, &
                     lp, mnproc, nproc, jpr, nreg, itdm, i0, &
                     xctilr, xcstop
   use mod_grid, only: scp2i
   use mod_state, only: dp, temp, saln, uflx, vflx, utflx, vtflx, &
                        usflx, vsflx, p, cau, cav, pbu, pbv
   use mod_tracers, only: trc, ntr

   implicit none

   private

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      cppm_limiting = 'non_oscillatory' ! CPPM limiting. Valid methods:
                                        ! 'monotonic', 'non_oscillatory'.

   ! Options derived from string options.
   integer :: &
      cppm_limiting_tag

   integer, parameter :: &
      cppm_limiting_monotonic       = 1, & ! Monotonic limiting.
      cppm_limiting_non_oscillatory = 2    ! Non-oscillatory limiting.

   real(r8), parameter :: &
      c0 = 0._r8, c1 = 1._r8, c2 = 2._r8, c3 = 3._r8, c4 = 4._r8, c6 = 6._r8, &
      c7 = 7._r8, c1_2 =1._r8/2._r8, c3_2 = 3._r8/2._r8, c1_3 = 1._r8/3._r8, &
      c2_3 = 2._r8/3._r8, c1_4 = 1._r8/4._r8, c1_5 = 1._r8/5._r8, &
      c1_6 = 1._r8/6._r8, c5_6 = 5._r8/6._r8, c1_12 = 1._r8/12._r8, &
      c7_12 = 7._r8/12._r8, &
      dpeps = 1.e-12_r8

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      evc1i, evc2i, evc3i, evc4i, slci, srci, scci, d2mi
   real(r8), dimension(1-nbdy:jdm+nbdy,1-nbdy:idm+nbdy) :: &
      evc1j, evc2j, evc3j, evc4j, slcj, srcj, sccj, d2mj

   integer :: ntr_loc

   public :: cppm_limiting, init_cppm, cppm

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   pure subroutine edge_value_coefficients(stencil_mask, evc1, evc2, evc3, evc4)
   ! ---------------------------------------------------------------------------
   ! Set coefficients for edge value reconstruction, accounting for the presence
   ! of eligible cell mean values in the 4-cell stencil.
   ! ---------------------------------------------------------------------------

      integer, dimension(4), intent(in) :: stencil_mask
      real(r8), intent(out) :: evc1, evc2, evc3, evc4

      if     (all(stencil_mask == [1, 1, 1, 1])) then
         evc1 = - c1_12
         evc2 =   c7_12
         evc3 =   c7_12
         evc4 = - c1_12
      elseif (all(stencil_mask == [1, 1, 1, 0])) then
         evc1 = - c1_6
         evc2 =   c5_6
         evc3 =   c1_3
         evc4 =   c0
      elseif (all(stencil_mask == [0, 1, 1, 1])) then
         evc1 =   c0
         evc2 =   c1_3
         evc3 =   c5_6
         evc4 = - c1_6
      elseif (all(stencil_mask == [0, 1, 1, 0])) then
         evc1 =   c0
         evc2 =   c1_2
         evc3 =   c1_2
         evc4 =   c0
      elseif (all(stencil_mask(1:2) == [1, 1])) then
         evc1 = - c1_2
         evc2 =   c3_2
         evc3 =   c0
         evc4 =   c0
      elseif (all(stencil_mask(3:4) == [1, 1])) then
         evc1 =   c0
         evc2 =   c0
         evc3 =   c3_2
         evc4 = - c1_2
      elseif (stencil_mask(2) == 1) then
         evc1 =   c0
         evc2 =   c1
         evc3 =   c0
         evc4 =   c0
      elseif (stencil_mask(3) == 1) then
         evc1 =   c0
         evc2 =   c0
         evc3 =   c1
         evc4 =   c0
      else
         evc1 =   c0
         evc2 =   c0
         evc3 =   c0
         evc4 =   c0
      endif

   end subroutine edge_value_coefficients

   pure subroutine slope_value_coefficients(stencil_mask, slc, src, scc)
   ! ---------------------------------------------------------------------------
   ! Set coefficients for left and right one-sided slopes and centered slopes,
   ! accounting for the presence of eligible cell mean values in the 3-cell
   ! stencil.
   ! ---------------------------------------------------------------------------

      integer, dimension(3), intent(in) :: stencil_mask
      real(r8), intent(out) :: slc, src, scc

      if     (any(stencil_mask == [0, 0, 0])) then
         slc = c0
         src = c0
         scc = c0
      else
         slc = c2
         src = c2
         scc = c1_2
      endif

   end subroutine slope_value_coefficients

   pure subroutine d2_mask(stencil_mask, d2m)
   ! ---------------------------------------------------------------------------
   ! Set mask for the value proportional to the second derivative of unlimited
   ! parabolas, accounting for the presence of eligible cell mean values in the
   ! 3-cell stencil.
   ! ---------------------------------------------------------------------------

      integer, dimension(3), intent(in) :: stencil_mask
      real(r8), intent(out) :: d2m

      if     (any(stencil_mask == [0, 0, 0])) then
         d2m = c0
      else
         d2m = c1
      endif

   end subroutine d2_mask

   pure subroutine edge_reconstruction(ijdm, ijs, ije, evc1, evc2, evc3, evc4, &
                                       hm, tm, hel, her, tel, ter)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge values using a 4th order scheme.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         evc1, evc2, evc3, evc4, hm
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: tel, ter

      real(r8) :: he, te
      integer :: i, nt

      do i = ijs, ije
         he = evc1(i)*hm(i-2) + evc2(i)*hm(i-1) &
            + evc3(i)*hm(i  ) + evc4(i)*hm(i+1)
         hel(i  ) = he
         her(i-1) = he      
         do nt = 1, ntr_loc
            te = evc1(i)*tm(nt,i-2) + evc2(i)*tm(nt,i-1) &
               + evc3(i)*tm(nt,i  ) + evc4(i)*tm(nt,i+1)
            tel(nt,i  ) = te
            ter(nt,i-1) = te
         enddo
      enddo

   end subroutine edge_reconstruction

   pure subroutine limit_non_oscillatory(ijdm, ijs, ije, slc, src, scc, d2m, &
                                         hm, tm, hel, her, tel, ter)
   ! ---------------------------------------------------------------------------
   ! Apply limiting of edge values to prevent a oscillatory reconstruction of
   ! piecewise parabolas.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         slc, src, scc, d2m, hm
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(inout) :: hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(inout) :: tel, ter

      real(r8), dimension(1-nbdy:ijdm+nbdy) :: d2h
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy) :: d2t

      real(r8) :: sl, sr, sc, d, q, r, a2
      integer :: i, nt

      do i = ijs-1, ije+1
         d2h(i) = d2m(i)*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            d2t(nt,i) = d2m(i)*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
         enddo
      enddo

      do i = ijs, ije

         if (d2h(i-1)*d2h(i  ) <= c0 .or. &
             d2h(i  )*d2h(i+1) <= c0) then
            sl = slc(i)*(hm(i  ) - hm(i-1))
            sr = src(i)*(hm(i+1) - hm(i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(hm(i+1) - hm(i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
                  hel(i) = hm(i) &
                         - sign(min(c1_2*abs(sc), &
                                    abs(hel(i) - hm(i))), sc)
               if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
                  her(i) = hm(i) &
                         + sign(min(c1_2*abs(sc), &
                                    abs(her(i) - hm(i))), sc)
               d = her(i) - hel(i)
               q = d*(c2*hm(i) - hel(i) - her(i))
               r = c1_3*d*d
               if     (  q > r) then
                  hel(i) = c3*hm(i) - c2*her(i)
               elseif (- r > q) then
                  her(i) = c3*hm(i) - c2*hel(i)
               endif
            else
               hel(i) = hm(i)
               her(i) = hm(i)
            endif
         endif

         hel(i) = max(hel(i), dpeps)
         her(i) = max(her(i), dpeps)
         sl = c2*(c3*hm(i) - c2*hel(i) - her(i))
         a2 = c3*(hel(i) - c2*hm(i) + her(i))
         sr = sl + c2*a2
         if (sl < c0 .and. sr > c0) then
            if (a2*hel(i) - c1_4*sl*sl < a2*dpeps) then
               q = c3*hm(i)/(c3*sl*sr + c4*a2*a2)
               hel(i) = sl*sl*q
               her(i) = sr*sr*q
            endif
         endif

         do nt = 1, ntr_loc
            if (d2t(nt,i-1)*d2t(nt,i  ) <= c0 .or. &
                d2t(nt,i  )*d2t(nt,i+1) <= c0) then
               sl = slc(i)*(tm(nt,i  ) - tm(nt,i-1))
               sr = src(i)*(tm(nt,i+1) - tm(nt,i  ))
               if (sl*sr > c0) then
                  sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
                  sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
                  if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                     tel(nt,i) = tm(nt,i) &
                               - sign(min(c1_2*abs(sc), &
                                          abs(tel(nt,i) - tm(nt,i))), sc)
                  if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                     ter(nt,i) = tm(nt,i) &
                               + sign(min(c1_2*abs(sc), &
                                          abs(ter(nt,i) - tm(nt,i))), sc)
                  d = ter(nt,i) - tel(nt,i)
                  q = d*(c2*tm(nt,i) - tel(nt,i) - ter(nt,i))
                  r = c1_3*d*d
                  if     (  q > r) then
                     tel(nt,i) = c3*tm(nt,i) - c2*ter(nt,i)
                  elseif (- r > q) then
                     ter(nt,i) = c3*tm(nt,i) - c2*tel(nt,i)
                  endif
               else
                  tel(nt,i) = tm(nt,i)
                  ter(nt,i) = tm(nt,i)
               endif
            endif
         enddo

         do nt = 2, ntr_loc
            tel(nt,i) = max(tel(nt,i), c0)
            ter(nt,i) = max(ter(nt,i), c0)
            sl = c2*(c3*tm(nt,i) - c2*tel(nt,i) - ter(nt,i))
            a2 = c3*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
            sr = sl + c2*a2
            if (sl < c0 .and. sr > c0) then
               if (a2*tel(nt,i) - c1_4*sl*sl < c0) then
                  q = c3*tm(nt,i)/(c3*sl*sr + c4*a2*a2)
                  tel(nt,i) = sl*sl*q
                  ter(nt,i) = sr*sr*q
               endif
            endif
         enddo

      enddo

   end subroutine limit_non_oscillatory

   pure subroutine limit_monotonic(ijdm, ijs, ije, slc, src, scc, &
                                   hm, tm, hel, her, tel, ter)
   ! ---------------------------------------------------------------------------
   ! Apply limiting of edge values to ensure a monotonic reconstruction of
   ! piecewise parabolas.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: slc, src, scc, hm
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(inout) :: hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(inout) :: tel, ter

      real(r8) :: sl, sr, sc, d, q, r
      integer :: i, nt

      do i = ijs, ije

         sl = slc(i)*(hm(i  ) - hm(i-1))
         sr = src(i)*(hm(i+1) - hm(i  ))
         if (sl*sr > c0) then
            sc = scc(i)*(hm(i+1) - hm(i-1))
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
            if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
               hel(i) = hm(i) &
                      - sign(min(c1_2*abs(sc), &
                                 abs(hel(i) - hm(i))), sc)
            if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
               her(i) = hm(i) &
                      + sign(min(c1_2*abs(sc), &
                                 abs(her(i) - hm(i))), sc)
            d = her(i) - hel(i)
            q = d*(c2*hm(i) - hel(i) - her(i))
            r = c1_3*d*d
            if     (  q > r) then
               hel(i) = c3*hm(i) - c2*her(i)
            elseif (- r > q) then
               her(i) = c3*hm(i) - c2*hel(i)
            endif
         else
            hel(i) = hm(i)
            her(i) = hm(i)
         endif

         do nt = 1, ntr_loc
            sl = slc(i)*(tm(nt,i  ) - tm(nt,i-1))
            sr = src(i)*(tm(nt,i+1) - tm(nt,i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                  tel(nt,i) = tm(nt,i) &
                            - sign(min(c1_2*abs(sc), &
                                       abs(tel(nt,i) - tm(nt,i))), sc)
               if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                  ter(nt,i) = tm(nt,i) &
                            + sign(min(c1_2*abs(sc), &
                                       abs(ter(nt,i) - tm(nt,i))), sc)
               d = ter(nt,i) - tel(nt,i)
               q = d*(c2*tm(nt,i) - tel(nt,i) - ter(nt,i))
               r = c1_3*d*d
               if     (  q > r) then
                  tel(nt,i) = c3*tm(nt,i) - c2*ter(nt,i)
               elseif (- r > q) then
                  ter(nt,i) = c3*tm(nt,i) - c2*tel(nt,i)
               endif
            else
               tel(nt,i) = tm(nt,i)
               ter(nt,i) = tm(nt,i)
            endif
         enddo

      enddo

   end subroutine limit_monotonic

   pure subroutine parabola_coeffs(ijdm, ijs, ije, &
                                   hm, tm, hel, her, tel, ter, &
                                   hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)
   ! ---------------------------------------------------------------------------
   ! Compute the coefficients of the piecewise parabolas.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: hm, hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm, tel, ter
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: &
         tpc0, tpc1, tpc2

      integer :: i, nt

      do i = ijs, ije
         hpc0(i) = hel(i)
         hpc1(i) = c6*hm(i) - c4*hel(i) - c2*her(i)
         hpc2(i) = c3*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            tpc0(nt,i) = tel(nt,i)
            tpc1(nt,i) = c6*tm(nt,i) - c4*tel(nt,i) - c2*ter(nt,i)
            tpc2(nt,i) = c3*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
         enddo
      enddo

   end subroutine parabola_coeffs

   pure subroutine flux_integration(ijdm, ijs, ije, &
                                    ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                    tpc0, tpc1, tpc2, &
                                    hf, htf)
   ! ---------------------------------------------------------------------------
   ! To obtain the fluxes, integrate the upstream piecewise parabolas over the
   ! flux area. For tracer fluxes, integrate the product of thickness and tracer
   ! parabolas. If the lower interface of the upstream cell (dl) is below the
   ! lowest edge interface (db), use a piecewise constant reconstruction of the
   ! upstream thickness (hb) set to the thickness of the upstream cell above db.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         ca, ai, db, du, dl, hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: &
         tpc0, tpc1, tpc2
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hf
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: htf

      real(r8) :: c, hb, p0, p1, p2, q1, q2, q3, q4
      integer :: i, nt

      do i = ijs, ije

         if (ca(i) < c0) then

           c = ca(i)*ai(i)

           if (dl(i) > db(i)) then

              hb = max(c0, db(i) - du(i))

              hf(i) = hb*ca(i)

              p0 =        hb
              p1 = - c1_2*hb*c
              p2 =   c1_3*hb*c*c

           else

              hf(i) = (hpc0(i) - (c1_2*hpc1(i) - c1_3*hpc2(i)*c)*c)*ca(i)

              p0 =         hpc0(i) - (c1_2*hpc1(i) - c1_3*hpc2(i)*c)*c 
              p1 = - (c1_2*hpc0(i) - (c1_3*hpc1(i) - c1_4*hpc2(i)*c)*c)*c
              p2 =   (c1_3*hpc0(i) - (c1_4*hpc1(i) - c1_5*hpc2(i)*c)*c)*c*c

           endif

           do nt = 1, ntr_loc
              htf(nt,i) = ( p0*tpc0(nt,i) &
                          + p1*tpc1(nt,i) &
                          + p2*tpc2(nt,i))*ca(i)
           enddo


         else

           c = ca(i)*ai(i-1)

           q1 = c1 - c1_2*c
           q2 = c1 - (c1 - c1_3*c)*c

           if (dl(i-1) > db(i)) then

              hb = max(c0, db(i) - du(i-1))

              hf(i) = hb*ca(i)

              p0 =    hb
              p1 = q1*hb
              p2 = q2*hb

           else

              hf(i) = (hpc0(i-1) + q1*hpc1(i-1) + q2*hpc2(i-1))*ca(i)

              q3 = c1_4*(c1 + c3*(c1 - c)*q2)
              q4 = c1_5*(c1 + c4*(c1 - c)*q3)
              p0 =    hpc0(i-1) + q1*hpc1(i-1) + q2*hpc2(i-1)
              p1 = q1*hpc0(i-1) + q2*hpc1(i-1) + q3*hpc2(i-1)
              p2 = q2*hpc0(i-1) + q3*hpc1(i-1) + q4*hpc2(i-1)

           endif

           do nt = 1, ntr_loc
              htf(nt,i) = ( p0*tpc0(nt,i-1) &
                          + p1*tpc1(nt,i-1) &
                          + p2*tpc2(nt,i-1))*ca(i)
           enddo

         endif

      enddo

   end subroutine flux_integration

   subroutine cppm_i_non_oscillatory(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: &
         tm, tel, ter, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 4, 0, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -3, idm + 4
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -3, idm + 4
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct edge values, apply limiting, define parabolas and
            ! obtain fluxes by integrating over the flux area.

            call edge_reconstruction(idm, -1, idm + 3, &
                                     evc1i(:,j), evc2i(:,j), &
                                     evc3i(:,j), evc4i(:,j), &
                                     hm, tm, hel, her, tel, ter)

            call limit_non_oscillatory(idm, 0, idm + 1, &
                                       slci(:,j), srci(:,j), scci(:,j), &
                                       d2mi(:,j), &
                                       hm, tm, hel, her, tel, ter)

            call parabola_coeffs(idm, 0, idm + 1, &
                                 hm, tm, hel, her, tel, ter, &
                                 hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_i_non_oscillatory

   subroutine cppm_j_non_oscillatory(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: &
         tm, tel, ter, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 0, 4, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -3, jdm + 4
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -3, jdm + 4
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct edge values, apply limiting, define parabolas and
            ! obtain fluxes by integrating over the flux area.

            call edge_reconstruction(jdm, -1, jdm + 3, &
                                     evc1j(:,i), evc2j(:,i), &
                                     evc3j(:,i), evc4j(:,i), &
                                     hm, tm, hel, her, tel, ter)

            call limit_non_oscillatory(jdm, 0, jdm + 1, &
                                       slcj(:,i), srcj(:,i), sccj(:,i), &
                                       d2mj(:,i), &
                                       hm, tm, hel, her, tel, ter)

            call parabola_coeffs(jdm, 0, jdm + 1, &
                                 hm, tm, hel, her, tel, ter, &
                                 hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_j_non_oscillatory

   subroutine cppm_i_monotonic(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying monotonic limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: &
         tm, tel, ter, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 3, 3, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -2, idm + 3
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -2, idm + 3
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct edge values, apply limiting, define parabolas and
            ! obtain fluxes by integrating over the flux area.

            call edge_reconstruction(idm, 0, idm + 2, &
                                     evc1i(:,j), evc2i(:,j), &
                                     evc3i(:,j), evc4i(:,j), &
                                     hm, tm, hel, her, tel, ter)

            call limit_monotonic(idm, 0, idm + 1, &
                                 slci(:,j), srci(:,j), scci(:,j), &
                                 hm, tm, hel, her, tel, ter)

            call parabola_coeffs(idm, 0, idm + 1, &
                                 hm, tm, hel, her, tel, ter, &
                                 hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_i_monotonic

   subroutine cppm_j_monotonic(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying monotonic limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: &
         tm, tel, ter, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 3, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 3, 3, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -2, jdm + 3
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -2, jdm + 3
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct edge values, apply limiting, define parabolas and
            ! obtain fluxes by integrating over the flux area.

            call edge_reconstruction(jdm, 0, jdm + 2, &
                                     evc1j(:,i), evc2j(:,i), &
                                     evc3j(:,i), evc4j(:,i), &
                                     hm, tm, hel, her, tel, ter)

            call limit_monotonic(jdm, 0, jdm + 1, &
                                 slcj(:,i), srcj(:,i), sccj(:,i), &
                                 hm, tm, hel, her, tel, ter)

            call parabola_coeffs(jdm, 0, jdm + 1, &
                                 hm, tm, hel, her, tel, ter, &
                                 hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_j_monotonic

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine init_cppm
   ! ---------------------------------------------------------------------------
   ! Resolve namelist options and set various coefficients and masks for CPPM.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         evc1j_perm, evc2j_perm, evc3j_perm, evc4j_perm, &
         slcj_perm, srcj_perm, sccj_perm, d2mj_perm

      integer, dimension(4) :: sm4
      integer, dimension(3) :: sm3
      real(r8) :: evc_tmp
      integer :: i, j

      ! Resolve namelist options.
      select case (trim(cppm_limiting))
         case ('monotonic')
            cppm_limiting_tag = cppm_limiting_monotonic
         case ('non_oscillatory')
            cppm_limiting_tag = cppm_limiting_non_oscillatory
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' init_cppm: cppm_limiting = ', trim(cppm_limiting), &
                  ' is unsupported!'
            call xcstop('(init_cppm)')
            stop '(init_cppm)'
      end select

      ! Set coefficients for edge value and slope reconstructions and mask for
      ! the value proportional to the second derivative of unlimited parabolas.

      evc1i(:,:) = c0
      evc2i(:,:) = c0
      evc3i(:,:) = c0
      evc4i(:,:) = c0
      slci(:,:) = c0
      srci(:,:) = c0
      scci(:,:) = c0
      d2mi(:,:) = c0
      evc1j_perm(:,:) = c0
      evc2j_perm(:,:) = c0
      evc3j_perm(:,:) = c0
      evc4j_perm(:,:) = c0
      slcj_perm(:,:) = c0
      srcj_perm(:,:) = c0
      sccj_perm(:,:) = c0
      d2mj_perm(:,:) = c0

      do j = 1, jj
         do i = 1, ii
            sm4 = ip(i-2:i+1,j)
            call edge_value_coefficients(sm4, &
                                         evc1i(i,j), evc2i(i,j), &
                                         evc3i(i,j), evc4i(i,j))
            sm3 = ip(i-1:i+1,j)
            call slope_value_coefficients(sm3, &
                                          slci(i,j), srci(i,j), scci(i,j))
            call d2_mask(sm3, d2mi(i,j))
            sm4 = ip(i,j-2:j+1)
            call edge_value_coefficients(sm4, &
                                         evc1j_perm(i,j), evc2j_perm(i,j), &
                                         evc3j_perm(i,j), evc4j_perm(i,j))
            sm3 = ip(i,j-1:j+1)
            call slope_value_coefficients(sm3, &
                                          slcj_perm(i,j), srcj_perm(i,j), &
                                          sccj_perm(i,j))
            call d2_mask(sm3, d2mj_perm(i,j))
         enddo
      enddo

      call xctilr(evc1i, 1, 1, nbdy, 0, halo_us)
      call xctilr(evc2i, 1, 1, nbdy, 0, halo_us)
      call xctilr(evc3i, 1, 1, nbdy, 0, halo_us)
      call xctilr(evc4i, 1, 1, nbdy, 0, halo_us)
      call xctilr(slci, 1, 1, nbdy, 0, halo_ps)
      call xctilr(srci, 1, 1, nbdy, 0, halo_ps)
      call xctilr(scci, 1, 1, nbdy, 0, halo_ps)
      call xctilr(d2mi, 1, 1, nbdy, 0, halo_ps)
      call xctilr(evc1j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(evc2j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(evc3j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(evc4j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(slcj_perm, 1, 1, 0, nbdy, halo_ps)
      call xctilr(srcj_perm, 1, 1, 0, nbdy, halo_ps)
      call xctilr(sccj_perm, 1, 1, 0, nbdy, halo_ps)
      call xctilr(d2mj_perm, 1, 1, 0, nbdy, halo_ps)

      ! With arctic patch, swap the order of of some the edge value
      ! coefficients.
      if (nreg == 2 .and. nproc == jpr) then
         j = jj
         do i = 1 - nbdy, ii + nbdy
            evc_tmp = evc1i(i,j)
            evc1i(i,j) = evc4i(i,j)
            evc4i(i,j) = evc_tmp
            evc_tmp = evc2i(i,j)
            evc2i(i,j) = evc3i(i,j)
            evc3i(i,j) = evc_tmp
         enddo
         do i = max(1, itdm/2 - i0 + 1), ii
            evc_tmp = evc1j_perm(i,j)
            evc1j_perm(i,j) = evc4j_perm(i,j)
            evc4j_perm(i,j) = evc_tmp
            evc_tmp = evc2j_perm(i,j)
            evc2j_perm(i,j) = evc3j_perm(i,j)
            evc3j_perm(i,j) = evc_tmp
         enddo
         do j = jj + 1, jj + nbdy
            do i = 1, ii
               evc_tmp = evc1j_perm(i,j)
               evc1j_perm(i,j) = evc4j_perm(i,j)
               evc4j_perm(i,j) = evc_tmp
               evc_tmp = evc2j_perm(i,j)
               evc2j_perm(i,j) = evc3j_perm(i,j)
               evc3j_perm(i,j) = evc_tmp
            enddo
         enddo
      endif

      do j = 1 - nbdy, jdm + nbdy
         do i = 1 - nbdy, idm + nbdy
            evc1j(j,i) = evc1j_perm(i,j)
            evc2j(j,i) = evc2j_perm(i,j)
            evc3j(j,i) = evc3j_perm(i,j)
            evc4j(j,i) = evc4j_perm(i,j)
            slcj(j,i) = slcj_perm(i,j)
            srcj(j,i) = srcj_perm(i,j)
            sccj(j,i) = sccj_perm(i,j)
            d2mj(j,i) = d2mj_perm(i,j)
         enddo
      enddo

      ! Local number of tracers where temperature and salinity is added to the
      ! ntr parameter.
      ntr_loc = 2 + ntr

   end subroutine init_cppm

   subroutine cppm(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! To complete a 2D transport step using CPPM, call 1D transport operators in
   ! alternate directional order (Strang splitting).
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      if (cppm_limiting_tag == cppm_limiting_non_oscillatory) then

         call xctilr(cau, 1, kdm, 4, 4, halo_uv)
         call xctilr(cav, 1, kdm, 4, 4, halo_vv)

         if (mod(nstep,2) == 1) then
            call cppm_i_non_oscillatory(m, n, mm, nn, k1m, k1n, &
                                        second_pass = .false.)
            call cppm_j_non_oscillatory(m, n, mm, nn, k1m, k1n, &
                                        second_pass = .true. )
         else
            call cppm_j_non_oscillatory(m, n, mm, nn, k1m, k1n, &
                                        second_pass = .false.)
            call cppm_i_non_oscillatory(m, n, mm, nn, k1m, k1n, &
                                        second_pass = .true. )
         endif

      else

         call xctilr(cau, 1, kdm, 3, 3, halo_uv)
         call xctilr(cav, 1, kdm, 3, 3, halo_vv)

         if (mod(nstep,2) == 1) then
            call cppm_i_monotonic(m, n, mm, nn, k1m, k1n, &
                                  second_pass = .false.)
            call cppm_j_monotonic(m, n, mm, nn, k1m, k1n, &
                                  second_pass = .true. )
         else
            call cppm_j_monotonic(m, n, mm, nn, k1m, k1n, &
                                  second_pass = .false.)
            call cppm_i_monotonic(m, n, mm, nn, k1m, k1n, &
                                  second_pass = .true. )
         endif

      endif

   end subroutine cppm

end module mod_cppm
