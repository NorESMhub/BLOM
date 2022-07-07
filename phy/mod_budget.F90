! ------------------------------------------------------------------------------
! Copyright (C) 2007-2022 Mats Bentsen
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

module mod_budget
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to budget computations.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, spcifh
   use mod_time, only: nstep, nstep1, delt1
   use mod_xc
   use mod_vcoord, only: vcoord_type_tag, isopyc_bulkml
   use mod_grid, only: scp2
   use mod_state, only: pb, dp, temp, saln
   use mod_forcing, only: surflx, surrlx, salflx, salrlx
   use mod_utility, only: util1, util2, util3, util4
#ifdef TRC
   use mod_tracers, only: itrtke, itrgls, trc, trflx
#endif

   implicit none

   private

   ! Constants.
   integer, parameter :: &
      ncalls = 7       ! Number of calls after which budgets are computed.
   logical :: &
      cnsvdi = .true.  ! Flag that indicates whether conservation diagnostics
                       ! are written.

   real(r8), dimension(ncalls, 2) :: &
      sdp, &           ! Global mass weighted sum of salinity.
      tdp              ! Global mass weighted sum of potential temperature.
#if defined(TRC) && defined(TKE)
   real(r8), dimension(ncalls, 2) :: &
      tkedp            ! Global mass weighted sum of TKE.
#  ifdef GLS
   real(r8), dimension(ncalls, 2) :: &
      glsdp            ! Global mass weighted sum of GLS.
#  endif
#endif
#ifdef TRC
   real(r8), dimension(ncalls, 2) :: &
      trdp             ! Global mass weighted sum of 1. tracer.
#endif

   real(r8) :: &
      mass0, &         ! Global sum of mass.
      sf, &            ! Global area weighted sum of surface salinity flux.
      tf               ! Global area weighted sum of surface salinity flux.
#ifdef TRC
   real(r8) :: &
      trf              ! Global area weighted sum of surface 1. tracer flux.
#endif

   public :: cnsvdi, budget_init, budget_sums, budget_output

contains

   subroutine budget_init
   ! ---------------------------------------------------------------------------
   ! Compute initial global mass sum.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

      if (.not. cnsvdi) return

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            util1(i, j) = pb(i, j, 1)*scp2(i, j)
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xcsum(mass0, util1, ips)

   end subroutine budget_init

   subroutine budget_sums(ncall, n, nn)
   ! ---------------------------------------------------------------------------
   ! Compute global mass weighted sums.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ncall, n, nn

      integer :: i, j, k, l, kn
      real(r8) :: q

      if (.not. cnsvdi) return

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            util1(i, j) = 0._r8
            util2(i, j) = 0._r8
#if defined(TRC) && defined(TKE)
            util3(i, j) = 0._r8
#  ifdef GLS
            util4(i, j) = 0._r8
#  endif
#endif
         enddo
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(k, kn, l, i, q)
      do j = 1,jj
         do k = 1,kk
            kn = k + nn
            do l = 1,isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               q = dp(i, j, kn)*scp2(i, j)
               util1(i, j) = util1(i, j) + saln(i, j, kn)*q
               util2(i, j) = util2(i, j) + temp(i, j, kn)*q
#if defined(TRC) && defined(TKE)
               util3(i, j) = util3(i, j) + trc(i, j, kn, itrtke)*q
#  ifdef GLS
               util4(i, j) = util4(i, j) + trc(i, j, kn, itrgls)*q
#  endif
#endif
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      call xcsum(sdp  (ncall, n), util1, ips)
      call xcsum(tdp  (ncall, n), util2, ips)
#if defined(TRC) && defined(TKE)
      call xcsum(tkedp(ncall, n), util3, ips)
#  ifdef GLS
      call xcsum(glsdp(ncall, n), util4, ips)
#  endif
#endif

#ifdef TRC
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            util1(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(k, kn, l, i, q)
      do j = 1, jj
         do k = 1, kk
            kn = k + nn
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               q = dp(i, j, kn)*scp2(i, j)
               util1(i, j) = util1(i, j) + trc(i, j, kn, 1)*q
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      call xcsum(trdp(ncall, n), util1, ips)
#endif

   end subroutine budget_sums

   subroutine budget_output(m)
   ! ---------------------------------------------------------------------------
   ! Output budgets.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m

      integer :: i, j, l

      if (.not.cnsvdi) return

      if (mnproc == 1 .and. nstep > nstep1 + 1) then

         if (vcoord_type_tag == isopyc_bulkml) then

            open (unit = nfu, file = 'salbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (sdp(2, m) - sdp(1, m))/mass0, &
               (sdp(3, m) - sdp(2, m))/mass0, &
               (sdp(4, m) - sdp(3, m))/mass0, &
               (sdp(5, m) - sdp(4, m) + sf*g)/mass0, &
               (sdp(6, m) - sdp(5, m))/mass0, &
               (sdp(7, m) - sdp(6, m))/mass0
            close (nfu)
            open (unit = nfu, file = 'tembud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (tdp(2, m) - tdp(1, m))/mass0, &
               (tdp(3, m) - tdp(2, m))/mass0, &
               (tdp(4, m) - tdp(3, m))/mass0, &
               (tdp(5, m) - tdp(4, m) + tf*g/spcifh)/mass0, &
               (tdp(6, m) - tdp(5, m))/mass0, &
               (tdp(7, m) - tdp(6, m))/mass0
            close (nfu)
#ifdef TRC
#  ifdef TKE
            open (unit = nfu, file = 'tkebud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (tkedp(2, m) - tkedp(1, m))/mass0, &
               (tkedp(3, m) - tkedp(2, m))/mass0, &
               (tkedp(4, m) - tkedp(3, m))/mass0, &
               (tkedp(5, m) - tkedp(4, m))/mass0, &
               (tkedp(6, m) - tkedp(5, m))/mass0, &
               (tkedp(7, m) - tkedp(6, m))/mass0
            close (nfu)
#    ifdef GLS
            open (unit = nfu, file = 'glsbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (glsdp(2, m) - glsdp(1, m))/mass0, &
               (glsdp(3, m) - glsdp(2, m))/mass0, &
               (glsdp(4, m) - glsdp(3, m))/mass0, &
               (glsdp(5, m) - glsdp(4, m))/mass0, &
               (glsdp(6, m) - glsdp(5, m))/mass0, &
               (glsdp(7, m) - glsdp(6, m))/mass0
            close (nfu)
#    endif
#  endif
            open (unit = nfu, file = 'trcbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (trdp(2, m) - trdp(1, m))/mass0, &
               (trdp(3, m) - trdp(2, m))/mass0, &
               (trdp(4, m) - trdp(3, m))/mass0, &
               (trdp(5, m) - trdp(4, m) + trf*g)/mass0, &
               (trdp(6, m) - trdp(5, m))/mass0, &
               (trdp(7, m) - trdp(6, m))/mass0
            close (nfu)
            open (unit = nfu, file = 'trcbudtot', position = 'append')
            write (nfu, '(i8,7e18.10)') nstep - 1, &
               trdp(1, m)/mass0, trdp(2, m)/mass0, trdp(3, m)/mass0, &
               trdp(4, m)/mass0, trdp(5, m)/mass0, trdp(6, m)/mass0, &
               trdp(7, m)/mass0
            close (nfu)
#endif

         else

            open (unit = nfu, file = 'salbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (sdp(2, m) - sdp(1, m))/mass0, &
               (sdp(3, m) - sdp(2, m))/mass0, &
               (sdp(4, m) - sdp(3, m) + sf*g)/mass0, &
               (sdp(5, m) - sdp(4, m))/mass0, &
               (sdp(6, m) - sdp(5, m))/mass0, &
               (sdp(7, m) - sdp(6, m))/mass0
            close (nfu)
            open (unit = nfu, file = 'tembud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (tdp(2, m) - tdp(1, m))/mass0, &
               (tdp(3, m) - tdp(2, m))/mass0, &
               (tdp(4, m) - tdp(3, m) + tf*g/spcifh)/mass0, &
               (tdp(5, m) - tdp(4, m))/mass0, &
               (tdp(6, m) - tdp(5, m))/mass0, &
               (tdp(7, m) - tdp(6, m))/mass0
            close (nfu)
#ifdef TRC
#  ifdef TKE
            open (unit = nfu, file = 'tkebud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (tkedp(2, m) - tkedp(1, m))/mass0, &
               (tkedp(3, m) - tkedp(2, m))/mass0, &
               (tkedp(4, m) - tkedp(3, m))/mass0, &
               (tkedp(5, m) - tkedp(4, m))/mass0, &
               (tkedp(6, m) - tkedp(5, m))/mass0, &
               (tkedp(7, m) - tkedp(6, m))/mass0
            close (nfu)
#    ifdef GLS
            open (unit = nfu, file = 'glsbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (glsdp(2, m) - glsdp(1, m))/mass0, &
               (glsdp(3, m) - glsdp(2, m))/mass0, &
               (glsdp(4, m) - glsdp(3, m))/mass0, &
               (glsdp(5, m) - glsdp(4, m))/mass0, &
               (glsdp(6, m) - glsdp(5, m))/mass0, &
               (glsdp(7, m) - glsdp(6, m))/mass0
            close (nfu)
#    endif
#  endif
            open (unit = nfu, file = 'trcbud', position = 'append')
            write (nfu, '(i8,6e12.4)') nstep - 1, &
               (trdp(2, m) - trdp(1, m))/mass0, &
               (trdp(3, m) - trdp(2, m))/mass0, &
               (trdp(4, m) - trdp(3, m) + trf*g)/mass0, &
               (trdp(5, m) - trdp(4, m))/mass0, &
               (trdp(6, m) - trdp(5, m))/mass0, &
               (trdp(7, m) - trdp(6, m))/mass0
            close (nfu)
            open (unit = nfu, file = 'trcbudtot', position = 'append')
            write (nfu, '(i8,7e18.10)') nstep - 1, &
               trdp(1, m)/mass0, trdp(2, m)/mass0, trdp(3, m)/mass0, &
               trdp(4, m)/mass0, trdp(5, m)/mass0, trdp(6, m)/mass0, &
               trdp(7, m)/mass0
            close (nfu)
#endif

         endif

      endif

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            util1(i, j) = (salflx(i, j) + salrlx(i, j))*scp2(i, j)*delt1
            util2(i, j) = (surflx(i, j) + surrlx(i, j))*scp2(i, j)*delt1
#ifdef TRC
            util3(i, j) = trflx(1, i, j)*scp2(i, j)*delt1
#endif
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xcsum(sf, util1, ips)
      call xcsum(tf, util2, ips)
#ifdef TRC
      call xcsum(trf, util3, ips)
#endif

   end subroutine budget_output

end module mod_budget
