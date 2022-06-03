! ------------------------------------------------------------------------------
! Copyright (C) 2021-2022 Mats Bentsen
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

subroutine cntiso_hybrid_forcing(m, n, mm, nn, k1m, k1n)
! ---------------------------------------------------------------------------
! Apply surface forcing to the water column.
! ---------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, spcifh, alpha0, onem, onemu
   use mod_xc
   use mod_eos, only: dsigdt0, dsigds0
   use mod_state, only: dp, temp, saln
   use mod_swabs, only: swbgal, swbgfc, swamxd
   use mod_forcing, only: surflx, sswflx, salflx, buoyfl, t_sw_nonloc
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   integer, intent(in) :: m, n, mm, nn, k1m, k1n

   real(r8) :: pres(kk+1)
   real(r8) :: cpi, pswamx, gaa, dsgdt, dsgds, lei, pswamxi, pswbot
   integer :: i, j, k, l, kswamx, kn

   ! Set some constants:
   cpi = 1._r8/spcifh      ! Multiplicative inverse of specific heat capacity.
   pswamx = swamxd*onem    ! Maximum pressure of shortwave absorption.
   gaa = g*alpha0*alpha0

!$omp parallel do private(l, i, dsgdt, dsgds, lei, pres, kswamx, k, kn, &
!$omp                     pswamxi, pswbot)
   do j = 1, jj
      do l = 1, isp(j)
      do i = max(1, ifp(j,l)), min(ii, ilp(j,l))

         ! Derivatives of potential density referenced at the surface.
         dsgdt = dsigdt0(temp(i,j,k1n), saln(i,j,k1n))
         dsgds = dsigds0(temp(i,j,k1n), saln(i,j,k1n))

         ! Compute surface buoyancy flux [cm2 s-3].
         buoyfl(i,j,1) = - (dsgdt*surflx(i,j)*cpi + dsgds*salflx(i,j))*gaa

         ! Compute shortwave penetration factors at layer interfaces.
         lei = 1._r8/(onem*swbgal(i,j))
         pres(1) = 0._r8
         kswamx = 1
         t_sw_nonloc(i,j,1) = 1._r8
         do k = 1, kk
            kn = k + nn
            pres(k+1) = pres(k) + dp(i,j,kn)
            if (dp(i,j,kn) > onemu) then
               t_sw_nonloc(i,j,k+1) = &
                  swbgfc(i,j)*exp( - lei*min(pswamx, pres(k+1)))
               kswamx = k
            else
               t_sw_nonloc(i,j,k+1) = t_sw_nonloc(i,j,k)
            endif
            if (pres(k+1) > pswamx) exit
         enddo

         ! Compute buoyancy flux at subsurface layer interfaces. Penetration
         ! factors are modified so that shortwave radiation destined to
         ! penetrate below the lowest model layer is evenly absorbed in the
         ! water column.
         pswamxi = 1._r8/min(pswamx, pres(kswamx+1))
         pswbot = t_sw_nonloc(i,j,kswamx+1)
         do k = kswamx+1, kk+1
            t_sw_nonloc(i,j,k) = 0._r8
            buoyfl(i,j,k) = 0._r8
         enddo
         do k = kswamx, 2, -1
            kn = k + nn
            if (dp(i,j,kn) > onemu) then
               t_sw_nonloc(i,j,k) = t_sw_nonloc(i,j,k) - pswbot*pres(k)*pswamxi
            else
               t_sw_nonloc(i,j,k) = t_sw_nonloc(i,j,k+1)
            endif
            buoyfl(i,j,k) = - dsgdt*t_sw_nonloc(i,j,k)*sswflx(i,j)*cpi*gaa
         enddo

      enddo
      enddo
   enddo
!$omp end parallel do

   if (csdiag) then
      if (mnproc == 1) then
         write (lp,*) 'cntiso_hybrid_forcing:'
      endif
      call chksummsk(buoyfl, ip, kk+1, 'buoyfl')
      call chksummsk(t_sw_nonloc, ip, kk+1, 't_sw_nonloc')
   endif

end subroutine cntiso_hybrid_forcing
