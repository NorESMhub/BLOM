! ------------------------------------------------------------------------------
! Copyright (C) 2004-2020 Mats Bentsen
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
   use mod_time, only: delt1
   use mod_xc
   use mod_eos, only: sig, dsigdt0, dsigds0
   use mod_state, only: dp, temp, saln, sigma
   use mod_swabs, only: swbgal, swbgfc, swamxd
   use mod_forcing, only: surflx, surrlx, sswflx, salflx, salrlx, buoyfl
   use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
   use mod_tracers, only: ntr, trc, trflx
#endif

   implicit none

   integer, intent(in) :: m, n, mm, nn, k1m, k1n

   real(r8) :: pradd, lei, pres, pswbas, pswup, pswlo, q
   integer :: i, j, k, l, kfmax, kn
#ifdef TRC
   integer :: nt
#endif

   ! Maximum pressure of shortwave radiation penetration.
   pradd = swamxd*onem

!$omp parallel do private(l, i, lei, pres, pswbas, pswup, kfmax, k, kn, pswlo, &
!$omp                     q  &
#ifdef TRC
!$omp                   , nt &
#endif
!$omp )
   do j = 1, jj
      do l = 1, isp(j)
      do i = max(1, ifp(j, l)), min(ii, ilp(j, l))

         ! Compute total buoyancy flux [cm2 s-3]
         buoyfl(i,j) = &
           - ( dsigdt0(temp(i, j, k1n),saln(i, j, k1n))*surflx(i,j)/spcifh &
             + dsigds0(temp(i, j, k1n),saln(i, j, k1n))*salflx(i,j)) &
             *g*alpha0*alpha0

         ! Modify temperature below top layer due to penetrating short-wave
         ! flux.
         lei = 1._r8/(onem*swbgal(i, j))
         pres = dp(i, j, k1n)
         pswbas = swbgfc(i, j)*exp( - lei*pres)
         pswup = pswbas
         kfmax = 1
         k = 2
         do while (k <= kk)
            kn = k + nn
            pres = pres + dp(i, j, kn)
            if (dp(i, j, kn) > onemu) then
               pswlo = swbgfc(i,j)*exp( - lei*min(pradd, pres))
               temp(i, j, kn) = temp(i, j, kn) &
                              - (pswup - pswlo)*sswflx(i, j)*delt1*g &
                                /(spcifh*dp(i, j, kn))
               pswup = pswlo
               kfmax = k
            endif
            k = k + 1
            if (pres > pradd) exit
         enddo

         ! Modify temperature and salinity in top layer due to surface heat and
         ! salt fluxes.
         q = delt1*g/dp(i, j, k1n)
         temp(i, j, k1n) = temp(i, j, k1n) &
                         - ( surflx(i, j ) - (pswbas - pswup)*sswflx(i, j) &
                           + surrlx(i, j))*q/spcifh
         saln(i, j, k1n) = saln(i, j, k1n) - (salflx(i,j) + salrlx(i,j))*q

#ifdef TRC
         ! Modify tracer content in top layer due to surface fluxes.
         do nt = 1, ntr
            trc(i, j, k1n, nt) = trc(i, j, k1n, nt) - trflx(nt, i, j)*q
         enddo
#endif

         ! Update potential density in modified layers.
         do k = 1, kfmax
            kn = k + nn
            sigma(i, j, kn) = sig(temp(i, j, kn), saln(i, j, kn))
         enddo

      enddo
      enddo
   enddo
!$omp end parallel do

   if (csdiag) then
      if (mnproc == 1) then
         write (lp,*) 'mxlayr:'
      endif
      call chksummsk(dp, ip, 2*kk, 'dp')
      call chksummsk(temp, ip, 2*kk, 'temp')
      call chksummsk(saln, ip, 2*kk, 'saln')
      call chksummsk(sigma, ip, 2*kk, 'sigma')
#ifdef TRC
      do nt=1,ntr
         call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
      enddo
#endif
   endif

end subroutine cntiso_hybrid_forcing
