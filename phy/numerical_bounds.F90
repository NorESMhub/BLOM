! ------------------------------------------------------------------------------
! Copyright (C) 2020 Mats Bentsen
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

subroutine numerical_bounds
! ---------------------------------------------------------------------------
! Set various numerical bounds.
! ---------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, spval, L_mks2cgs
   use mod_time, only: baclin
   use mod_xc
   use mod_grid, only: scqx, scqy, scpx, scpy, scuy, scvx, scp2, depths
   use mod_diffusion, only: difmxp, difmxq
   use mod_utility, only: umax, vmax
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   real(r8) :: dx2, dy2, btdtmx, umaxmin, vmaxmin, umaxmax, vmaxmax
   integer :: i, j, l

   ! Determine upper bound of lateral diffusivity based on numerical stability
   ! concerns.
!$omp parallel do private(i, dx2, dy2)
   do j = 1 - nbdy, jj + nbdy
      do i = 1 - nbdy, ii + nbdy
         dx2 = scpx(i, j)*scpx(i, j)
         dy2 = scpy(i, j)*scpy(i, j)
         difmxp(i, j) = .9_r8*.5_r8*dx2*dy2 &
                        /max(1._r8,(dx2 + dy2)*(baclin + baclin))
         dx2 = scqx(i, j)*scqx(i, j)
         dy2 = scqy(i, j)*scqy(i, j)
         difmxq(i, j) = .9_r8*.5_r8*dx2*dy2 &
                        /max(1._r8,(dx2 + dy2)*(baclin + baclin))
      enddo
   enddo
!$omp end parallel do

   ! Estimate maximum barotropic time step.
   btdtmx = 86400._r8
!$omp parallel do private(l, i) reduction(min:btdtmx)
   do j = 1, jj
      do l = 1, isp(j)
      do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
         btdtmx = min(btdtmx, &
                      scpx(i, j)*scpy(i, j) &
                      /sqrt(g*depths(i, j)*L_mks2cgs*( scpx(i, j)*scpx(i, j) &
                                                   + scpy(i, j)*scpy(i, j))))
      enddo
      enddo
   enddo
!$omp end parallel do
   call xcmin(btdtmx)
   if (mnproc == 1) then
      write (lp, *) 'estimated max. barotropic time step:', btdtmx/sqrt(2._r8)
      call flush(lp)
   endif


   ! Set maximum velocities allowable ensuring stability of the upwind scheme.

   umaxmin = spval
   vmaxmin = spval
   umaxmax = 0._r8
   vmaxmax = 0._r8
!$omp parallel do private(l, i) &
!$omp reduction(min:umaxmin, vmaxmin) reduction(max:umaxmax, vmaxmax)
   do j = 1, jj
     do l = 1, isu(j)
     do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
       umax(i, j) = .9_r8*.125_r8*min(scp2(i - 1, j), scp2(i, j)) &
                    /(scuy(i, j)*baclin)
       umaxmin = min(umaxmin, umax(i, j))
       umaxmax = max(umaxmax, umax(i, j))
     enddo
     enddo
     do l = 1, isv(j)
     do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
       vmax(i, j) = .9_r8*.125_r8*min(scp2(i, j - 1), scp2(i, j)) &
                    /(scvx(i, j)*baclin)
       vmaxmin = min(vmaxmin, vmax(i, j))
       vmaxmax = max(vmaxmax, vmax(i, j))
     enddo
     enddo
   enddo
!$omp end parallel do

   call xctilr(umax, 1, 1, nbdy, nbdy, halo_us)
   call xctilr(vmax, 1, 1, nbdy, nbdy, halo_vs)

   call xcmin(umaxmin)
   call xcmax(umaxmax)
   call xcmin(vmaxmin)
   call xcmax(vmaxmax)

   if (mnproc == 1) then
      write (lp, *) 'min/max umax:', umaxmin, umaxmax
      write (lp, *) 'min/max vmax:', vmaxmin, vmaxmax
      call flush(lp)
   endif

   if (csdiag) then
      if (mnproc == 1) then
         write (lp, *) 'numerical_bounds:'
      endif
      call chksummsk(difmxp, ip, 1,'difmxp')
      call chksummsk(difmxq, iq, 1,'difmxq')
      call chksummsk(umax, iu, 1,'umax')
      call chksummsk(vmax, iv, 1,'vmax')
   endif

end subroutine numerical_bounds
