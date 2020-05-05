! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen
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

      subroutine export_mct(o2x_o, lsize, perm, jjcpl, nsend, sbuff,
     .                      tlast_coupled)

      ! Uses modules
      
      use mct_mod
      use shr_const_mod, only: SHR_CONST_TKFRZ
      use types, only : r8
      use blom_cpl_indices
      use mod_xc

      implicit none

      ! Input/output arguments

      type (mct_aVect)         , intent(inout) :: o2x_o
      integer                  , intent(in)    :: lsize
      integer, dimension(lsize), intent(in)    :: perm
      integer                  , intent(in)    :: jjcpl
      integer                  , intent(in)    :: nsend
      real (r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nsend),
     .   intent(inout) :: sbuff
      real (r8)                , intent(inout) :: tlast_coupled

      ! Include files

#include "common_geo.h"
#include "common_forc.h"

      ! Local variables

      integer i, j, n
      real (r8) :: tfac, utmp, vtmp

      tfac = 1._r8/tlast_coupled

      ! ----------------------------------------------------------------
      ! Interpolate onto scalar points, rotate, and pack surface
      ! velocity (m/s) and surface gradient (m/m)
      ! ----------------------------------------------------------------

      call xctilr(sbuff(1-nbdy,1-nbdy,index_o2x_So_u), 
     .            1,1, 1,1, halo_uv)
      call xctilr(sbuff(1-nbdy,1-nbdy,index_o2x_So_v),
     .            1,1, 1,1, halo_vv)
      call xctilr(sbuff(1-nbdy,1-nbdy,index_o2x_So_dhdx),
     .            1,1, 1,1, halo_uv)
      call xctilr(sbuff(1-nbdy,1-nbdy,index_o2x_So_dhdy),
     .            1,1, 1,1, halo_vv)

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            utmp = .5_r8*( sbuff(i  ,j,index_o2x_So_u)
     .                   + sbuff(i+1,j,index_o2x_So_u))
            vtmp = .5_r8*( sbuff(i,j  ,index_o2x_So_v)
     .                   + sbuff(i,j+1,index_o2x_So_v))
            o2x_o%rattr(index_o2x_So_u,n) =
     .         (utmp*cosang(i,j) - vtmp*sinang(i,j))*tfac*1.e-2_r8
            o2x_o%rattr(index_o2x_So_v,n) =
     .         (utmp*sinang(i,j) + vtmp*cosang(i,j))*tfac*1.e-2_r8
            utmp = ( sbuff(i  ,j,index_o2x_So_dhdx)*iu(i  ,j)
     .             + sbuff(i+1,j,index_o2x_So_dhdx)*iu(i+1,j))
     .             /max(1,iu(i,j) + iu(i+1,j))
            vtmp = ( sbuff(i,j  ,index_o2x_So_dhdy)*iv(i,j  )
     .             + sbuff(i,j+1,index_o2x_So_dhdy)*iv(i,j+1))
     .             /max(1,iv(i,j) + iv(i,j+1))
            o2x_o%rAttr(index_o2x_So_dhdx,n) =
     .         (utmp*cosang(i,j) - vtmp*sinang(i,j))*tfac
            o2x_o%rAttr(index_o2x_So_dhdy,n) =
     .         (utmp*sinang(i,j) + vtmp*cosang(i,j))*tfac
         enddo
      enddo

      ! ----------------------------------------------------------------
      ! Pack temperature (K), salinity (psu), freezing/melting potential
      ! (W/m~2)
      ! ----------------------------------------------------------------

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            o2x_o%rAttr(index_o2x_So_t,n) =
     .         sbuff(i,j,index_o2x_So_t)*tfac + SHR_CONST_TKFRZ
            o2x_o%rAttr(index_o2x_So_s,n) =
     .         sbuff(i,j,index_o2x_So_s)*tfac
            if (sbuff(i,j,index_o2x_Fioo_q) > 0._r8) then
              o2x_o%rAttr(index_o2x_Fioo_q,n) =
     .           sbuff(i,j,index_o2x_Fioo_q)*tfac
            else
              o2x_o%rAttr(index_o2x_Fioo_q,n) =
     .           mltpot(i,j)*tfac
            endif
         enddo
      enddo

      ! ----------------------------------------------------------------
      ! Pack dms flux (kmol DMS/m^2/s), if requested
      ! ----------------------------------------------------------------

      if (index_o2x_Faoo_fdms_ocn > 0) then
        n = 0
         do j = 1, jjcpl
            do i = 1, ii
               n = n + 1
               o2x_o%rAttr(index_o2x_Faoo_fdms_ocn,n) =
     .            sbuff(i,j,index_o2x_Faoo_fdms_ocn)*tfac
            enddo
         enddo
      else
         if (mnproc.eq.1)
     .      write (lp,*) 'export_mct: dms flux not sent to coupler'
      endif

      ! ----------------------------------------------------------------
      ! Pack co2 flux (kg CO2/m^2/s), if requested
      ! ----------------------------------------------------------------

      if (index_o2x_Faoo_fco2_ocn > 0) then
        n = 0
         do j = 1, jjcpl
            do i = 1, ii
               n = n + 1
               o2x_o%rAttr(index_o2x_Faoo_fco2_ocn,n) = 
     .            sbuff(i,j,index_o2x_Faoo_fco2_ocn)*tfac
            enddo
         enddo
      else
         if (mnproc.eq.1)
     .      write (lp,*) 'export_mct: co2 flux not sent to coupler'
      endif

      tlast_coupled = 0._r8

      !-----------------------------------------------------------------
      ! permute in-place before returning the output state
      !-----------------------------------------------------------------

      call mct_aVect_permute(o2x_o, perm)

      end subroutine export_mct
