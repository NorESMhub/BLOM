! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen, Alok Kumar Gupta
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

      subroutine import_mct(x2o_o, lsize, perm, jjcpl)

      ! Uses modules
      
      use mct_mod
      use shr_const_mod, only: SHR_CONST_RHOSW, SHR_CONST_LATICE
      use types, only : r8
      use blom_cpl_indices
      use mod_xc

      implicit none

      ! Input/output arguments

      type (mct_aVect)         , intent(inout) :: x2o_o
      integer                  , intent(in)    :: lsize
      integer, dimension(lsize), intent(in)    :: perm
      integer                  , intent(in)    :: jjcpl

      ! Include files

#include "common_blocks.h"
#include "common_geo.h"
#include "common_cesm.h"

      ! Local parameters

      real (r8), parameter ::
     .   mval = -1.e12_r8,
     .   fval = -1.e13_r8

      ! Local variables

      integer i, j, l, n
      real (r8) :: utmp, vtmp

      !-----------------------------------------------------------------
      ! unpermute in-place before unpacking
      !-----------------------------------------------------------------

      call mct_aVect_unpermute(x2o_o, perm)

      ! Update time level indices

      if (l1ci == 1 .and. l2ci == 1) then
         l1ci = 2
         l2ci = 2
      else
         l1ci = l2ci
         l2ci = 3 - l2ci
      endif

      ! Unpack

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            if     (ip(i,j) == 0) then
               util1(i,j) = mval
               util2(i,j) = mval
               ustarw_da(i,j,l2ci) = mval
            elseif (cplmsk(i,j) == 0) then
               util1(i,j) = fval
               util2(i,j) = fval
               ustarw_da(i,j,l2ci) = fval
            else
               utmp = x2o_o%rAttr(index_x2o_Foxx_taux,n)
               vtmp = x2o_o%rAttr(index_x2o_Foxx_tauy,n)
               util1(i,j) =   utmp*cosang(i,j) + vtmp*sinang(i,j)
               util2(i,j) = - utmp*sinang(i,j) + vtmp*cosang(i,j)

               ! Friction velocity (m/s)
               ustarw_da(i,j,l2ci) = sqrt(sqrt(utmp*utmp+vtmp*vtmp)
     .                              /SHR_CONST_RHOSW)
            endif
         enddo
      enddo

      call fill_global(mval, fval, halo_pv, util1)
      call fill_global(mval, fval, halo_pv, util2)
      call fill_global(mval, fval, halo_ps,
     .                 ustarw_da(1-nbdy,1-nbdy,l2ci))

      call xctilr(util1, 1,1, 1,1, halo_pv)
      call xctilr(util2, 1,1, 1,1, halo_pv)

      do j=1,jj
         do l=1,isu(j)
         do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            ! x-component of momentum flux (kg/m/s^2)
            ztx_da(i,j,l2ci) = .5_r8*(util1(i-1,j) + util1(i,j))
         enddo
         enddo
         do l=1,isv(j)
         do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            ! y-component of momentum flux (kg/m/s^2)
            mty_da(i,j,l2ci) = .5_r8*(util2(i,j-1) + util2(i,j))
         enddo
         enddo
      enddo
      
      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            if     (ip(i,j) == 0) then
               lip_da(i,j,l2ci) = mval
               sop_da(i,j,l2ci) = mval
               eva_da(i,j,l2ci) = mval
               rnf_da(i,j,l2ci) = mval
               rfi_da(i,j,l2ci) = mval
               fmltfz_da(i,j,l2ci) = mval
               sfl_da(i,j,l2ci) = mval
               swa_da(i,j,l2ci) = mval
               nsf_da(i,j,l2ci) = mval
               hmlt_da(i,j,l2ci) = mval
               slp_da(i,j,l2ci) = mval
               ficem_da(i,j,l2ci) = mval
               abswnd_da(i,j,l2ci) = mval
            elseif (cplmsk(i,j) == 0) then
               lip_da(i,j,l2ci) = 0._r8
               sop_da(i,j,l2ci) = 0._r8
               eva_da(i,j,l2ci) = 0._r8
               rnf_da(i,j,l2ci) = 0._r8
               rfi_da(i,j,l2ci) = 0._r8
               fmltfz_da(i,j,l2ci) = 0._r8
               sfl_da(i,j,l2ci) = 0._r8
               swa_da(i,j,l2ci) = 0._r8
               nsf_da(i,j,l2ci) = 0._r8
               hmlt_da(i,j,l2ci) = 0._r8
               slp_da(i,j,l2ci) = fval
               ficem_da(i,j,l2ci) = fval
               abswnd_da(i,j,l2ci) = fval
            else

               ! Liquid water flux, positive downwards [kg/m^2/s]
               lip_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Faxa_rain,n)

               ! Solid precipitation, positive downwards [kg/m^2/s]
               sop_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Faxa_snow,n)

               ! Evaporation, positive downwards [kg/m^2/s]
               eva_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Foxx_evap,n)

               ! Liquid runoff, positive downwards [kg/m^2/s]
               rnf_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Foxx_rofl,n)

               ! Frozen runoff, positive downwards [kg/m^2/s]
               rfi_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Foxx_rofi,n)

               ! Fresh water due to melting/freezing, positive downwards
               ! [kg/m^2/s]
               fmltfz_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Fioi_meltw,n)

               ! Salt flux, positive downwards [kg/m^2/s]
               sfl_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Fioi_salt,n)

               ! Shortwave heat flux, positive downwards [W/m^2]
               swa_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Foxx_swnet,n)

               ! Non-solar heat flux, positive downwards [W/m^2]
               nsf_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Foxx_lat,n)
     .                          + x2o_o%rAttr(index_x2o_Foxx_sen,n)
     .                          + x2o_o%rAttr(index_x2o_Foxx_lwup,n)
     .                          + x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
     .                          - ( x2o_o%rAttr(index_x2o_Faxa_snow,n)
     .                            + x2o_o%rAttr(index_x2o_Foxx_rofi,n))
     .                            *SHR_CONST_LATICE

               ! Heat flux due to melting, positive downwards [W/m^2]
               hmlt_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Fioi_melth,n)

               ! Sea level pressure [kg/m/s^2]
               slp_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Sa_pslv,n)

               ! Ice fraction []
               ficem_da(i,j,l2ci) = x2o_o%rAttr(index_x2o_Si_ifrac,n)

               ! 10m wind speed [m/s]
               abswnd_da(i,j,l2ci) =
     .            sqrt(x2o_o%rAttr(index_x2o_So_duu10n,n))

            endif

         enddo
      enddo

      if (nreg == 2) then
         call xctilr(lip_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(sop_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(eva_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(rnf_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(rfi_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(fmltfz_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(sfl_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(swa_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(nsf_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
         call xctilr(hmlt_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      endif

      call fill_global(mval, fval, halo_ps,
     .                 slp_da(1-nbdy,1-nbdy,l2ci))
      call fill_global(mval, fval, halo_ps,
     .                 ficem_da(1-nbdy,1-nbdy,l2ci))
      call fill_global(mval, fval, halo_ps,
     .                 abswnd_da(1-nbdy,1-nbdy,l2ci))

#ifdef PROGCO2
      if (index_x2o_Sa_co2prog > 0) then
         n = 0
         do j = 1, jjcpl
            do i = 1, ii
               n = n + 1
               if     (ip(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = mval
               elseif (cplmsk(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = fval
               else
                  ! Atmospheric co2 concentration [ppmv?]
                  atmco2_da(i,j,l2ci) =
     .               x2o_o%rAttr(index_x2o_Sa_co2prog,n)
               endif
            enddo
         enddo
         call fill_global(mval, fval, halo_ps,
     .                    atmco2_da(1-nbdy,1-nbdy,l2ci))
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: prog. atmospheric co2 read'
      else
         do j = 1, jj
            do i = 1, ii
               if (ip(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = mval
               else
                  atmco2_da(i,j,l2ci) = -1
               endif
            enddo
         enddo
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: prog. atmospheric co2 not read'
      endif

#elif defined(DIAGCO2)
      if (index_x2o_Sa_co2diag > 0) then
         n = 0
         do j = 1, jjcpl
            do i = 1, ii
               n = n + 1
               if     (ip(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = mval
               elseif (cplmsk(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = fval
               else
                  ! Atmospheric co2 concentration [ppmv?]
                  atmco2_da(i,j,l2ci) =
     .               x2o_o%rAttr(index_x2o_Sa_co2diag,n)
               endif
            enddo
         enddo
         call fill_global(mval, fval, halo_ps,
     .                    atmco2_da(1-nbdy,1-nbdy,l2ci))
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: diag. atmospheric co2 read'
      else
         do j = 1, jj
            do i = 1, ii
               if (ip(i,j) == 0) then
                  atmco2_da(i,j,l2ci) = mval
               else
                  atmco2_da(i,j,l2ci) = -1
               endif
            enddo
         enddo
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: diag. atmospheric co2 not read'
      endif
#else
      do j = 1, jj
         do i = 1, ii
            if (ip(i,j) == 0) then
               atmco2_da(i,j,l2ci) = mval
            else
               atmco2_da(i,j,l2ci) = -1
            endif
         enddo
      enddo
      if (mnproc.eq.1)
     .   write (lp,*) 'import_mct: atmospheric co2 not read'
#endif      

      if (csdiag) then
         if (mnproc.eq.1) then
            write (lp,*) 'import_mct:'
         endif
         call chksummsk(ustarw_da(1-nbdy,1-nbdy,l2ci),ip,1,'ustarw')
         call chksummsk(ztx_da(1-nbdy,1-nbdy,l2ci),iu,1,'ztx')
         call chksummsk(mty_da(1-nbdy,1-nbdy,l2ci),iv,1,'mty')
         call chksummsk(lip_da(1-nbdy,1-nbdy,l2ci),ip,1,'lip')
         call chksummsk(sop_da(1-nbdy,1-nbdy,l2ci),ip,1,'sop')
         call chksummsk(eva_da(1-nbdy,1-nbdy,l2ci),ip,1,'eva')
         call chksummsk(rnf_da(1-nbdy,1-nbdy,l2ci),ip,1,'rnf')
         call chksummsk(rfi_da(1-nbdy,1-nbdy,l2ci),ip,1,'rfi')
         call chksummsk(fmltfz_da(1-nbdy,1-nbdy,l2ci),ip,1,'fmltfz')
         call chksummsk(sfl_da(1-nbdy,1-nbdy,l2ci),ip,1,'sfl')
         call chksummsk(swa_da(1-nbdy,1-nbdy,l2ci),ip,1,'swa')
         call chksummsk(nsf_da(1-nbdy,1-nbdy,l2ci),ip,1,'nsf')
         call chksummsk(hmlt_da(1-nbdy,1-nbdy,l2ci),ip,1,'hmlt')
         call chksummsk(slp_da(1-nbdy,1-nbdy,l2ci),ip,1,'slp')
         call chksummsk(ficem_da(1-nbdy,1-nbdy,l2ci),ip,1,'ficem')
         call chksummsk(abswnd_da(1-nbdy,1-nbdy,l2ci),ip,1,'abswnd')
         call chksummsk(atmco2_da(1-nbdy,1-nbdy,l2ci),ip,1,'abswnd')
      endif

      end subroutine import_mct
