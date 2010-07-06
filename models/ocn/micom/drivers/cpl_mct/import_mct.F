      subroutine import_mct(x2o_o, perm)

      ! Uses modules
      
      use mct_mod
!     use seq_flds_mod
      use seq_flds_indices
      use shr_const_mod, only: SHR_CONST_RHOSW, SHR_CONST_LATICE
      use types, only : r8
      use mod_xc

      implicit none

      ! Input/output arguments

      type (mct_aVect), intent(inout) :: x2o_o
      integer, dimension(ii*jj), intent(in)    :: perm

      ! Include files

#include "common_blocks.h"
#include "common_geo.h"
#include "common_forc.h"
#include "common_ice.h"

      ! Local parameters

      real (r8), parameter ::
     .  mval = -1.e12_r8,
     .  fval = -1.e13_r8

      ! Local variables

      integer i, j, l, n
      real (r8) :: utmp, vtmp

      !-----------------------------------------------------------------
      ! unpermute in-place before unpacking
      !-----------------------------------------------------------------

      call mct_aVect_unpermute(x2o_o, perm)

      ! Unpack

      n = 0
      do j = 1, jj
         do i = 1, ii
            n = n + 1
            if     (ip(i,j) == 0) then
               util1(i,j) = mval
               util2(i,j) = mval
               ustarw(i,j) = mval
            elseif (cplmsk(i,j) == 0) then
               util1(i,j) = fval
               util2(i,j) = fval
               ustarw(i,j) = fval
            else
               utmp = x2o_o%rAttr(index_x2o_Foxx_taux,n)
               vtmp = x2o_o%rAttr(index_x2o_Foxx_tauy,n)
               util1(i,j) =   utmp*cos(angle(i,j))
     .                      + vtmp*sin(angle(i,j))
               util2(i,j) = - utmp*sin(angle(i,j))
     .                      + vtmp*cos(angle(i,j))

               ! Friction velocity (m/s)
               ustarw(i,j) = sqrt(sqrt(utmp*utmp+vtmp*vtmp)
     .                      /SHR_CONST_RHOSW)
            endif
         enddo
      enddo

      call fill_global(mval, fval, halo_pv, util1)
      call fill_global(mval, fval, halo_pv, util2)
      call fill_global(mval, fval, halo_ps, ustarw)

      call xctilr(util1, 1,1, 1,1, halo_pv)
      call xctilr(util2, 1,1, 1,1, halo_pv)

      do j=1,jj
         do l=1,isu(j)
         do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            ! x-component of momentum flux (kg/m/s^2)
            ztx(i,j) = .5_r8*(util1(i-1,j) + util1(i,j))
         enddo
         enddo
         do l=1,isv(j)
         do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            ! y-component of momentum flux (kg/m/s^2)
            mty(i,j) = .5_r8*(util2(i,j-1) + util2(i,j))
         enddo
         enddo
      enddo
      
      n = 0
      do j = 1, jj
         do i = 1, ii
            n = n + 1
            if     (ip(i,j) == 0) then
               lip(i,j) = mval
               sop(i,j) = mval
               eva(i,j) = mval
               rnf(i,j) = mval
               rfi(i,j) = mval
               fmltfz(i,j) = mval
               sfl(i,j) = mval
               swa(i,j) = mval
               nsf(i,j) = mval
               hmlt(i,j) = mval
               slp(i,j) = mval
               ficem(i,j) = mval
               abswnd(i,j) = mval
            elseif (cplmsk(i,j) == 0) then
               lip(i,j) = 0._r8
               sop(i,j) = 0._r8
               eva(i,j) = 0._r8
               rnf(i,j) = 0._r8
               rfi(i,j) = 0._r8
               fmltfz(i,j) = 0._r8
               sfl(i,j) = 0._r8
               swa(i,j) = 0._r8
               nsf(i,j) = 0._r8
               hmlt(i,j) = 0._r8
               slp(i,j) = fval
               ficem(i,j) = fval
               abswnd(i,j) = fval
            else

               ! Liquid water flux, positive downwards [kg/m^2/s]
               lip(i,j) = x2o_o%rAttr(index_x2o_Foxx_rain,n)

               ! Solid precipitation, positive downwards [kg/m^2/s]
               sop(i,j) = x2o_o%rAttr(index_x2o_Foxx_snow,n)

               ! Evaporation, positive downwards [kg/m^2/s]
               eva(i,j) = x2o_o%rAttr(index_x2o_Foxx_evap,n)

               ! Liquid runoff, positive downwards [kg/m^2/s]
               rnf(i,j) = x2o_o%rAttr(index_x2o_Forr_roff,n)

               ! Frozen runoff, positive downwards [kg/m^2/s]
               rfi(i,j) = x2o_o%rAttr(index_x2o_Forr_ioff,n)

               ! Fresh water due to melting/freezing, positive downwards
               ! [kg/m^2/s]
               fmltfz(i,j) = x2o_o%rAttr(index_x2o_Foxx_meltw,n)

               ! Salt flux, positive downwards [kg/m^2/s]
               sfl(i,j) = x2o_o%rAttr(index_x2o_Foxx_salt,n)

               ! Shortwave heat flux, positive downwards [W/m^2]
               swa(i,j) = x2o_o%rAttr(index_x2o_Foxx_swnet,n)

               ! Non-solar heat flux, positive downwards [W/m^2]
               nsf(i,j) = x2o_o%rAttr(index_x2o_Foxx_lat,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_sen,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_lwup,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_lwdn,n)
     .                  - ( x2o_o%rAttr(index_x2o_Foxx_snow,n)
     .                    + x2o_o%rAttr(index_x2o_Forr_ioff,n))
     .                    *SHR_CONST_LATICE

               ! Heat flux due to melting, positive downwards [W/m^2]
               hmlt(i,j) = x2o_o%rAttr(index_x2o_Foxx_melth,n)

               ! Sea level pressure [kg/m/s^2]
               slp(i,j) = x2o_o%rAttr(index_x2o_Sa_pslv,n)

               ! Ice fraction []
               ficem(i,j) = x2o_o%rAttr(index_x2o_Si_ifrac,n)

               ! 10m wind speed [m/s]
               abswnd(i,j) = sqrt(x2o_o%rAttr(index_x2o_Sx_duu10n,n))

            endif

         enddo
      enddo

      call fill_global(mval, fval, halo_ps, slp)
      call fill_global(mval, fval, halo_ps, ficem)
      call fill_global(mval, fval, halo_ps, abswnd)

      if (index_x2o_Sa_co2prog > 0) then
         n = 0
         do j = 1, jj
            do i = 1, ii
               n = n + 1
               if     (ip(i,j) == 0) then
                  atmco2(i,j) = mval
               elseif (cplmsk(i,j) == 0) then
                  atmco2(i,j) = fval
               else
                  ! Atmospheric co2 concentration [ppmv?]
                  atmco2(i,j) = x2o_o%rAttr(index_x2o_Sa_co2prog,n)
               endif
            enddo
         enddo
         call fill_global(mval, fval, halo_ps, atmco2)
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: atmospheric co2 read'
      else
         do j = 1, jj
            do i = 1, ii
               if (ip(i,j) == 0) then
                  atmco2(i,j) = mval
               else
                  atmco2(i,j) = -1
               endif
            enddo
         enddo
         if (mnproc.eq.1)
     .      write (lp,*) 'import_mct: atmospheric co2 not read'
      endif
      
      if (csdiag) then
         if (mnproc.eq.1) then
            write (lp,*) 'import_mct:'
         endif
         call chksummsk(ustarw,ip,1,'ustarw')
         call chksummsk(ztx,iu,1,'ztx')
         call chksummsk(mty,iv,1,'mty')
         call chksummsk(lip,ip,1,'lip')
         call chksummsk(sop,ip,1,'sop')
         call chksummsk(eva,ip,1,'eva')
         call chksummsk(rnf,ip,1,'rnf')
         call chksummsk(rfi,ip,1,'rfi')
         call chksummsk(fmltfz,ip,1,'fmltfz')
         call chksummsk(sfl,ip,1,'sfl')
         call chksummsk(swa,ip,1,'swa')
         call chksummsk(nsf,ip,1,'nsf')
         call chksummsk(hmlt,ip,1,'hmlt')
         call chksummsk(slp,ip,1,'slp')
         call chksummsk(ficem,ip,1,'ficem')
         call chksummsk(abswnd,ip,1,'abswnd')
      endif

      end subroutine import_mct
