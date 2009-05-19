      subroutine import_mct(x2o_o, perm)

      ! Uses modules
      
      use mct_mod
!     use seq_flds_mod
      use seq_flds_indices
      use shr_const_mod, only:
     .   SHR_CONST_RHOSW, SHR_CONST_LATVAP, SHR_CONST_LATICE
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
               sfl(i,j) = mval
               swa(i,j) = mval
               nsf(i,j) = mval
               slp(i,j) = mval
               ficem(i,j) = mval
            elseif (cplmsk(i,j) == 0) then
               lip(i,j) = 0._r8
               sop(i,j) = 0._r8
               eva(i,j) = 0._r8
               rnf(i,j) = 0._r8
               sfl(i,j) = 0._r8
               swa(i,j) = 0._r8
               nsf(i,j) = 0._r8
               slp(i,j) = fval
               ficem(i,j) = fval
            else

               ! Liquid water flux, positive downwards (kg/m^2/s)
               lip(i,j) = x2o_o%rAttr(index_x2o_Foxx_rain,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_meltw,n)

               ! Solid precipitation, positive downwards (kg/m^2/s)
               sop(i,j) = x2o_o%rAttr(index_x2o_Foxx_snow,n)

               ! Evaporation, positive downwards (kg/m^2/s)
               eva(i,j) = x2o_o%rAttr(index_x2o_Foxx_evap,n)

               ! Runoff, positive downwards (kg/m^2/s)
               rnf(i,j) = x2o_o%rAttr(index_x2o_Forr_roff,n)

               ! Salt flux, positive downwards (kg/m^2/s)
               sfl(i,j) = x2o_o%rAttr(index_x2o_Foxx_salt,n)

               ! Shortwave heat flux, positive downwards (W/m^2)
               swa(i,j) = x2o_o%rAttr(index_x2o_Foxx_swnet,n)

               ! Total heat flux, positive downwards (W/m^2)
               nsf(i,j) = x2o_o%rAttr(index_x2o_Foxx_evap,n)
     .                    *SHR_CONST_LATVAP
     .                  + x2o_o%rAttr(index_x2o_Foxx_sen,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_lwup,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_lwdn,n)
     .                  + x2o_o%rAttr(index_x2o_Foxx_melth,n)
     .                  - x2o_o%rAttr(index_x2o_Foxx_snow,n)
     .                    *SHR_CONST_LATICE

               ! Sea level pressure (kg/m/s^2)
               slp(i,j) = x2o_o%rAttr(index_x2o_Sa_pslv,n)

               ! Fields that are not used
               ficem(i,j) = x2o_o%rAttr(index_x2o_Si_ifrac,n)
               util1(i,j) = x2o_o%rAttr(index_x2o_Sx_duu10n,n)

               ! Diagnostics
               dfl(i,j) = 0.

            endif

         enddo
      enddo

      call fill_global(mval, fval, halo_ps, slp)
      call fill_global(mval, fval, halo_ps, ficem)

      if (index_x2o_Sa_co2prog > 0) then
         n = 0
         do j = 1, jj
            do i = 1, ii
               n = n + 1
               ! Only for testing!!!
               dfl(i,j) = x2o_o%rAttr(index_x2o_Sa_co2prog,n)
            enddo
         enddo
      endif
      
      do j=1,jj
         do i=1,ii
            surrlx(i,j)=0.
            salrlx(i,j)=0.
         enddo
      enddo

      end subroutine import_mct
