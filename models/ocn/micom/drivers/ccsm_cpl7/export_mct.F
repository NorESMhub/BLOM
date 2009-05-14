      subroutine export_mct(o2x_o, perm, nsend, sbuff,
     .                      tlast_coupled)

      ! Uses modules
      
      use mct_mod
!     use seq_flds_mod
      use seq_flds_indices
      use shr_const_mod, only: SHR_CONST_TKFRZ
      use types, only : r8
      use mod_xc

      implicit none

      ! Input/output arguments

      type (mct_aVect)         , intent(inout) :: o2x_o
      integer, dimension(ii*jj), intent(in)    :: perm
      integer                  , intent(in)    :: nsend
      real (r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nsend),
     .   intent(in) :: sbuff
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
      do j = 1, jj
         do i = 1, ii
            n = n + 1
            utmp = .5_r8*( sbuff(i  ,j,index_o2x_So_u)
     .                   + sbuff(i+1,j,index_o2x_So_u))
            vtmp = .5_r8*( sbuff(i,j  ,index_o2x_So_v)
     .                   + sbuff(i,j+1,index_o2x_So_v))
            o2x_o%rattr(index_o2x_So_u,n) =
     .         ( utmp*cos(angle(i,j))
     .         - vtmp*sin(angle(i,j)))*tfac*1.e-2_r8
            o2x_o%rattr(index_o2x_So_v,n) =
     .         ( utmp*sin(angle(i,j))
     .         + vtmp*cos(angle(i,j)))*tfac*1.e-2_r8
            utmp = ( sbuff(i  ,j,index_o2x_So_dhdx)*iu(i  ,j)
     .             + sbuff(i+1,j,index_o2x_So_dhdx)*iu(i+1,j))
     .             /max(1,iu(i,j) + iu(i+1,j))
            vtmp = ( sbuff(i,j  ,index_o2x_So_dhdy)*iv(i,j  )
     .             + sbuff(i,j+1,index_o2x_So_dhdy)*iv(i,j+1))
     .             /max(1,iv(i,j) + iv(i,j+1))
            o2x_o%rAttr(index_o2x_So_dhdx,n) =
     .         ( utmp*cos(angle(i,j))
     .         - vtmp*sin(angle(i,j)))*tfac
            o2x_o%rAttr(index_o2x_So_dhdy,n) =
     .         ( utmp*sin(angle(i,j))
     .         + vtmp*cos(angle(i,j)))*tfac
         enddo
      enddo

      ! ----------------------------------------------------------------
      ! Pack temperature (K), salinity (psu), freezing/melting potential
      ! (W/m~2)
      ! ----------------------------------------------------------------

      n = 0
      do j = 1, jj
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
      ! Pack co2 flux (kg CO2/m^2/s), if requested
      ! ----------------------------------------------------------------

      if (index_o2x_Faoo_fco2 > 0) then
        n = 0
         do j = 1, jj
            do i = 1, ii
               n = n + 1
               o2x_o%rAttr(index_o2x_Faoo_fco2,n) = 
     .            sbuff(i,j,index_o2x_Faoo_fco2)*tfac
            enddo
         enddo
      endif

      tlast_coupled = 0._r8

      !-----------------------------------------------------------------
      ! permute in-place before returning the output state
      !-----------------------------------------------------------------

      call mct_aVect_permute(o2x_o,perm)

      end subroutine export_mct
