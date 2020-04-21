! ------------------------------------------------------------------------------
! Copyright (C) 2008-2012 Mats Bentsen
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

      subroutine domain_mct(gsMap_ocn, dom_ocn, lsize, perm, jjcpl)

      ! Uses modules
      
      use mct_mod
      use seq_flds_mod
      use shr_const_mod, only: SHR_CONST_REARTH
      use types, only : r8
      use mod_xc

      implicit none

      ! Input/output arguments

      type (mct_gsMap)         , intent(in)    :: gsMap_ocn
      type (mct_ggrid)         , intent(inout) :: dom_ocn
      integer                  , intent(in)    :: lsize
      integer, dimension(lsize), intent(in)    :: perm
      integer                  , intent(in)    :: jjcpl

      ! Include files

#include "common_geo.h"
#include "common_blocks.h"

      ! Local variables

      integer, pointer :: idata(:)
      real (r8), pointer :: rdata(:)
      integer i, j, n
      real (r8) :: radius

      ! ----------------------------------------------------------------
      ! Initialize mct domain type
      ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 
      ! 0 (non-ocean)
      ! ----------------------------------------------------------------

      call mct_gGrid_init(GGrid = dom_ocn,
     .                    CoordChars = trim(seq_flds_dom_coord),
     .                    OtherChars = trim(seq_flds_dom_other),
     .                    lsize = lsize)
      allocate(rdata(lsize))

      ! ----------------------------------------------------------------
      ! Determine global gridpoint number attribute, GlobGridNum, which
      ! is set automatically by MCT
      ! ----------------------------------------------------------------

      call mct_gsMap_orderedPoints(gsMap_ocn, mnproc - 1, idata)
      call mct_gGrid_importIAttr(dom_ocn, 'GlobGridNum', idata, lsize)

      ! ----------------------------------------------------------------
      ! Determine domain (numbering scheme is: West to East and South to
      ! North to South pole)
      ! Initialize attribute vector with special value
      ! ----------------------------------------------------------------

      rdata(:) = -9999.0_r8 
      call mct_gGrid_importRAttr(dom_ocn, "lat"  , rdata, lsize) 
      call mct_gGrid_importRAttr(dom_ocn, "lon"  , rdata, lsize) 
      call mct_gGrid_importRAttr(dom_ocn, "area" , rdata, lsize) 
      call mct_gGrid_importRAttr(dom_ocn, "aream", rdata, lsize) 
      rdata(:) = 0.0_r8     
      call mct_gGrid_importRAttr(dom_ocn, "mask", rdata, lsize) 
      call mct_gGrid_importRAttr(dom_ocn, "frac", rdata, lsize) 

      ! ----------------------------------------------------------------
      ! Fill in correct values for domain components
      ! ----------------------------------------------------------------

      ! A correction for north pole mapping of velocity fields in the
      ! coupler requires longitudes in the range [0, 360) degrees to
      ! work.
      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            rdata(n) = modulo(plon(i,j), 360._r8)
         enddo
      enddo
      call mct_gGrid_importRattr(dom_ocn, "lon", rdata, lsize) 

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            rdata(n) = plat(i,j)
         enddo
      enddo
      call mct_gGrid_importRattr(dom_ocn, "lat", rdata, lsize) 

      radius = SHR_CONST_REARTH*1.e2_r8 ! Earth's radius in cm
      
      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            rdata(n) = scp2(i,j)/(radius*radius)
         enddo
      enddo
      call mct_gGrid_importRattr(dom_ocn, "area", rdata, lsize) 

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            rdata(n) = real(cplmsk(i,j), kind = r8)
         enddo
      enddo
      call mct_gGrid_importRattr(dom_ocn, "mask", rdata, lsize) 
      call mct_gGrid_importRattr(dom_ocn, "frac", rdata, lsize) 

      !-----------------------------------------------------------------
      ! Permute dom_ocn to have ascending order
      !-----------------------------------------------------------------

      call mct_gGrid_permute(dom_ocn, perm)

      deallocate(rdata)

      end subroutine domain_mct
