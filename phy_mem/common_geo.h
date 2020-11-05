! ------------------------------------------------------------------------------
! Copyright (C) 2004-2015 Mats Bentsen
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

c --- ------------------------------------------------------------------
c --- common blocks related to the model grids geographical coordinates
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     .  qclat,qclon,   ! grid coordinates q-cell corners
     .  pclat,pclon,   ! grid coordinates p-cell corners
     .  uclat,uclon,   ! grid coordinates u-cell corners
     .  vclat,vclon    ! grid coordinates v-cell corners
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  qlat,qlon,     ! grid coordinates q-points
     .  plat,plon,     ! grid coordinates p-points
     .  ulat,ulon,     ! grid coordinates u-points
     .  vlat,vlon,     ! grid coordinates v-points
     .  angle,         ! local angle of i-direction and eastward
                       ! direction at p-points
     .  cosang,sinang  ! cosine and sine of local angle of i-direction
                       ! and eastward direction at p-points

      integer ::
     .  nwp            ! number of wet grid cells
c
      common /geo/ qclat,qclon,pclat,pclon,uclat,uclon,vclat,vclon,
     .             qlat,qlon,plat,plon,ulat,ulon,vlat,vlon,angle,
     .             cosang,sinang,nwp
