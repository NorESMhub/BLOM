! ------------------------------------------------------------------------------
! Copyright (C) 2004-2017 Mats Bentsen
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
c --- common block related to model calendar and time management
c --- ------------------------------------------------------------------
c
      real time0
      integer, dimension(12) :: nd_in_m
      integer nday1,nday2,nday,nmonth,nyear,nday0,nmonth0,nyear0,
     .        nday_in_year,nday_of_year,nstep_in_day
      character*19 calendar
c
      common /clndr/ time0,nd_in_m,nday1,nday2,nday,nmonth,nyear,nday0,
     .               nmonth0,nyear0,nday_in_year,nday_of_year,
     .               nstep_in_day,calendar
