! ------------------------------------------------------------------------------
! Copyright (C) 2009 Mats Bentsen
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

      logical function isnan(a) 
c
c --- check if a real number is NaN
c
      implicit none
c
      real a 
c
      if (a.ne.a) then 
        isnan=.true. 
      else 
        isnan=.false. 
      endif 
c
      return 
      end 
