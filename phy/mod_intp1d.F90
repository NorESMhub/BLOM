! ------------------------------------------------------------------------------
! Copyright (C) 2002-2015 Mats Bentsen
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

      real function intp1d(d1,d2,d3,d4,d5,x)
c
c --- ------------------------------------------------------------------
c --- Interpolate data in one dimension.
c --- ------------------------------------------------------------------
c
      implicit none
c
      real d1,d2,d3,d4,d5,x
c
      real a1,a2,a3,b1,b2,b3,b4,b5,c1,c2
      parameter(a1=-3./7.,a2=-15./7.,a3= 3./2.,
     .          b1= 4./7.,b2=-16./7.,b3=15./7.,b4=-5./7.,b5=2./7.,
     .          c1=-1./7.,c2=  9./14)
      real a,b,c
c
      a=a1*(d1+d5)+a2*d3+a3*(d2+d4)
      b=b1*d1+b2*d2+b3*d3+b4*d4+b5*d5
      c=c1*(d1+d4)+c2*(d2+d3)
c
      intp1d=(a*x+b)*x+c
c
      return
      end function intp1d
