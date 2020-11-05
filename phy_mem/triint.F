! ------------------------------------------------------------------------------
! Copyright (C) 2006 Mats Bentsen
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

      subroutine triint(ac,x1,y1,x2,y2,x3,y3,a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
c
      implicit none
c
      real ac,x1,y1,x2,y2,x3,y3,a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
      real axxx,ayyy,axxy,axyy
#endif
c
      real r1_3,r1_6,r1_12
      parameter (r1_3=1./3.,r1_6=1./6.,r1_12=1./12.)
#if defined(TRC) && defined(ATRC)
      real r1_10,r1_30
      parameter (r1_10=.1,r1_30=1./30.)
#endif
c
      real xx,yy,xy1,xy2,xy3,xy
c
      xx=x1*x2+x2*x3+x1*x3
      yy=y1*y2+y2*y3+y1*y3
      xy1=x1*y1
      xy2=x2*y2
      xy3=x3*y3
      xy=xy1+xy2+xy3
c
      a=.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac
c
      ax=r1_3*(x1+x2+x3)
      ay=r1_3*(y1+y2+y3)
      axx=r1_6*(9.*ax*ax-xx)
      ayy=r1_6*(9.*ay*ay-yy)
      axy=r1_12*(9.*ax*ay+xy)
#if defined(TRC) && defined(ATRC)
      axxx=r1_10*((18.*axx-3.*xx)*ax+x1*x2*x3)
      ayyy=r1_10*((18.*ayy-3.*yy)*ay+y1*y2*y3)
      axxy=r1_30*(18.*axx*ay+3.*ax*xy+x1*xy1+x2*xy2+x3*xy3)
      axyy=r1_30*(18.*ayy*ax+3.*ay*xy+y1*xy1+y2*xy2+y3*xy3)
#endif
c
      ax=ax*a
      ay=ay*a
      axx=axx*a
      ayy=ayy*a
      axy=axy*a
#if defined(TRC) && defined(ATRC)
      axxx=axxx*a
      ayyy=ayyy*a
      axxy=axxy*a
      axyy=axyy*a
#endif
c
      return
      end
