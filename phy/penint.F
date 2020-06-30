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

      subroutine penint(ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,
     .                  a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
c
      implicit none
c
      real ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,a,ax,ay,axx,ayy,axy
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
      real xx123,yy123,xx135,yy135,xx345,yy345,xy1,xy2,xy3,xy4,xy5,
     .     xy123,xy135,xy345,a123,a135,a345,ax123,ax135,ax345,
     .     ay123,ay135,ay345,axx123,axx135,axx345,ayy123,ayy135,ayy345,
     .     axy123,axy135,axy345
#if defined(TRC) && defined(ATRC)
      real axxx123,axxx135,axxx345,ayyy123,ayyy135,ayyy345,
     .     axxy123,axxy135,axxy345,axyy123,axyy135,axyy345
#endif
c
      xx123=x1*x2+x2*x3+x1*x3
      yy123=y1*y2+y2*y3+y1*y3
      xx135=x1*x3+x3*x5+x1*x5
      yy135=y1*y3+y3*y5+y1*y5
      xx345=x3*x4+x4*x5+x3*x5
      yy345=y3*y4+y4*y5+y3*y5
      xy1=x1*y1
      xy2=x2*y2
      xy3=x3*y3
      xy4=x4*y4
      xy5=x5*y5
      xy123=xy1+xy2+xy3
      xy135=xy1+xy3+xy5
      xy345=xy3+xy4+xy5
c
      a123=.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac
      a135=.5*((x3-x1)*(y5-y1)-(y3-y1)*(x5-x1))*ac
      a345=.5*((x4-x3)*(y5-y3)-(y4-y3)*(x5-x3))*ac
c
      ax123=r1_3*(x1+x2+x3)
      ay123=r1_3*(y1+y2+y3)
      ax135=r1_3*(x1+x3+x5)
      ay135=r1_3*(y1+y3+y5)
      ax345=r1_3*(x3+x4+x5)
      ay345=r1_3*(y3+y4+y5)
      axx123=r1_6*(9.*ax123*ax123-xx123)
      ayy123=r1_6*(9.*ay123*ay123-yy123)
      axy123=r1_12*(9.*ax123*ay123+xy123)
      axx135=r1_6*(9.*ax135*ax135-xx135)
      ayy135=r1_6*(9.*ay135*ay135-yy135)
      axy135=r1_12*(9.*ax135*ay135+xy135)
      axx345=r1_6*(9.*ax345*ax345-xx345)
      ayy345=r1_6*(9.*ay345*ay345-yy345)
      axy345=r1_12*(9.*ax345*ay345+xy345)
#if defined(TRC) && defined(ATRC)
      axxx123=r1_10*((18.*axx123-3.*xx123)*ax123+x1*x2*x3)
      ayyy123=r1_10*((18.*ayy123-3.*yy123)*ay123+y1*y2*y3)
      axxy123=r1_30*(18.*axx123*ay123+3.*ax123*xy123
     .              +x1*xy1+x2*xy2+x3*xy3)
      axyy123=r1_30*(18.*ayy123*ax123+3.*ay123*xy123
     .              +y1*xy1+y2*xy2+y3*xy3)
      axxx135=r1_10*((18.*axx135-3.*xx135)*ax135+x1*x3*x5)
      ayyy135=r1_10*((18.*ayy135-3.*yy135)*ay135+y1*y3*y5)
      axxy135=r1_30*(18.*axx135*ay135+3.*ax135*xy135
     .              +x1*xy1+x3*xy3+x5*xy5)
      axyy135=r1_30*(18.*ayy135*ax135+3.*ay135*xy135
     .              +y1*xy1+y3*xy3+y5*xy5)
      axxx345=r1_10*((18.*axx345-3.*xx345)*ax345+x3*x4*x5)
      ayyy345=r1_10*((18.*ayy345-3.*yy345)*ay345+y3*y4*y5)
      axxy345=r1_30*(18.*axx345*ay345+3.*ax345*xy345
     .              +x3*xy3+x4*xy4+x5*xy5)
      axyy345=r1_30*(18.*ayy345*ax345+3.*ay345*xy345
     .              +y3*xy3+y4*xy4+y5*xy5)
#endif
c
      a=a123+a135+a345
      ax=ax123*a123+ax135*a135+ax345*a345
      ay=ay123*a123+ay135*a135+ay345*a345
      axx=axx123*a123+axx135*a135+axx345*a345
      ayy=ayy123*a123+ayy135*a135+ayy345*a345
      axy=axy123*a123+axy135*a135+axy345*a345
#if defined(TRC) && defined(ATRC)
      axxx=axxx123*a123+axxx135*a135+axxx345*a345
      ayyy=ayyy123*a123+ayyy135*a135+ayyy345*a345
      axxy=axxy123*a123+axxy135*a135+axxy345*a345
      axyy=axyy123*a123+axyy135*a135+axyy345*a345
#endif
c
      return
      end
