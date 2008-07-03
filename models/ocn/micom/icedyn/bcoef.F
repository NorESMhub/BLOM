C***********************************************************************
      subroutine bcoef (lrhs)
C***********************************************************************
C=======================================================================
C     THIS SUBROUTINE CALCULATES THE CONTRIBUTIONS FROM THE LOCAL GRID
C     POINT, IE (I,J), TO THE FINITE DIFFERENCE APPROXIMATION TO THE
C     MOMENTUM EQUATIONS FOR THE VISCOUS AND ADVECTION TERMS.  THOSE
C     CONTRIBUTIONS TO FROM THE LOCAL RATE OF CHANGE AND FROM THE
C     OCEANIC DRAG TERM WILL BE ADDED AFTERWARDS IN SUBROUTINE  RELCON.
C
C     NOTE: THE FOLLOWING CODE IS BASICALLY A VERSION
C     OF EQUATION 5 IN HIBLER (1979) IN WHICH WE HAVE EXPLICITLY DONE
C     THE DIFFERENTIATION AND THEN FOUND THOSE TERMS THAT
C     ARE PROPORTIONAL TO U(I,J) OR V(I,J)
C=======================================================================
      implicit none
c
      integer lrhs
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      integer i,j,jm1
c-----------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(jm1)
      do 10 j=0,m
        jm1=mod(j-1+mp,mp)
        do 10 i=2,l
c-----------------------------------------------------------------------
c     first do the bulk viscosity terms
c     add in -d/dx( zeta *(du/dx) )
c-----------------------------------------------------------------------
          bu(i,j) = sx2*
     1            ((zeta(i,j)
     3            +zeta(i,jm1))/(pn(i,j)+pn(i,jm1))
     4            *(pm(i,j)+pm(i,jm1))
     5            +(zeta(i-1,j)
     7            +zeta(i-1,jm1))/(pn(i-1,j)+pn(i-1,jm1))
     8            *(pm(i-1,j)+pm(i-1,jm1)))
c
c-----------------------------------------------------------------------
c     add -d/dy ( zeta * (dv/dy) )  to v equation
c-----------------------------------------------------------------------
          bv(i,j) = sy2*
     1            ((zeta(i,j)
     3            +zeta(i-1,j))/(pm(i,j)+pm(i-1,j))
     4            *(pn(i,j)+pn(i-1,j))
     5            +(zeta(i,jm1)
     7            +zeta(i-1,jm1))/(pm(i,jm1)+pm(i-1,jm1))
     8            *(pn(i,jm1)+pn(i-1,jm1)))
c
c-----------------------------------------------------------------------
c     next add in the shear viscosity terms
c     add -d/dx ( eta * (du/dx) ) to u equation
c-----------------------------------------------------------------------
          bu(i,j) = bu(i,j) + sx2*
     1            ((eta(i,j)
     3            +eta(i,jm1))/(pn(i,j)+pn(i,jm1))
     4            *(pm(i,j)+pm(i,jm1))
     5            +(eta(i-1,j)
     7            +eta(i-1,jm1))/(pn(i-1,j)+pn(i-1,jm1))
     8            *(pm(i-1,j)+pm(i-1,jm1)))
c
c-----------------------------------------------------------------------
c     add -d/dy ( eta * (dv/dy) ) to v equation
c-----------------------------------------------------------------------
          bv(i,j) = bv(i,j) + sy2*
     1            ((eta(i,j)
     3            +eta(i-1,j))/(pm(i,j)+pm(i-1,j))
     4            *(pn(i,j)+pn(i-1,j))
     5            +(eta(i,jm1)
     7            +eta(i-1,jm1))/(pm(i,jm1)+pm(i-1,jm1))
     8            *(pn(i,jm1)+pn(i-1,jm1)))
c
c-----------------------------------------------------------------------
c     add -d/dy ( eta * (du/dy) ) to u equation
c-----------------------------------------------------------------------
          bu(i,j) = bu(i,j) + sy2*
     1            ((eta(i,j)
     3            +eta(i-1,j))/(pm(i,j)+pm(i-1,j))
     4            *(pn(i,j)+pn(i-1,j))
     2            +(eta(i,jm1)
     7            +eta(i-1,jm1))/(pm(i,jm1)+pm(i-1,jm1))
     8            *(pn(i,jm1)+pn(i-1,jm1)))
c
c-----------------------------------------------------------------------
c     add -d/dx ( eta * (dv/dx) ) to v equation
c-----------------------------------------------------------------------
          bv(i,j) = bv(i,j) + sx2*
     1            ((eta(i,j)
     3            +eta(i,jm1))/(pn(i,j)+pn(i,jm1))
     4            *(pm(i,j)+pm(i,jm1))
     5            +(eta(i-1,j)
     7            +eta(i-1,jm1))/(pn(i-1,j)+pn(i-1,jm1))
     8            *(pm(i-1,j)+pm(i-1,jm1)))		! corrected by m.s.
c
c-----------------------------------------------------------------------
c     add in horizontal advection terms:
c     note - this term enters because of the othonormal coordinates
c
c advection of momentum is neglected because it is the smallest term. 
c M. Harder, 25-jan-95
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     add in m * u * d/dx (u)
c-----------------------------------------------------------------------
cmh          bu(i,j) = bu(i,j)+amas(i,j)*u(i,j,lrhs)*
cmh     1     (1./(pn(i,j)+pn(i,j-1))
cmh     2     -1./(pn(i-1,j)+pn(i-1,j-1)))/dx
c
c-----------------------------------------------------------------------
c     add in m * v * d/dy (u)
c-----------------------------------------------------------------------
cmh          bu(i,j) = bu(i,j)+amas(i,j)*v(i,j,lrhs)*
cmh     3     (1./(pm(i,j)+pm(i-1,j))
cmh     4     -1./(pm(i,j-1)+pm(i-1,j-1)))/dy
c
c-----------------------------------------------------------------------
c     add in m * u d/dx (v)
c-----------------------------------------------------------------------
cmh          bv(i,j) = bv(i,j)+amas(i,j)*u(i,j,lrhs)*
cmh     1     (1./(pn(i,j)+pn(i,j-1))
cmh     2     -1./(pn(i-1,j)+pn(i-1,j-1)))/dx
c
c-----------------------------------------------------------------------
c     add in m * v d/dy (v)
c-----------------------------------------------------------------------
cmh         bv(i,j) = bv(i,j)+amas(i,j)*v(i,j,lrhs)*
cmh     3     (1./(pm(i,j)+pm(i-1,j))
cmh     4     -1./(pm(i,j-1)+pm(i-1,j-1)))/dy

   10 continue
c$OMP END PARALLEL DO
c
      return
      end
