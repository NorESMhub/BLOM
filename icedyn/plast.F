c***********************************************************************
      subroutine plast (lrhs)
c***********************************************************************
c Calculate viscosities at scalar grid points.
c
c Variables
c      zeta      out      modified      bulk  viscosity
c      eta      out      modified      shear viscosity
c      delta      out      modified      a function of strain rates
c
c Constants
c      ecm2            1/(e**2)
c      gmin            lower bound for delta
c      om            outflow mask
c
c Modified by Markus Harder, January 1995:
c Ice pressure -> 0 for delta -> 0.
      implicit none

      integer lrhs

#include "param_icedyn.h"
#include "common_icedyn.h"

      real pi,Pplast,e11,e22,e12,delta
      integer i,j,jp1
c-----------------------------------------------------------------------
      pi =4.*atan(1.)
c-----------------------------------------------------------------------
c*** Ice pressure and viscosities
c-----------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(Pplast,delta)
      do j=0,M
        jp1=mod(j+1,MP)
        do i=1,LM
c-----------------------------------------------------------------------
c*** Calculate ice strength for plastic flow
c-----------------------------------------------------------------------
          Pplast = Pstar(i,j)*h(i,j)
     &              *exp(-Cstar*(1.0-A(i,j)))*om(i,j)
c-----------------------------------------------------------------------
c*** Calculate strain rates e11, e12, e22
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c e11 = d(u)/dx
c-----------------------------------------------------------------------
          e11 = 0.5*(pm(i,j)/dx)*
     1      ( (u(i+1,j,lrhs)+u(i+1,jp1,lrhs))
     2     -  (u(i  ,j,lrhs)+u(i  ,jp1,lrhs)) )
c-----------------------------------------------------------------------
c e22 = d(v)/dy
c-----------------------------------------------------------------------
          e22 = 0.5*(pn(i,j)/dy)*
     1      ( (v(i,jp1,lrhs)+v(i+1,jp1,lrhs))
     2     -  (v(i,j  ,lrhs)+v(i+1,j  ,lrhs)) )
c-----------------------------------------------------------------------
c e12 = ( d(u)/dy + d(v)/dx ) / 2
c-----------------------------------------------------------------------
          e12 = (pm(i,j)/dx)*
     1      ( (v(i+1,j,lrhs)+v(i+1,jp1,lrhs))
     2     -  (v(i  ,j,lrhs)+v(i  ,jp1,lrhs)) )

          e12 = e12+(pn(i,j)/dy)*
     1      ( (u(i,jp1,lrhs)+u(i+1,jp1,lrhs))
     2     -  (u(i,j  ,lrhs)+u(i+1,j  ,lrhs)) )

          e12 = 0.25*e12
c-----------------------------------------------------------------------
c ELLIPTICAL YIELD CURVE
c-----------------------------------------------------------------------
c*** Calculate delta
c In case of small strain rates, delta becomes very small,
c yielding extremly large viscosities.
c To deal with this case, a lower limit is imposed on delta.
c-----------------------------------------------------------------------
          delta = sqrt((e11**2+e22**2)*(1.0+ecm2)
     1                 +4.0*ecm2*e12**2
     2                 +2.0*e11*e22*(1.0-ecm2)     )
c-----------------------------------------------------------------------
c Modify ice pressure for very small strain rates
c-----------------------------------------------------------------------
          P(i,j) = Pplast*delta/(delta+gmin)
c-----------------------------------------------------------------------
c Calculate viscosities
c-----------------------------------------------------------------------
          zeta(i,j) = 0.5*Pplast/(delta+gmin)
          eta (i,j) = ecm2*zeta(i,j)
c-----------------------------------------------------------------------
c Set viscosities to zero at outflow points
c (In previous codes, this was done by subroutine outbcs)
c-----------------------------------------------------------------------
          zeta(i,j) = zeta(i,j)*om(i,j)
          eta (i,j) = eta (i,j)*om(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
      return
      end            ! subroutine plast
