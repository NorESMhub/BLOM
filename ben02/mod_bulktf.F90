! ------------------------------------------------------------------------------
! Copyright (C) 2002-2020 Mats Bentsen
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

      module mod_bulktf
c
c --- ------------------------------------------------------------------
c --- This module contains routine and functions for computing transfer
c --- coefficients of momentum, sensible heat, and latent heat.
c --- ------------------------------------------------------------------
c
      implicit none
c
      private
c
      public :: bulktf
c
      contains
c
c --- ------------------------------------------------------------------
c
      real function psiu(zeta)
c
c --- Monin-Obukhov similarity velocity profile function
c
      real zeta
c
      real athird,pi,sqrt3,sqrt3i
      parameter(athird=1./3.,pi=3.141592653589793,
     .          sqrt3=1.732050807568877,sqrt3i=.5773502691896258)
c
      real x,psik,y,psic,f
c
      if (zeta.eq.0.) then
        psiu=0.
      elseif (zeta.gt.0.) then
        psiu=-4.7*zeta
      else
        x=(1.-16.*zeta)**.25
        psik=2.*log((1.+x)*.5)+log((1.+x*x)*.5)-2.*atan(x)+pi*.5
        y=(1.-12.87*zeta)**athird
        psic=1.5*log((y*y+y+1.)*athird)-sqrt3*atan((2.*y+1.)*sqrt3i)
     .      +pi*sqrt3i
        f=1./(1.+zeta*zeta)
        psiu=f*psik+(1.-f)*psic
      endif
c
      return
      end function psiu
c
c --- ------------------------------------------------------------------
c
      real function psitq(zeta)
c
      implicit none
c
c --- Monin-Obukhov similarity temperature and humidity profile function
c
      real zeta
c
      real athird,pi,sqrt3,sqrt3i
      parameter(athird=1./3.,pi=3.141592653589793,
     .          sqrt3=1.732050807568877,sqrt3i=.5773502691896258)
c
      real x,psik,y,psic,f
c
      if (zeta.eq.0.) then
        psitq=0.
      elseif (zeta.gt.0.) then
        psitq=-4.7*zeta
      else
        x=(1.-16.*zeta)**.25
        psik=2.*log((1+x*x)*.5)
        y=(1.-12.87*zeta)**athird
        psic=1.5*log((y*y+y+1.)*athird)-sqrt3*atan((2.*y+1.)*sqrt3i)
     .      +pi*sqrt3i
        f=1./(1.+zeta*zeta)
        psitq=f*psik+(1.-f)*psic
      endif
c
      return
      end function psitq
c
c --- ------------------------------------------------------------------
c
      subroutine lkb(reu,ret,req)
c
      implicit none
c
c --- computes the roughness Reynolds for temperature and humidity
c --- following Liu, Katsaros and Businger (1979), J. Atmos.  Sci., 36,
c --- 1722-1735.
c
      real reu,ret,req
c
      real, dimension(8,2) :: a,b
      real, dimension(8) :: re
      integer i
c
      data a/0.177,1.376,1.026,1.625,4.661,34.904,1667.19,5.88e5,
     .       0.292,1.808,1.393,1.956,4.994,30.709,1448.68,2.98e5/
      data b/0.,0.929,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
     .       0.,0.826,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
      data re/0.11,0.825,3.0,10.0,30.0,100.,300.,1000./
c
      i=1
 10   if (reu.gt.re(i).and.i.lt.8) then
        i=i+1
        goto 10
      endif
c
      ret=a(i,1)*reu**b(i,1)
      req=a(i,2)*reu**b(i,2)
c
      return
      end subroutine lkb
c
c --- ------------------------------------------------------------------
c
      subroutine bulktf(du,zu,ta,zt,qa,zq,ts,qs,icec,cd,ch,ce,wg2)
c
c --- this routine computes momentum, sensible heat, and latent heat
c --- transfer coefficients by one iteration of a bulk flux algorithm
c --- (Fairall et al. 1996)
c
      implicit none
c
c --- input variables:
c ---   du     - velocity difference (wind speed - current) [m/s]
c ---   zu     - measurement height of wind speed [m]
c ---   ta     - air temperature [K]
c ---   zt     - measurement height of air temperature [m]
c ---   qa     - specific humidity of air [kg/kg]
c ---   zq     - measurement height of specific humidity [m]
c ---   ts     - sea surface temperature [K]
c ---   qs     - specific humidity of air at the surface [kg/kg]
c ---   icec   - ice concentration (0 open water, 1 sea ice cover)
c ---   cd     - momentum transfer coefficient
c ---   ch     - sensible heat transfer coefficient
c ---   ce     - latent heat transfer coefficient
c ---   wg2    - gustiness squared [m^2/s^2]
c
      real du,zu,ta,zt,qa,zq,ts,qs,icec,cd,ch,ce,wg2
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c ---   t0     - freezing point of water [K]
c ---   zi     - inversion height [m]
c ---   g      - acceleration due to gravity [m/s^2] 
c ---   beta   - empirical constant relating gustiness to convective
c ---            scaling velocity
c ---   alpha  - Charnock constant
c ---   k      - von Karman constant
c
      real eps,cv,t0,zi,g,beta,athird,alpha,gi,k,ki
      parameter(eps=.62197,cv=1./eps-1.,t0=273.15,zi=600.,g=9.8,
     .          beta=1.2,athird=1./3.,alpha=.011,gi=1./g,k=.4,ki=2.5)
c
      real tv,tac,visca,dt,dq,du1,du2,s,ustar2,ustar,fac,tstar,qstar,
     .     tvstar,li,w3,wg,zetau,zetat,zetaq,z0,cd2,reu,ret,req,
     .     z0t,z0q,ct2,cq2
c
c --- virtual temperature
      tv=ta*(1.+cv*qa)
c
c --- kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
c --- [m^2/s]
      tac=ta-t0
      visca=1.326e-5*(1.+tac*(6.542e-3+tac*(8.301e-6-tac*4.84e-9)))
c
c --- temperature difference
      dt=ta-ts+.0098*zt
c
c --- humidity difference
      dq=qa-qs
c
c --- initial values
      du1=max(du,1.e-2)
      du2=du1*du1
      s=sqrt(du2+wg2)
      ustar2=cd*s*du1
      ustar=sqrt(ustar2)
      fac=ustar/(cd*du1)
      tstar=fac*ch*dt
      qstar=fac*ce*dq
c
c --- inverse Monin-Obukov lenght
      tvstar=tstar*(1+cv*qa)+cv*ta*qstar
      li=min(3./zu,g*k*tvstar/(ustar2*tv))
c
c --- gustiness
      w3=-zi*g*ustar*tvstar/ta
      wg=max(.1,beta*max(0.,w3)**athird)
c
c --- wind speed included gustiness
      s=sqrt(du2+wg*wg)
c
c --- nondimensional heights
      zetau=zu*li
      zetat=zt*li
      zetaq=zq*li
c
c --- roughness lenght - choose roughness lenght for sea ice
c --- corresponding to smooth multiyear ice, Guest and Davidson (1991)
c --- JGR 96, 4709-4721
      z0=icec*2.e-3+(1.-icec)*(0.11*visca/ustar+alpha*ustar2*gi)
c
c --- square root of the momentum transfer coefficient
      cd2=k/max(7.,log(zu/z0)-psiu(zetau))
c
c --- update friction velocity
      ustar=cd2*sqrt(s*du1)
c
c --- roughness Reynolds numbers
      reu=ustar*z0/visca
      call lkb(reu,ret,req)
c
c --- temperature and humidity roughness scales
      fac=visca/ustar
      z0t=fac*ret
      z0q=fac*req
c
c --- transfer coefficient components for temperature and humidity
      ct2=k/max(7.,log(zt/z0t)-psitq(zetat))
      cq2=k/max(7.,log(zq/z0q)-psitq(zetaq))
c
c --- momentum transfer coefficient
      cd=cd2*cd2
c
c --- heat transfer coefficients
      ch=cd2*ct2
      ce=cd2*cq2
c
c --- gustiness squared
      wg2=wg*wg
c
      return
      end subroutine bulktf
c
c --- ------------------------------------------------------------------
c
      end module mod_bulktf
