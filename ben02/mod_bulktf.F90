! ------------------------------------------------------------------------------
! Copyright (C) 2002-2020 Mats Bentsen

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_bulktf

  ! --- ------------------------------------------------------------------
  ! --- This module contains routine and functions for computing transfer
  ! --- coefficients of momentum, sensible heat, and latent heat.
  ! --- ------------------------------------------------------------------

  implicit none
  private

  public :: bulktf

  private :: psiu
  private :: psitq
  private :: lkb

contains

  ! --- ------------------------------------------------------------------

  real function psiu(zeta)

    ! --- Monin-Obukhov similarity velocity profile function

    ! Arguments
    real, intent(in) :: zeta

    real :: athird,pi,sqrt3,sqrt3i
    parameter(athird=1./3.,pi = 3.141592653589793, &
         sqrt3=1.732050807568877,sqrt3i = .5773502691896258)

    real :: x,psik,y,psic,f

    if (zeta == 0.) then
      psiu = 0.
    else if (zeta > 0.) then
      psiu = -4.7*zeta
    else
      x = (1.-16.*zeta)**.25
      psik = 2.*log((1.+x)*.5)+log((1.+x*x)*.5)-2.*atan(x)+pi*.5
      y = (1.-12.87*zeta)**athird
      psic = 1.5*log((y*y+y+1.)*athird)-sqrt3*atan((2.*y+1.)*sqrt3i) &
           +pi*sqrt3i
      f = 1./(1.+zeta*zeta)
      psiu = f*psik+(1.-f)*psic
    end if

  end function psiu

  ! --- ------------------------------------------------------------------

  real function psitq(zeta)

    ! --- Monin-Obukhov similarity temperature and humidity profile function

    ! Arguments
    real, intent(in) :: zeta

    real :: athird,pi,sqrt3,sqrt3i
    parameter(athird=1./3.,pi = 3.141592653589793, &
         sqrt3=1.732050807568877,sqrt3i = .5773502691896258)
    real :: x,psik,y,psic,f

    if (zeta == 0.) then
      psitq = 0.
    else if (zeta > 0.) then
      psitq = -4.7*zeta
    else
      x = (1.-16.*zeta)**.25
      psik = 2.*log((1+x*x)*.5)
      y = (1.-12.87*zeta)**athird
      psic = 1.5*log((y*y+y+1.)*athird)-sqrt3*atan((2.*y+1.)*sqrt3i) &
           +pi*sqrt3i
      f = 1./(1.+zeta*zeta)
      psitq = f*psik+(1.-f)*psic
    end if

  end function psitq

  ! --- ------------------------------------------------------------------

  subroutine lkb(reu,ret,req)

    ! --- computes the roughness Reynolds for temperature and humidity
    ! --- following Liu, Katsaros and Businger (1979), J. Atmos.  Sci., 36,
    ! --- 1722-1735.

    ! Arguments
    real, intent(in)  :: reu
    real, intent(out) :: ret
    real, intent(out) :: req

    real, dimension(8,2) :: a,b
    real, dimension(8) :: re
    integer :: i

    data a/0.177,1.376,1.026,1.625,4.661,34.904,1667.19,5.88e5, &
         0.292,1.808,1.393,1.956,4.994,30.709,1448.68,2.98e5/
    data b/0.,0.929,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935, &
         0.,0.826,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
    data re/0.11,0.825,3.0,10.0,30.0,100.,300.,1000./

    i = 1
10  if (reu > re(i).and.i < 8) then
      i = i+1
      goto 10
    end if

    ret = a(i,1)*reu**b(i,1)
    req = a(i,2)*reu**b(i,2)

  end subroutine lkb

  ! --- ------------------------------------------------------------------

  subroutine bulktf(du,zu,ta,zt,qa,zq,ts,qs,icec,cd,ch,ce,wg2)

    ! --- this routine computes momentum, sensible heat, and latent heat
    ! --- transfer coefficients by one iteration of a bulk flux algorithm
    ! --- (Fairall et al. 1996)

    ! Arguments
    real, intent(in)  :: du   ! velocity difference (wind speed - current) [m/s]
    real, intent(in)  :: zu   ! measurement height of wind speed [m]
    real, intent(in)  :: ta   ! air temperature [K]
    real, intent(in)  :: zt   ! measurement height of air temperature [m]
    real, intent(in)  :: qa   ! specific humidity of air [kg/kg]
    real, intent(in)  :: zq   ! measurement height of specific humidity [m]
    real, intent(in)  :: ts   ! sea surface temperature [K]
    real, intent(in)  :: qs   ! specific humidity of air at the surface [kg/kg]
    real, intent(in)  :: icec ! ice concentration (0 open water, 1 sea ice cover)
    real, intent(out) :: cd   ! momentum transfer coefficient
    real, intent(out) :: ch   ! sensible heat transfer coefficient
    real, intent(out) :: ce   ! latent heat transfer coefficient
    real, intent(out) :: wg2  ! gustiness squared [m^2/s^2]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour
    ! ---   t0     - freezing point of water [K]
    ! ---   zi     - inversion height [m]
    ! ---   g      - acceleration due to gravity [m/s^2]
    ! ---   beta   - empirical constant relating gustiness to convective
    ! ---            scaling velocity
    ! ---   alpha  - Charnock constant
    ! ---   k      - von Karman constant

    real :: eps,cv,t0,zi,g,beta,athird,alpha,gi,k,ki
    parameter(eps=.62197,cv=1./eps-1.,t0=273.15,zi=600.,g = 9.8, &
         beta=1.2,athird=1./3.,alpha=.011,gi=1./g,k=.4,ki = 2.5)

    ! Local variables
    real :: tv,tac,visca,dt,dq,du1,du2,s,ustar2,ustar,fac,tstar,qstar, &
         tvstar,li,w3,wg,zetau,zetat,zetaq,z0,cd2,reu,ret,req, &
         z0t,z0q,ct2,cq2

    ! --- virtual temperature
    tv = ta*(1.+cv*qa)

    ! --- kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
    ! --- [m^2/s]
    tac = ta-t0
    visca = 1.326e-5*(1.+tac*(6.542e-3+tac*(8.301e-6-tac*4.84e-9)))

    ! --- temperature difference
    dt = ta-ts+.0098*zt

    ! --- humidity difference
    dq = qa-qs

    ! --- initial values
    du1 = max(du,1.e-2)
    du2 = du1*du1
    s = sqrt(du2+wg2)
    ustar2 = cd*s*du1
    ustar = sqrt(ustar2)
    fac = ustar/(cd*du1)
    tstar = fac*ch*dt
    qstar = fac*ce*dq

    ! --- inverse Monin-Obukov lenght
    tvstar = tstar*(1+cv*qa)+cv*ta*qstar
    li = min(3./zu,g*k*tvstar/(ustar2*tv))

    ! --- gustiness
    w3 = -zi*g*ustar*tvstar/ta
    wg = max(.1,beta*max(0.,w3)**athird)

    ! --- wind speed included gustiness
    s = sqrt(du2+wg*wg)

    ! --- nondimensional heights
    zetau = zu*li
    zetat = zt*li
    zetaq = zq*li

    ! --- roughness lenght - choose roughness lenght for sea ice
    ! --- corresponding to smooth multiyear ice, Guest and Davidson (1991)
    ! --- JGR 96, 4709-4721
    z0 = icec*2.e-3+(1.-icec)*(0.11*visca/ustar+alpha*ustar2*gi)

    ! --- square root of the momentum transfer coefficient
    cd2 = k/max(7.,log(zu/z0)-psiu(zetau))

    ! --- update friction velocity
    ustar = cd2*sqrt(s*du1)

    ! --- roughness Reynolds numbers
    reu = ustar*z0/visca
    call lkb(reu,ret,req)

    ! --- temperature and humidity roughness scales
    fac = visca/ustar
    z0t = fac*ret
    z0q = fac*req

    ! --- transfer coefficient components for temperature and humidity
    ct2 = k/max(7.,log(zt/z0t)-psitq(zetat))
    cq2 = k/max(7.,log(zq/z0q)-psitq(zetaq))

    ! --- momentum transfer coefficient
    cd = cd2*cd2

    ! --- heat transfer coefficients
    ch = cd2*ct2
    ce = cd2*cq2

    ! --- gustiness squared
    wg2 = wg*wg

  end subroutine bulktf

  ! --- ------------------------------------------------------------------

end module mod_bulktf
