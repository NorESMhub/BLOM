! ------------------------------------------------------------------------------
! Copyright (C) 2002-2015 Mats Bentsen

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

module mod_ben02func

  ! --- ------------------------------------------------------------------
  ! --- This module contains supporting functions for generating surface
  ! --- fluxes from NCEP or ERA40 reanalysis following Bentsen and Drange
  ! --- (2002).
  ! ---  -----------------------------------------------------------------

  implicit none
  private

  public :: spherdist,qsatw,dqsatw,qsati,dqsati,rhoair

contains

  ! --- ------------------------------------------------------------------

  real function spherdist(r,lon1,lat1,lon2,lat2)

    ! --- ------------------------------------------------------------------
    ! --- Computes the distance between geo. pos. lon1,lat1 and lon2,lat2 on
    ! --- a sphere with radius r
    ! --- ------------------------------------------------------------------

    real, intent(in) :: r,lon1,lat1,lon2,lat2

    real :: rad
    parameter(rad = 1.74532925199432958e-02)
    real :: lambda,phi,x1,y1,z1,x2,y2,z2

    phi = lon1*rad
    lambda = lat1*rad
    x1 = cos(lambda)*cos(phi)
    y1 = cos(lambda)*sin(phi)
    z1 = sin(lambda)

    phi = lon2*rad
    lambda = lat2*rad
    x2 = cos(lambda)*cos(phi)
    y2 = cos(lambda)*sin(phi)
    z2 = sin(lambda)

    spherdist = r*acos(min(1.,max(-1.,x1*x2+y1*y2+z1*z2)))

  end function spherdist

  ! --- ------------------------------------------------------------------

  real function qsatw(t,p)

    ! --- computes saturation specific humidity [kg/kg] over open water,
    ! --- Buck (1981) JAM 20, 1527-1532

    ! Arguments
    real, intent(in) :: t ! air temperature [K]
    real, intent(in) :: p ! atmospheric pressure [Pa]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour
    real :: eps,c1,c2,c3,c4,t0,t1
    parameter(eps=0.62197,c1=611.21,c2=1.0007,c3=3.46e-8,c4 = 17.502, &
         t0=273.15,t1 = 32.19)

    real :: tl,e

    tl = max(150.,t)
    e = c1*(c2+c3*p)*exp(c4*(tl-t0)/(tl-t1))
    qsatw = eps*e/(p-(1.-eps)*e)

  end function qsatw

  ! --- ------------------------------------------------------------------

  real function dqsatw(t,p)

    ! --- computes the derivative of saturation specific humidity [kg/kg/K]
    ! --- over open water with respect to temperature,
    ! --- Buck (1981) JAM 20, 1527-1532

    ! Arguments
    real, intent(in) :: t ! air temperature [K]
    real, intent(in) :: p ! atmospheric pressure [Pa]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour

    real :: eps,c1,c2,c3,c4,t0,t1
    parameter(eps=0.62197,c1=611.21,c2=1.0007,c3=3.46e-8,c4 = 17.502, &
         t0=273.15,t1 = 32.19)

    real :: tl,e,dedt

    tl = max(150.,t)
    e = c1*(c2+c3*p)*exp(c4*(tl-t0)/(tl-t1))
    dedt = e*c4*(t0-t1)/(tl-t1)**2
    dqsatw = dedt*eps*p/(p-(1.-eps)*e)**2

  end function dqsatw

  ! --- ------------------------------------------------------------------

  real function qsati(t,p)

    ! --- computes saturation specific humidity [kg/kg] over sea ice,
    ! --- Parkinson and Washington (1979) JGR 84, 311-337

    ! Arguments
    real, intent(in) :: t ! air temperature [K]
    real, intent(in) :: p ! atmospheric pressure [Pa]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour

    real :: eps,c1,c2,t0,t1
    parameter(eps=0.62197,c1=611.,c2=9.5,t0=273.15,t1 = 7.66)

    real :: tl,e

    tl = max(150.,t)
    e = c1*10.**(c2*(tl-t0)/(tl-t1))
    qsati = eps*e/(p-(1.-eps)*e)

  end function qsati

  ! --- ------------------------------------------------------------------

  real function dqsati(t,p)

    ! --- computes the derivative of saturation specific humidity [kg/kg/k]
    ! --- over sea ice with respect to temperature,
    ! --- Parkinson and Washington (1979) JGR 84, 311-337

    ! Arguments
    real, intent(in) :: t ! air temperature [K]
    real, intent(in) :: p ! atmospheric pressure [Pa]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour

    real :: eps,c1,c2,t0,t1
    parameter(eps=0.62197,c1=611.,c2=9.5,t0=273.15,t1 = 7.66)

    real :: tl,e,dedt

    tl = max(150.,t)
    e = c1*10.**(c2*(tl-t0)/(tl-t1))
    dedt = e*c2*(t0-t1)*log(10.)/(tl-t1)**2
    dqsati = dedt*eps*p/(p-(1.-eps)*e)**2

  end function dqsati

  ! --- ------------------------------------------------------------------

  real function rhoair(t,q,p)

    ! --- computes air density [kg/m^3]

    ! Arguments
    real, intent(in) :: t ! air temperature [K]
    real, intent(in) :: p ! atmospheric pressure [Pa]
    real, intent(in) :: q ! specific humidity of air [kg/kg]

    ! --- parameters:
    ! ---   eps    - molecular weight ratio of dry air and water vapour
    ! ---   rgas   - gas constant for dry air [J/(kg*K)]

    real :: eps,cv,rgas
    parameter(eps=.62197,cv=1./eps-1.,rgas = 287.04)

    real :: tv

    ! --- virtual temperature
    tv = t*(1.+cv*q)

    ! --- moist air density [kg/m^3]
    rhoair = p/(rgas*tv)

  end function rhoair

  ! --- ------------------------------------------------------------------

end module mod_ben02func
