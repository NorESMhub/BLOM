! ------------------------------------------------------------------------------
! Copyright (C) 2013-2015 Mehmet Ilicak, Mats Bentsen
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
c --- common blocks related to second order turbulence closure
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: 
     .  Prod,          ! shear production
     .  Buoy,          ! buoyancy production
     .  Shear2,        ! square of the shear frequency
     .  L_scale        ! dissipative length scale
c
      common /tke1/ Prod,Buoy,Shear2,L_scale
c
c --- various coefficients
      real gls_s0,gls_s1,gls_s2,gls_s3,gls_s4,gls_s5,gls_s6,
     .     gls_b0,gls_b1,gls_b2,gls_b3,gls_b4,gls_b5,
     .     sqrt2,cmu_fac1,cmu_fac2,cmu_fac3,tke_exp1,gls_exp1,gls_fac6
c
      common /tke2/ gls_s0,gls_s1,gls_s2,gls_s3,gls_s4,gls_s5,gls_s6, 
     .              gls_b0,gls_b1,gls_b2,gls_b3,gls_b4,gls_b5,
     .              sqrt2,cmu_fac1,cmu_fac2,cmu_fac3,tke_exp1,
     .              gls_exp1,gls_fac6
