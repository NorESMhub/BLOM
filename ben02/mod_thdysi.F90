! ------------------------------------------------------------------------------
! Copyright (C) 2015 Mats Bentsen
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

      module mod_thdysi
c
c --- ------------------------------------------------------------------
c --- This module contains parameters and variables for the
c --- thermodynamic sea ice.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      private
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  tsrfm,ticem
c
      real ::
     .  albi_f,        ! max albedo over ice []
     .  albi_m,        ! max albedo over melting ice []
     .  albs_f,        ! albedo over snow []
     .  albs_m,        ! albedo over melting snow []
     .  rhoice,        ! density of ice [kg m-3]
     .  rhosnw,        ! density of snow [kg m-3]
     .  rkice,         ! ice conductivity [W m-1 K-1]
     .  rksnw,         ! snow conductivity [W m-1 K-1]
     .  fusi,          ! heat of fusion of ice [J kg-1]
     .  fuss,          ! heat of fusion of snow [J kg-1]
     .  fice_max,      ! maximum fractional ice cover []
     .  tice_m,        ! melting point of ice [K]
     .  tsnw_m,        ! melting point of snow [K]
     .  hice_nhmn,     ! min. ice thickness northern hemi. [m]
     .  hice_shmn,     ! min. ice thickness southern hemi. [m]
     .  sagets,        ! snow aging timescale [s-1]
     .  sice,          ! salinity of seaice [g kg-1]
     .  cwi,           ! ice-ocean heat transfer coeff. []
     .  cuc            ! const. for heat flux associated with
                       ! under-cooled water, resulting in a temp.
                       ! adjustment of a 20 m mixed layer towards
                       ! freezing point with an e-folding timescale of
                       ! approx. one day. [W m-2 K-1]
c
      data albi_f   /.70/,
     .     albi_m   /.60/,
     .     albs_f   /.85/,
     .     albs_m   /.75/,
     .     rhoice   /906./,
     .     rhosnw   /330./,
     .     rkice    /2.04/,
     .     rksnw    /.31/,
     .     fusi     /3.02e8/,
     .     fuss     /1.10e8/,
     .     fice_max /.995/,
     .     tice_m   /273.05/,
     .     tsnw_m   /273.15/,
     .     hice_nhmn/.50/,
     .     hice_shmn/.30/,
     .     sagets   /2.e-7/,
     .     sice     /6./,
     .     cwi      /0.006/,
     .     cuc      /1.e3/
c
      public :: tsrfm,ticem,albi_f,albi_m,albs_f,albs_m,rhoice,rhosnw,
     .          rkice,rksnw,fusi,fuss,fice_max,tice_m,tsnw_m,hice_nhmn,
     .          hice_shmn,sagets,sice,cwi,cuc
c
      contains
c
      end module mod_thdysi
