! ------------------------------------------------------------------------------
! Copyright (C) 2007-2020 Mats Bentsen
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

module mod_eos
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to the equation of
! state.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_config, only: expcnf
   use mod_xc, only: mnproc, lp, xcstop

   implicit none

   private

   ! Coefficients for the functional fit of in situ density.
   real(r8), parameter :: &
      a11 =  9.9985372432159340e-01_r8, &
      a12 =  1.0380621928183473e-02_r8, &
      a13 =  1.7073577195684715e-03_r8, &
      a14 = -3.6570490496333680e-05_r8, &
      a15 = -7.3677944503527477e-06_r8, &
      a16 = -3.5529175999643348e-06_r8, &
      b11 =  1.7083494994335439e-10_r8, &
      b12 =  7.1567921402953455e-13_r8, &
      b13 =  1.2821026080049485e-13_r8, &
      a21 =  1.0_r8                   , &
      a22 =  1.0316374535350838e-02_r8, &
      a23 =  8.9521792365142522e-04_r8, &
      a24 = -2.8438341552142710e-05_r8, &
      a25 = -1.1887778959461776e-05_r8, &
      a26 = -4.0163964812921489e-06_r8, &
      b21 =  1.1995545126831476e-10_r8, &
      b22 =  5.5234008384648383e-13_r8, &
      b23 =  8.4310335919950873e-14_r8

   ! Reference pressure [g cm-1 s-2].
   real(r8) :: pref

   ! Coefficients for potential density in sigma units with reference pressure
   ! at pref.
   real(r8) :: ap11, ap12, ap13, ap14, ap15, ap16, &
               ap21, ap22, ap23, ap24, ap25, ap26

   ! Coefficients for potential density in sigma units with reference pressure
   ! at the surface.
   real(r8) :: ap110, ap120, ap130, ap140, ap150, ap160, &
               ap210, ap220, ap230, ap240, ap250, ap260

   ! Coefficients for freezing temperature
   real(r8) :: atf, btf, ctf

   public :: pref, &
             ap11, ap12, ap13, ap14, ap15, ap16, &
             ap21, ap22, ap23, ap24, ap25, ap26, &
             atf, btf, ctf, &
             inieos, rho, alp, sig, sig0, dsigdt, dsigdt0, dsigds, dsigds0, &
             tofsig, sofsig, p_alpha, p_p_alpha, delphi

contains

   subroutine inieos
   ! ---------------------------------------------------------------------------
   ! Initialize coefficients related to the equation of state.
   ! ---------------------------------------------------------------------------

      ! In situ density [kg m-3] as a function of pressure, potential
      ! temperature and salinity is approximated by the functional form
      !    rho(p,th,s) = P1(p,th,s)/P2(p,th,s)
      ! where
      !    P1(p,th,s) = a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s
      !               + (b11 + b12*th + b13*s)*p
      ! and
      !    P2(p,th,s) = a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s
      !               + (b21 + b22*th + b23*s)*p
      ! Here we compute the coefficients needed for an expression for potential
      ! density [g cm-3] in sigma units of the form
      !    sig(th,s) = R1(th,s)/R2(th,s) 
      ! where
      !    R1(p,th,s) = ap11 + (ap12 + ap14*th + ap15*s)*th + (ap13 + ap16*s)*s
      ! and
      !    R2(p,th,s) = ap21 + (ap22 + ap24*th + ap25*s)*th + (ap23 + ap26*s)*s

      ap21 = a21 + b21*pref
      ap22 = a22 + b22*pref
      ap23 = a23 + b23*pref
      ap24 = a24
      ap25 = a25
      ap26 = a26
      ap11 = a11 + b11*pref - ap21
      ap12 = a12 + b12*pref - ap22
      ap13 = a13 + b13*pref - ap23
      ap14 = a14 - ap24
      ap15 = a15 - ap25
      ap16 = a16 - ap26

      ap210 = a21
      ap220 = a22
      ap230 = a23
      ap240 = a24
      ap250 = a25
      ap260 = a26
      ap110 = a11 - ap210
      ap120 = a12 - ap220
      ap130 = a13 - ap230
      ap140 = a14 - ap240
      ap150 = a15 - ap250
      ap160 = a16 - ap260

      ! Coefficients for freezing temperature.
      select case (trim(expcnf))
         case ('cesm')
            atf =  0._r8
            btf = -1.8_r8
            ctf =  0._r8
         case ('ben02clim', 'ben02syn', 'fuk95', 'single_column','channel')
            atf = -0.0547_r8
            btf =  0._r8
            ctf =  0._r8
         case ('isomip1', 'isomip2')
            atf = -5.7846e-2_r8
            btf =  1.0307e-1_r8
            ctf = -7.7961e-9_r8
         case default
            if (mnproc == 1) then
               write (lp,'(3a)') ' inieos: expcnf = ', trim(expcnf), &
                                 ' is unsupported!'
            endif
            call xcstop('(inieos)')
                   stop '(inieos)'
      end select

   end subroutine inieos

   pure real(r8) function rho(p, th, s)
   ! ---------------------------------------------------------------------------
   ! In situ density [g cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p,    & ! Pressure [g cm-1 s-2].
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      rho =  ( a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s &
             + (b11 + b12*th + b13*s)*p) &
            /( a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s &
             + (b21 + b22*th + b23*s)*p)

   end function rho

   pure real(r8) function alp(p, th, s)
   ! ---------------------------------------------------------------------------
   ! Specific volume [cm3 g-1].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p,    & ! Pressure [g cm-1 s-2].
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      alp =  ( a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s &
             + (b21 + b22*th + b23*s)*p) &
            /( a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s &
             + (b11 + b12*th + b13*s)*p)

   end function alp

   pure real(r8) function sig(th, s)
   ! ---------------------------------------------------------------------------
   ! Potential density in sigma units [g cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      sig =  (ap11 + (ap12 + ap14*th + ap15*s)*th + (ap13 + ap16*s)*s) &
            /(ap21 + (ap22 + ap24*th + ap25*s)*th + (ap23 + ap26*s)*s)

   end function sig

   pure real(r8) function sig0(th, s)
   ! ---------------------------------------------------------------------------
   ! Potential density in sigma units with reference pressure at the surface
   ! [g cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      sig0 =  (ap110 + (ap120 + ap140*th + ap150*s)*th + (ap130 + ap160*s)*s) &
             /(ap210 + (ap220 + ap240*th + ap250*s)*th + (ap230 + ap260*s)*s)

   end function sig0

   pure real(r8) function dsigdt(th, s)
   ! ---------------------------------------------------------------------------
   ! Derivative of potential density with respect to potential temperature
   ! [g cm-3 K-1].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8) :: r1, r2i

      r1 = ap11 + (ap12 + ap14*th + ap15*s)*th + (ap13 + ap16*s)*s
      r2i = 1._r8/(ap21 + (ap22 + ap24*th + ap25*s)*th + (ap23 + ap26*s)*s)

      dsigdt = ( ap12 + 2._r8*ap14*th + ap15*s &
               - (ap22 + 2._r8*ap24*th + ap25*s)*r1*r2i)*r2i

   end function dsigdt

   pure real(r8) function dsigdt0(th, s)
   ! ---------------------------------------------------------------------------
   ! Derivative of potential density referenced at the surface with respect to
   ! potential temperature [g cm-3 K-1].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8) :: r1, r2i

      r1 = ap110 + (ap120 + ap140*th + ap150*s)*th + (ap130 + ap160*s)*s
      r2i = 1._r8/( ap210 + (ap220 + ap240*th + ap250*s)*th &
                  + (ap230 + ap260*s)*s)

      dsigdt0 = ( ap120 + 2._r8*ap140*th + ap150*s &
                - (ap220 + 2._r8*ap240*th + ap250*s)*r1*r2i)*r2i

   end function dsigdt0

   pure real(r8) function dsigds(th, s)
   ! ---------------------------------------------------------------------------
   ! Derivative of potential density with respect to salinity [kg cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8) :: r1, r2i

      r1 = ap11 + (ap12 + ap14*th + ap15*s)*th + (ap13 + ap16*s)*s
      r2i = 1._r8/(ap21 + (ap22 + ap24*th + ap25*s)*th + (ap23 + ap26*s)*s)

      dsigds = ( ap13 + ap15*th + 2._r8*ap16*s &
               - (ap23 + ap25*th + 2._r8*ap26*s)*r1*r2i)*r2i

   end function dsigds

   pure real(r8) function dsigds0(th, s)
   ! ---------------------------------------------------------------------------
   ! Derivative of potential density referenced at the surface with respect to
   ! salinity [kg cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8) :: r1, r2i

      r1 = ap110 + (ap120 + ap140*th + ap150*s)*th + (ap130 + ap160*s)*s
      r2i = 1._r8/( ap210 + (ap220 + ap240*th + ap250*s)*th &
                  + (ap230 + ap260*s)*s)

      dsigds0 = ( ap130 + ap150*th + 2._r8*ap160*s &
                - (ap230 + ap250*th + 2._r8*ap260*s)*r1*r2i)*r2i

   end function dsigds0

   pure real(r8) function tofsig(sg, s)
   ! ---------------------------------------------------------------------------
   ! Potential temperature as function of potential density and salinity
   ! [deg C].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         sg,   & ! Potental density [g cm-3].
         s       ! Salinity [g kg-1].

      real(r8) :: a, b, c

      a = ap14 - ap24*sg
      b = ap12 - ap22*sg + (ap15 - ap25*sg)*s
      c = ap11 - ap21*sg + (ap13 - ap23*sg + (ap16 - ap26*sg)*s)*s

      tofsig = ( - b - sqrt(b*b - 4._r8*a*c))/(2._r8*a)

   end function tofsig

   pure real(r8) function sofsig(sg, th)
   ! ---------------------------------------------------------------------------
   ! Salinity as function of potential density and potential temperature
   ! [g kg-1].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         sg,   & ! Potental density [g cm-3].
         th      ! Potential temperature [deg C].

      real(r8) :: a, b, c

      a = ap16 - ap26*sg
      b = ap13 - ap23*sg + (ap15 - ap25*sg)*th
      c = ap11 - ap21*sg + (ap12 - ap22*sg + (ap14 - ap24*sg)*th)*th

      sofsig = ( - b + sqrt(b*b - 4._r8*a*c))/(2._r8*a)

   end function sofsig

   pure real(r8) function p_alpha(p1, p2, th, s)
   ! ---------------------------------------------------------------------------
   ! Integral of specific volume with respect to pressure [cm2 s-2].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p1,   & ! Lower pressure bound [g cm-1 s-2].
         p2,   & ! Upper pressure bound [g cm-1 s-2].
         th,   & ! Potential temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8), parameter :: &
         r1_3 = 1._r8/3._r8, &
         r1_5 = 1._r8/5._r8, &
         r1_7 = 1._r8/7._r8, &
         r1_9 = 1._r8/9._r8

      real(r8) :: a1, a2, b1, b2, pm, r, q, qq

      a1 = a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s
      a2 = a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s
      b1 = b11 + b12*th + b13*s
      b2 = b21 + b22*th + b23*s

      ! The analytic solution of the integral is
      !    p_alpha = ( b2*(p2 - p1)
      !              + (a2 - a1*b2/b1)*log((a1 + b1*p2)/(a1 + b1*p1)))/b1
      ! A truncated series expansion of the integral is used that provides
      ! better computational efficiency and accuarcy for most relevant
      ! parameters.

      pm = .5_r8*(p2 + p1)
      r = .5_r8*(p2 - p1)/(a1 + b1*pm)
      q = b1*r
      qq = q*q

      p_alpha = 2._r8*r*( a2 + b2*pm &
                        + (a2 - a1*b2/b1)*qq*( r1_3 &
                                             + qq*( r1_5 &
                                                  + qq*( r1_7 &
                                                       + qq*  r1_9))))

   end function p_alpha

   pure real(r8) function p_p_alpha(p1, p2, th, s)
   ! ---------------------------------------------------------------------------
   ! Double integral of specific volume with respect to pressure [cm g s-4].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p1,   & ! Lower pressure bound [g cm-1 s-2].
         p2,   & ! Upper pressure bound [g cm-1 s-2].
         th,   & ! Potential temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8), parameter :: &
         r1_3  = 1._r8/3._r8, &
         r1_5  = 1._r8/5._r8, &
         r1_7  = 1._r8/7._r8, &
         r1_9  = 1._r8/9._r8, &
         r1_10 = 1._r8/10._r8

      real(r8) :: a1, a2, b1, b2, pm, dp, r, q

      a1 = a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s
      a2 = a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s
      b1 = b11 + b12*th + b13*s
      b2 = b21 + b22*th + b23*s

      ! The analytic solution of the integral is
      !    p_p_alpha = 
      !       ( .5*b2*(p2 - p1)**2
      !       + (a2 - a1*b2/b1)*( (a1/b1 + p2)*log((a1 + b1*p2)/(a1 + b1*p1))
      !                         - (p2 - p1)))/b1
      ! A truncated series expansion of the integral is used that provide
      ! better computational efficiency and accuarcy for most relevant
      ! parameters.

      pm = .5_r8*(p2 + p1)
      dp = .5_r8*(p2-p1)
      r = dp/(a1 + b1*pm)
      q = b1*r

      p_p_alpha = 2._r8*dp*r*( a2 + b2*pm &
                             + (a2-a1*b2/b1)*q*( r1_3 + q*(r1_3 &
                                               + q*( r1_5 + q*(r1_5 &
                                                   + q*( r1_7 + q*(r1_7 &
                                                       + q*( r1_9 + q*(r1_9 &
                                                           + q* r1_10)))))))))

   end function p_p_alpha

   subroutine delphi(p1, p2, th, s, dphi, alp1, alp2)
   ! ---------------------------------------------------------------------------
   ! Integrate specific volume with respect to pressure to find the difference
   ! in geopotential between two pressure levels
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p1,   & ! Lower pressure bound [g cm-1 s-2].
         p2,   & ! Upper pressure bound [g cm-1 s-2].
         th,   & ! Potential temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8), intent(out) :: &
         dphi, & ! Geopotential difference [cm2 s-2].
         alp1, & ! Specific volume at lower pressure bound [cm3 g-1].
         alp2    ! Specific volume at upper pressure bound [cm3 g-1].
      
      real(r8), parameter :: &
         r1_3 = 1._r8/3._r8, &
         r1_5 = 1._r8/5._r8, &
         r1_7 = 1._r8/7._r8, &
         r1_9 = 1._r8/9._r8

      real(r8) :: a1, a2, b1, b2, pm, r, q, qq

      a1 = a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s
      a2 = a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s
      b1 = b11 + b12*th + b13*s
      b2 = b21 + b22*th + b23*s

      ! The analytic solution of the integral is
      !    dphi = - ( b2*(p2 - p1)
      !             + (a2 - a1*b2/b1)*log((a1 + b1*p2)/(a1 + b1*p1)))/b1
      ! A truncated series expansion of the integral is used that provides
      ! better computational efficiency and accuarcy for most relevant
      ! parameters.

      pm = .5_r8*(p2 + p1)
      r = .5_r8*(p2 - p1)/(a1 + b1*pm)
      q = b1*r
      qq = q*q

      dphi = - 2._r8*r*( a2 + b2*pm &
                       + (a2 - a1*b2/b1)*qq*( r1_3 &
                                            + qq*( r1_5 &
                                                 + qq*( r1_7 &
                                                      + qq*  r1_9))))

      alp1 = (a2 + b2*p1)/(a1 + b1*p1)
      alp2 = (a2 + b2*p2)/(a1 + b1*p2)

   end subroutine delphi

end module mod_eos
