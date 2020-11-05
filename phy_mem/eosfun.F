! ------------------------------------------------------------------------------
! Copyright (C) 2007-2018 Mats Bentsen
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
c
c --- coefficients for the functional fit of in situ density
c
      real ::  
     . a11= 9.9985372432159340e-01,
     . a12= 1.0380621928183473e-02,
     . a13= 1.7073577195684715e-03,
     . a14=-3.6570490496333680e-05,
     . a15=-7.3677944503527477e-06,
     . a16=-3.5529175999643348e-06,
     . b11= 1.7083494994335439e-10,
     . b12= 7.1567921402953455e-13,
     . b13= 1.2821026080049485e-13,
     . a21= 1.0                   ,
     . a22= 1.0316374535350838e-02,
     . a23= 8.9521792365142522e-04,
     . a24=-2.8438341552142710e-05,
     . a25=-1.1887778959461776e-05,
     . a26=-4.0163964812921489e-06,
     . b21= 1.1995545126831476e-10,
     . b22= 5.5234008384648383e-13,
     . b23= 8.4310335919950873e-14,
c
c --- reference pressure (dyn/cm^2)
c
     . pref=2000.e5,

     . ap11,ap12,ap13,ap14,ap15,ap16,
     . ap21,ap22,ap23,ap24,ap25,ap26,
     . ap110,ap120,ap130,ap140,ap150,ap160,
     . ap210,ap220,ap230,ap240,ap250,ap260,
     . atf,btf,ctf    

      
c$acc declare copyin(a11,a12,a13,a14,a15,a16,b11,b12,b13
c$acc.    ,a21,a22,a23,a24,a25,a26,b21,b22,b23
c$acc.    ,pref) create(ap11,ap12,ap13,ap14,ap15,ap16
c$acc.    ,ap21,ap22,ap23,ap24,ap25,ap26
c$acc.    ,ap110,ap120,ap130,ap140,ap150,ap160
c$acc.    ,ap210,ap220,ap230,ap240,ap250,ap260
c$acc.    ,atf,btf,ctf)
      
      end module



      
      subroutine inieos
c
c --- ------------------------------------------------------------------
c --- Initialize coefficients related to the equation of state
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_eos
c
      implicit none
c
#include "common_blocks.h"
cinclude "common_eos.h"
c
c --- In situ density [kg/m^3] as a function of pressure, potential
c --- temperature and salinity is approximated by the functional form
c ---   rho(p,th,s)=P1(p,th,s)/P2(p,th,s)
c --- where
c ---   P1(p,th,s)=a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s+(b11+b12*th+b13*s)*p
c --- and
c ---   P2(p,th,s)=a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s+(b21+b22*th+b23*s)*p
c --- Here we compute the coefficients needed for an expression for
c --- potential density [g/cm^3] in sigma units of the form
c ---   sig(th,s)=R1(th,s)/R2(th,s) 
c --- where
c ---   R1(p,th,s)=ap11+(ap12+ap14*th+ap15*s)*th+(ap13+ap16*s)*s
c --- and
c ---   R2(p,th,s)=ap21+(ap22+ap24*th+ap25*s)*th+(ap23+ap26*s)*s
c
      ap21=a21+b21*pref
      ap22=a22+b22*pref
      ap23=a23+b23*pref
      ap24=a24
      ap25=a25
      ap26=a26
      ap11=a11+b11*pref-ap21
      ap12=a12+b12*pref-ap22
      ap13=a13+b13*pref-ap23
      ap14=a14-ap24
      ap15=a15-ap25
      ap16=a16-ap26
c
      ap210=a21
      ap220=a22
      ap230=a23
      ap240=a24
      ap250=a25
      ap260=a26
      ap110=a11-ap210
      ap120=a12-ap220
      ap130=a13-ap230
      ap140=a14-ap240
      ap150=a15-ap250
      ap160=a16-ap260
c
c --- Coefficients for freezing temperature
      if     (expcnf.eq.'cesm') then
        atf=0.
        btf=-1.8
        ctf=0.
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        atf=-0.0547
        btf=0.
        ctf=0.
      elseif (expcnf.eq.'isomip1'.or.expcnf.eq.'isomip2') then
        atf=-5.7846e-2
        btf= 1.0307e-1
        ctf=-7.7961e-9
      elseif (expcnf.eq.'channel') then
        atf=-0.0547
        btf=0.
        ctf=0.
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(inieos)')
               stop '(inieos)'
      endif

c$acc update device(
c$acc.    ap11,ap12,ap13,ap14,ap15,ap16
c$acc.    ,ap21,ap22,ap23,ap24,ap25,ap26
c$acc.    ,ap110,ap120,ap130,ap140,ap150,ap160
c$acc.    ,ap210,ap220,ap230,ap240,ap250,ap260
c$acc.    ,atf,btf,ctf)

      return
      end
c
c --- ------------------------------------------------------------------
c
      real function rho(p,th,s)
c
!$acc routine(rho) seq
      use mod_eos
      implicit none
c
      real p,th,s
c
cinclude "common_eos.h"
c
      rho=(a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s+(b11+b12*th+b13*s)*p)
     .   /(a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s+(b21+b22*th+b23*s)*p)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function alp(p,th,s)
      use mod_eos
      implicit none
c
      real p,th,s
c
cinclude "common_eos.h"
c
      alp=(a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s+(b21+b22*th+b23*s)*p)
     .   /(a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s+(b11+b12*th+b13*s)*p)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function sig(th,s)
!$acc routine(sig) seq

c --- potential density [g/cm^3] in sigma units as a function of
c --- potential temperature and salinity
c
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      sig=(ap11+(ap12+ap14*th+ap15*s)*th+(ap13+ap16*s)*s)
     .   /(ap21+(ap22+ap24*th+ap25*s)*th+(ap23+ap26*s)*s)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function sig0(th,s)
c
c --- potential density [g/cm^3] with reference pressure at the surface
c --- in sigma units as a function of potential temperature and salinity
c
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      sig0=(ap110+(ap120+ap140*th+ap150*s)*th+(ap130+ap160*s)*s)
     .    /(ap210+(ap220+ap240*th+ap250*s)*th+(ap230+ap260*s)*s)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dsigdt(th,s)
c$acc routine(dsigdt) seq
c
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      real r1,r2i
c
      r1=ap11+(ap12+ap14*th+ap15*s)*th+(ap13+ap16*s)*s
      r2i=1./(ap21+(ap22+ap24*th+ap25*s)*th+(ap23+ap26*s)*s)
c
      dsigdt=(ap12+2.*ap14*th+ap15*s
     .       -(ap22+2.*ap24*th+ap25*s)*r1*r2i)*r2i
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dsigdt0(th,s)
c
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      real r1,r2i
c
      r1=ap110+(ap120+ap140*th+ap150*s)*th+(ap130+ap160*s)*s
      r2i=1./(ap210+(ap220+ap240*th+ap250*s)*th+(ap230+ap260*s)*s)
c
      dsigdt0=(ap120+2.*ap140*th+ap150*s
     .        -(ap220+2.*ap240*th+ap250*s)*r1*r2i)*r2i
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dsigds(th,s)
c$acc routine(dsigds) seq
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      real r1,r2i
c
      r1=ap11+(ap12+ap14*th+ap15*s)*th+(ap13+ap16*s)*s
      r2i=1./(ap21+(ap22+ap24*th+ap25*s)*th+(ap23+ap26*s)*s)
c
      dsigds=(ap13+ap15*th+2.*ap16*s
     .       -(ap23+ap25*th+2.*ap26*s)*r1*r2i)*r2i
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dsigds0(th,s)
c
      use mod_eos
      implicit none
c
      real th,s
c
cinclude "common_eos.h"
c
      real r1,r2i
c
      r1=ap110+(ap120+ap140*th+ap150*s)*th+(ap130+ap160*s)*s
      r2i=1./(ap210+(ap220+ap240*th+ap250*s)*th+(ap230+ap260*s)*s)
c
      dsigds0=(ap130+ap150*th+2.*ap160*s
     .        -(ap230+ap250*th+2.*ap260*s)*r1*r2i)*r2i
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function tofsig(sg,s)
      use mod_eos
      implicit none
c
      real sg,s
c
cinclude "common_eos.h"
c
      real a,b,c
c
      a=ap14-ap24*sg
      b=ap12-ap22*sg+(ap15-ap25*sg)*s
      c=ap11-ap21*sg+(ap13-ap23*sg+(ap16-ap26*sg)*s)*s
c
      tofsig=(-b-sqrt(b*b-4.*a*c))/(2.*a)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function sofsig(sg,th)
!$acc routine(sofsig) seq      
c
      use mod_eos
      implicit none
c
      real sg,th
c
cinclude "common_eos.h"
c
      real a,b,c
c
      a=ap16-ap26*sg
      b=ap13-ap23*sg+(ap15-ap25*sg)*th
      c=ap11-ap21*sg+(ap12-ap22*sg+(ap14-ap24*sg)*th)*th
c
      sofsig=(-b+sqrt(b*b-4.*a*c))/(2.*a)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function p_alpha(p1,p2,th,s)
c
c --- integrate specific volume with respect to pressure
c
      use mod_eos
      implicit none
c
      real p1,p2,th,s
c
cinclude "common_eos.h"
c
      real r1_3,r1_5,r1_7,r1_9
      parameter (r1_3=1./3.,r1_5=1./5.,r1_7=1./7.,r1_9=1./9.)
c
      real a1,a2,b1,b2,pm,r,q,qq
c
      a1=a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s
      a2=a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s
      b1=b11+b12*th+b13*s
      b2=b21+b22*th+b23*s
c
c --- the analytic solution of the integral is
c       p_alpha=(b2*(p2-p1)
c               +(a2-a1*b2/b1)*log((a1+b1*p2)/(a1+b1*p1)))/b1
c --- a truncated series expansion of the integral is used that provide
c --- better computational efficiency and accuarcy for most relevant
c --- parameters
c
      pm=.5*(p2+p1)
      r=.5*(p2-p1)/(a1+b1*pm)
      q=b1*r
      qq=q*q
c
      p_alpha=2.*r*(a2+b2*pm
     .             +(a2-a1*b2/b1)*qq*(r1_3+qq*(r1_5+qq*(r1_7+qq*r1_9))))
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function p_p_alpha(p1,p2,th,s)
c
c --- double integration of specific volume with respect to pressure
c
      use mod_eos
      implicit none
c
      real p1,p2,th,s
c
cinclude "common_eos.h"
c
      real r1_3,r1_5,r1_7,r1_9,r1_10
      parameter (r1_3=1./3.,r1_5=1./5.,r1_7=1./7.,r1_9=1./9.,
     .           r1_10=1./10.)
c
      real a1,a2,b1,b2,pm,dp,r,q
c
      a1=a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s
      a2=a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s
      b1=b11+b12*th+b13*s
      b2=b21+b22*th+b23*s
c
c --- the analytic solution of the integral is
c       p_p_alpha=(.5*b2*(p2-p1)**2
c                 +(a2-a1*b2/b1)*((a1/b1+p2)*log((a1+b1*p2)/(a1+b1*p1))
c                                -(p2-p1)))/b1
c --- a truncated series expansion of the integral is used that provide
c --- better computational efficiency and accuarcy for most relevant
c --- parameters
c
      pm=.5*(p2+p1)
      dp=.5*(p2-p1)
      r=dp/(a1+b1*pm)
      q=b1*r
c
      p_p_alpha=2.*dp*r*(a2+b2*pm
     .                  +(a2-a1*b2/b1)*q*(r1_3+q*(r1_3+
     .                                 q*(r1_5+q*(r1_5+
     .                                 q*(r1_7+q*(r1_7+
     .                                 q*(r1_9+q*(r1_9+
     .                                 q* r1_10)))))))))
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine delphi(p1,p2,th,s,dphi,alp1,alp2)
c
c --- integrate specific volume with respect to pressure to find the
c --- difference in geopotential between two pressure levels
c
      use mod_eos
      implicit none
c
      real p1,p2,th,s,dphi,alp1,alp2
c
cinclude "common_eos.h"
c
      real r1_3,r1_5,r1_7,r1_9
      parameter (r1_3=1./3.,r1_5=1./5.,r1_7=1./7.,r1_9=1./9.)
c
      real a1,a2,b1,b2,pm,r,q,qq
c
      a1=a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s
      a2=a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s
      b1=b11+b12*th+b13*s
      b2=b21+b22*th+b23*s
c
c --- the analytic solution of the integral is
c       dphi=-(b2*(p2-p1)
c             +(a2-a1*b2/b1)*log((a1+b1*p2)/(a1+b1*p1)))/b1
c --- a truncated series expansion of the integral is used that provide
c --- better computational efficiency and accuarcy for most relevant
c --- parameters
c
      pm=.5*(p2+p1)
      r=.5*(p2-p1)/(a1+b1*pm)
      q=b1*r
      qq=q*q
c
      dphi=-2.*r*(a2+b2*pm
     .           +(a2-a1*b2/b1)*qq*(r1_3+qq*(r1_5+qq*(r1_7+qq*r1_9))))
c
      alp1=(a2+b2*p1)/(a1+b1*p1)
      alp2=(a2+b2*p2)/(a1+b1*p2)
c
      return
      end
