      subroutine eosini
c
      implicit none
c
#include "common_eos.h"
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
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function rho(p,th,s)
c
      implicit none
c
      real p,th,s
c
#include "common_eos.h"
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
c
      implicit none
c
      real p,th,s
c
#include "common_eos.h"
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
c
c --- potential density [g/cm^3] in sigma units as a function of
c --- potential temperature and salinity
c
      implicit none
c
      real th,s
c
#include "common_eos.h"
c
      sig=(ap11+(ap12+ap14*th+ap15*s)*th+(ap13+ap16*s)*s)
     .   /(ap21+(ap22+ap24*th+ap25*s)*th+(ap23+ap26*s)*s)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dsigdt(th,s)
c
      implicit none
c
      real th,s
c
#include "common_eos.h"
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
      real function dsigds(th,s)
c
      implicit none
c
      real th,s
c
#include "common_eos.h"
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
      real function tofsig(sg,s)
c
      implicit none
c
      real sg,s
c
#include "common_eos.h"
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
c
      implicit none
c
      real sg,th
c
#include "common_eos.h"
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
      subroutine delphi(p1,p2,th,s,dphi,alp1,alp2)
c
c --- integrate specific volume with respect to pressure to find the
c --- difference in geopotential between two pressure levels
c
      implicit none
c
      real p1,p2,th,s,dphi,alp1,alp2
c
#include "common_eos.h"
c
      real a1,a2,b1,b2,b1i,a1b1i
c
      a1=a11+(a12+a14*th+a15*s)*th+(a13+a16*s)*s
      a2=a21+(a22+a24*th+a25*s)*th+(a23+a26*s)*s
      b1=b11+b12*th+b13*s
      b2=b21+b22*th+b23*s

      b1i=1./b1
      a1b1i=a1*b1i
c
      dphi=-(b2*(p2-p1)+(a2-a1b1i*b2)*log((a1b1i+p2)/(a1b1i+p1)))*b1i
c
      alp1=(a2+b2*p1)/(a1+b1*p1)
      alp2=(a2+b2*p2)/(a1+b1*p2)
c
      return
      end
