      subroutine geoenv_test
c
c --- ------------------------------------------------------------------
c --- Define bathymetry, grid specification and Coriolis parameter for
c --- test case
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
c
      integer i,j,k,l
c
      depths=0.
      nwp=0
      qlon=0.
      qlat=0.
      plon=0.
      plat=0.
      ulon=0.
      ulat=0.
      vlon=0.
      vlat=0.
      qclon=0.
      qclat=0.
      pclon=0.
      pclat=0.
      uclon=0.
      uclat=0.
      vclon=0.
      vclat=0.
      scqx=0.
      scqy=0.
      scpx=0.
      scpy=0.
      scux=0.
      scuy=0.
      scvx=0.
      scvy=0.
      scq2=0.
      scp2=0.
      scu2=0.
      scv2=0.
      angle=0.
      corioq=0.
      coriop=0.
      betafp=0.
      cosang=0.
      sinang=0.
c
      return
      end
