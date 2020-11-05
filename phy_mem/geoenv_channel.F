! ------------------------------------------------------------------------------
! Copyright (C) 2020 Aleksi Nummelin, Mats Bentsen
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

      subroutine geoenv_channel
c
c --- ------------------------------------------------------------------
c --- Define bathymetry, grid specification and Coriolis parameter for
c --- test case [cm s]
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
c
      intrinsic random_seed, random_number,tanh,sin

      integer i,j,l,ios
      integer, parameter :: ncorru=10
      real :: r
      real :: sldepth,sfdepth,rdepth,cwidth,swidth,scxy,corio0,beta0,
     .        d_corru
      real, dimension(ncorru) :: acorru=0.0,wlcorru=0.0
      logical :: fexist
c      real, parameter :: PI0=4.*atan(1.0)
c      real, parameter :: cwidth=75.E5 !slope width
c      real, parameter :: swidth=150.E5 !slope mid-point distance 
c
      namelist /idlgeo/ sldepth,sfdepth,rdepth,acorru,wlcorru,cwidth,
     .                  swidth,scxy,corio0,beta0

      inquire(file='limits',exist=fexist)
      if (fexist) then
          open (unit=nfu,file='limits',status='old',action='read')
      else
          write (lp,*) 'geoenv_channel: could not find namelist file!'
          call xchalt('(geoenv_channel)')
          stop '(geoenv_channel)'
      endif
c
      read (unit=nfu,nml=idlgeo,iostat=ios)
      close (unit=nfu)
c
c --- broadcast the namelist variables
      call xcbcst(sldepth)
      call xcbcst(sfdepth)
      call xcbcst(rdepth)
      call xcbcst(acorru)
      call xcbcst(wlcorru)
      call xcbcst(cwidth)
      call xcbcst(swidth)
      call xcbcst(scxy)
      call xcbcst(corio0)
      call xcbcst(beta0)
c
c --- initialize the random numbe generator
      call random_seed()

c --- initialize depth to 0 everywhere     
      depths=0.
c
c --- Number of wet points (southern and northern most rows are land)
      nwp=jtdm*itdm-2.*itdm !jj*ii-2.*ii
c --- the latitudes need to be set, but actually the values don't matter it seems      
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
c --- these are the grid dimensions (in cm!)
      scqx=scxy !2.e5
      scqy=scxy !2.e5
      scpx=scxy !2.e5
      scpy=scxy !2.e5
      scux=scxy !2.e5
      scuy=scxy !2.e5
      scvx=scxy !2.e5
      scvy=scxy !2.e5
c --- square of the grid dimensions
      scq2=scqx*scqy
      scp2=scpx*scpy
      scu2=scux*scuy
      scv2=scvx*scvy
c --- here we just set a f-plane (beta=0) 
c      angle=0.
      corioq=corio0 !1.e-4
      coriop=corio0 !1.e-4
      betafp=beta0  !0.
      angle=0. !we leave these hard coded - grid will be cartesian
      cosang=1.
      sinang=0.
c
c --- this is now a tanh slope that will have the same shape
c --- independent of the grid size - i.e. depends on grid
c --- dimensions (no interpolation done though)
c --- 
c$OMP PARALLEL DO PRIVATE(i,r)
      do j=1,jj
        if (j0+j.gt.1) then
        if (j0+j.lt.jtdm) then
          do i=1,ii
             call random_number(r)
             if ((scpy(i,j)*(j0+j)).lt.(swidth+cwidth)) then
               l=1
               d_corru=0.0
               do while (acorru(l).gt.0.0)
                 d_corru=d_corru
     .           +acorru(l)*sin(2.*pi*scpx(i,j)*(i0+i)/wlcorru(l))
                 l=l+1
               enddo
               depths(i,j)=
     .         sfdepth+rdepth*r+0.5*sldepth*(1.+tanh(pi*
     .         (scpy(i,j)*(j0+j)-swidth
     .         -d_corru)/cwidth))
c     .         -acorru*sin(2.*pi*scpx(i,j)*(i0+i)/wlcorru))/cwidth))
c     .         250.+10.*r+1.e3*(1+tanh(pi*scpy(i,j)*(j0+j)/
c     .         cwidth-(swidth/cwidth)*pi))
             elseif ((jtdm-(j0+j))*scpy(i,j).lt.(swidth+cwidth)) then
               l=1
               d_corru=0.0
               do while (acorru(l).gt.0.0)
                 d_corru=d_corru
     .           +acorru(l)*sin(2.*pi*scpx(i,j)*(i0+i)/wlcorru(l))
                 l=l+1
               enddo
               depths(i,j)=
     .         sfdepth+rdepth*r+0.5*sldepth*(1.+tanh(pi*
     .         (scpy(i,j)*(jtdm-(j0+j))-swidth
     .         -d_corru)/cwidth))
c     .         -acorru*sin(2.*pi*scpx(i,j)*(i0+i)/wlcorru))/cwidth))
c     .         250.+10.*r+1.e3*(1+tanh(pi*scpy(i,j)*(jtdm-(j0+j))/
c     .         cwidth-(swidth/cwidth)*pi))
             else
c               depths(i,j)=250.+2.e3+10.*r
               depths(i,j)=sfdepth+rdepth*r+sldepth
             endif
          enddo
        endif
        endif
      enddo
c$OMP END PARALLEL DO

      return
      end
