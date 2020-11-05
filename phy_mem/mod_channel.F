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

      module mod_channel
c
c --- ------------------------------------------------------------------
c --- This module contains routines for generating surface fluxes
c --- from NCEP or ERA40 reanalysis following Bentsen and Drange (2002).
c --- ------------------------------------------------------------------
c
      use mod_xc
c      use mod_ben02func
c      use mod_bulktf
c
      implicit none
c
c      private
c
      public :: getfrc_channel
c
      contains
c
c
c
c
      subroutine getfrc_channel
c
c --- ------------------------------------------------------------------
c --- Get climatological forcing
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_forc.h"
#include "common_blocks.h"
c
c      intrinsic tanh
c
      integer i,j,l,k,ios
c      real, parameter :: cwidth=75.e5
c      real, parameter :: ztx0=-0.05
      real :: ztx0,mty0
      logical :: fexist
c      
      namelist /idlfor/ ztx0,mty0
c
c --- READ NAMELIST FILE
      inquire(file='limits',exist=fexist)
      if (fexist) then
          open (unit=nfu,file='limits',status='old',action='read')
      else
          write (lp,*) 'getfrc_channel: could not find namelist file!'
          call xchalt('(getfrc_channel)')
          stop '(getfrc_channel)'
      endif
c
c --- READ AND BROADCAST
      read (unit=nfu,nml=idlfor,iostat=ios)
      close (unit=nfu)
      call xcbcst(ztx0)
      call xcbcst(mty0)

c      call forfun_channel !not needed

c --- THESE WILL BE SET TO ZERO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          swa(i,j)=0.    ! shortwave
          nsf(i,j)=0.    ! non-solar
          hmltfz(i,j)=0. ! heat flux due to melting/freezing
          hmlt(i,j)=0.   ! heat flux due to melting
          dfl(i,j)=0.    ! derivate of non-solar in respect to T
c
          alb(i,j)=0. ! albedo
          eva(i,j)=0. ! evaporation
          lip(i,j)=0. ! liquid precip
          sop(i,j)=0. ! solid precip 
c
          rnf(i,j)=0.       ! runoff
          rfi(i,j)=0.       ! runoff ice
          fmltfz(i,j)=0.    ! fresh water flux due to melting/freezing
          sfl(i,j)=0.       ! salt flux
          ustarw(i,j)=0.005 ! friction velocity for open water
          slp(i,j)=1000.    ! sea level pressure
          abswnd(i,j)=0.    ! wind speed at measurement height -zu-
          albw(i,j)=0.      ! daily mean open water albedo
          frzpot(i,j)=0.    ! freezing potential
          mltpot(i,j)=0.    ! melting potential
          atmco2(i,j)=400.  ! atmospheric co2 concentration
          flxco2(i,j)=0.    ! air-sea co2 flux
          flxdms(i,j)=0.    ! sea-air dms flux
c
c --- --------------------------------------- ---
c --- DEFINE SST CLIMATOLOGY AND WIND STRESS
c --- THESE SHOULD BE CONSISTENT WITH INITIAL
C --- CONDITIONS
c --- --------------------------------------- ---
c
          do k=1,12
            sstclm(i,j,k)=10.+3.
            sssclm(i,j,k)=35.
          enddo
        enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
C ---   CASE 12
          ztx(i,j)=ztx0
c
c ---   CASES 8-11
c ---   THIS SETUP IS CONSISTENT WITH SLOPE DEFINITION
c ---   Divergence/Convergence peaks at the slope midpoint 
c ---   leads to upwelling over the southern slope and
c ---   downwelling over the northern slope 
c
c           if ((scpy(i,j)*(j0+j)).lt.(2.*cwidth)) then
c               ztx(i,j)=
c     .         0.5*ztx0*(1+tanh(pi*scpy(i,j)*(j0+j)/
c     .         (0.5*cwidth)-pi))
c           elseif ((jtdm-(j0+j))*scpy(i,j).lt.(2.*cwidth)) then
c               ztx(i,j)=
c     .         0.5*ztx0*(1+tanh(pi*scpy(i,j)*(jtdm-(j0+j))/
c     .         (0.5*cwidth)-pi))
c           else
c               ztx(i,j)=ztx0
c           endif
c ---   --------------------------------------------------------
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          mty(i,j)=mty0 !0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end subroutine getfrc_channel
c
c --- ------------------------------------------------------------------
c

      end module mod_channel
