! ------------------------------------------------------------------------------
! Copyright (C) 2021 Aleksi Nummelin, Mats Bentsen
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
! ------------------------------------------------------------------------------
! This module contains variables and routines related to a periodic channel
! experiment with 'continental shelves and slopes' on the southern and northern
! sides of the domain. This setup was originally created for the KeyCLIM project
! to study eddy-topography interactions.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, rearth, pi, radian
   use mod_xc
   use mod_grid, only: sigmar, &
                       qclon, qclat, pclon, pclat, uclon, uclat, vclon, vclat, &
                       scqx, scqy, scpx, scpy, scux, scuy, scvx, scvy, &
                       scq2, scp2, scu2, scv2, &
                       qlon, qlat, plon, plat, ulon, ulat, vlon, vlat, &
                       depths, corioq, coriop, betafp, angle, cosang, sinang, &
                       nwp

   use mod_eos, only: tofsig
   use mod_ben02, only: ntda, alb, albw, dfl
   use mod_forcing, only: surflx, surrlx, sswflx, salflx, brnflx, salrlx, &
                          taux, tauy, ustarw, slp, swa, atmco2, nsf, mty, &
                          ztx, hmltfz, eva, lip, sop, rnf, rfi, & 
                          fmltfz, sfl, abswnd, flxco2, flxdms, sstclm, sssclm
   !use mod_mxlayr, only: mltmin
   use mod_state, only:  temp, saln, sigma, phi
   use mod_checksum, only: csdiag, chksummsk
   
   implicit none
   
   private
   
   public :: geoenv_channel, inifrc_channel, ictsz_channel
   
contains
   
   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------
   
   subroutine geoenv_channel
   ! ---------------------------------------------------------------------------
   ! Define bathymetry, grid specification and Coriolis parameter for the
   ! channel configuration
   ! ---------------------------------------------------------------------------
      intrinsic random_seed, random_number, tanh, sin
 
      integer, parameter :: ncorru=10
      real(r8), dimension(ncorru) :: acorru, wlcorru
      real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: r0
      real(r8), dimension(itdm,jtdm) :: rtmp
      real(r8) :: sldepth,sfdepth,rdepth,cwidth,swidth,scxy, &
                  corio0, beta0, d_corru, r
      integer :: i,j,l,ios
      integer, dimension(:), allocatable :: seed
      logical :: fexist
      
      ! Read parameters from the namelist
      namelist /idlgeo/ sldepth,sfdepth,rdepth,acorru,wlcorru, &
                        cwidth,swidth,scxy,corio0,beta0
      acorru(:)=0._r8
      wlcorru(:)=0._r8
      inquire(file='limits', exist=fexist)
      if (fexist) then
         open (unit=nfu,file='limits',status='old',action='read')
      else
         write (lp,*) 'geoenv_channel: could not find namelist file!'
         call xchalt('(geoenv_channel)')
         stop '(geoenv_channel)'
      endif
      
      read (unit=nfu,nml=idlgeo,iostat=ios)
      close (unit=nfu)
      ! Done reading
      
      ! broadcast the namelist variables
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
      !
      ! Number of wet points (southern and northern most rows are land)
      nwp=jtdm*itdm-2*itdm
      !
      if (mnproc == 1) then
        call random_seed(size = i)
        allocate(seed(i))
        seed = 1144153914 !hard-coded seed
        call random_seed(put = seed)
        call random_number(rtmp)
      endif
      call xcaput(rtmp, r0, 1)
      !$omp parallel do private(i)
         do j = 1, jj
            do i = 1, ii
               ! Set lon/lat information to 0. We will use cartesian
               ! coordinates and control coriolis and beta instead
               qlon(i,j)=0._r8
               qlat(i,j)=0._r8
               plon(i,j)=0._r8
               plat(i,j)=0._r8
               ulon(i,j)=0._r8
               ulat(i,j)=0._r8
               vlon(i,j)=0._r8
               vlat(i,j)=0._r8
               qclon(i,j,:)=0._r8
               qclat(i,j,:)=0._r8
               pclon(i,j,:)=0._r8
               pclat(i,j,:)=0._r8
               uclon(i,j,:)=0._r8
               uclat(i,j,:)=0._r8
               vclon(i,j,:)=0._r8
               vclat(i,j,:)=0._r8
               ! Grid dimensions (in cm!)
               scqx(i,j)=scxy
               scqy(i,j)=scxy
               scpx(i,j)=scxy
               scpy(i,j)=scxy
               scux(i,j)=scxy
               scuy(i,j)=scxy
               scvx(i,j)=scxy
               scvy(i,j)=scxy
               ! Square of the grid dimensions
               scq2(i,j)=scqx(i,j)*scqy(i,j)
               scp2(i,j)=scpx(i,j)*scpy(i,j)
               scu2(i,j)=scux(i,j)*scuy(i,j)
               scv2(i,j)=scvx(i,j)*scvy(i,j)
               ! Namelist controls coriolis and beta
               corioq=corio0
               coriop=corio0
               betafp=beta0
               ! Cartesian grid hard coded
               angle=0._r8
               cosang=1._r8
               sinang=0._r8
               ! initialize depth to 0
               depths(i,j)=0._r8               
            enddo
         enddo
      !$omp end parallel do
      
      ! Set the bottom topography to be a tanh function.
      ! The resulting slope will have the same shape independent of the 
      ! grid size (no interpolation done though).
      !$omp parallel do private(i,r,l,d_corru)
         do j=1,jj
            if (j0+j.gt.1) then
            if (j0+j.lt.jtdm) then
               do i=1,ii
                  !r=r0(i,j)-0.5_r8
                  if ((scpy(i,j)*(j0+j)).lt.(swidth+cwidth)) then
                     l=1
                     d_corru=0._r8
                     do while (acorru(l).gt.0._r8)
                        d_corru=d_corru &
                        +acorru(l)*sin(2._r8*pi*scpx(i,j)*(i0+i)/wlcorru(l))
                        l=l+1
                     enddo
                     depths(i,j) = sfdepth+rdepth*r0(i,j)+.5_r8*sldepth* &
                                  (1._r8+tanh(pi*(scpy(i,j)*(j0+j)- &
                                  swidth-d_corru)/cwidth))
                  elseif ((jtdm-(j0+j))*scpy(i,j).lt.(swidth+cwidth)) then
                     l=1
                     d_corru=0._r8
                     do while (acorru(l).gt.0._r8)
                        d_corru=d_corru &
                                +acorru(l)*sin(2._r8*pi*scpx(i,j)*(i0+i)/ &
                                wlcorru(l))
                        l=l+1
                     enddo
                     depths(i,j) = sfdepth+rdepth*r0(i,j)+.5_r8*sldepth* &
                                   (1._r8+tanh(pi*(scpy(i,j)*(jtdm-(j0+j)) &
                                   -swidth-d_corru)/cwidth))
                  else
                     depths(i,j)=sfdepth+rdepth*r0(i,j)+sldepth
                  endif
               enddo
            endif
            endif
         enddo
      !$omp end parallel do
      
      end subroutine geoenv_channel
      
      subroutine ictsz_channel
      ! define layer temperature and salinity, and geopotential at layer interfaces.
         real(r8), dimension(1 - nbdy:idm + nbdy, &
                             1 - nbdy:jdm + nbdy, kdm + 1) :: z
         real(r8), dimension(1 - nbdy:idm + nbdy, &
                             1 - nbdy:jdm + nbdy, kdm) :: dz
         real(r8), dimension(kdm) :: sigmr0, dz0
         real(r8) :: S0,sig0,sig0dz,sigdz,sigscl,dztop,dzmax,dzscl
         integer i,j,k,l,ios
         logical :: fexist
         
         namelist /idlini/ S0,sig0,sig0dz,sigdz,sigscl,dztop,dzmax,dzscl
         
         if (mnproc.eq.1) then
            write (lp,'(1a)') ' idealized initial conditions'
            call flush(lp)
         endif
         
         ! 
         inquire(file='limits',exist=fexist)
         if (fexist) then
            open (unit=nfu,file='limits',status='old',action='read')
         else
            write (lp,*) 'ictsz_channel: could not find namelist file!'
            call xchalt('(ictsz_channel)')
            stop '(ictsz_channel)'
         endif
         !
         read (unit=nfu,nml=idlini,iostat=ios)
         close (unit=nfu)
         call xcbcst(ios)
         !
         ! Idealized functional form for the sigma profile.
         ! Parameters are defined by the user namelist
         sigmr0(1)=sig0
         sigmr0(2)=sig0
         dz0(1)=dztop
         dz0(2)=dztop
         do k=3,kk
            sigmr0(k) = sigmr0(k-1) + sig0dz + sigdz* &
            (1.-tanh(sigscl*pi*(k-1)/kk))
            dz0(k) = dzmax*tanh(dzscl*pi*(k-1)/kk)
         enddo
      !$omp parallel do private(k, l, i)
         do j=1,jj
            do k=1,kk
               do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                     dz(i,j,k) = dz0(k)
                     saln(i,j,k) = S0
                     sigmar(i,j,k) = sigmr0(k)*1.e-3_r8 !convert units to g/cm^3
                     temp(i,j,k) = tofsig(sigmar(i,j,k),saln(i,j,k))
                  enddo
               enddo
            enddo
         enddo
      !$omp end parallel do
      !
      !$omp parallel do private(l, i)
         do j=1,jj
            do l=1,isp(j)
               do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  z(i,j,1)=0._r8
               enddo
            enddo
         enddo
      !$omp end parallel do
      !
      !$omp parallel do private(k, l, i)
         do j=1,jj
            do k=1,kk
            do l=1,isp(j)
               do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  z(i,j,k+1)=min(depths(i,j)*1.e2_r8,z(i,j,k)+dz(i,j,k)*1.e2_r8)
               enddo
            enddo
            enddo
         enddo
      !$omp end parallel do
      !
      !$omp parallel do private(k, l, i)
         do j=1,jj
            do k=2,kk
               do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                     if ((z(i,j,kk+1)-z(i,j,k)).lt.1.e-4_r8) then
                        z(i,j,k)=depths(i,j)*1.e2_r8
                     endif
                  enddo
               enddo
            enddo
            do l=1,isp(j)
               do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  z(i,j,kk+1)=depths(i,j)*1.e2_r8
               enddo
            enddo
         enddo
      !$omp end parallel do
      !
      ! Compute layer interface geopotential.
      !$omp parallel do private(k, l, i)
         do j = 1, jj
            do k = 1, kk + 1
               do l = 1, isp(j)
                  do i = max(1, ifp(j, l)),min(ii, ilp(j, l))
                     phi(i, j, k) = - g*z(i, j, k)
                  enddo
               enddo
            enddo
         enddo
      !$omp end parallel do
      
      end subroutine ictsz_channel
     
      subroutine inifrc_channel
      ! ---------------------------------------------------
      ! Define forcing - This could be used as a getfrc
      ! function as well if time varying forcing is desired
      ! (time dependency would need to be added).
      ! ---------------------------------------------------
      
         intrinsic tanh
         
         real(r8) :: ztx0,mty0,sst0,sss0
         integer i,j,l,k,ios
         logical :: fexist
         
         namelist /idlfor/ ztx0,mty0,sst0,sss0
         
         ! --- READ NAMELIST FILE
         inquire(file='limits',exist=fexist)
         if (fexist) then
            open (unit=nfu,file='limits',status='old',action='read')
         else
            write (lp,*) 'inifrc_channel: could not find namelist file!'
            call xchalt('(inifrc_channel)')
            stop '(inifrc_channel)'
         endif
         !
         ! --- READ AND BROADCAST
         read (unit=nfu,nml=idlfor,iostat=ios)
         close (unit=nfu)
         call xcbcst(ztx0)
         call xcbcst(mty0)
         call xcbcst(sst0)
         call xcbcst(sss0)
         !
         ! Most variables will be set to 0, but it is useful to keep them here
         ! to facilitate future studies with more complex forcing.
         ! try removing the omp loop, this is done once anyway
      !$omp parallel do private(l, i, k)
            do j = 1, jj
              do l = 1, isp(j)
              do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                 ustarw(i,j) = 0.005_r8 ! friction velocity for open water
                 slp(i,j)    = 1000._r8 ! sea level pressure
                 atmco2(i,j) = 400._r8  ! atmospheric co2 concentration
                 swa(i,j)    = 0._r8    ! shortwave
                 nsf(i,j)    = 0._r8    ! non-solar
                 hmltfz(i,j) = 0._r8    ! heat flux due to melting/freezing
                 !dfl(i,j)    = 0._r8    ! derivate of non-solar in respect to T
                 !alb(i,j)    = 0._r8    ! albedo
                 eva(i,j)    = 0._r8    ! evaporation
                 lip(i,j)    = 0._r8    ! liquid precip
                 sop(i,j)    = 0._r8    ! solid precip 
                 rnf(i,j)    = 0._r8    ! runoff
                 rfi(i,j)    = 0._r8    ! runoff ice
                 fmltfz(i,j) = 0._r8    ! fresh water flux due to melting/freezing
                 sfl(i,j)    = 0._r8    ! salt flux
                 abswnd(i,j) = 0._r8    ! wind speed at measurement height -zu-
                 !albw(i,j)   = 0._r8    ! daily mean open water albedo
                 flxco2(i,j) = 0._r8    ! air-sea co2 flux
                 flxdms(i,j) = 0._r8    ! sea-air dms flux
                 !
                 ! -------------------------------------------------
                 ! SST and SSS climatologies are set here (to values
                 ! defined in the namelist), but their usage is controlled 
                 ! by namelist timescales that are 0. by default. 
                 ! -------------------------------------------------
                 !
                 do k=1,12
                    sstclm(i,j,k) = sst0
                    sssclm(i,j,k) = sss0
                 enddo
              enddo
              enddo
              do l = 1, isu(j)
              do i = max(1,ifu(j,l)), min(ii,ilu(j,l))
                 taux(i,j) = ztx0*10._r8
              enddo
              enddo
              do l = 1, isv(j)
              do i = max(1,ifv(j,l)), min(ii,ilv(j,l))
                 tauy(i,j) = mty0*10._r8
              enddo
              enddo
            enddo
      !$omp end parallel do
      
      
      end subroutine inifrc_channel
      
      subroutine icaux_channel
      ! This is mostly an empty subroutine, 
      ! but could be used to setup sea ice
      
         ntda=0
      
      end subroutine icaux_channel
      
      subroutine sfcstr_channel(m,n,mm,nn,k1m,k1n)
      !
      ! -----------------------------------------------------------------
      ! Compute the surface stress. This subroutine is kept here as a place
      ! holder in case of sea ice coupling. For the moment taux and tauy
      ! are set in the inifrc_channel.
      ! ------------------------------------------------------------------

      integer m,n,mm,nn,k1m,k1n

      integer i,j,l

      !$omp parallel do private(l,i)
      do j=1,jj
         do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
               taux(i,j)=10._r8*ztx(i,j)
            enddo
         enddo
         do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
               tauy(i,j)=10._r8*mty(i,j)
            enddo
         enddo
      enddo
      !$omp end parallel do
      !
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'sfcstr:'
        endif
        call chksummsk(taux,iu,1,'taux')
        call chksummsk(tauy,iv,1,'tauy')
      endif
      !
      end subroutine sfcstr_channel
      
end module mod_channel
