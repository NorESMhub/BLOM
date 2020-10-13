      subroutine ictsz_channel
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from file to define layer temperature and
c --- salinity and geopotential at layer interfaces.
c --- ------------------------------------------------------------------
c
      use mod_xc
c      use netcdf
c
      implicit none
c
#include "common_blocks.h"
c
      intrinsic random_seed, random_number, atan, sin, tanh, cos
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dz
      real, dimension(kdm) :: sigmr0
      real, dimension(kdm) :: dz0
      real :: S0,sig0,sig0dz,sigdz,sigscl,dztop,dzmax,dzscl
c      real, dimension(kdm) :: temp0
c      
c      real :: r
c      real, parameter :: PI0=4.*atan(1.0)
c      real, parameter :: S0=35.
      integer i,j,k,l,ios
      logical :: fexist
      
      real sofsig,tofsig,sig,rho
      external sofsig,tofsig,sig,rho
c
      namelist /idlini/ S0,sig0,sig0dz,sigdz,sigscl,dztop,dzmax,dzscl
c
      if (mnproc.eq.1) then
        write (lp,'(1a)') ' idealized initial conditions'
        call flush(lp)
      endif
c
c --- DEFINE IDEALIZED PROFILES OF sigmar and DZ
c --- salinity is constant (S0) and temperature is then inverted below
c
c --- read parameters from the namelist     
      inquire(file='limits',exist=fexist)
      if (fexist) then
          open (unit=nfu,file='limits',status='old',action='read')
      else
          write (lp,*) 'ictsz_channel: could not find namelist file!'
          call xchalt('(ictsz_channel)')
          stop '(ictsz_channel)'
      endif
c --- 
      read (unit=nfu,nml=idlini,iostat=ios)
      close (unit=nfu)
      call xcbcst(ios)
c     
c ---     
      sigmr0(1)=sig0 !33.2133 !34.59 cases _01-_03
      sigmr0(2)=sig0 !33.2133 !34.59
      dz0(1)=dztop !10.
      dz0(2)=dztop !10.
      do k=3,kk
         sigmr0(k)=sigmr0(k-1) + sig0dz + sigdz*
     .   (1.-tanh(sigscl*pi*(k-1)/kk))
         dz0(k)=dzmax*tanh(dzscl*pi*(k-1)/kk)
c ---    define a sigma profie and invert temperature
c ---    cases _01-_03 
c         sigmr0(k)=sigmr0(k-1) + 0.05*k*(1.+cos(pi*(k-1.)/kk))/kk
c        cases _04
c         sigmr0(k)=sigmr0(k-1) + 5.e-6+0.5*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    cases _05-07
c         sigmr0(k)=sigmr0(k-1) + 5.e-4+0.5*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    cases _08
c         sigmr0(k)=sigmr0(k-1) + 5.e-4+0.4*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    cases _09
c         sigmr0(k)=sigmr0(k-1) + 5.e-4+0.3*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    cases _10
c         sigmr0(k)=sigmr0(k-1) + 5.e-4+0.2*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    cases _11
c         sigmr0(k)=sigmr0(k-1) + 5.e-4+0.1*(1.-tanh(1.5*pi*(k-1)/kk))
c ---    define a temperature profile and invert sigma
c         temp0(k)=temp0(k-1) + 0.2*(k/kk)*(1.+cos(PI_8*(k-1)/kk))
c         sigmr0(k)=sig(temp0(k),S0)
c ---    cases _01-_03
c         dz0(k)=10.+0.3*k*exp(k/(0.5*kk))
c         dz0(k)=50.*tanh(1.5*pi*(k-1)/kk)
      enddo
c
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c                call random_number(r) #add r to the list of privates if used
                dz(i,j,k)=dz0(k) !min(k*10.,100.)
                saln(i,j,k)=S0 !-(j0+j)/(0.5*jtdm*exp(k/2.))
                sigmar(i,j,k)=sigmr0(k)*1.e-3 !convert units to g/cm^3
                temp(i,j,k)=tofsig(sigmar(i,j,k),saln(i,j,k)) 
c                temp(i,j,k)=10.+(j0+j)/(0.125*jtdm*exp(k/2.))
c                temp(i,j,k)=10.+0.1*r+(j0+j)/(0.125*jtdm*exp(k/10.))
c                temp(i,j,k)=10.+(10.+j0+j)/(0.0625*jtdm*exp(k/10.))
C                temp(i,j,k)=10.+(10.+j0+j)/
C     .          ((j0+j)/(jtdm/20.)+0.0625*jtdm*exp(k/10.))
c                temp(i,j,k)=10.+(5.+10.*sin((j0+j)*PI_8/jtdm))/
c     .          (exp((k-1)/10.))
c --- cases LargeBasin5-7
c                temp(i,j,k)=10.+(2.5+5.*sin((j0+j)*PI_8/jtdm))/
c     .          (exp((k-1)/10.))
c --- cases LargeBasin8-9
c                 temp(i,j,k)=10.+(2.5+2.5*sin((j0+j)*PI_8/jtdm))/
c     .           (exp((k-1)/10.))
c --- cases LargeBasin10-11
c                 temp(i,j,k)=10.+(1.+2.5*sin((j0+j)*PI_8/jtdm))/
c     .                 (exp((k-1)/10.))
c --- case LargeBasin12 and LargeBasin13
c                temp(i,j,k)=10.+(2.5/exp((k-1)/10.))
c --- case LargeBasin14
c                temp(i,j,k)=10.+(5./exp((k-1)/10.))
c --- case LargeBasin15
c                temp(i,j,k)=10.+(5./exp((k-1)/5.))
c --- case LargeBasin16
c                temp(i,j,k)=10.+(3.0/exp((k-1)/10.))
c --- case LargeBasin17
c                temp(i,j,k)=10.+(3.0/exp((k-1)/5.))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO

  
c      sigmar(1-nbdy,1-nbdy,k)
c
c
c --- Construct interface depths [cm] from layer thicknesses [m] and
c --- convert unit of reference potential density from [kg/m^3] to
c --- [g/cm^3]
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,1)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      do k=1,kk
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,k+1)=min(depths(i,j)*1.e2,z(i,j,k)+dz(i,j,k)*1.e2)
c            sigmar(i,j,k)=sig(temp(i,j,k),saln(i,j,k)) !sigmar should not vary spatially!!
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=2,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (z(i,j,kk+1)-z(i,j,k).lt.1.e-4)
     .        z(i,j,k)=depths(i,j)*1.e2
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,kk+1)=depths(i,j)*1.e2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute layer interface geopotential
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk+1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            phi(i,j,k)=-g*z(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
