      subroutine iniicedyn(scpx,scpy,corio)
c
      use mod_xc, only: nbdy,idm,jdm,itdm,jtdm,ip,mnproc,ii,jj,xcmaxr,
     .                  xchalt,xcaget 
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  scpx,scpy,corio
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      real, dimension(itdm,jtdm) :: tmp2d
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util1
      integer, dimension(itdm,jtdm) :: itp
      integer, dimension(0:l,0:m) :: ones
      integer i,j,i1,j1,jm1
      real pi,rad,swin,swat
c
c----------------------------------------------------------------------
c dummy mask of ones
c----------------------------------------------------------------------
      if (mnproc.eq.1) then
c$OMP PARALLEL DO
        do j=0,m
          do i=0,l
            ones(i,j) = 1
          enddo
        enddo
c$OMP END PARALLEL DO
c
c----------------------------------------------------------------------
c check grid dimension consistency between ice model and ocean model
c----------------------------------------------------------------------
        if (itdm.ne.LP.or.jtdm.ne.MP) then
          write (*,*)
     .      'Grid dimension inconsistency between ice and ocean model!'
          write (*,*)
     .      'Check dimensions.h and param_icedyn.h!'
          call xchalt('(iniicedyn)')
        endif
c
c----------------------------------------------------------------------
c spherical parameters and user supplied parameters (use mks units)
c----------------------------------------------------------------------
        pi =4.*atan(1.)
        rad=pi/180.
      endif

      dx = 1.
      dy = 1.
c$OMP PARALLEL DO REDUCTION(max:dx,dy)
      do 101 j=1,jj
      do 101 i=1,ii
      dx=max(dx,scpx(i,j))
  101 dy=max(dy,scpy(i,j))
c$OMP END PARALLEL DO
      call xcmaxr(dx)
      call xcmaxr(dy)

c----------------------------------------------------------------------
c set coordinate mapping factors
c cartesian coordinates: pm=1, pn=1 and dx=dx at mean latitude of
c the domain.
c spherical coordinates: pm=1/cosphi, pn=1 and dx=dx at phi=0 deg.
c where phi is the latitude of the corresponding scalar grid point
c----------------------------------------------------------------------
      call xcaget(tmp2d,scpx,1)
      if (mnproc.eq.1) then
c$OMP PARALLEL DO
        do j=1,jtdm
          do i=1,itdm
            tmp2d(i,j) = dx/max(1.,tmp2d(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
        call intp_mp2ip(itdm,jtdm,ones,tmp2d,pm)
      endif
      call xcaget(tmp2d,scpy,1)
      if (mnproc.eq.1) then
c$OMP PARALLEL DO
        do j=1,jtdm
          do i=1,itdm
            tmp2d(i,j) = dy/max(1.,tmp2d(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
        call intp_mp2ip(itdm,jtdm,ones,tmp2d,pn)
      endif

      if (mnproc.eq.1) then
        dx=dx*1.e-2
        dy=dy*1.e-2
c----------------------------------------------------------------------
c Derived constants
c----------------------------------------------------------------------
        dxsq=dx*dx
        dysq=dy*dy
        sx2 = 0.5/dxsq
        sy2 = 0.5/dysq
        sxy = 0.25/(dx*dy)

c-----------------------------------------------------------------------
c Grid parameters and constants for subroutines RELCON and RELAX
c-----------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(jm1)
        do 103 j=0,M
        jm1=mod(j-1+MP,MP)
        do 103 i=1,L
        pmpnv(i,j) = 0.25*(pm(i,j)*pn(i,j)
     .                    +pm(i,jm1)*pn(i,jm1)
     .                    +pm(i-1,j)*pn(i-1,j)
     .                    +pm(i-1,jm1)*pn(i-1,jm1))

        sx2p(i,j)=sx2*(pm(i  ,j)+pm(i  ,jm1))/(pn(i  ,j)+pn(i  ,jm1))
        sx2m(i,j)=sx2*(pm(i-1,j)+pm(i-1,jm1))/(pn(i-1,j)+pn(i-1,jm1))

        sy2p(i,j)=sy2*(pn(i,j)  +pn(i-1,j)  )/(pm(i,j)  +pm(i-1,j)  )
  103   sy2m(i,j)=sy2*(pn(i,jm1)+pn(i-1,jm1))/(pm(i,jm1)+pm(i-1,jm1))
c$OMP END PARALLEL DO
      endif

c-----------------------------------------------------------------------
c Coriolis parameters at vector grid points
c   f is a function of the geographical latitude
c-----------------------------------------------------------------------
      call xcaget(tmp2d,corio,1)
      if (mnproc.eq.1) then
        call intp_mq2iq(itdm,jtdm,tmp2d,f)
      endif

c-----------------------------------------------------------------------
c Get mask for the full ocean domain
c-----------------------------------------------------------------------
      do j=1,jj
        do i=1,ii
          util1(i,j)=ip(i,j)
        enddo
      enddo
      call xcaget(tmp2d,util1,1)
      if (mnproc.eq.1) then
        do j=1,jtdm
          do i=1,itdm
            itp(i,j)=tmp2d(i,j)
          enddo
        enddo
      endif

      if (mnproc.eq.1) then

c----------------------------------------------------------------------
c Gravity constant [m/s**2]
c----------------------------------------------------------------------
        g=9.806
c----------------------------------------------------------------------
c Time step [s], 1 day = 86400 sec
c----------------------------------------------------------------------
        dt=86400.
c-----------------------------------------------------------------------
c Drag coefficients
c-----------------------------------------------------------------------
        cdwat=5.5E-3      ! do not forget E-3   ! SIOM standard
c-----------------------------------------------------------------------
c Turning angle of wind
c For surface winds (ECMWF), select an angle of 0 degree.
c For geostrophic winds, an angle of about 20 .. 25 degree is OK
c-----------------------------------------------------------------------
        angwin=.0
        swin  =sin(angwin*rad)
        coswin=cos(angwin*rad)
c-----------------------------------------------------------------------
c Turning angle of water
c-----------------------------------------------------------------------
        angwat=.0 !25.0
        swat  =sin(angwat*rad)
        coswat=cos(angwat*rad)
c-----------------------------------------------------------------------
c Southern hemisphere: Ekman turning in the opposite direction
c-----------------------------------------------------------------------
c$OMP PARALLEL DO
        do j=0,M
           do i=1,L
              sinwin(i,j)=sign(1.,f(i,j))*swin
              sinwat(i,j)=sign(1.,f(i,j))*swat
           enddo
        enddo
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
c Densities            in [kg / m**3]
c-----------------------------------------------------------------------
        rhowat=1024.         ! sea water
        rhoice= 900.         ! sea ice
c----------------------------------------------------------------------
c Rheology parameters
c ===================
c Ice pressure parameter P*
c   P* has the unit [N/m**2]
c       P* must be set according to the time scale of the forcing
c   High P* for strong ice, hard to deform
c   Low  P* for weak   ice, easy to deform

c For the Weddell Sea, reasonable values are
c   P* = 20 000   for daily forcing
c   P* = 10 000   for monthly forcing (climatology)

c Ice pressure parameter C*
c   The higher C*, the higher is the minimum ice concentration
c   for significant internal forces.
c   Hibler (1979) suggests C* = 20
c   For C* = 20, there is deformation for A >= 85%

c Eccentricity e of the elliptic yield curve for internal forces
c   Hibler (1979) suggests e = 2

c gmin: minimum value of DELTA, used in subroutine PLAST

c               Markus Harder, April 1994
c----------------------------------------------------------------------
c   pstar=0.   ! Arctic Test
c..   pstar=25000.   ! Arctic


c$OMP PARALLEL DO
        do j=0,M
          do i=0,L
#if defined(SYN) || defined(BCM)
            pstar(i,j)=.5*((1.+sign(1.,f(i,j)))*20000.
     .                    +(1.-sign(1.,f(i,j)))*20000.)
#else
            pstar(i,j)=.5*((1.+sign(1.,f(i,j)))*7500.
     .                    +(1.-sign(1.,f(i,j)))*7500.)
#endif
          enddo
        enddo
c$OMP END PARALLEL DO

        cstar=20.
c-----------------------------------------------------------------------
c Eccentricity of yield curve
c-----------------------------------------------------------------------
        eccen=2.0
        ecm2=1.0/eccen**2
c-----------------------------------------------------------------------
c Regime parameter
c-----------------------------------------------------------------------
        gmin=5.0E-9      ! Neuer Standard
c-----------------------------------------------------------------------
c parameters for the relax subroutine
c-----------------------------------------------------------------------
        mmax=   10000      ! maximum number of iterations
        mwt=   100      ! iterations with overrelaxation
c       mmax=   1      ! maximum number of iterations
c       mwt=   1      ! iterations with overrelaxation
        vrmax=   1.E-04      ! error margin
        wt=   1.45      ! overrelaxation factor
        vrwt=   2.      ! threshold to stop overrelaxation
        amasmin=rhoice*1.E-6   ! minimum ice mass for relax
c-----------------------------------------------------------------------
c Define HM and OM mask for scalar grid points
c-----------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(j1,i1)
        do j=1,jtdm
          j1=j-1
          do i=1,itdm
            i1=i-1
            if (itp(i,j).eq.1) then
              hm(i1,j1)=1
              om(i1,j1)=1
            else
              hm(i1,j1)=0
              om(i1,j1)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
c Define VM mask for vector grid points
c-----------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(jm1)
        do j=0,M
          jm1=mod(j-1+MP,MP)
          do i=1,L
            if (hm(i-1,jm1)+hm(i-1,j)+hm(i,jm1)+hm(i,j).eq.4) then
              vm(i,j)=1
            else
              vm(i,j)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
C----------------------------------------------------------------------
c Initial conditions for vector points
C----------------------------------------------------------------------
c$OMP PARALLEL DO
         do 105 j=0,M
          do 105 i=1,L
c----------------------------------------------------------------------
c Dummy values for RELAX variables
c----------------------------------------------------------------------
            ru (i,j) = 0.
            rv (i,j) = 0.
            bu (i,j) = 0.
            bv (i,j) = 0.
            fx (i,j) = 0.
            fy (i,j) = 0.
            asy (i,j)= 0.
            den (i,j)= 1.
            amas (i,j)=0.

  105    continue
c$OMP END PARALLEL DO
C----------------------------------------------------------------------
c Initial conditions for scalar points
C----------------------------------------------------------------------
c$OMP PARALLEL DO
        do 106 j=0,M
          do 106 i=0,L
c-----------------------------------------------------------------------
c Dummy values for rheology parameters
c----------------------------------------------------------------------
             zeta(i,j)=0.
             eta(i,j) =0.
  106    continue
c$OMP END PARALLEL DO
 
       endif

       return
       end
