      subroutine ice2mod(ustari,tauxice,tauyice,uicem,vicem,hicem,ficem,
     .                   itp)
c
      use mod_xc, only: nbdy,idm,jdm,itdm,jtdm,mnproc,xcaput,xctilr,
     .                  halo_ps,ii,jj,xcaget
c
      implicit none
c
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ustari,tauxice,tauyice,uicem,vicem,hicem,ficem
      integer, dimension(itdm,jtdm) :: itp
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
c --- if ice volume is less than -vol_min- on both sides of a velocity
c --- point, set the ice velocity to zero
      real vol_min
      parameter (vol_min=.1)
c
      real, dimension(L,0:M) :: ustrq,tauxq,tauyq
      real, dimension(itdm,jtdm) :: tmp2da,tmp2db
      integer, dimension(itdm,jtdm) :: itu,itv,iuvt
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util1
      real du2,q
      integer i,j,i1,j1,jm1
c
c --- ------------------------------------------------------------------
c --- compute friction velocity under ice  ustari  [m/s] and ice-water
c --- stress  tauxice,tauyice  [N/m^2] felt by the ocean
c ---  ustari is defined on the  p  point in a c-grid
c ---  tauxice,tauyice  are defined on the  u,v  points in a c-grid
c ---  cdwat     , is the ice-ocean drag coefficient       [ - ]
c ---  rhowat    , is the mean surface density of seawater [kg/m^3]
c ---  u,v       , is the ice velocity on b-grid           [m/s]
c ---  uwat,vwat , is the ml velocity components on b-grid [m/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
c$OMP PARALLEL DO PRIVATE(du2,q)
        do j=0,M
          do i=1,L
            du2=(u(i,j,lold)-uwat(i,j))**2+(v(i,j,lold)-vwat(i,j))**2
            ustrq(i,j)=sqrt(cdwat*du2)*vm(i,j)
            q=cdwat*rhowat*sqrt(du2)
            tauxq(i,j)=q*(u(i,j,lold)-uwat(i,j))*vm(i,j)
            tauyq(i,j)=q*(v(i,j,lold)-vwat(i,j))*vm(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do j=1,jtdm
          jm1=mod(j-2+jtdm,jtdm)+1
          do i=1,itdm
            if (itp(i,j).eq.1.and.itp(max(1,i-1),j).eq.1) then
              itu(i,j)=1
            else
              itu(i,j)=0
            endif
            if (itp(i,j).eq.1.and.itp(i,jm1).eq.1) then
              itv(i,j)=1
            else
              itv(i,j)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call intp_iq2mp(itdm,jtdm,itp,ustrq,tmp2da)
      endif
      call xcaput(tmp2da,ustari,1)
      if (mnproc.eq.1) then
        call intp_iq2mu(itdm,jtdm,itu,tauxq,tmp2da)
      endif
      call xcaput(tmp2da,tauxice,1)
      if (mnproc.eq.1) then
        call intp_iq2mv(itdm,jtdm,itv,tauyq,tmp2da)
      endif
      call xcaput(tmp2da,tauyice,1)
c
c --- ------------------------------------------------------------------
c --- get advective ice velocities [m/s]
c --- ------------------------------------------------------------------
c
      call xctilr(hicem, 1,1, 1,1, halo_ps)
      call xctilr(ficem, 1,1, 1,1, halo_ps)
      do j=1,jj
        do i=1,ii
          util1(i,j)=hicem(i,j)*ficem(i,j)
        enddo
      enddo
      call xcaget(tmp2db,util1,1)
      if (mnproc.eq.1) then
c$OMP PARALLEL DO
        do j=1,jtdm
          do i=1,itdm
            iuvt(i,j)=0
            if (itu(i,j).eq.1) then
              if (tmp2db(i-1,j).gt.vol_min.or.
     .            tmp2db(i  ,j).gt.vol_min) then
                iuvt(i,j)=1
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call intp_iq2mu(itdm,jtdm,iuvt,u(1,0,lold),tmp2da)
      endif
      call xcaput(tmp2da,uicem,1)
      if (mnproc.eq.1) then
c$OMP PARALLEL DO PRIVATE(jm1)
        do j=1,jtdm
          jm1=mod(j-2+jtdm,jtdm)+1
          do i=1,itdm
            iuvt(i,j)=0
            if (itv(i,j).eq.1) then
              if (tmp2db(i,jm1).gt.vol_min.or.
     .            tmp2db(i,j  ).gt.vol_min) then
                iuvt(i,j)=1
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call intp_iq2mv(itdm,jtdm,iuvt,v(1,0,lold),tmp2da)
      endif
      call xcaput(tmp2da,vicem,1)
c
      return
      end
