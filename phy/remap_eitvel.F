! ------------------------------------------------------------------------------
! Copyright (C) 2006-2020 Mats Bentsen, Mehmet Ilicak
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

      subroutine remap_eitvel(scuy,scvx,scp2i,scp2,pbmin,pbu,pbv,plo,
     .                        u,v,dt,mrg,dp,temp,saln,uflx,vflx,
     .                        utflx,vtflx,usflx,vsflx
#ifdef TRC
     .                       ,k,trc
#endif
     .                       )
c
c --- ------------------------------------------------------------------
c --- Advection of layer pressure thickness and tracers by incremental
c --- remapping.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
c --- Argument variables:
c ---   scuy   - length of cell boundary with u-point as midpoint.
c ---   scvx   - length of cell boundary with v-point as midpoint.
c ---   scp2i  - inverse of grid cell area.
c ---   scp2   - grid cell area.
c ---   pbmin  - minimum bottom pressure of a grid cell and its
c ---            neighbors.
c ---   pbu    - bottom pressure at u-point, defined as
c ---            min(pb(i-1,j),pb(i,j)).
c ---   pbv    - bottom pressure at v-point, defined as
c ---            min(pb(i,j-1),pb(i,j)).
c ---   plo    - lower interface pressure of layer pressure thickness.
c ---   u      - u-component of velocity.
c ---   v      - v-component of velocity.
c ---   dt     - time step.
c ---   dp     - layer pressure thickness.
c ---   temp   - temperature.
c ---   saln   - salinity.
c ---   uflx   - u-component of mass flux.
c ---   vflx   - v-component of mass flux.
c ---   utflx  - u-component of heat flux.
c ---   vtflx  - v-component of heat flux.
c ---   usflx  - u-component of salt flux.
c ---   vsflx  - v-component of salt flux.
c ---   mrg    - margin of halo that must be valid upon return.
c ---   k      - layer index
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  scuy,scvx,scp2i,scp2,pbmin,pbu,pbv,plo,u,v,dp,temp,saln,
     .  uflx,vflx,utflx,vtflx,usflx,vsflx
      real dt
      integer mrg
#ifdef TRC
#  include "param_trc.h"
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm,ntr) :: trc
      integer k
#endif
c
c --- Local variables.
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  pup,dx,dy,xd,yd,tx,ty,td,sx,sy,sd,cu,cv,cuc,cvc,fdu,fdv,ftu,ftv,
     .  fsu,fsv
      real dxi,dyi,dpw,dpe,dps,dpn,dpsw,dpse,dpc,dpnw,dpne,
     .     dgmx,dfmx,dfmn,q,q1,q2,q3,q4,tgmx,tgmn,tfmx,tfmn,
     .     sgmx,sgmn,sfmx,sfmn,xm,ym,xc0,xc1,yc0,yc1,x2,y2,x4,y4,
     .     a,ax,ay,axx,ayy,axy,dl,fd,qx,qy
      integer i,j,l,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne,nw
#ifdef TRC
#  ifdef ATRC
      real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  trx,try,trd,ftru,ftrv
      real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ag,agx,agy,agd,fagu,fagv
      real xdt,ydt,axxx,ayyy,axxy,axyy,qxx,qyy,qxy,fdt
      integer nt,nat
#  else
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  trx,try,trd,ftru,ftrv
      integer nt
#  endif
#endif
c
c --- Parameters:
c ---   dpeps - small layer pressure thickness (equivalent to approx.
c ---           10^(-16) m).
c ---   treps - small tracer concentration
c
      real dpeps
      parameter (dpeps=1.e-11)
#if defined(TRC) && defined(ATRC)
c
      real treps
      parameter (treps=1.e-14)
#endif
c
c --- ------------------------------------------------------------------
c --- General information:
c ---   Logical arrangment of variables is as follows: Layer pressure
c ---   thickness dp(i,j), as an example of a scalar variable, is the
c ---   mean layer pressure thickness of grid cell (i,j). Velocity
c ---   component u(i,j) is located at the midpoint of the cell boundary
c ---   separating grid cells (i-1,j) and (i,j). Velocity component
c ---   v(i,j) is located at the midpoint of the cell boundary
c ---   separating grid cells (i,j-1) and (i,j). A corner variable with
c ---   index (i,j) is located at the common grid cell corner of grid
c ---   cells (i-1,j-1), (i,j-1), (i-1,j), and (i,j).
c ---
c ---   The divergence of the velocity field is defined as follows:
c ---     (u(i+1,j)*scuy(i+1,j)-u(i,j)*scuy(i,j)
c ---     +v(i,j+1)*scvy(i,j+1)-v(i,j)*scvx(i,j))*scp2i(i,j)
c ---   By construction, the "fluxing areas" used in obtaining fluxes
c ---   trough cell boundaries containing u(i,j) and v(i,j), are equal
c ---   to u(i,j)*scuy(i,j)*dt and v(i,j)*scvx(i,j)*dt, respectively. 
c --- ------------------------------------------------------------------
c
c --- ------------------------------------------------------------------
c --- Add small number to density field and initialize some variables.
c --- ------------------------------------------------------------------
c
      do j=1-mrg-2,jj+mrg+2
        do l=1,isp(j)
        do i=max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
          dp(i,j)=max(0.,dp(i,j))+dpeps
          pup(i,j)=plo(i,j)-dp(i,j)
        enddo
        enddo
        do i=1-mrg-1,ii+mrg+1
          fdu(i,j)=0.
          fdv(i,j)=0.
          ftu(i,j)=0.
          ftv(i,j)=0.
          fsu(i,j)=0.
          fsv(i,j)=0.
#ifdef TRC
#  ifdef ATRC
          do nt=1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            ftru(nt,i,j)=0.
            ftrv(nt,i,j)=0.
          enddo
          do nt=1,natr
            fagu(nt,i,j)=0.
            fagv(nt,i,j)=0.
          enddo
#  else
          do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            ftru(nt,i,j)=0.
            ftrv(nt,i,j)=0.
          enddo
#  endif
#endif
          cu(i,j)=0.
          cv(i,j)=0.
        enddo
      enddo
#if defined(TRC) && defined(ATRC)
c
      do nt=1,natr
        nat=ntr-natr+nt
        do j=1-mrg-2,jj+mrg+2
          do l=1,isp(j)
          do i=max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
            trc(i,j,k,nt)=max(0.,trc(i,j,k,nt))+treps
            ag(nt,i,j)=trc(i,j,k,nat)/trc(i,j,k,nt)
          enddo
          enddo
        enddo
      enddo
#endif
c
c --- ------------------------------------------------------------------
c --- Compute limited gradients, center of mass coordinates, and
c --- non-dimensional velocities.
c --- ------------------------------------------------------------------
c
      do j=1-mrg-1,jj+mrg+1
c
        do l=1,isp(j)
        do i=max(1-mrg-1,ifp(j,l)),min(ii+mrg+1,ilp(j,l))
c
c --- --- Define indices for grid cell neighbors, ensuring that only wet
c --- --- points are used.
          iw=i-iu(i  ,j)
          ie=i+iu(i+1,j)
          js=j-iv(i,j  )
          jn=j+iv(i,j+1)
          isw=i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw=j*(1-ip(iw,js))+js*ip(iw,js)
          ise=i*(1-ip(ie,js))+ie*ip(ie,js)
          jse=j*(1-ip(ie,js))+js*ip(ie,js)
          inw=i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw=j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine=i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne=j*(1-ip(ie,jn))+jn*ip(ie,jn)
c
          dxi=1./max(1,ie-iw)
          dyi=1./max(1,jn-js)
c
c --- --- Compute limited gradient for layer pressure thickness and
c --- --- center of mass coordinate.
          dpsw=max(dpeps,min(pbmin(i,j)-pup(isw,jsw),dp(isw,jsw)))
          dps =max(dpeps,min(pbmin(i,j)-pup(i  ,js ),dp(i  ,js )))
          dpse=max(dpeps,min(pbmin(i,j)-pup(ise,jse),dp(ise,jse)))
          dpw =max(dpeps,min(pbmin(i,j)-pup(iw ,j  ),dp(iw ,j  )))
          dpc =max(dpeps,min(pbmin(i,j)-pup(i  ,j  ),dp(i  ,j  )))
          dpe =max(dpeps,min(pbmin(i,j)-pup(ie ,j  ),dp(ie ,j  )))
          dpnw=max(dpeps,min(pbmin(i,j)-pup(inw,jnw),dp(inw,jnw)))
          dpn =max(dpeps,min(pbmin(i,j)-pup(i  ,jn ),dp(i  ,jn )))
          dpne=max(dpeps,min(pbmin(i,j)-pup(ine,jne),dp(ine,jne)))
          dx(i,j)=(dpe-dpw)*dxi
          dy(i,j)=(dpn-dps)*dyi
          dgmx=.5*(abs(dx(i,j))+abs(dy(i,j)))
          dfmx=max(0.,max(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          dfmn=min(0.,min(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          if (dfmx.gt.0..and.dfmn.lt.0.) then
            q=min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
            dx(i,j)=dx(i,j)*q
            dy(i,j)=dy(i,j)*q
            xd(i,j)=dx(i,j)/(12.*dp(i,j))
            yd(i,j)=dy(i,j)/(12.*dp(i,j))
          else
            dx(i,j)=0.
            dy(i,j)=0.
            xd(i,j)=0.
            yd(i,j)=0.
          endif
c
c --- --- Compute limited gradients for temperature, salinity, and
c --- --- density
          tx(i,j)=(temp(ie,j)-temp(iw,j))*dxi
          ty(i,j)=(temp(i,jn)-temp(i,js))*dyi
          q1=tx(i,j)*(-.5-xd(i,j))
          q2=tx(i,j)*( .5-xd(i,j))
          q3=ty(i,j)*(-.5-yd(i,j))
          q4=ty(i,j)*( .5-yd(i,j))
          tgmx=max(q1,q2)+max(q3,q4)
          tgmn=min(q1,q2)+min(q3,q4)
          tfmx=max(0.,max(temp(isw,jsw),temp(i  ,js ),
     .                    temp(ise,jse),temp(iw ,j  ),
     .                    temp(ie ,j  ),temp(inw,jnw),
     .                    temp(i  ,jn ),temp(ine,jne))
     .               -temp(i,j))
          tfmn=min(0.,min(temp(isw,jsw),temp(i  ,js ),
     .                    temp(ise,jse),temp(iw ,j  ),
     .                    temp(ie ,j  ),temp(inw,jnw),
     .                    temp(i  ,jn ),temp(ine,jne))
     .               -temp(i,j))
          if (tfmx.gt.0..and.tfmn.lt.0.) then
            q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
            tx(i,j)=tx(i,j)*q
            ty(i,j)=ty(i,j)*q
            td(i,j)=temp(i,j)-tx(i,j)*xd(i,j)-ty(i,j)*yd(i,j)
          else
            tx(i,j)=0.
            ty(i,j)=0.
            td(i,j)=temp(i,j)
          endif
c
          sx(i,j)=(saln(ie,j)-saln(iw,j))*dxi
          sy(i,j)=(saln(i,jn)-saln(i,js))*dyi
          q1=sx(i,j)*(-.5-xd(i,j))
          q2=sx(i,j)*( .5-xd(i,j))
          q3=sy(i,j)*(-.5-yd(i,j))
          q4=sy(i,j)*( .5-yd(i,j))
          sgmx=max(q1,q2)+max(q3,q4)
          sgmn=min(q1,q2)+min(q3,q4)
          sfmx=max(0.,max(saln(isw,jsw),saln(i  ,js ),
     .                    saln(ise,jse),saln(iw ,j  ),
     .                    saln(ie ,j  ),saln(inw,jnw),
     .                    saln(i  ,jn ),saln(ine,jne))
     .               -saln(i,j))
          sfmn=min(0.,min(saln(isw,jsw),saln(i  ,js ),
     .                    saln(ise,jse),saln(iw ,j  ),
     .                    saln(ie ,j  ),saln(inw,jnw),
     .                    saln(i  ,jn ),saln(ine,jne))
     .               -saln(i,j))
          if (sfmx.gt.0..and.sfmn.lt.0.) then
            q=min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
            sx(i,j)=sx(i,j)*q
            sy(i,j)=sy(i,j)*q
            sd(i,j)=saln(i,j)-sx(i,j)*xd(i,j)-sy(i,j)*yd(i,j)
          else
            sx(i,j)=0.
            sy(i,j)=0.
            sd(i,j)=saln(i,j)
          endif
#ifdef TRC
#  ifdef ATRC
c
c --- --- Compute limited gradient for tracers.
          do nt=1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            trx(nt,i,j)=(trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
            try(nt,i,j)=(trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
            q1=trx(nt,i,j)*(-.5-xd(i,j))
            q2=trx(nt,i,j)*( .5-xd(i,j))
            q3=try(nt,i,j)*(-.5-yd(i,j))
            q4=try(nt,i,j)*( .5-yd(i,j))
            tgmx=max(q1,q2)+max(q3,q4)
            tgmn=min(q1,q2)+min(q3,q4)
            tfmx=max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt),
     .                      trc(ise,jse,k,nt),trc(iw ,j  ,k,nt),
     .                      trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt),
     .                      trc(i  ,jn ,k,nt),trc(ine,jne,k,nt))
     .                 -trc(i,j,k,nt))
            tfmn=min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt),
     .                      trc(ise,jse,k,nt),trc(iw ,j  ,k,nt),
     .                      trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt),
     .                      trc(i  ,jn ,k,nt),trc(ine,jne,k,nt))
     .                 -trc(i,j,k,nt))
            if (tfmx.gt.0..and.tfmn.lt.0.) then
              q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
              trx(nt,i,j)=trx(nt,i,j)*q
              try(nt,i,j)=try(nt,i,j)*q
              trd(nt,i,j)=trc(i,j,k,nt)
     .                   -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
            else
              trx(nt,i,j)=0.
              try(nt,i,j)=0.
              trd(nt,i,j)=trc(i,j,k,nt)
            endif
          enddo
c
c --- --- Compute limited gradient for age tracers.
          do nt=1,natr
            nat=ntr-natr+nt
            agx(nt,i,j)=(ag(nt,ie,j)-ag(nt,iw,j))*dxi
            agy(nt,i,j)=(ag(nt,i,jn)-ag(nt,i,js))*dyi
            q=1./(12.*trc(i,j,k,nt))
            xdt=(12.*xd(i,j)*trd(nt,i,j)+trx(nt,i,j))*q
            ydt=(12.*yd(i,j)*trd(nt,i,j)+try(nt,i,j))*q
            q1=agx(nt,i,j)*(-.5-xdt)
            q2=agx(nt,i,j)*( .5-xdt)
            q3=agy(nt,i,j)*(-.5-ydt)
            q4=agy(nt,i,j)*( .5-ydt)
            tgmx=max(q1,q2)+max(q3,q4)
            tgmn=min(q1,q2)+min(q3,q4)
            tfmx=max(0.,max(ag(nt,isw,jsw),ag(nt,i  ,js ),
     .                      ag(nt,ise,jse),ag(nt,iw ,j  ),
     .                      ag(nt,ie ,j  ),ag(nt,inw,jnw),
     .                      ag(nt,i  ,jn ),ag(nt,ine,jne))
     .                 -ag(nt,i,j))
            tfmn=min(0.,min(ag(nt,isw,jsw),ag(nt,i  ,js ),
     .                      ag(nt,ise,jse),ag(nt,iw ,j  ),
     .                      ag(nt,ie ,j  ),ag(nt,inw,jnw),
     .                      ag(nt,i  ,jn ),ag(nt,ine,jne))
     .                 -ag(nt,i,j))
            if (tfmx.gt.0..and.tfmn.lt.0.) then
              q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
              agx(nt,i,j)=agx(nt,i,j)*q
              agy(nt,i,j)=agy(nt,i,j)*q
              agd(nt,i,j)=ag(nt,i,j)-agx(nt,i,j)*xdt-agy(nt,i,j)*ydt
            else
              agx(nt,i,j)=0.
              agy(nt,i,j)=0.
              agd(nt,i,j)=ag(nt,i,j)
            endif
          enddo
#  else
c
c --- --- Compute limited gradient for tracers.
          do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            trx(nt,i,j)=(trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
            try(nt,i,j)=(trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
            q1=trx(nt,i,j)*(-.5-xd(i,j))
            q2=trx(nt,i,j)*( .5-xd(i,j))
            q3=try(nt,i,j)*(-.5-yd(i,j))
            q4=try(nt,i,j)*( .5-yd(i,j))
            tgmx=max(q1,q2)+max(q3,q4)
            tgmn=min(q1,q2)+min(q3,q4)
            tfmx=max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt),
     .                      trc(ise,jse,k,nt),trc(iw ,j  ,k,nt),
     .                      trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt),
     .                      trc(i  ,jn ,k,nt),trc(ine,jne,k,nt))
     .                 -trc(i,j,k,nt))
            tfmn=min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt),
     .                      trc(ise,jse,k,nt),trc(iw ,j  ,k,nt),
     .                      trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt),
     .                      trc(i  ,jn ,k,nt),trc(ine,jne,k,nt))
     .                 -trc(i,j,k,nt))
            if (tfmx.gt.0..and.tfmn.lt.0.) then
              q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
              trx(nt,i,j)=trx(nt,i,j)*q
              try(nt,i,j)=try(nt,i,j)*q
              trd(nt,i,j)=trc(i,j,k,nt)
     .                   -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
            else
              trx(nt,i,j)=0.
              try(nt,i,j)=0.
              trd(nt,i,j)=trc(i,j,k,nt)
            endif
          enddo
#  endif
#endif
c
        enddo
        enddo
c
      enddo
c
c --- Compute non-dimensional velocities.
c
      do j=1-mrg-1,jj+mrg+1
        do l=1,isu(j)
        do i=max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))
          if (u(i,j).gt.0.) then
            cu(i,j)=u(i,j)*dt*scuy(i,j)*scp2i(i-1,j)
          else
            cu(i,j)=u(i,j)*dt*scuy(i,j)*scp2i(i  ,j)
          endif
        enddo
        enddo
      enddo
c
      do j=1-mrg,jj+mrg+1
        do l=1,isv(j)
        do i=max(1-mrg-1,ifv(j,l)),min(ii+mrg+1,ilv(j,l))
          if (v(i,j).gt.0.) then
            cv(i,j)=v(i,j)*dt*scvx(i,j)*scp2i(i,j-1)
          else
            cv(i,j)=v(i,j)*dt*scvx(i,j)*scp2i(i,j  )
          endif
        enddo
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- Compute corner velocities. The velocity components are computed as
c --- the harmonic mean of the nearest C-grid velocity components with
c --- the following exeptions: The corner velocity component is set to
c --- zero if the nearest C-grid components have different sign, or one
c --- or tree of the neighboring grid cells are wet, or two neighbors
c --- are wet and are arranged diagonally. This construction of corner
c --- velocities will ensure that the entire fluxing area is located
c --- upwind of the cell boundary.
c --- ------------------------------------------------------------------
c
      do j=1-mrg,jj+mrg+1
        do i=1-mrg,ii+mrg+1
          nw=ip(i-1,j-1)+ip(i,j-1)+ip(i-1,j)+ip(i,j)
          if     (nw.eq.4) then
            if (cu(i,j-1)*cu(i,j).le.0.) then
              cuc(i,j)=0.
            else
              cuc(i,j)=2.*cu(i,j-1)*cu(i,j)/(cu(i,j-1)+cu(i,j))
            endif
            if (cv(i-1,j)*cv(i,j).le.0.) then
              cvc(i,j)=0.
            else
              cvc(i,j)=2.*cv(i-1,j)*cv(i,j)/(cv(i-1,j)+cv(i,j))
            endif
          elseif (nw.eq.2) then
            if     (ip(i-1,j-1)+ip(i  ,j-1).eq.2) then
              cuc(i,j)=cu(i,j-1)
              cvc(i,j)=0.
            elseif (ip(i-1,j  )+ip(i  ,j  ).eq.2) then
              cuc(i,j)=cu(i,j  )
              cvc(i,j)=0.
            elseif (ip(i-1,j-1)+ip(i-1,j  ).eq.2) then
              cuc(i,j)=0.
              cvc(i,j)=cv(i-1,j)
            elseif (ip(i  ,j-1)+ip(i  ,j  ).eq.2) then
              cuc(i,j)=0.
              cvc(i,j)=cv(i  ,j)
            else
              cuc(i,j)=0.
              cvc(i,j)=0.
            endif
          else
            cuc(i,j)=0.
            cvc(i,j)=0.
          endif
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- Compute cell boundary fluxes.
c --- ------------------------------------------------------------------
c
c --- - u-components of fluxes.
c
      do j=1-mrg,jj+mrg
c
        do l=1,isu(j)
        do i=max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))
c
c --- --- Assuming coordinate [0,0] at the u-point, the non-dimensional
c --- --- fluxing area is defined as the area of a polygon with vertices
c --- --- [0,1/2], [-cuc(i,j+1),-cvc(i,j+1)+1/2], [xm,ym],
c --- --- [-cuc(i,j),-cvc(i,j)-1/2], and [0,-1/2]. The vertex [xm,ym] is
c --- --- defined so that the polygon area is equal to cu(i,j).
c
          ym=-.5*(cvc(i,j)+cvc(i,j+1))
          xm=((ym+.5)*cuc(i,j)-(ym-.5)*cuc(i,j+1)-2.*cu(i,j))
     .       /(1.+cvc(i,j)-cvc(i,j+1))
c
          if (cu(i,j).gt.0.) then
c
            if (cvc(i,j).gt.0.) then
c
c --- ------- Add contributions from grid cell (i-1,j-1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [xc1+1/2,1/2], [-cuc(i,j)+1/2,-cvc(i,j)+1/2], and
c --- ------- [1/2,1/2].
c
              xc0=(xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1=xc0*scp2(i-1,j)*scp2i(i-1,j-1)
              x4=xc0+.5
              y4=-.5
              call triint(scp2(i-1,j-1),
     .          xc1+.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5,.5,.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i-1,j-1),max(0.,pbu(i,j)-pup(i-1,j-1)))
              fd=a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdu(i,j)=fdu(i,j)+fd
              qx=ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy=ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftu(i,j)=ftu(i,j)+fd*td(i-1,j-1)
     .                +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsu(i,j)=fsu(i,j)+fd*sd(i-1,j-1)
     .                +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
              qyy=ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
              qxy=axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
              do nt=1,natr
                fdt=fd*trd(nt,i-1,j-1)
     .             +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                ftru(nt,i,j)=ftru(nt,i,j)+fdt
                fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j-1)
     .                      +(qx *trd(nt,i-1,j-1)
     .                       +qxx*trx(nt,i-1,j-1)
     .                       +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1)
     .                      +(qy *trd(nt,i-1,j-1)
     .                       +qxy*trx(nt,i-1,j-1)
     .                       +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j-1)
     .                      +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j-1)
     .                      +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
              enddo
#  endif
#endif
            else
              x4=-cuc(i,j)+.5
              y4=-cvc(i,j)-.5
            endif
c
            if (cvc(i,j+1).lt.0.) then
c
c --- ------- Add contributions from grid cell (i-1,j+1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [xc1+1/2,-1/2], [1/2,-1/2], and
c --- ------- [-cuc(i,j+1)+1/2,-cvc(i,j+1)-1/2].
c
              xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1=xc0*scp2(i-1,j)*scp2i(i-1,j+1)
              x2=xc0+.5
              y2=.5
              call triint(scp2(i-1,j+1),
     .          xc1+.5,-.5,.5,-.5,-cuc(i,j+1)+.5,-cvc(i,j+1)-.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i-1,j+1),max(0.,pbu(i,j)-pup(i-1,j+1)))
              fd=a*dl+ax*dx(i-1,j+1)+ay*dy(i-1,j+1)
              fdu(i,j)=fdu(i,j)+fd
              qx=ax*dl+axx*dx(i-1,j+1)+axy*dy(i-1,j+1)
              qy=ay*dl+axy*dx(i-1,j+1)+ayy*dy(i-1,j+1)
              ftu(i,j)=ftu(i,j)+fd*td(i-1,j+1)
     .                +qx*tx(i-1,j+1)+qy*ty(i-1,j+1)
              fsu(i,j)=fsu(i,j)+fd*sd(i-1,j+1)
     .                +qx*sx(i-1,j+1)+qy*sy(i-1,j+1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i-1,j+1)+axxy*dy(i-1,j+1)
              qyy=ayy*dl+axyy*dx(i-1,j+1)+ayyy*dy(i-1,j+1)
              qxy=axy*dl+axxy*dx(i-1,j+1)+axyy*dy(i-1,j+1)
              do nt=1,natr
                fdt=fd*trd(nt,i-1,j+1)
     .             +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                ftru(nt,i,j)=ftru(nt,i,j)+fdt
                fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j+1)
     .                      +(qx *trd(nt,i-1,j+1)
     .                       +qxx*trx(nt,i-1,j+1)
     .                       +qxy*try(nt,i-1,j+1))*agx(nt,i-1,j+1)
     .                      +(qy *trd(nt,i-1,j+1)
     .                       +qxy*trx(nt,i-1,j+1)
     .                       +qyy*try(nt,i-1,j+1))*agy(nt,i-1,j+1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j+1)
     .                      +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j+1)
     .                      +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
              enddo
#  endif
#endif
            else
              x2=-cuc(i,j+1)+.5
              y2=-cvc(i,j+1)+.5
            endif
c
c --- ----- Add contributions from grid cell (i-1,j). Assuming
c --- ----- coordinate [0,0] at the cell center, the contributions are
c --- ----- flux integrals over the pentagon with vertices [1/2,1/2],
c --- ----- [x2,y2], [xm+1/2,ym], [x4,y4], and [1/2,-1/2].
c
            call penint(scp2(i-1,j),
     .        .5,.5,x2,y2,xm+.5,ym,x4,y4,.5,-.5,
     .                  a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
            dl=min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
            fd=a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
            fdu(i,j)=fdu(i,j)+fd
            qx=ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
            qy=ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
            ftu(i,j)=ftu(i,j)+fd*td(i-1,j)
     .              +qx*tx(i-1,j)+qy*ty(i-1,j)
            fsu(i,j)=fsu(i,j)+fd*sd(i-1,j)
     .              +qx*sx(i-1,j)+qy*sy(i-1,j)
#ifdef TRC
#  ifdef ATRC
            qxx=axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
            qyy=ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
            qxy=axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
            do nt=1,natr
              fdt=fd*trd(nt,i-1,j)
     .           +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
              ftru(nt,i,j)=ftru(nt,i,j)+fdt
              fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j)
     .                    +(qx *trd(nt,i-1,j)
     .                     +qxx*trx(nt,i-1,j)
     .                     +qxy*try(nt,i-1,j))*agx(nt,i-1,j)
     .                    +(qy *trd(nt,i-1,j)
     .                     +qxy*trx(nt,i-1,j)
     .                     +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
            enddo
            do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j)
     .                    +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
            enddo
#  else
            do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j)
     .                    +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
            enddo
#  endif
#endif
c
          else
c
            if (cvc(i,j).gt.0.) then
c
c --- ------- Add contributions from grid cell (i,j-1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [xc1-1/2,1/2], [-cuc(i,j)-1/2,-cvc(i,j)+1/2], and
c --- ------- [-1/2,1/2].
c
              xc0=(xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1=xc0*scp2(i,j)*scp2i(i,j-1)
              x4=xc0-.5
              y4=-.5
              call triint(scp2(i,j-1),
     .          xc1-.5,.5,-cuc(i,j)-.5,-cvc(i,j)+.5,-.5,.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i,j-1),max(0.,pbu(i,j)-pup(i,j-1)))
              fd=a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
              fdu(i,j)=fdu(i,j)+fd
              qx=ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
              qy=ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
              ftu(i,j)=ftu(i,j)+fd*td(i,j-1)
     .                +qx*tx(i,j-1)+qy*ty(i,j-1)
              fsu(i,j)=fsu(i,j)+fd*sd(i,j-1)
     .                +qx*sx(i,j-1)+qy*sy(i,j-1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
              qyy=ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
              qxy=axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
              do nt=1,natr
                fdt=fd*trd(nt,i,j-1)
     .             +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                ftru(nt,i,j)=ftru(nt,i,j)+fdt
                fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j-1)
     .                      +(qx *trd(nt,i,j-1)
     .                       +qxx*trx(nt,i,j-1)
     .                       +qxy*try(nt,i,j-1))*agx(nt,i,j-1)
     .                      +(qy *trd(nt,i,j-1)
     .                       +qxy*trx(nt,i,j-1)
     .                       +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j-1)
     .                      +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j-1)
     .                      +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
              enddo
#  endif
#endif
            else
              x4=-cuc(i,j)-.5
              y4=-cvc(i,j)-.5
            endif
c
            if (cvc(i,j+1).lt.0.) then
c
c --- ------- Add contributions from grid cell (i,j+1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [xc1-1/2,-1/2], [-1/2,-1/2], and
c --- ------- [-cuc(i,j+1)-1/2,-cvc(i,j+1)-1/2].
c
              xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1=xc0*scp2(i,j)*scp2i(i,j+1)
              x2=xc0-.5
              y2=.5
              call triint(scp2(i,j+1),
     .          xc1-.5,-.5,-.5,-.5,-cuc(i,j+1)-.5,-cvc(i,j+1)-.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i,j+1),max(0.,pbu(i,j)-pup(i,j+1)))
              fd=a*dl+ax*dx(i,j+1)+ay*dy(i,j+1)
              fdu(i,j)=fdu(i,j)+fd
              qx=ax*dl+axx*dx(i,j+1)+axy*dy(i,j+1)
              qy=ay*dl+axy*dx(i,j+1)+ayy*dy(i,j+1)
              ftu(i,j)=ftu(i,j)+fd*td(i,j+1)
     .                +qx*tx(i,j+1)+qy*ty(i,j+1)
              fsu(i,j)=fsu(i,j)+fd*sd(i,j+1)
     .                +qx*sx(i,j+1)+qy*sy(i,j+1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i,j+1)+axxy*dy(i,j+1)
              qyy=ayy*dl+axyy*dx(i,j+1)+ayyy*dy(i,j+1)
              qxy=axy*dl+axxy*dx(i,j+1)+axyy*dy(i,j+1)
              do nt=1,natr
                fdt=fd*trd(nt,i,j+1)
     .             +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                ftru(nt,i,j)=ftru(nt,i,j)+fdt
                fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j+1)
     .                      +(qx *trd(nt,i,j+1)
     .                       +qxx*trx(nt,i,j+1)
     .                       +qxy*try(nt,i,j+1))*agx(nt,i,j+1)
     .                      +(qy *trd(nt,i,j+1)
     .                       +qxy*trx(nt,i,j+1)
     .                       +qyy*try(nt,i,j+1))*agy(nt,i,j+1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j+1)
     .                      +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j+1)
     .                      +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
              enddo
#  endif
#endif
            else
              x2=-cuc(i,j+1)-.5
              y2=-cvc(i,j+1)+.5
            endif
c
c --- ----- Add contributions from grid cell (i,j). Assuming
c --- ----- coordinate [0,0] at the cell center, the contributions are
c --- ----- flux integrals over the pentagon with vertices [-1/2,1/2],
c --- ----- [x2,y2], [xm-1/2,ym], [x4,y4], and [-1/2,-1/2].
c
            call penint(scp2(i,j),
     .        -.5,.5,x2,y2,xm-.5,ym,x4,y4,-.5,-.5,
     .                  a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
            dl=min(dp(i,j),max(0.,pbu(i,j)-pup(i,j)))
            fd=a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdu(i,j)=fdu(i,j)+fd
            qx=ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy=ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftu(i,j)=ftu(i,j)+fd*td(i,j)
     .              +qx*tx(i,j)+qy*ty(i,j)
            fsu(i,j)=fsu(i,j)+fd*sd(i,j)
     .              +qx*sx(i,j)+qy*sy(i,j)
#ifdef TRC
#  ifdef ATRC
            qxx=axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
            qyy=ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
            qxy=axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
            do nt=1,natr
              fdt=fd*trd(nt,i,j)
     .           +qx*trx(nt,i,j)+qy*try(nt,i,j)
              ftru(nt,i,j)=ftru(nt,i,j)+fdt
              fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j)
     .                    +(qx *trd(nt,i,j)
     .                     +qxx*trx(nt,i,j)
     .                     +qxy*try(nt,i,j))*agx(nt,i,j)
     .                    +(qy *trd(nt,i,j)
     .                     +qxy*trx(nt,i,j)
     .                     +qyy*try(nt,i,j))*agy(nt,i,j)
            enddo
            do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j)
     .                    +qx*trx(nt,i,j)+qy*try(nt,i,j)
            enddo
#  else
            do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j)
     .                    +qx*trx(nt,i,j)+qy*try(nt,i,j)
            enddo
#  endif
#endif
c
          endif
c
c --- --- u-component of mass, heat and salt flux.
          uflx(i,j)=fdu(i,j)
          utflx(i,j)=ftu(i,j)
          usflx(i,j)=fsu(i,j)
c
        enddo
        enddo
c
      enddo
c
c --- v-components of fluxes.
c
      do j=1-mrg,jj+mrg+1
c
        do l=1,isv(j)
        do i=max(1-mrg,ifv(j,l)),min(ii+mrg,ilv(j,l))
c
c --- --- Assuming coordinate [0,0] at the v-point, the non-dimensional
c --- --- fluxing area is defined as the area of a polygon with vertices
c --- --- [-1/2,0], [-cuc(i,j)-1/2,-cvc(i,j)], [xm,ym],
c --- --- [-cuc(i+1,j)+1/2,-cvc(i+1,j)], and [1/2,0]. The vertex [xm,ym]
c --- --- is defined so that the polygon area is equal to cv(i,j).
c
          xm=-.5*(cuc(i,j)+cuc(i+1,j))
          ym=((xm+.5)*cvc(i,j)-(xm-.5)*cvc(i+1,j)-2.*cv(i,j))
     .       /(1.+cuc(i,j)-cuc(i+1,j))
c
          if (cv(i,j).gt.0) then
c
            if (cuc(i,j).gt.0.) then
c
c --- ------- Add contributions from grid cell (i-1,j-1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [1/2,yc1+1/2], [1/2,1/2], and
c --- ------- [-cuc(i,j)+1/2,-cvc(i,j)+1/2].
c
              yc0=(ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1=yc0*scp2(i,j-1)*scp2i(i-1,j-1)
              x2=-.5
              y2=yc0+.5
              call triint(scp2(i-1,j-1),
     .          .5,yc1+.5,.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i-1,j-1),max(0.,pbv(i,j)-pup(i-1,j-1)))
              fd=a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdv(i,j)=fdv(i,j)+fd
              qx=ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy=ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftv(i,j)=ftv(i,j)+fd*td(i-1,j-1)
     .                +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsv(i,j)=fsv(i,j)+fd*sd(i-1,j-1)
     .                +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
              qyy=ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
              qxy=axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
              do nt=1,natr
                fdt=fd*trd(nt,i-1,j-1)
     .             +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
                fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i-1,j-1)
     .                      +(qx *trd(nt,i-1,j-1)
     .                       +qxx*trx(nt,i-1,j-1)
     .                       +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1)
     .                      +(qy *trd(nt,i-1,j-1)
     .                       +qxy*trx(nt,i-1,j-1)
     .                       +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i-1,j-1)
     .                      +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i-1,j-1)
     .                      +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
              enddo
#  endif
#endif
            else
              x2=-cuc(i,j)-.5
              y2=-cvc(i,j)+.5
            endif
c
            if (cuc(i+1,j).lt.0.) then
c
c --- ------- Add contributions from grid cell (i+1,j-1). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [-1/2,yc1+1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)+1/2], and
c --- ------- [-1/2,1/2].
c
              yc0=(ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1=yc0*scp2(i,j-1)*scp2i(i+1,j-1)
              x4=.5
              y4=yc0+.5
              call triint(scp2(i+1,j-1),
     .          -.5,yc1+.5,-cuc(i+1,j)-.5,-cvc(i+1,j)+.5,-.5,.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i+1,j-1),max(0.,pbv(i,j)-pup(i+1,j-1)))
              fd=a*dl+ax*dx(i+1,j-1)+ay*dy(i+1,j-1)
              fdv(i,j)=fdv(i,j)+fd
              qx=ax*dl+axx*dx(i+1,j-1)+axy*dy(i+1,j-1)
              qy=ay*dl+axy*dx(i+1,j-1)+ayy*dy(i+1,j-1)
              ftv(i,j)=ftv(i,j)+fd*td(i+1,j-1)
     .                +qx*tx(i+1,j-1)+qy*ty(i+1,j-1)
              fsv(i,j)=fsv(i,j)+fd*sd(i+1,j-1)
     .                +qx*sx(i+1,j-1)+qy*sy(i+1,j-1)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i+1,j-1)+axxy*dy(i+1,j-1)
              qyy=ayy*dl+axyy*dx(i+1,j-1)+ayyy*dy(i+1,j-1)
              qxy=axy*dl+axxy*dx(i+1,j-1)+axyy*dy(i+1,j-1)
              do nt=1,natr
                fdt=fd*trd(nt,i+1,j-1)
     .             +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
                fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i+1,j-1)
     .                      +(qx *trd(nt,i+1,j-1)
     .                       +qxx*trx(nt,i+1,j-1)
     .                       +qxy*try(nt,i+1,j-1))*agx(nt,i+1,j-1)
     .                      +(qy *trd(nt,i+1,j-1)
     .                       +qxy*trx(nt,i+1,j-1)
     .                       +qyy*try(nt,i+1,j-1))*agy(nt,i+1,j-1)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i+1,j-1)
     .                      +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i+1,j-1)
     .                      +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
              enddo
#  endif
#endif
            else
              x4=-cuc(i+1,j)+.5
              y4=-cvc(i+1,j)+.5
            endif
c
c --- ----- Add contributions from grid cell (i,j-1). Assuming
c --- ----- coordinate [0,0] at the cell center, the contributions are
c --- ----- flux integrals over the pentagon with vertices [-1/2,1/2],
c --- ----- [x2,y2], [xm,ym+1/2], [x4,y4], and [1/2,1/2].
c
            call penint(scp2(i,j-1),
     .        -.5,.5,x2,y2,xm,ym+.5,x4,y4,.5,.5,
     .                  a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
            dl=min(dp(i,j-1),max(0.,pbv(i,j)-pup(i,j-1)))
            fd=a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
            fdv(i,j)=fdv(i,j)+fd
            qx=ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
            qy=ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
            ftv(i,j)=ftv(i,j)+fd*td(i,j-1)
     .              +qx*tx(i,j-1)+qy*ty(i,j-1)
            fsv(i,j)=fsv(i,j)+fd*sd(i,j-1)
     .              +qx*sx(i,j-1)+qy*sy(i,j-1)
#ifdef TRC
#  ifdef ATRC
            qxx=axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
            qyy=ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
            qxy=axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
            do nt=1,natr
              fdt=fd*trd(nt,i,j-1)
     .           +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
              ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
              fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i,j-1)
     .                    +(qx *trd(nt,i,j-1)
     .                     +qxx*trx(nt,i,j-1)
     .                     +qxy*try(nt,i,j-1))*agx(nt,i,j-1)
     .                    +(qy *trd(nt,i,j-1)
     .                     +qxy*trx(nt,i,j-1)
     .                     +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
            enddo
            do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i,j-1)
     .                    +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
            enddo
#  else
            do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i,j-1)
     .                    +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
            enddo
#  endif
#endif
c
          else
c
            if (cuc(i,j).gt.0.) then
c
c --- ------- Add contributions from grid cell (i-1,j). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [1/2,yc1-1/2], [1/2,-1/2], and 
c --- ------- [-cuc(i,j)+1/2,-cvc(i,j)-1/2].
c
              yc0=(ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1=yc0*scp2(i,j)*scp2i(i-1,j)
              x2=-.5
              y2=yc0-.5
              call triint(scp2(i-1,j),
     .          .5,yc1-.5,.5,-.5,-cuc(i,j)+.5,-cvc(i,j)-.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i-1,j),max(0.,pbv(i,j)-pup(i-1,j)))
              fd=a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
              fdv(i,j)=fdv(i,j)+fd
              qx=ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
              qy=ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
              ftv(i,j)=ftv(i,j)+fd*td(i-1,j)
     .                +qx*tx(i-1,j)+qy*ty(i-1,j)
              fsv(i,j)=fsv(i,j)+fd*sd(i-1,j)
     .                +qx*sx(i-1,j)+qy*sy(i-1,j)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
              qyy=ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
              qxy=axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
              do nt=1,natr
                fdt=fd*trd(nt,i-1,j)
     .             +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
                fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i-1,j)
     .                      +(qx *trd(nt,i-1,j)
     .                       +qxx*trx(nt,i-1,j)
     .                       +qxy*try(nt,i-1,j))*agx(nt,i-1,j)
     .                      +(qy *trd(nt,i-1,j)
     .                       +qxy*trx(nt,i-1,j)
     .                       +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i-1,j)
     .                      +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i-1,j)
     .                      +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
              enddo
#  endif
#endif
            else
              x2=-cuc(i,j)-.5
              y2=-cvc(i,j)-.5
            endif
c
            if (cuc(i+1,j).lt.0.) then
c
c --- ------- Add contributions from grid cell (i+1,j). Assuming
c --- ------- coordinate [0,0] at the cell center, the contributions are
c --- ------- flux integrals over the triangle with vertices
c --- ------- [-1/2,yc1-1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)-1/2], and
c --- ------- [-1/2,-1/2].
c
              yc0=(ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1=yc0*scp2(i,j)*scp2i(i+1,j)
              x4=.5
              y4=yc0-.5
              call triint(scp2(i+1,j),
     .          -.5,yc1-.5,-cuc(i+1,j)-.5,-cvc(i+1,j)-.5,-.5,-.5,
     .                    a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                   ,axxx,ayyy,axxy,axyy
#endif
     .                   )
              dl=min(dp(i+1,j),max(0.,pbv(i,j)-pup(i+1,j)))
              fd=a*dl+ax*dx(i+1,j)+ay*dy(i+1,j)
              fdv(i,j)=fdv(i,j)+fd
              qx=ax*dl+axx*dx(i+1,j)+axy*dy(i+1,j)
              qy=ay*dl+axy*dx(i+1,j)+ayy*dy(i+1,j)
              ftv(i,j)=ftv(i,j)+fd*td(i+1,j)
     .                +qx*tx(i+1,j)+qy*ty(i+1,j)
              fsv(i,j)=fsv(i,j)+fd*sd(i+1,j)
     .                +qx*sx(i+1,j)+qy*sy(i+1,j)
#ifdef TRC
#  ifdef ATRC
              qxx=axx*dl+axxx*dx(i+1,j)+axxy*dy(i+1,j)
              qyy=ayy*dl+axyy*dx(i+1,j)+ayyy*dy(i+1,j)
              qxy=axy*dl+axxy*dx(i+1,j)+axyy*dy(i+1,j)
              do nt=1,natr
                fdt=fd*trd(nt,i+1,j)
     .             +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
                fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i+1,j)
     .                      +(qx *trd(nt,i+1,j)
     .                       +qxx*trx(nt,i+1,j)
     .                       +qxy*try(nt,i+1,j))*agx(nt,i+1,j)
     .                      +(qy *trd(nt,i+1,j)
     .                       +qxy*trx(nt,i+1,j)
     .                       +qyy*try(nt,i+1,j))*agy(nt,i+1,j)
              enddo
              do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i+1,j)
     .                      +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
              enddo
#  else
              do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i+1,j)
     .                      +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
              enddo
#  endif
#endif
            else
              x4=-cuc(i+1,j)+.5
              y4=-cvc(i+1,j)-.5
            endif
c
c --- ----- Add contributions from grid cell (i,j). Assuming
c --- ----- coordinate [0,0] at the cell center, the contributions are
c --- ----- flux integrals over the pentagon with vertices [-1/2,-1/2],
c --- ----- [x2,y2], [xm,ym-1/2], [x4,y4], and [1/2,-1/2].
c
            call penint(scp2(i,j),
     .        -.5,-.5,x2,y2,xm,ym-.5,x4,y4,.5,-.5,
     .                  a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
     .                 ,axxx,ayyy,axxy,axyy
#endif
     .                 )
            dl=min(dp(i,j),max(0.,pbv(i,j)-pup(i,j)))
            fd=a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdv(i,j)=fdv(i,j)+fd
            qx=ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy=ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftv(i,j)=ftv(i,j)+fd*td(i,j)
     .              +qx*tx(i,j)+qy*ty(i,j)
            fsv(i,j)=fsv(i,j)+fd*sd(i,j)
     .              +qx*sx(i,j)+qy*sy(i,j)
#ifdef TRC
#  ifdef ATRC
            qxx=axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
            qyy=ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
            qxy=axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
            do nt=1,natr
              fdt=fd*trd(nt,i,j)
     .           +qx*trx(nt,i,j)+qy*try(nt,i,j)
              ftrv(nt,i,j)=ftrv(nt,i,j)+fdt
              fagv(nt,i,j)=fagv(nt,i,j)+fdt*agd(nt,i,j)
     .                    +(qx *trd(nt,i,j)
     .                     +qxx*trx(nt,i,j)
     .                     +qxy*try(nt,i,j))*agx(nt,i,j)
     .                    +(qy *trd(nt,i,j)
     .                     +qxy*trx(nt,i,j)
     .                     +qyy*try(nt,i,j))*agy(nt,i,j)
            enddo
            do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i,j)
     .                    +qx*trx(nt,i,j)+qy*try(nt,i,j)
            enddo
#  else
            do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftrv(nt,i,j)=ftrv(nt,i,j)+fd*trd(nt,i,j)
     .                    +qx*trx(nt,i,j)+qy*try(nt,i,j)
            enddo
#  endif
#endif
c
          endif
c
c --- --- v-component of mass, heat and salt flux.
          vflx(i,j)=fdv(i,j)
          vtflx(i,j)=ftv(i,j)
          vsflx(i,j)=fsv(i,j)
c
        enddo
        enddo
c
      enddo
c
c --- ------------------------------------------------------------------
c --- Update fields
c --- ------------------------------------------------------------------
c
      do j=1-mrg,jj+mrg
        do l=1,isp(j)
        do i=max(1-mrg,ifp(j,l)),min(ii+mrg,ilp(j,l))
          q=dp(i,j)
          dp(i,j)=q-(fdu(i+1,j)-fdu(i,j)
     .              +fdv(i,j+1)-fdv(i,j))*scp2i(i,j)
          temp(i,j)=(q*temp(i,j)
     .              -(ftu(i+1,j)-ftu(i,j)
     .               +ftv(i,j+1)-ftv(i,j))*scp2i(i,j))
     .              /dp(i,j)
          saln(i,j)=(q*saln(i,j)
     .              -(fsu(i+1,j)-fsu(i,j)
     .               +fsv(i,j+1)-fsv(i,j))*scp2i(i,j))
     .              /dp(i,j)
#ifdef TRC
#  ifdef ATRC
          do nt=1,natr
            nat=ntr-natr+nt
            trc(i,j,k,nt)=
     .        max(0.,(q*trc(i,j,k,nt)
     .               -(ftru(nt,i+1,j)-ftru(nt,i,j)
     .                +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j))
     .               /dp(i,j)-treps)
            trc(i,j,k,nat)=(q*trc(i,j,k,nat)
     .                     -(fagu(nt,i+1,j)-fagu(nt,i,j)
     .                      +fagv(nt,i,j+1)-fagv(nt,i,j))*scp2i(i,j))
     .                     /dp(i,j)
          enddo
          do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            trc(i,j,k,nt)=(q*trc(i,j,k,nt)
     .                    -(ftru(nt,i+1,j)-ftru(nt,i,j)
     .                     +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j))
     .                    /dp(i,j)
          enddo
#  else
          do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
            if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
            trc(i,j,k,nt)=(q*trc(i,j,k,nt)
     .                    -(ftru(nt,i+1,j)-ftru(nt,i,j)
     .                     +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j))
     .                    /dp(i,j)
          enddo
#  endif
#endif
          dp(i,j)=max(0.,dp(i,j)-dpeps)
        enddo
        enddo
      enddo
c
      return
      end
