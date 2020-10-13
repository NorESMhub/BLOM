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

      subroutine thermf_channel(m,n,mm,nn,k1m,k1n)
c
c --- thermf for the channel config
c
      use mod_xc
      use mod_ben02, only: ntda
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_seaice.h"
#include "common_forc.h"
#include "common_geo.h"
#include "common_clndr.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vrtsfl
c
      integer i,j,k,l,m1,m2,m3,m4,m5
      real dt,cpsw,rnf_fac,sag_fac,y,
     .     dpotl,hotl,totl,sotl,dpmxl,hmxl,tmxl,smxl,tice_f,hice_min,
     .     fice,hice,hsnw,tsrf,fice0,hice0,hsnw0,qsww,qnsw,tice,albi,
     .     tsmlt,albi_h,qswi,dh,qsnwf,fcond,qdamp,qsmlt,qo2i,qbot,swfac,
     .     dtml,q,volice,df,dvi,dvs,fwflx,sstc,rice,trxflx,sssc,srxflx,
     .     totsfl,totwfl,sflxc,totsrp,totsrn
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ttrsf,ttrav
      real tottrsf,tottrav,trflxc
#  ifdef TKE
#    include "param_tke.h"
#  endif
#endif
c
      real intp1d
      external intp1d
c
c --- Set various constants
      dt=baclin                         ! Time step
      cpsw=spcifh*1.e3                  ! Specific heat of seawater
      rnf_fac=baclin/real(nrfets*86400) ! Runoff reservoar detrainment rate
c
c --- Set parameters for time interpolation when applying diagnosed heat
c --- and salt relaxation fluxes
      y=(nday_of_year-1+mod(nstep,nstep_in_day)/real(nstep_in_day))*48.
     .  /real(nday_in_year)
      m3=int(y)+1
      y=y-real(m3-1)
      m1=mod(m3+45,48)+1
      m2=mod(m3+46,48)+1
      m4=mod(m3   ,48)+1
      m5=mod(m3+ 1,48)+1
c
c --- Time level for diagnosing heat and salt relaxation fluxes
      k=m3
c
      if (ditflx.or.disflx) nflxdi(k)=nflxdi(k)+1
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,dpotl,hotl,totl,sotl,dpmxl,hmxl,tmxl,smxl,tice_f,hice_min,
c$OMP+ fice,hice,hsnw,tsrf,fice0,hice0,hsnw0,qsww,qnsw,tice,albi,tsmlt,
c$OMP+ albi_h,qswi,dh,qsnwf,fcond,qdamp,qsmlt,qo2i,qbot,swfac,dtml,q,
c$OMP+ volice,df,dvi,dvs,fwflx,sstc,rice,trxflx,sssc,srxflx
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ ) SHARED(xmi,y)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Initialize variables describing the state of the ocean top
c --- --- layer, the mixed layer and ice/snow fraction
          dpotl=dp(i,j,k1n)
          hotl=dpotl/onem
          totl=temp(i,j,k1n)+t0deg
          sotl=saln(i,j,k1n)
c
          dpmxl=dp(i,j,1+nn)+dp(i,j,2+nn)
          hmxl=dpmxl/onem
          tmxl=(temp(i,j,1+nn)*dp(i,j,1+nn)
     .         +temp(i,j,2+nn)*dp(i,j,2+nn))/dpmxl+t0deg
          smxl=(saln(i,j,1+nn)*dp(i,j,1+nn)
     .         +saln(i,j,2+nn)*dp(i,j,2+nn))/dpmxl
c
c
c --- --- Accumulate the runoff in a reservoar to delay the discharge
c --- --- into the ocean (by nrfets days approximately 1/e of runoff
cc --- --- added will by discharged).
          rnfres(i,j)=rnfres(i,j)+rnf(i,j)
          rnfflx(i,j)=rnfres(i,j)*rnf_fac
          rnfres(i,j)=rnfres(i,j)*(1.-rnf_fac)
          rfiflx(i,j)=rfi(i,j)
c
c --- --- Fresh water flux due to melting/freezing [kg m-2 s-1]
c --- --- (positive downwards)
          fmltfz(i,j)=0.0 !-(dvi*rhoice+dvs*rhosnw)/dt
c
c --- --- Fresh water flux [kg m-2 s-1] (positive downwards)
          fwflx=eva(i,j)+lip(i,j)+sop(i,j)+rnfflx(i,j)+rfiflx(i,j)
     .         +fmltfz(i,j)
c
c --- --- Salt flux [kg m-2 s-1] (positive downwards)
          sfl(i,j)=0.0 !-sice*dvi*rhoice/dt*1.e-3
c
c --- --- Salt flux due to brine rejection of freezing sea
c --- --- ice [kg m-2 m-1] (positive downwards)
          brnflx(i,j)=0.0 !max(0.,-sotl*fmltfz(i,j)*1.e-3+sfl(i,j))
c
c --- --- Virtual salt flux [kg m-2 s-1] (positive downwards)
          vrtsfl(i,j)=-sotl*fwflx*1.e-3
c
c --- --- Store area weighted virtual salt flux and fresh water flux
          util1(i,j)=vrtsfl(i,j)*scp2(i,j)
          util2(i,j)=fwflx*scp2(i,j)
c
c --- --- Heat flux due to melting/freezing [W m-2] (positive downwards)
          hmltfz(i,j)=0.0 !(dvi*fusi+dvs*fuss)/dt
c
c --- --- Total heat flux in MICOM units [W cm-2] (positive upwards)
          surflx(i,j)=-(swa(i,j)+nsf(i,j)+hmltfz(i,j))*1.e-4
c
c --- --- Short-wave heat flux in MICOM units [W cm-2] (positive
c --- --- upwards)
          sswflx(i,j)=0.0 !-qsww*(1.-fice0)*1.e-4
c
#ifdef TRC
c --- ------------------------------------------------------------------
c --- --- Tracer fluxes (positive downwards)
c --- ------------------------------------------------------------------
c
          do nt=1,ntr
#  ifdef TKE
            if (nt.eq.itrtke) then
              trflx(nt,i,j)=0.
              ttrsf(nt,i,j)=0.
              ttrav(nt,i,j)=0.
              cycle
            endif
#    ifdef GLS
            if (nt.eq.itrgls) then
              trflx(nt,i,j)=-gls_n*difdia(i,j,1)*(gls_cmu0**gls_p)
     .                       *(trc(i,j,k1n,itrtke)**gls_m)
     .                       *(vonKar**gls_n)*Zos**(gls_n-1.)
              ttrsf(nt,i,j)=0.
              ttrav(nt,i,j)=0.
              cycle
            endif
#    else
            if (nt.eq.itrgls) then
              trflx(nt,i,j)=0.
              ttrsf(nt,i,j)=0.
              ttrav(nt,i,j)=0.
              cycle
            endif
#    endif
#  endif
            trflx(nt,i,j)=-trc(i,j,k1n,nt)*fwflx*1.e-3
            ttrsf(nt,i,j)=trflx(nt,i,j)*scp2(i,j)
            ttrav(nt,i,j)=trc(i,j,k1n,nt)*scp2(i,j)
          enddo
#endif
c
c --- ------------------------------------------------------------------
c --- --- Relaxation fluxes
c --- ------------------------------------------------------------------
c
          surrlx(i,j)=0.
c
c --- --- If  trxday>0 , apply relaxation towards observed sst
          if (trxday.gt.epsil) then
            sstc=intp1d(sstclm(i,j,l1mi),sstclm(i,j,l2mi),
     .                  sstclm(i,j,l3mi),sstclm(i,j,l4mi),
     .                  sstclm(i,j,l5mi),xmi)
            rice=intp1d(ricclm(i,j,l1mi),ricclm(i,j,l2mi),
     .                  ricclm(i,j,l3mi),ricclm(i,j,l4mi),
     .                  ricclm(i,j,l5mi),xmi)
            sstc=(1.-rice)*sstc !max(sstc,tice_f)+rice*tice_f
            trxflx=spcifh*100.*min(hmxl,trxdpt)/(trxday*86400.)
     .                        *min(trxlim,max(-trxlim,sstc-tmxl))
            surrlx(i,j)=-trxflx
          else
            trxflx=0.
          endif
c
c --- --- If aptflx=.true., apply diagnosed relaxation flux
          if (aptflx) then
            surrlx(i,j)=surrlx(i,j)
     .        -intp1d(tflxap(i,j,m1),tflxap(i,j,m2),tflxap(i,j,m3),
     .                tflxap(i,j,m4),tflxap(i,j,m5),y)
          endif
c
c --- --- If ditflx=.true., diagnose relaxation flux by accumulating the
c --- --- relaxation flux
          if (ditflx) then
            tflxdi(i,j,k)=tflxdi(i,j,k)+trxflx
          endif
c
          salrlx(i,j)=0.
c
c --- --- if  srxday>0 , apply relaxation towards observed sss
          if (srxday.gt.epsil) then
            sssc=intp1d(sssclm(i,j,l1mi),sssclm(i,j,l2mi),
     .                  sssclm(i,j,l3mi),sssclm(i,j,l4mi),
     .                  sssclm(i,j,l5mi),xmi)
            srxflx=100.*min(hmxl,srxdpt)/(srxday*86400.)
     .                 *min(srxlim,max(-srxlim,sssc-smxl))
            salrlx(i,j)=-srxflx
            util3(i,j)=max(0.,salrlx(i,j))*scp2(i,j)
            util4(i,j)=min(0.,salrlx(i,j))*scp2(i,j)
          else
            srxflx=0.
          endif
c
c --- --- If apsflx=.true., apply diagnosed relaxation flux
          if (apsflx) then
            salrlx(i,j)=salrlx(i,j)
     .        -intp1d(sflxap(i,j,m1),sflxap(i,j,m2),sflxap(i,j,m3),
     .                sflxap(i,j,m4),sflxap(i,j,m5),y)
          endif
c
c --- --- If disflx=.true., diagnose relaxation flux by accumulating the
c --- --- relaxation flux
          if (disflx) then
            sflxdi(i,j,k)=sflxdi(i,j,k)+srxflx
          endif
c
c --- -------------------------------------------------------------------
c --- --- Compute friction velocity (cm/s)
c --- -------------------------------------------------------------------
c
          ustar(i,j)=ustarw(i,j)*1.e2
c          ustar(i,j)=(min(ustari(i,j),.8e-2)*fice0
c     .               +ustarw(i,j)*(1.-fice0))*1.e2
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute correction to the virtual salt flux so it is globally
c --- consistent with a salt flux based on some reference salinity.
c --- Also combine virtual and true salt flux and convert salt fluxes
c --- used later to unit [10e-3 g cm-2 s-1] and positive upwards.
c --- ------------------------------------------------------------------
c
      call xcsum(totsfl,util1,ips)
      call xcsum(totwfl,util2,ips)
c
c --- Correction for the virtual salt flux [kg m-2 s-1]
      sflxc=(-sref*totwfl*1.e-3-totsfl)/area
      if (mnproc.eq.1) then
        write (lp,*) 'thermf_channel: totsfl/area,sflxc',
     .               totsfl/area,sflxc
      endif
c
c --- Apply the virtual salt flux correction and the compute the total
c --- salt flux by combining the virtual and true salt flux
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          salflx(i,j)=-(vrtsfl(i,j)+sflxc+sfl(i,j))*1.e2
          brnflx(i,j)=-brnflx(i,j)*1.e2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- if  srxday>0  and  srxbal=.true. , balance the sss relaxation flux
c --- so the net input of salt in grid cells connected to the world
c --- ocean is zero
      if (srxday.gt.epsil.and.srxbal) then
        call xcsum(totsrp,util3,ipwocn)
        call xcsum(totsrn,util4,ipwocn)
        if (abs(totsrp).gt.abs(totsrn)) then
          q=-totsrn/totsrp
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (salrlx(i,j).gt.0..and.ipwocn(i,j).eq.1) then
                salrlx(i,j)=q*salrlx(i,j)
              endif
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          q=-totsrp/totsrn
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (salrlx(i,j).lt.0..and.ipwocn(i,j).eq.1) then
                salrlx(i,j)=q*salrlx(i,j)
              endif
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
      endif
c
#ifdef TRC
      do nt=1,ntr
c
#  ifdef TKE
        if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
        trflxc=0.
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trflx(nt,i,j)=-(trflx(nt,i,j)+trflxc)*1.e2
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c 
      enddo
#endif
c
c --- ------------------------------------------------------------------
c --- number of accumulated fields for flux calculations
c --- ------------------------------------------------------------------
c
      ntda=ntda+1
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'thermf_channel:'
        endif
        call chksummsk(rnfres,ip,1,'rnfres')
        call chksummsk(surflx,ip,1,'surflx')
        call chksummsk(sswflx,ip,1,'sswflx')
        call chksummsk(salflx,ip,1,'salflx')
        call chksummsk(surrlx,ip,1,'surrlx')
        call chksummsk(salrlx,ip,1,'salrlx')
        call chksummsk(ustar,ip,1,'ustar')
      endif
c
      return
      end
