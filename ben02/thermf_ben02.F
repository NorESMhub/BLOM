! ------------------------------------------------------------------------------
! Copyright (C) 2002-2018 Mats Bentsen
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

      subroutine thermf_ben02(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. 
c
      use mod_xc
      use mod_swtfrz, only: swtfrz
      use mod_swabs, only: swbgal,swbgfc
      use mod_ben02, only: tsi_tda,tml_tda,sml_tda,alb_tda,fice_tda,tsi,
     .                     ntda
      use mod_thdysi, only: tsrfm,ticem,albi_f,albi_m,albs_f,albs_m,
     .                      rhoice,rhosnw,rkice,rksnw,fusi,fuss,
     .                      fice_max,tice_m,tsnw_m,hice_nhmn,hice_shmn,
     .                      sagets,sice,cwi,cuc
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
c --- Due to conservation, the ratio of ice and snow density must be
c --- equal to the ratio of ice and snow heat of fusion
      if (abs(fuss/fusi-rhosnw/rhoice).gt.epsil) then
        if (mnproc.eq.1) then
          write (lp,*)
     .      'thermf: check consistency between snow/ice densities'
          write (lp,*)
     .      'and heat of fusion!'
          stop
        endif
      endif
c
c --- Set various constants
      dt=baclin                         ! Time step
      cpsw=spcifh*1.e3                  ! Specific heat of seawater
      rnf_fac=baclin/real(nrfets*86400) ! Runoff reservoar detrainment rate
      sag_fac=exp(-sagets*dt)            ! Snow aging rate
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
          fice=ficem(i,j)
          hice=hicem(i,j)
          hsnw=hsnwm(i,j)
          tsrf=tsrfm(i,j)
c
          fice0=fice
          hice0=hice
          hsnw0=hsnw
c
c --- --- Freezing point of sea water (in k)
          tice_f=swtfrz(p(i,j,1),sotl)+t0deg
c
c --- --- Minmimum ice thickness
          if (plat(i,j).gt.0.) then
            hice_min=hice_nhmn
          else
            hice_min=hice_shmn
          endif
c
          if     (fice*hice.lt.1.e-5) then
c
c --- ------------------------------------------------------------------
c --- ----- At most, a small amount of ice is present in the grid cell.
c --- ----- Melt the remainders of ice and snow.
c --- ------------------------------------------------------------------
c
            hice=0.
            hsnw=0.
            fice=0.
c
c --- ----- Mean albedo of grid cell
            alb(i,j)=albw(i,j)
c
c --- ----- Solar heat flux that enters the ocean
            qsww=swa(i,j)
c
c --- ----- Non solar heat flux that enters the ocean
            qnsw=nsf(i,j)
c
c --- ----- Set surface temperature and ice surface temperature
            tsrf=totl
            tice=totl
c
          else
c
c --- ------------------------------------------------------------------
c --- ----- Do thermodynamics for an ice slab
c --- ------------------------------------------------------------------
c
            if (fice*hsnw.gt.1.e-3) then
c
c --- ------- Set various variables in the case of a snow layer
c
c --- ------- Albedo
              if (tsrf.gt.tsnw_m-.1) then
                albi=albs_m
              else
                albi=albs_f
              endif
c
c --- ------- Surface melting temperature
              tsmlt=tsnw_m
c
            else
c
c --- ------- Set various variables in the case a thin or non existent
c --- ------- snow layer
c
c --- ------- Albedo
              albi_h=.065+.44*hice**.28
              if (tsrf.gt.tice_m-.1) then
                albi=min(albi_m,albi_h)
              else
                albi=min(albi_f,albi_h)
              endif
c
c --- ------- Surface melting temperature
              tsmlt=tice_m
c
            endif
c
c --- ----- Mean albedo of the grid cell
            alb(i,j)=albi*fice+albw(i,j)*(1.-fice)
c
c --- ----- Short wave radiation trough the ice covered fraction
            qswi=swa(i,j)*(1.-albi)/(1.-alb(i,j))
c
c --- ----- Solar heat flux trough the open water fraction
            qsww=swa(i,j)*(1.-albw(i,j))/(1.-alb(i,j))
c
c --- ----- Update snow thickness due to precipitation
            dh=sop(i,j)*dt/rhosnw
            hsnw=hsnw+dh
c
c --- ----- Heat flux from snow to ice to balance the latent heat of
c --- ----- snow fall
            qsnwf=dh*fuss/dt
c
c --- ----- Conductive factor in snow and ice layer
            fcond=rkice*rksnw/(rksnw*hice+rkice*hsnw)
c
c --- ----- Find the snow surface temperature and the non solar heat
c --- ----- flux that enters the open water fraction
            if (abs(fcond-dfl(i,j)*(2.-fice)).lt.1.e-3) then
              tsrf=tice_f+(qswi+nsf(i,j))/fcond
              qnsw=nsf(i,j)
              qdamp=0.
            else
              tsrf=(qswi+nsf(i,j)-dfl(i,j)*(tsi(i,j)+(1.-fice)*totl)
     .             +fcond*tice_f)
     .             /(fcond-dfl(i,j)*(2.-fice))
              qnsw=nsf(i,j)+dfl(i,j)*fice*(totl-min(tsrf,tsmlt))
              qdamp=dfl(i,j)*(min(tsrf,tsmlt)-tsi(i,j))
            endif
c
c --- ----- If the new surface temperature is above the snow melting
c --- ----- temperature, determine the heat that goes to melting at the
c --- ----- surface
            if (tsrf.gt.tsmlt) then
              tsrf=tsmlt
              qsmlt=qswi+nsf(i,j)
     .             +dfl(i,j)*((1.-fice)*(tsrf-totl)+tsrf-tsi(i,j))
     .             +fcond*(tice_f-tsrf)
            else
              qsmlt=0.
            endif
c
c --- ----- Set ice surface temperature
            tice=tice_f-fcond*(tice_f-tsrf)*hice/rkice
c
c --- ----- Heat flux from ocean to ice (Maykut and McPhee 1995)
            qo2i=rhowat*cpsw*cwi*max(ustari(i,j),.2e-2)
     .           *min(tice_f-totl,0.)+cuc*max(tice_f-totl,0.)
c
c --- ----- Heat budget at bottom of ice
            qbot=-fcond*(tice_f-tsrf)-qo2i-qdamp+qsnwf
c
c --- ----- Update snow thickness due to melting
            dh=-qsmlt*dt/fuss
            if (hsnw+dh.lt.0.) then
              qsmlt=qsmlt-hsnw*fuss/dt
              hsnw=0.
            else
              qsmlt=0.
              hsnw=hsnw+dh
            endif
c
c --- ----- Update ice thickness due to melting/freezing
            hice=max(0.,hice-(qbot+qsmlt)*dt/fusi)
c
c --- ----- Convert snow to ice due to aging
            hice=hice+hsnw*(1.-sag_fac)*rhosnw/rhoice
            hsnw=hsnw*sag_fac
c
c --- ----- Convert snow to ice if snow load is larger than the updrift
c --- ----- of ice
            dh=(hsnw*rhosnw-hice*(rhowat-rhoice))/rhowat
            if (dh.gt.0.) then
              hice=hice+dh
              hsnw=hsnw-dh*rhoice/rhosnw
            endif
c
          endif
c
c --- ------------------------------------------------------------------
c --- --- Do thermodynamics for open water fraction of the grid cell
c --- ------------------------------------------------------------------
c
c --- --- Predict temperature change in mixed layer after a leapfrog
c --- --- time step due to heat fluxes
          swfac=1.-swbgfc(i,j)*exp(-hotl/swbgal(i,j))
          dtml=(swfac*qsww+qnsw)*2.*dt/(cpsw*rhowat*hotl)
c
          if     (totl+dtml.lt.tice_f) then
c
c --- ----- Heat flux required to change the mixed layer temperature
c --- ----- to the freezing point after a leapfrog time step
            q=.5*(tice_f-totl)*cpsw*rhowat*hotl/dt
c
c --- ----- Ice volume that has to freeze to balance the heat budget
            volice=-(qsww+qnsw-q)*(1.-fice)*dt/fusi
c
            if (volice.gt.epsil) then
c
c --- ------- New ice in the lead is formed with a specified thickness.
c --- ------- Estimate the change in ice fraction
              df=volice/hice_min
c
c --- ------- Redistribute ice and snow over an updated ice fraction
              hice=(hice*fice+volice)/min(fice_max,fice+df)
              hsnw=hsnw*fice/min(fice_max,fice+df)
              fice=min(fice_max,fice+df)
c
            endif
c
          elseif (swfac*qsww+qnsw.gt.0.) then
c
c --- ----- If the lead is warming, let the fraction  (1 - fice)  go to
c --- ----- warm the lead, and the fraction  fice  to melt ice laterally
            fice=fice-(swfac*qsww+qnsw)*fice*dt
     .                /max(hice*fusi+hsnw*fuss,epsil)
            if (fice.lt.0.) then
              fice=0.
              hice=0.
              hsnw=0.
            endif
c
          endif
c
c --- ------------------------------------------------------------------
c --- --- Store the updated ice/snow state variables
c --- ------------------------------------------------------------------
c
          ficem(i,j)=fice
          hicem(i,j)=hice
          hsnwm(i,j)=hsnw
          tsrfm(i,j)=tsrf
          ticem(i,j)=tice
c
c --- ------------------------------------------------------------------
c --- --- Accumutate variables to produce averages in flux calculations
c --- ------------------------------------------------------------------
c
          alb_tda(i,j)=alb_tda(i,j)+alb(i,j)
          tml_tda(i,j)=tml_tda(i,j)+totl
          sml_tda(i,j)=sml_tda(i,j)+sotl
          fice_tda(i,j)=fice_tda(i,j)+fice
          tsi_tda(i,j)=tsi_tda(i,j)+tsrf
c
c --- ------------------------------------------------------------------
c --- --- Compute fluxes of heat and salt to the ocean
c --- ------------------------------------------------------------------
c
c --- --- Ice volume change
          dvi=hice*fice-hice0*fice0
c
c --- --- Snow volume change
          dvs=hsnw*fice-hsnw0*fice0
c
c --- --- Accumulate the runoff in a reservoar to delay the discharge
c --- --- into the ocean (by nrfets days approximately 1/e of runoff
c --- --- added will by discharged).
          rnfres(i,j)=rnfres(i,j)+rnf(i,j)
          rnfflx(i,j)=rnfres(i,j)*rnf_fac
          rnfres(i,j)=rnfres(i,j)*(1.-rnf_fac)
          rfiflx(i,j)=rfi(i,j)
c
c --- --- Fresh water flux due to melting/freezing [kg m-2 s-1]
c --- --- (positive downwards)
          fmltfz(i,j)=-(dvi*rhoice+dvs*rhosnw)/dt
c
c --- --- Fresh water flux [kg m-2 s-1] (positive downwards)
          fwflx=eva(i,j)+lip(i,j)+sop(i,j)+rnfflx(i,j)+rfiflx(i,j)
     .         +fmltfz(i,j)
c
c --- --- Salt flux [kg m-2 s-1] (positive downwards)
          sfl(i,j)=-sice*dvi*rhoice/dt*1.e-3
c
c --- --- Salt flux due to brine rejection of freezing sea
c --- --- ice [kg m-2 m-1] (positive downwards)
          brnflx(i,j)=max(0.,-sotl*fmltfz(i,j)*1.e-3+sfl(i,j))
c
c --- --- Virtual salt flux [kg m-2 s-1] (positive downwards)
          vrtsfl(i,j)=-sotl*fwflx*1.e-3
c
c --- --- Store area weighted virtual salt flux and fresh water flux
          util1(i,j)=vrtsfl(i,j)*scp2(i,j)
          util2(i,j)=fwflx*scp2(i,j)
c
c --- --- Heat flux due to melting/freezing [W m-2] (positive downwards)
          hmltfz(i,j)=(dvi*fusi+dvs*fuss)/dt
c
c --- --- Total heat flux in BLOM units [W cm-2] (positive upwards)
          surflx(i,j)=-(swa(i,j)+nsf(i,j)+hmltfz(i,j))*1.e-4
c
c --- --- Short-wave heat flux in BLOM units [W cm-2] (positive
c --- --- upwards)
          sswflx(i,j)=-qsww*(1.-fice0)*1.e-4
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
            sstc=(1.-rice)*max(sstc,tice_f)+rice*tice_f
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
c --- ------------------------------------------------------------------
c --- --- Update age of ice
c --- ------------------------------------------------------------------
c
          if (fice*hice.lt.1.e-5) then
            iagem(i,j)=0.
          else
            iagem(i,j)=(iagem(i,j)+dt/86400.)
     .                 *(1.-max(0.,dvi)/(fice*hice))
          endif
c
c --- -------------------------------------------------------------------
c --- --- Compute friction velocity (cm/s)
c --- -------------------------------------------------------------------
c
          ustar(i,j)=(min(ustari(i,j),.8e-2)*fice0
     .               +ustarw(i,j)*(1.-fice0))*1.e2
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
        write (lp,*) 'thermf: totsfl/area,sflxc',totsfl/area,sflxc
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
cc$OMP PARALLEL DO PRIVATE(l,i)
c        do j=1,jj
c          do l=1,isp(j)
c          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c            util1(i,j)=ttrsf(nt,i,j)
c            util2(i,j)=ttrav(nt,i,j)
c          enddo
c          enddo
c        enddo
cc$OMP END PARALLEL DO
cc
c        call xcsum(tottrsf,util1,ips)
c        call xcsum(tottrav,util2,ips)
cc
c        tottrav=tottrav/area
cc
c        trflxc=(-tottrsf)/area
c        trflxc=(-tottrav*totwfl*1.e-3-tottrsf)/area
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
          write (lp,*) 'thermf_ben02:'
        endif
        call chksummsk(alb,ip,1,'alb')
        call chksummsk(ficem,ip,1,'ficem')
        call chksummsk(hicem,ip,1,'hicem')
        call chksummsk(hsnwm,ip,1,'hsnwm')
        call chksummsk(tsrfm,ip,1,'tsrfm')
        call chksummsk(ticem,ip,1,'ticem')
        call chksummsk(alb_tda,ip,1,'alb_tda')
        call chksummsk(tml_tda,ip,1,'tml_tda')
        call chksummsk(sml_tda,ip,1,'sml_tda')
        call chksummsk(fice_tda,ip,1,'fice_tda')
        call chksummsk(tsi_tda,ip,1,'tsi_tda')
        call chksummsk(rnfres,ip,1,'rnfres')
        call chksummsk(surflx,ip,1,'surflx')
        call chksummsk(sswflx,ip,1,'sswflx')
        call chksummsk(salflx,ip,1,'salflx')
        call chksummsk(surrlx,ip,1,'surrlx')
        call chksummsk(salrlx,ip,1,'salrlx')
        call chksummsk(iagem,ip,1,'iagem')
        call chksummsk(ustar,ip,1,'ustar')
      endif
c
      return
      end
