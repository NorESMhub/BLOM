      subroutine thermf(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. 
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_ice.h"
#include "common_forc.h"
#include "common_geo.h"
#if defined(SYN)
#  include "common_clndr.h"
#endif
c
      integer i,j,k,l,m1,m2,m3,m4,m5
      real dt,h2o_snw,fwf_fac,rnf_fac,sag_fac,y,
     .     tml,sml,hml,tice_f,hice_min,fice,hice,hsnw,tsrf,fice0,hice0,
     .     hsnw0,qsww,qnsw,tice,albi,tsmlt,albi_h,qswi,dh,qsnwf,fcond,
     .     qdamp,qsmlt,qo2i,qbot,swfac,dtml,q,volice,df,dvi,dvs,fwflx,
     .     hflx,sstc,rice,trxflx,sssc,srxflx,
     .     totsfl,totdvi,totwfl,sflxc
      real swtfrz,intpol
      external swtfrz,intpol
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
      dt=baclin                          ! Time step
      h2o_snw=1000./rhosnw               ! Volume converting factor
                                         ! between water and snow
      fwf_fac=1.e-3/86400.               ! Conversion factor mm/day -> m/s
      rnf_fac=baclin/real(nrfets*86400) ! Runoff reservoar detrainment rate
      sag_fac=exp(-gamma*dt)             ! Snow aging rate
c
c --- Set parameters for time interpolation when applying diagnosed heat
c --- and salt relaxation fluxes
#if   defined(CLIM)
      y=mod(time,360.)/7.5
#elif defined(SYN)
      y=(real(nday_of_year-1)+mod(time,1.))*48./real(nday_in_year)
#endif
      m3=int(y)+1
      y=y-real(m3-1)
      m1=mod(m3+45,48)+1
      m2=mod(m3+46,48)+1
      m4=mod(m3   ,48)+1
      m5=mod(m3+ 1,48)+1
c
c --- Find time level for diagnosing heat and salt relaxation fluxes
#if   defined(CLIM)
      k=int(mod(time,360.)/7.5)+1
#elif defined(SYN)
      k=int((real(nday_of_year-1)+mod(time,1.))*48.
     .      /real(nday_in_year))+1
#endif
c
      if (ditflx.or.disflx) nflxdi(k)=nflxdi(k)+1
c
c$OMP PARALLEL DO PRIVATE(tml,sml,hml,tice_f,hice_min,fice,hice,hsnw,
c$OMP+  tsrf,fice0,hice0,hsnw0,qsww,qnsw,tice,albi,tsmlt,albi_h,qswi,dh,
c$OMP+  qsnwf,fcond,qdamp,qsmlt,qo2i,qbot,swfac,dtml,q,volice,df,dvi,
c$OMP+  dvs,fwflx,hflx,sstc,rice,trxflx,sssc,srxflx) SHARED(x,y)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Initialize variables describing the state of the mixed layer
c --- --- and ice/snow fraction
          tml=temp(i,j,k1n)+t0deg
          sml=saln(i,j,k1n)
          hml=dp(i,j,k1n)/onem
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
          tice_f=swtfrz(sml)+t0deg
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
            tsrf=tml
            tice=tml
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
            dh=sop(i,j)*fwf_fac*h2o_snw*dt
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
              tsrf=(qswi+nsf(i,j)-dfl(i,j)*(tsi(i,j)+(1.-fice)*tml)
     .             +fcond*tice_f)
     .             /(fcond-dfl(i,j)*(2.-fice))
              qnsw=nsf(i,j)+dfl(i,j)*fice*(tml-min(tsrf,tsmlt))
              qdamp=dfl(i,j)*(min(tsrf,tsmlt)-tsi(i,j))
            endif
c
c --- ----- If the new surface temperature is above the snow melting
c --- ----- temperature, determine the heat that goes to melting at the
c --- ----- surface
            if (tsrf.gt.tsmlt) then
              tsrf=tsmlt
              qsmlt=qswi+nsf(i,j)
     .             +dfl(i,j)*((1.-fice)*(tsrf-tml)+tsrf-tsi(i,j))
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
     .           *min(tice_f-tml,0.)+cuc*max(tice_f-tml,0.)
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
          swfac=1.-(1.-redfac)*exp(-hml/betabl)
          dtml=(swfac*qsww+qnsw)*2.*dt/(cpsw*rhowat*hml)
c
          if     (tml+dtml.lt.tice_f) then
c
c --- ----- Heat flux required to change the mixed layer temperature
c --- ----- to the freezing point after a leapfrog time step
            q=.5*(tice_f-tml)*cpsw*rhowat*hml/dt
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
          tml_tda(i,j)=tml_tda(i,j)+tml
          sml_tda(i,j)=sml_tda(i,j)+sml
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
c
c --- --- Fresh water flux to the ocean
          fwflx=(eva(i,j)+lip(i,j)+sop(i,j)+rnfflx(i,j))*fwf_fac
     .         -dvs/(h2o_snw*dt)
c
c --- --- Heat flux to the ocean
          hflx=swa(i,j)+nsf(i,j)+(dvi*fusi+dvs*fuss)/dt
c
c --- --- Convert heat flux to MICOM type [W cm-2] (positive upwards)
          surflx(i,j)=-hflx*spcifh*100./(cpsw*rhowat)     
c
c --- --- Short-wave heat flux from the ocean [W cm-2] (positive
c --- --- upwards)
          sswflx(i,j)=-qsww*(1.-fice0)*spcifh*100./(cpsw*rhowat)
c
c --- --- Salt flux from ocean [kg m-2 s-1] (positive upwards)
          salflx(i,j)=sml*fwflx-max(0.,sml-sice)*dvi/dt
c
c --- --- Accumulate salt flux, ice volume change, fresh water flux,
c --- --- and area
          util1(i,j)=salflx(i,j)*scp2(i,j)
          util2(i,j)=dvi*scp2(i,j)
          util3(i,j)=fwflx*scp2(i,j)
c
c --- ------------------------------------------------------------------
c --- --- Relaxation fluxes
c --- ------------------------------------------------------------------
c
          surrlx(i,j)=0.
c
c --- --- If  trxday>0 , apply relaxation towards observed sst
          if (trxday.gt.epsil) then
            sstc=intpol(sstclm(i,j,l1),sstclm(i,j,l2),sstclm(i,j,l3),
     .                  sstclm(i,j,l4),sstclm(i,j,l5),x)
            rice=intpol(ricclm(i,j,l1),ricclm(i,j,l2),ricclm(i,j,l3),
     .                  ricclm(i,j,l4),ricclm(i,j,l5),x)
            sstc=(1.-rice)*max(sstc,tice_f)+rice*tice_f
            trxflx=spcifh*100.*hml/(trxday*86400.)
     .             *min(1.,50./hml)*min(1.5,max(-1.5,sstc-tml))
            surrlx(i,j)=-trxflx
          else
            trxflx=0.
          endif
c
c --- --- If aptflx=.true., apply diagnosed relaxation flux
          if (aptflx) then
            surrlx(i,j)=surrlx(i,j)
     .        -intpol(tflxap(i,j,m1),tflxap(i,j,m2),tflxap(i,j,m3),
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
            sssc=intpol(sssclm(i,j,l1),sssclm(i,j,l2),sssclm(i,j,l3),
     .                  sssclm(i,j,l4),sssclm(i,j,l5),x)
            srxflx=100.*hml/(srxday*86400.)
     .             *min(1.,50./hml)*min( .5,max( -.5,sssc-sml))
            salrlx(i,j)=-srxflx
          else
            srxflx=0.
          endif
c
c --- --- If apsflx=.true., apply diagnosed relaxation flux
          if (apsflx) then
            salrlx(i,j)=salrlx(i,j)
     .        -intpol(sflxap(i,j,m1),sflxap(i,j,m2),sflxap(i,j,m3),
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
c --- --- compute friction velocity (cm/s)
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
c --- Compute correction to the salt flux so it is globally consistent
c --- with a salt flux based on some reference salinities. Also convert
c --- salt flux to MICOM units [10e-3 g cm-2 s-1]
c --- ------------------------------------------------------------------
c
#ifdef ARCTIC
      do j=jtdm-j0,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
        enddo
        enddo
      enddo
#endif
      call xcsum(totsfl,util1,ip)
      call xcsum(totdvi,util2,ip)
      call xcsum(totwfl,util3,ip)
c
c --- Correction for salt flux
      sflxc=(srefg*totwfl-(srefi-sice)*totdvi/dt-totsfl)/area
cdiag if (mnproc.eq.1) then
cdiag   write (lp,*) 'thermf: totsfl/area,sflxc',totsfl/area,sflxc
cdiag endif
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          salflx(i,j)=(salflx(i,j)+sflxc)*1.e2/alpha0
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- number of accumulated fields for flux calculations
c --- ------------------------------------------------------------------
c
      ntda=ntda+1
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'thermf:'
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
