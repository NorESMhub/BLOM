      subroutine thermf_ccsm(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. To be used when coupled to CCSM
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
#include "common_clndr.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vrtsfl
c
      integer i,j,k,l,m1,m2,m3,m4,m5
      real y,dpotl,hotl,totl,sotl,dpmxl,hmxl,tmxl,smxl,tice_f,fwflx,
     .     sstc,rice,trxflx,sssc,srxflx,totsfl,totwfl,sflxc,totsrp,
     .     totsrn,q
c
#ifdef TKE
#  include "param_tke.h"
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ttrsf,ttrav
      real tottrsf,tottrav,trflxc
#endif
c
      real swtfrz,intpol
      external swtfrz,intpol
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
c$OMP+ dpotl,hotl,totl,sotl,dpmxl,hmxl,tmxl,smxl,tice_f,fwflx,sstc,rice,
c$OMP+ trxflx,sssc,srxflx)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Set some quantities
c --- ------------------------------------------------------------------
c
c --- --- ocean top layer quantities
          dpotl=dp(i,j,k1n)
          hotl=dpotl/onem
          totl=temp(i,j,k1n)+t0deg
          sotl=saln(i,j,k1n)
c
c --- --- ocean mixed layer quantities
          dpmxl=dp(i,j,1+nn)+dp(i,j,2+nn)
          hmxl=dpmxl/onem
          tmxl=(temp(i,j,1+nn)*dp(i,j,1+nn)
     .         +temp(i,j,2+nn)*dp(i,j,2+nn))/dpmxl+t0deg
          smxl=(saln(i,j,1+nn)*dp(i,j,1+nn)
     .         +saln(i,j,2+nn)*dp(i,j,2+nn))/dpmxl
c
          tice_f=swtfrz(sotl)+t0deg
c
c --- ------------------------------------------------------------------
c --- --- Fresh water and salt fluxes
c --- ------------------------------------------------------------------
c
c --- --- Copy runoff flux to rnfflx for consistency with uncoupled
c --- --- configuration
          rnfflx(i,j)=rnf(i,j)
          rfiflx(i,j)=rfi(i,j)
c
c --- --- Fresh water flux [kg m-2 s-1] (positive downwards)
          fwflx=eva(i,j)+lip(i,j)+sop(i,j)+rnfflx(i,j)+rfiflx(i,j)
     .         +fmltfz(i,j)
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
c --- ------------------------------------------------------------------
c --- --- Heat fluxes
c --- ------------------------------------------------------------------
c
c --- --- Freezing/melting potential [J m-2]. A positive flux means the ocean
c --- --- surface has a temperature below freezing temperature and must
c --- --- be heated. Note the freezing potential is multiplied by 1/2
c --- --- due to the leap-frog time stepping. The melting potential uses
c --- --- time averaged quantities since it is not accumulated.
          frzpot(i,j)=max(0.,tice_f-totl)*spcifh*dpotl/(2.*g)*1.e4
          mltpot(i,j)=min(0.,swtfrz(.5*(saln(i,j,k1m)+saln(i,j,k1n)))
     .                      -.5*(temp(i,j,k1m)+temp(i,j,k1n)))
     .                *spcifh*.5*(dp(i,j,k1m)+dp(i,j,k1n))/g*1.e4
c
c --- --- Heat flux due to melting/freezing [W m-2] (positive downwards)
          hmltfz(i,j)=hmlt(i,j)+frzpot(i,j)/baclin
c
c --- --- Total heat flux in MICOM units [W cm-2] (positive upwards)
          surflx(i,j)=-(swa(i,j)+nsf(i,j)+hmltfz(i,j))*1.e-4
c
c --- --- Short-wave heat flux in MICOM units [W cm-2] (positive
c --- --- upwards)
          sswflx(i,j)=-swa(i,j)*1.e-4
c
#ifdef TKE
          tkeflx(i,j)=0.
#  ifdef GLS
          gls_psiflx(i,j)=-gls_n*difdia(i,j,1)*(gls_cmu0**gls_p)*
     .                     (tke(i,j,k1n)**gls_m)*(vonKar**gls_n)*
     .                     Zos**(gls_n-1.)
#  endif
#endif
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
#ifdef TRC
c --- ------------------------------------------------------------------
c --- --- Virtual tracer fluxes (positive downwards)
c --- ------------------------------------------------------------------
c
          do nt=1,ntrc
            trflx(nt,i,j)=-trc(i,j,1,nt)*fwflx*1.e-3
            ttrsf(nt,i,j)=trflx(nt,i,j)*scp2(i,j)
            ttrav(nt,i,j)=trc(i,j,1,nt)*scp2(i,j)
          enddo
#endif
c
c --- -------------------------------------------------------------------
c --- --- Friction velocity (cm/s)
c --- -------------------------------------------------------------------
c
          ustar(i,j)=ustarw(i,j)*1.e2
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
      do nt=1,ntrc
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=ttrsf(nt,i,j)
            util2(i,j)=ttrav(nt,i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcsum(tottrsf,util1,ips)
        call xcsum(tottrav,util2,ips)
c
        tottrav=tottrav/area
c
        trflxc=(-tottrsf)/area
c        trflxc=(-tottrav*totwfl*1.e-3-tottrsf)/area
c        trflxc=0.
c
c$OMP PARALLEL DO
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
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'thermf_ccsm:'
        endif
        call chksummsk(surflx,ip,1,'surflx')
        call chksummsk(sswflx,ip,1,'sswflx')
        call chksummsk(salflx,ip,1,'salflx')
        call chksummsk(brnflx,ip,1,'brnflx')
        call chksummsk(surrlx,ip,1,'surrlx')
        call chksummsk(salrlx,ip,1,'salrlx')
        call chksummsk(rnfflx,ip,1,'rnfflx')
        call chksummsk(rfiflx,ip,1,'rfiflx')
        call chksummsk(hmltfz,ip,1,'hmltfz')
        call chksummsk(ustar,ip,1,'ustar')
        call chksummsk(frzpot,ip,1,'frzpot')
        call chksummsk(mltpot,ip,1,'mltpot')
      endif
c
c
      return
      end
