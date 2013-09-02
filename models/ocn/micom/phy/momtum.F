      subroutine momtum(m,n,mm,nn,k1m,k1n)
c
c --- version 2.8.1 -- cyclic and noncyclic b.c. combined
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  drag,ubrhs,vbrhs,stress,dpmx,vsc2,vsc4,vort,wgtia,wgtib,
     .  wgtja,wgtjb,dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,uflux1,vflux1
      real cutoff,thkbop,tsfac,dt1inv,pbotl,ptopl,ubot,vbot,ubbl,deform,
     .     dpxy,dpia,dpib,dpja,dpjb,vsc2a,vsc2b,vsc4a,vsc4b,q,botstr,pgf
      integer i,j,k,l,kn,km,kan
      character text*20
c
      real hfharm
      external hfharm
c
      cutoff=onem
      thkbop=thkbot*onem
      tsfac=dlt/delt1
      dt1inv = 1./delt1
c
c$OMP PARALLEL DO PRIVATE(km)
      do j=-1,jj+2
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- bottom drag (standard bulk formula)
c
c$OMP PARALLEL DO PRIVATE(kn,pbotl,ptopl,ubot,vbot,ubbl)
      do 803 j=0,jj
      do 803 l=1,isp(j)
c
      do 800 i=max(0,ifp(j,l)),min(ii,ilp(j,l))
      util1(i,j)=0.
 800  util2(i,j)=0.
c
      do 801 k=1,kk
      kn=k+nn
      do 801 i=max(0,ifp(j,l)),min(ii,ilp(j,l))
      pbotl=max(p(i,j,k+1),p(i,j,kk+1)-thkbop)
      ptopl=max(p(i,j,k  ),p(i,j,kk+1)-thkbop)
      util1(i,j)=util1(i,j)+(u(i,j,kn)+u(i+1,j,kn))*(pbotl-ptopl)
 801  util2(i,j)=util2(i,j)+(v(i,j,kn)+v(i,j+1,kn))*(pbotl-ptopl)
c
      do 803 i=max(0,ifp(j,l)),min(ii,ilp(j,l))
      ubot=(ubflxs_p(i  ,j,n)/max(epsil,pbu(i  ,j,n)*scuy(i  ,j))
     .     +ubflxs_p(i+1,j,n)/max(epsil,pbu(i+1,j,n)*scuy(i+1,j)))*tsfac
     .    +util1(i,j)/thkbop
      vbot=(vbflxs_p(i,j  ,n)/max(epsil,pbv(i,j  ,n)*scvx(i,j  ))
     .     +vbflxs_p(i,j+1,n)/max(epsil,pbv(i,j+1,n)*scvx(i,j+1)))*tsfac
     .    +util2(i,j)/thkbop
      ubbl=.5*sqrt(ubot*ubot+vbot*vbot)
      q=cb*(ubbl+cbar)
      drag(i,j)=q*g/(alpha0*thkbop)
 803  ustarb(i,j)=sqrt(q*ubbl)
c$OMP END PARALLEL DO
c
c --- store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c --- store wind forcing in -stresx,stresy-
c
c$OMP PARALLEL DO
      do 69 j=1,jj
c
      do 68 l=1,isu(j)
      do 68 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 68   ubrhs(i,j)=ubcors_p(i,j)*tsfac
c
      do 69 l=1,isv(j)
      do 69 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 69   vbrhs(i,j)=vbcors_p(i,j)*tsfac
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 804 j=0,jj+1
      do 804 i=0,ii+1
      dl2u(i,j)=0.
 804  dl2v(i,j)=0.
c$OMP END PARALLEL DO
c
      do 8 k=1,kk
      km=k+mm
c$OMP PARALLEL DO
      do 81 j=-1,jj+2
      do 82 l=1,isu(j)
      do 82 i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
 82   pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,km)
      do 81 l=1,isv(j)
      do 81 i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
 81   pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,km)
c$OMP END PARALLEL DO
 8    continue
c
      call xctilr(difwgt,  1,1, 2,2, halo_ps)
c
c$OMP PARALLEL DO
c$OMP+PRIVATE(km,kn,dpmx,wgtja,wgtjb,wgtia,wgtib,i,vort,potvor,j,
c$OMP+  defor1,dl2uja,dl2ujb,dl2via,dl2vib,deform,vsc2,vsc4,dpxy,
c$OMP+  dpja,dpjb,vsc2a,vsc2b,vsc4a,vsc4b,stress,ptopl,pbotl,uflux1,
c$OMP+  botstr,pgf,dpia,dpib,vflux1,q)
c$OMP+FIRSTPRIVATE(defor2,utotm,utotn,uflux,vtotm,vtotn,vflux,uja,ujb,
c$OMP+  dl2u,via,vib,dl2v,uflux2,uflux3,vflux2,vflux3)
      do 9 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
      do 700 j=0,jj+2
      do 700 i=0,ii+2
 700  dpmx(i,j)=8.*cutoff
c
      do 701 j=0,jj+2
      do 701 l=1,isu(j)
      do 701 i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
 701  dpmx(i,j  )=max(dpmx(i,j  ),dp(i,j,km)+dp(i-1,j,km))
c
      do 702 j=-1,jj+1
      do 702 l=1,isu(j)
      do 702 i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
 702  dpmx(i,j+1)=max(dpmx(i,j+1),dp(i,j,km)+dp(i-1,j,km))
c
      do 703 j=0,jj+2
      do 703 l=1,isv(j)
      do 703 i=max(0,ifv(j,l)),min(ii+2,ilv(j,l))
 703  dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,km)+dp(i,j-1,km))
c
      do 704 j=0,jj+2
      do 704 l=1,isv(j)
      do 704 i=max(-1,ifv(j,l)),min(ii+1,ilv(j,l))
 704  dpmx(i+1,j)=max(dpmx(i+1,j),dp(i,j,km)+dp(i,j-1,km))
c
      do 807 j=0,jj+1
      do 807 l=1,isu(j)
      do 807 i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
      utotm(i,j)=u(i,j,km)+ubflxs_p(i,j,m)*tsfac/(pbu(i,j,m)*scuy(i,j))
 807  uflux(i,j)=utotm(i,j)*max(dpu(i,j,km),cutoff)
c
      do 808 j=-1,jj+2
      do 808 l=1,isu(j)
      do 808 i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
 808  utotn(i,j)=u(i,j,kn)+ubflxs_p(i,j,n)*tsfac/(pbu(i,j,n)*scuy(i,j))
c
      do 809 j=0,jj+1
      do 809 l=1,isv(j)
      do 809 i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
      vtotm(i,j)=v(i,j,km)+vbflxs_p(i,j,m)*tsfac/(pbv(i,j,m)*scvx(i,j))
 809  vflux(i,j)=vtotm(i,j)*max(dpv(i,j,km),cutoff)
c
      do 810 j=-1,jj+2
      do 810 l=1,isv(j)
      do 810 i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
 810  vtotn(i,j)=v(i,j,kn)+vbflxs_p(i,j,n)*tsfac/(pbv(i,j,n)*scvx(i,j))
c
c --- define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- sidewall friction along near-vertical bottom slopes. wgtja,wgtjb,wgtia,
c --- wgtib indicate the extent to which a sidewall is present.
c
      do 805 j=-1,jj+2
      do 805 l=1,isu(j)
      do 805 i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
      wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-pbu(i,j-1,m))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
      wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-pbu(i,j+1,m))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
      uja(i,j)=(1.-wgtja(i,j))*utotn(i,j-1)+wgtja(i,j)*slip*utotn(i,j)
      ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,j+1)+wgtjb(i,j)*slip*utotn(i,j)
 805  dl2u(i,j)=utotn(i,j)
     .         -.25*(utotn(i+1,j)+utotn(i-1,j)+uja(i,j)+ujb(i,j))
c --- (to switch from biharmonic to laplacian friction, delete previous line)
c
      do 806 j=0,jj+2
      do 806 l=1,isv(j)
      do 806 i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
      wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-pbv(i-1,j,m))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
      wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-pbv(i+1,j,m))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
      via(i,j)=(1.-wgtia(i,j))*vtotn(i-1,j)+wgtia(i,j)*slip*vtotn(i,j)
      vib(i,j)=(1.-wgtib(i,j))*vtotn(i+1,j)+wgtib(i,j)*slip*vtotn(i,j)
 806  dl2v(i,j)=vtotn(i,j)
     .         -.25*(vtotn(i,j+1)+vtotn(i,j-1)+via(i,j)+vib(i,j))
c --- (to switch from biharmonic to laplacian friction, delete previous line)
c
c --- vorticity, pot.vort., defor. at lateral boundary points
      do 884 j=1,jj
      do 884 l=1,isv(j)
      i=ifv(j,l)
      if (i.ge.1.and.i.le.ii+1) then
        vort(i  ,j)= vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i  ,j)
        potvor(i  ,j)=(vort(i  ,j)+corioq(i  ,j))*8.
     .                /max(4.*(dp(i,j,km)+dp(i,j-1,km)),
     .                     dpmx(i,j),dpmx(i+1,j))
      endif
      i=ilv(j,l)
      if (i.ge.0.and.i.le.ii) then
        vort(i+1,j)=-vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i+1,j)
        potvor(i+1,j)=(vort(i+1,j)+corioq(i+1,j))*8.
     .                /max(4.*(dp(i,j,km)+dp(i,j-1,km)),
     .                     dpmx(i,j),dpmx(i+1,j))
      endif
 884  continue
c
      do 885 j=0,jj+2
      do 885 l=1,isv(j)
      i=ifv(j,l)
      if (i.ge.0) then
        defor2(i  ,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i  ,j)
      endif
      i=ilv(j,l)
      if (i.lt.ii+2) then
        defor2(i+1,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i+1,j)
      endif
 885  continue
c
      do 886 i=1,ii+1
      do 886 l=1,jsu(i)
      j=jfu(i,l)
      if (j.ge.1.and.j.le.jj+1) then
        vort(i,j  )=-utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j  )
        potvor(i,j  )=(vort(i,j  )+corioq(i,j  ))*8.
     .                /max(4.*(dp(i,j,km)+dp(i-1,j,km)),
     .                     dpmx(i,j),dpmx(i,j+1))
      endif
      j=jlu(i,l)
      if (j.ge.0.and.j.le.jj) then
        vort(i,j+1)= utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j+1)
        potvor(i,j+1)=(vort(i,j+1)+corioq(i,j+1))*8.
     .                /max(4.*(dp(i,j,km)+dp(i-1,j,km)),
     .                     dpmx(i,j),dpmx(i,j+1))
      endif
 886  continue
c
      do 887 i=0,ii+2
      do 887 l=1,jsu(i)
      j=jfu(i,l)
      if (j.ge.0) then
        defor2(i,j  )=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,j  )
      endif
      j=jlu(i,l)
      if (j.lt.jj+2) then
        defor2(i,j+1)=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,j+1)
      endif
 887  continue
c
c --- vorticity, pot.vort., defor. at interior points (incl. promontories).
c --- defor1 = du/dx-dv/dy at mass points, defor2 = dv/dx+du/dy at vort. points
c
      do 63 j=-1,jj+1
      do 63 l=1,isp(j)
      do 63 i=max(-1,ifp(j,l)),min(ii+1,ilp(j,l))
 63   defor1(i,j)=((utotn(i+1,j)*scuy(i+1,j)-utotn(i,j)*scuy(i,j))
     .            -(vtotn(i,j+1)*scvx(i,j+1)-vtotn(i,j)*scvx(i,j)))**2
     .            *scp2i(i,j)
c
      do 64 j=1,jj+1
      do 64 l=1,isq(j)
      do 64 i=max(1,ifq(j,l)),min(ii+1,ilq(j,l))
      vort(i,j)=(vtotm(i,j)*scvy(i,j)-vtotm(i-1,j)*scvy(i-1,j)
     .          -utotm(i,j)*scux(i,j)+utotm(i,j-1)*scux(i,j-1))
     .          *scq2i(i,j)
 64   potvor(i,j)=(vort(i,j)+corioq(i,j)) * 8.
     .   /max(2.*(dp(i,j  ,km)+dp(i-1,j  ,km)+
     .            dp(i,j-1,km)+dp(i-1,j-1,km))
     .   ,dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),dpmx(i,j-1),dpmx(i,j+1))
c
      do 65 j=0,jj+2
      do 65 l=1,isq(j)
      do 65 i=max(0,ifq(j,l)),min(ii+2,ilq(j,l))
 65   defor2(i,j)=(vib(i-1,j)*scvy(i,j)-via(i,j)*scvy(i-1,j)
     .            +ujb(i,j-1)*scux(i,j)-uja(i,j)*scux(i,j-1))**2
     .            *scq2i(i,j)
c
c --- define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to imple-
c --- ment biharmonic sidewall friction along near-vertical bottom slopes.
c
      do 906 j=1,jj
c
      do 905 l=1,isu(j)
      do 905 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,j-1)+wgtja(i,j)*slip*dl2u(i,j)
 905  dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,j+1)+wgtjb(i,j)*slip*dl2u(i,j)
c
      do 906 l=1,isv(j)
      do 906 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      dl2via(i,j)=(1.-wgtia(i,j))*dl2v(i-1,j)+wgtia(i,j)*slip*dl2v(i,j)
 906  dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(i+1,j)+wgtib(i,j)*slip*dl2v(i,j)
c
c --- ----------
c --- u equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
      do 37 j=0,jj+1
      do 37 l=1,isu(j)
      do 37 i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
      q=.5*(difwgt(i-1,j)+difwgt(i,j))
      deform=sqrt(.5*(defor1(i,j)+defor1(i-1,j)
     .               +defor2(i,j)+defor2(i,j+1)))
      vsc2(i,j)=max( q*mdv2hi+(1.-q)*mdv2lo,
     .              (q*vsc2hi+(1.-q)*vsc2lo)*deform)
 37   vsc4(i,j)=max( q*mdv4hi+(1.-q)*mdv4lo,
     .              (q*vsc4hi+(1.-q)*vsc4lo)*deform)
c
      do 822 j=1,jj
c
      do 820 l=1,isu(j)
      i=ifu(j,l)
      if (i.gt.0   ) then
        vsc2(i-1,j)=vsc2(i,j)
        vsc4(i-1,j)=vsc4(i,j)
      end if
      i=ilu(j,l)
      if (i.lt.ii+1) then
        vsc2(i+1,j)=vsc2(i,j)
        vsc4(i+1,j)=vsc4(i,j)
      end if
 820  continue
c
c --- longitudinal turb. momentum flux (at mass points)
c
      do 824 l=1,isp(j)
      do 824 i=max(0,ifp(j,l)),min(ii,ilp(j,l))
      if (iu(i,j)+iu(i+1,j).gt.0) then
       dpxy=max(dpu(i  ,j,km),onemm)
       dpib=max(dpu(i+1,j,km),onemm)
       uflux1(i,j)=min(difmxp(i,j),(vsc2(i,j)+vsc2(i+1,j))*scpy(i,j))
     .             *hfharm(dpxy,dpib)*(utotn(i,j)-utotn(i+1,j))
     .            +min(.125*difmxp(i,j),
     .                 (vsc4(i,j)+vsc4(i+1,j))*scpy(i,j))
     .             *hfharm(dpxy,dpib)*(dl2u(i,j)-dl2u(i+1,j))
      endif
 824  continue
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
      do 822 l=1,isu(j)
      do 822 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      dpxy=max(dpu(i,j  ,km),onemm)
      dpja=max(dpu(i,j-1,km),onemm)
      dpja=dpja+wgtja(i,j)*(dpxy-dpja)
      dpjb=max(dpu(i,j+1,km),onemm)
      dpjb=dpjb+wgtjb(i,j)*(dpxy-dpjb)
c
      if (iu(i,j-1).eq.0) then
        vsc2a=vsc2(i,j  )
        vsc4a=vsc4(i,j  )
      else
        vsc2a=vsc2(i,j-1)
        vsc4a=vsc4(i,j-1)
      end if
      if (iu(i,j+1).eq.0) then
        vsc2b=vsc2(i,j  )
        vsc4b=vsc4(i,j  )
      else
        vsc2b=vsc2(i,j+1)
        vsc4b=vsc4(i,j+1)
      end if
      uflux2(i,j)=min(difmxq(i,j  ),(vsc2(i,j)+vsc2a)*scqx(i,j  ))
     .            *hfharm(dpja,dpxy)*(uja(i,j)-utotn(i,j))
     .           +min(.125*difmxq(i,j  ),(vsc4(i,j)+vsc4a)*scqx(i,j  ))
     .            *hfharm(dpja,dpxy)*(dl2uja(i,j)-dl2u(i,j))
 822  uflux3(i,j)=min(difmxq(i,j+1),(vsc2(i,j)+vsc2b)*scqx(i,j+1))
     .            *hfharm(dpjb,dpxy)*(utotn(i,j)-ujb(i,j))
     .           +min(.125*difmxq(i,j+1),(vsc4(i,j)+vsc4b)*scqx(i,j+1))
     .            *hfharm(dpjb,dpxy)*(dl2u(i,j)-dl2ujb(i,j))
c
      if (k.eq.1) then
c
c --- store wind forcing in -stress-
c
      do 70 j=1,jj
      do 70 l=1,isu(j)
      do 70 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 70   stress(i,j)=-2.*taux(i,j)*g*scux(i,j)/(p(i,j,2)+p(i-1,j,2))
      else
      do 72 j=1,jj
      do 72 l=1,isu(j)
      do 72 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 72   stress(i,j)=0.
      endif
c
#ifdef NEST
      do 6 j=2,jj-2
      do 6 l=1,isu(j)
      do 6 i=max(3,ifu(j,l)),min(ii-2,ilu(j,l))
#else
      do 6 j=1,jj
      do 6 l=1,isu(j)
      do 6 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
#endif
c
      ptopl=min(pbu(i,j,m),.5*(p(i,j,k  )+p(i-1,j,k  )))
      pbotl=min(pbu(i,j,m),.5*(p(i,j,k+1)+p(i-1,j,k+1)))
c
c --- bottom boundary layer stress. stress profile is assumed linear and
c --- the drag is treated implicitly
      q=.5*(drag(i,j)+drag(i-1,j))*
     .  (max(pbu(i,j,m)-thkbop,          pbotl       )
     .  -max(pbu(i,j,m)-thkbop,min(ptopl,pbotl-onemm)))
     .  /max(dpu(i,j,km),onemm)
      botstr=-utotn(i,j)*q/(1.+delt1*q)
c
c --- time averaged pressure gradient term
      pgf=(1.-2.*wpgf)*pgfx(i,j,km)+wpgf*(pgfxo(i,j,k)+pgfx(i,j,kn))
c
c --- time smoothing of -u- field  (part 1)
      u(i,j,km)=u(i,j,km)*(wuv1*dpu(i,j,km)+onemm)
     .         +u(i,j,kn)* wuv2*dpuold(i,j,k)
c
 6    u(i,j,kn)=u(i,j,kn)+delt1*(-scuxi(i,j)*(-pgf+stress(i,j)
     .+.25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,j+1)**2
     .     -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,j+1)**2))
     .+.125*(vflux(i  ,j)+vflux(i  ,j+1)+vflux(i-1,j)+vflux(i-1,j+1))
     .     *(potvor(i,j)+potvor(i,j+1)) - ubrhs(i,j) + botstr
     .-(uflux1(i,j)-uflux1(i-1,j)
     . +uflux3(i,j)-uflux2(i  ,j))/(scu2(i,j)*max(dpu(i,j,km),onemm)))
c
c --- ----------
c --- v equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
      do 38 j=0,jj+1
      do 38 l=1,isv(j)
      do 38 i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
      q=.5*(difwgt(i,j-1)+difwgt(i,j))
      deform=sqrt(.5*(defor1(i,j)+defor1(i,j-1)
     .               +defor2(i,j)+defor2(i+1,j)))
      vsc2(i,j)=max( q*mdv2hi+(1.-q)*mdv2lo,
     .              (q*vsc2hi+(1.-q)*vsc2lo)*deform)
 38   vsc4(i,j)=max( q*mdv4hi+(1.-q)*mdv4lo,
     .              (q*vsc4hi+(1.-q)*vsc4lo)*deform)
c
      do 821 i=0,ii+1
      do 821 l=1,jsv(i)
      j=jfv(i,l)
      if (j.gt.0   ) then
        vsc2(i,j-1)=vsc2(i,j)
        vsc4(i,j-1)=vsc4(i,j)
      end if
      j=jlv(i,l)
      if (j.lt.jj+1) then
        vsc2(i,j+1)=vsc2(i,j)
        vsc4(i,j+1)=vsc4(i,j)
      end if
 821  continue
c
c --- longitudinal turb. momentum flux (at mass points)
c
      do 825 j=0,jj
      do 825 l=1,isp(j)
      do 825 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      if (iv(i,j)+iv(i,j+1).gt.0) then
       dpxy=max(dpv(i,j  ,km),onemm)
       dpjb=max(dpv(i,j+1,km),onemm)
       vflux1(i,j)=min(difmxp(i,j),(vsc2(i,j)+vsc2(i,j+1))*scpx(i,j))
     .             *hfharm(dpxy,dpjb)*(vtotn(i,j)-vtotn(i,j+1))
     .            +min(.125*difmxp(i,j),
     .                 (vsc4(i,j)+vsc4(i,j+1))*scpx(i,j))
     .             *hfharm(dpxy,dpjb)*(dl2v(i,j)-dl2v(i,j+1))
      endif
 825  continue
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
      do 823 j=1,jj
      do 823 l=1,isv(j)
      do 823 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c
      dpxy=max(dpv(i  ,j,km),onemm)
      dpia=max(dpv(i-1,j,km),onemm)
      dpia=dpia+wgtia(i,j)*(dpxy-dpia)
      dpib=max(dpv(i+1,j,km),onemm)
      dpib=dpib+wgtib(i,j)*(dpxy-dpib)
c
      if (iv(i-1,j).eq.0) then
        vsc2a=vsc2(i  ,j)
        vsc4a=vsc4(i  ,j)
      else
        vsc2a=vsc2(i-1,j)
        vsc4a=vsc4(i-1,j)
      end if
      if (iv(i+1,j).eq.0) then
        vsc2b=vsc2(i  ,j)
        vsc4b=vsc4(i  ,j)
      else
        vsc2b=vsc2(i+1,j)
        vsc4b=vsc4(i+1,j)
      end if
      vflux2(i,j)=min(difmxq(i  ,j),(vsc2(i,j)+vsc2a)*scqy(i  ,j))
     .            *hfharm(dpia,dpxy)*(via(i,j)-vtotn(i,j))
     .           +min(.125*difmxq(i  ,j),(vsc4(i,j)+vsc4a)*scqy(i  ,j))
     .            *hfharm(dpia,dpxy)*(dl2via(i,j)-dl2v(i,j))
 823  vflux3(i,j)=min(difmxq(i+1,j),(vsc2(i,j)+vsc2b)*scqy(i+1,j))
     .            *hfharm(dpib,dpxy)*(vtotn(i,j)-vib(i,j))
     .           +min(.125*difmxq(i+1,j),(vsc4(i,j)+vsc4b)*scqy(i+1,j))
     .            *hfharm(dpib,dpxy)*(dl2v(i,j)-dl2vib(i,j))
c
      if (k.eq.1) then
c
c --- store wind forcing in -stress-
c
      do 71 j=1,jj
      do 71 l=1,isv(j)
      do 71 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 71   stress(i,j)=-2.*tauy(i,j)*g*scvy(i,j)/(p(i,j,2)+p(i,j-1,2))
      else
      do 73 j=1,jj
      do 73 l=1,isv(j)
      do 73 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 73   stress(i,j)=0.
      endif
c
#ifdef NEST
      do 7 j=3,jj-2
      do 7 l=1,isv(j)
      do 7 i=max(2,ifv(j,l)),min(ii-2,ilv(j,l))
#else
      do 7 j=1,jj
      do 7 l=1,isv(j)
      do 7 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
#endif
c
      ptopl=min(pbv(i,j,m),.5*(p(i,j,k  )+p(i,j-1,k  )))
      pbotl=min(pbv(i,j,m),.5*(p(i,j,k+1)+p(i,j-1,k+1)))
c
c --- bottom boundary layer stress. stress profile is assumed linear and
c --- the drag is treated implicitly
      q=.5*(drag(i,j)+drag(i,j-1))*
     .  (max(pbv(i,j,m)-thkbop,          pbotl       )
     .  -max(pbv(i,j,m)-thkbop,min(ptopl,pbotl-onemm)))
     .  /max(dpv(i,j,km),onemm)
      botstr=-vtotn(i,j)*q/(1.+delt1*q)
c
c --- time averaged pressure gradient term
      pgf=(1.-2.*wpgf)*pgfy(i,j,km)+wpgf*(pgfyo(i,j,k)+pgfy(i,j,kn))
c
c --- time smoothing of -v- field  (part 1)
      v(i,j,km)=v(i,j,km)*(wuv1*dpv(i,j,km)+onemm)
     .         +v(i,j,kn)* wuv2*dpvold(i,j,k)
c
 7    v(i,j,kn)=v(i,j,kn)+delt1*(-scvyi(i,j)*(-pgf+stress(i,j)
     .+.25*(vtotm(i,j+1)**2+utotm(i,j  )**2+utotm(i+1,j  )**2
     .     -vtotm(i,j-1)**2-utotm(i,j-1)**2-utotm(i+1,j-1)**2))
     .-.125*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,j-1)+uflux(i+1,j-1))
     .     *(potvor(i,j)+potvor(i+1,j)) - vbrhs(i,j) + botstr
     .-(vflux1(i,j)-vflux1(i,j-1)
     . +vflux3(i,j)-vflux2(i,j  ))/(scv2(i,j)*max(dpv(i,j,km),onemm)))
c
 9    continue
c$OMP END PARALLEL DO
c
c --- substitute depth-weighted averages for (u,v) at massless grid points.
c --- (scan layers in top-down direction to save time.)
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
c$OMP PARALLEL DO PRIVATE(km,kn,kan,q)
      do j=1,jj
c
      do 31 l=1,isu(j)
c
      do 32 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 32   utotn(i,j)=0.
      do 33 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 33 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      q=min(dpu(i,j,km),dpu(i,j,kn),onem)
      u(i,j,kn)=(u(i,j,kn)*q+u(i,j,kan)*(onem-q))/onem
      u(i,j,kn)=max(-umax(i,j),min(umax(i,j),u(i,j,kn)+ub(i,j,m)))
     .         -ub(i,j,m)
 33   utotn(i,j)=utotn(i,j)+u(i,j,kn)*dpu(i,j,kn)
      do 31 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 31   utotn(i,j)=utotn(i,j)/pbu_p(i,j)
c
      do 30 l=1,isv(j)
c
      do 34 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 34   vtotn(i,j)=0.
      do 35 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 35 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      q=min(dpv(i,j,km),dpv(i,j,kn),onem)
      v(i,j,kn)=(v(i,j,kn)*q+v(i,j,kan)*(onem-q))/onem
      v(i,j,kn)=max(-vmax(i,j),min(vmax(i,j),v(i,j,kn)+vb(i,j,m)))
     .         -vb(i,j,m)
 35   vtotn(i,j)=vtotn(i,j)+v(i,j,kn)*dpv(i,j,kn)
      do 30 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 30   vtotn(i,j)=vtotn(i,j)/pbv_p(i,j)
      end do
c$OMP END PARALLEL DO
c
c --- time smoothing of -u,v- fields  (part 2)
c
      do 22 k=1,kk
      km=k+mm
      kn=k+nn
c
c$OMP PARALLEL DO
#ifdef NEST
      do 24 j=2,jj-2
      do 24 l=1,isu(j)
      do 24 i=max(3,ifu(j,l)),min(ii-2,ilu(j,l))
#else
      do 24 j=1,jj
      do 24 l=1,isu(j)
      do 24 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
#endif
      u(i,j,kn)=u(i,j,kn)-utotn(i,j)
 24   u(i,j,km)=(u(i,j,km)+u(i,j,kn)*wuv2*dpu(i,j,kn))/
     .(wuv1*dpu(i,j,km)+onemm+wuv2*(dpuold(i,j,k)+dpu(i,j,kn)))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
#ifdef NEST
      do 25 j=3,jj-2
      do 25 l=1,isv(j)
      do 25 i=max(2,ifv(j,l)),min(ii-2,ilv(j,l))
#else
      do 25 j=1,jj
      do 25 l=1,isv(j)
      do 25 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
#endif
      v(i,j,kn)=v(i,j,kn)-vtotn(i,j)
 25   v(i,j,km)=(v(i,j,km)+v(i,j,kn)*wuv2*dpv(i,j,kn))/
     .(wuv1*dpv(i,j,km)+onemm+wuv2*(dpvold(i,j,k)+dpv(i,j,kn)))
c$OMP END PARALLEL DO
c
 22   continue
c
c$OMP PARALLEL DO
      do 866 j=1,jj
c
      do 865 l=1,isu(j)
      do 865 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 865  utotn(i,j)=utotn(i,j)*dt1inv
c
      do 866 l=1,isv(j)
      do 866 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 866  vtotn(i,j)=vtotn(i,j)*dt1inv
c$OMP END PARALLEL DO
c
c --- store 'old' interface pressures in -pu,pv- (to be used later for
c --- momentum redistribution)
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
        do 16 k=1,kk
        kn=k+nn
c
        do 15 l=1,isu(j)
        do 15 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 15     pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
c
        do 16 l=1,isv(j)
        do 16 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 16     pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
      end do
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'momtum:'
        endif
        call chksummsk(drag,ip,1,'drag')
        call chksummsk(ubrhs,iu,1,'ubrhs')
        call chksummsk(vbrhs,iv,1,'vbrhs')
        call chksummsk(dpu,iu,2*kk,'dpu')
        call chksummsk(dpv,iv,2*kk,'dpv')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
        call chksummsk(utotn,iu,1,'utotn')
        call chksummsk(vtotn,iv,1,'vtotn')
      endif
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 1995 - changed min.depth in pot.vort. calculation from 1 mm to 1 m
c>             (loops 812,802,803)
c> July 1997 - transferred -visc- and -vort- from common to local
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - added some loops that used to be in barotp.f
c> Aug. 1997 - transferred -wgtia,wgtib,wgtja,wgtjb- from common to local
c> Mar. 1998 - added -th3d-
c> Sep. 1998 - switched to biharmonic friction
c> Apr. 2000 - changed i/j loop nesting to j/i
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Aug. 2000 - separated dpmx(i,j) and dpmx(i,jb) calculation in loop 807 
c>             to avoid coincidental execution in multi-threaded runs
