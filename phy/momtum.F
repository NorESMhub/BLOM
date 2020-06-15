! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2020 Mats Bentsen, Lars Inge Enstad, Mehmet Ilicak
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

      subroutine momtum(m,n,mm,nn,k1m,k1n)
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
c --- Parameters:
      real c1,c2,c3,slope
      parameter (c1=1.-1.5*.5,c2=1.-.5,c3=2.,slope=.5)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  drag,ubrhs,vbrhs,stress,dpmx,vsc2,vsc4,vort,wgtia,wgtib,
     .  wgtja,wgtjb,dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,ke,
     .  uh_min,uh_max,vh_min,vh_max,cau,cav,uflux1,vflux1
      real cutoff,thkbop,tsfac,dt1inv,pbotl,ptopl,ubot,vbot,ubbl,up,um,
     .     vp,vm,up2,um2,vp2,vm2,uhc,uhm,vhc,vhm,temp1,temp2,deform,
     .     dpxy,dpia,dpib,dpja,dpjb,vsc2a,vsc2b,vsc4a,vsc4b,q,botstr,pgf
      integer i,j,k,l,kn,km,kan
c
      real hfharm
      external hfharm
c
      cutoff=onem
      thkbop=thkbot*onem
      tsfac=dlt/delt1
      dt1inv=1./delt1
c
      if (mommth.eq.'enedis') then
        uh_min=0.
        uh_max=0.
        vh_min=0.
        vh_max=0.
      endif
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
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
c$OMP PARALLEL DO PRIVATE(l,i,k,kn,pbotl,ptopl,ubot,vbot,ubbl,q)
      do j=0,jj
        do l=1,isp(j)
c
        do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
        enddo
c
        do k=1,kk
          kn=k+nn
          do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
            pbotl=max(p(i,j,k+1),p(i,j,kk+1)-thkbop)
            ptopl=max(p(i,j,k  ),p(i,j,kk+1)-thkbop)
            util1(i,j)=util1(i,j)+(u(i,j,kn)+u(i+1,j,kn))*(pbotl-ptopl)
            util2(i,j)=util2(i,j)+(v(i,j,kn)+v(i,j+1,kn))*(pbotl-ptopl)
          enddo
        enddo
c
        do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
          ubot=(ubflxs_p(i  ,j,n)
     .          /max(epsil,pbu(i  ,j,n)*scuy(i  ,j))
     .         +ubflxs_p(i+1,j,n)
     .          /max(epsil,pbu(i+1,j,n)*scuy(i+1,j)))*tsfac
     .        +util1(i,j)/thkbop
          vbot=(vbflxs_p(i,j  ,n)
     .          /max(epsil,pbv(i,j  ,n)*scvx(i,j  ))
     .         +vbflxs_p(i,j+1,n)
     .          /max(epsil,pbv(i,j+1,n)*scvx(i,j+1)))*tsfac
     .        +util2(i,j)/thkbop
          ubbl=.5*sqrt(ubot*ubot+vbot*vbot)
          q=cb*(ubbl+cbar)
          drag(i,j)=q*g/(alpha0*thkbop)
          ustarb(i,j)=sqrt(q*ubbl)
        enddo
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c --- store wind forcing in -stresx,stresy-
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubrhs(i,j)=ubcors_p(i,j)*tsfac
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vbrhs(i,j)=vbcors_p(i,j)*tsfac
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=0,jj+1
        do i=0,ii+1
          dl2u(i,j)=0.
          dl2v(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        km=k+mm
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,km)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,km)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      call xctilr(difwgt,  1,1, 2,2, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ km,kn,j,i,l,dpmx,wgtja,wgtjb,wgtia,wgtib,vort,potvor,defor1,
c$OMP+ dl2uja,dl2ujb,dl2via,dl2vib,ke,up,um,vp,vm,up2,um2,vp2,vm2,uhc,
c$OMP+ uhm,vhc,vhm,cau,cav,temp1,temp2,q,deform,vsc2,vsc4,dpxy,dpja,
c$OMP+ dpjb,dpia,dpib,vsc2a,vsc2b,vsc4a,vsc4b,stress,ptopl,pbotl,uflux1,
c$OMP+ botstr,pgf,vflux1)
c$OMP+FIRSTPRIVATE(
c$OMP+ defor2,utotm,utotn,uflux,vtotm,vtotn,vflux,uja,ujb,dl2u,via,vib,
c$OMP+ dl2v,uh_min,uh_max,vh_min,vh_max,uflux2,uflux3,vflux2,vflux3)
      do k=1,kk
        km=k+mm
        kn=k+nn
c
c --- - store total (barotropic plus baroclinic) flow at old and mid
c --- - time in -utotn,vtotn- and -utotm,vtotm- respectively. store
c --- - minimum thickness values for use in pot.vort. calculation in
c --- - -dpmx-.
c
        do j=0,jj+2
          do i=0,ii+2
            dpmx(i,j)=8.*cutoff
          enddo
        enddo
c
        do j=0,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            dpmx(i,j  )=max(dpmx(i,j  ),dp(i,j,km)+dp(i-1,j,km))
          enddo
          enddo
        enddo
c
        do j=-1,jj+1
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            dpmx(i,j+1)=max(dpmx(i,j+1),dp(i,j,km)+dp(i-1,j,km))
          enddo
          enddo
        enddo
c
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(0,ifv(j,l)),min(ii+2,ilv(j,l))
            dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,km)+dp(i,j-1,km))
          enddo
          enddo
        enddo
c
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+1,ilv(j,l))
            dpmx(i+1,j)=max(dpmx(i+1,j),dp(i,j,km)+dp(i,j-1,km))
          enddo
          enddo
        enddo
c
        do j=0,jj+1
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=u(i,j,km)
     .                +ubflxs_p(i,j,m)*tsfac/(pbu(i,j,m)*scuy(i,j))
            uflux(i,j)=utotm(i,j)*max(dpu(i,j,km),cutoff)
          enddo
          enddo
        enddo
c
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            utotn(i,j)=u(i,j,kn)
     .                +ubflxs_p(i,j,n)*tsfac/(pbu(i,j,n)*scuy(i,j))
          enddo
          enddo
        enddo
c
        do j=0,jj+1
          do l=1,isv(j)
          do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
            vtotm(i,j)=v(i,j,km)
     .                +vbflxs_p(i,j,m)*tsfac/(pbv(i,j,m)*scvx(i,j))
            vflux(i,j)=vtotm(i,j)*max(dpv(i,j,km),cutoff)
          enddo
          enddo
        enddo
c
        do j=-1,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            vtotn(i,j)=v(i,j,kn)
     .                +vbflxs_p(i,j,n)*tsfac/(pbv(i,j,n)*scvx(i,j))
          enddo
          enddo
        enddo
c
c --- - define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- - sidewall friction along near-vertical bottom slopes.
c --- - wgtja,wgtjb,wgtia, wgtib indicate the extent to which a sidewall
c --- - is present.
c
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-pbu(i,j-1,m))
     .                 /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-pbu(i,j+1,m))
     .                 /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            uja(i,j)=(1.-wgtja(i,j))*utotn(i,j-1)
     .              +wgtja(i,j)*slip*utotn(i,j)
            ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,j+1)
     .              +wgtjb(i,j)*slip*utotn(i,j)
            dl2u(i,j)=utotn(i,j)
     .               -.25*(utotn(i+1,j)+utotn(i-1,j)+uja(i,j)+ujb(i,j))
          enddo
          enddo
        enddo
c --- - (to switch from biharmonic to laplacian friction, delete
c --- - previous line)
c
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-pbv(i-1,j,m))
     .                 /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-pbv(i+1,j,m))
     .                 /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            via(i,j)=(1.-wgtia(i,j))*vtotn(i-1,j)
     .              +wgtia(i,j)*slip*vtotn(i,j)
            vib(i,j)=(1.-wgtib(i,j))*vtotn(i+1,j)
     .              +wgtib(i,j)*slip*vtotn(i,j)
            dl2v(i,j)=vtotn(i,j)
     .               -.25*(vtotn(i,j+1)+vtotn(i,j-1)+via(i,j)+vib(i,j))
          enddo
          enddo
        enddo
c --- - (to switch from biharmonic to laplacian friction, delete
c --- - previous line)
c
c --- - vorticity, pot.vort., defor. at lateral boundary points
        do j=1,jj
          do l=1,isv(j)
            i=ifv(j,l)
            if (i.ge.1.and.i.le.ii+1) then
              vort(i  ,j)= vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i  ,j)
              absvor(i  ,j,k)=vort(i  ,j)+corioq(i  ,j)
              dpvor(i  ,j,k)=.125*max(4.*(dp(i,j,km)+dp(i,j-1,km)),
     .                                dpmx(i,j),dpmx(i+1,j))
              potvor(i  ,j)=absvor(i  ,j,k)/dpvor(i  ,j,k)
            endif
            i=ilv(j,l)
            if (i.ge.0.and.i.le.ii) then
              vort(i+1,j)=-vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i+1,j)
              absvor(i+1,j,k)=vort(i+1,j)+corioq(i+1,j)
              dpvor(i+1,j,k)=.125*max(4.*(dp(i,j,km)+dp(i,j-1,km)),
     .                                dpmx(i,j),dpmx(i+1,j))
              potvor(i+1,j)=absvor(i+1,j,k)/dpvor(i+1,j,k)
            endif
          enddo
        enddo
c
        do j=0,jj+2
          do l=1,isv(j)
            i=ifv(j,l)
            if (i.ge.0) then
              defor2(i  ,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2
     .                      *scq2i(i  ,j)
            endif
            i=ilv(j,l)
            if (i.lt.ii+2) then
              defor2(i+1,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2
     .                      *scq2i(i+1,j)
            endif
          enddo
        enddo
c
        do i=1,ii+1
          do l=1,jsu(i)
            j=jfu(i,l)
            if (j.ge.1.and.j.le.jj+1) then
              vort(i,j  )=-utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j  )
              absvor(i,j  ,k)=vort(i,j  )+corioq(i,j  )
              dpvor(i,j  ,k)=.125*max(4.*(dp(i,j,km)+dp(i-1,j,km)),
     .                                dpmx(i,j),dpmx(i,j+1))
              potvor(i,j  )=absvor(i,j  ,k)/dpvor(i,j  ,k)
            endif
            j=jlu(i,l)
            if (j.ge.0.and.j.le.jj) then
              vort(i,j+1)= utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j+1)
              absvor(i,j+1,k)=vort(i,j+1)+corioq(i,j+1)
              dpvor(i,j+1,k)=.125*max(4.*(dp(i,j,km)+dp(i-1,j,km)),
     .                                dpmx(i,j),dpmx(i,j+1))
              potvor(i,j+1)=absvor(i,j+1,k)/dpvor(i,j+1,k)
            endif
          enddo
        enddo
c
        do i=0,ii+2
          do l=1,jsu(i)
            j=jfu(i,l)
            if (j.ge.0) then
              defor2(i,j  )=(utotn(i,j)*(1.-slip)*scux(i,j))**2
     .                      *scq2i(i,j  )
            endif
            j=jlu(i,l)
            if (j.lt.jj+2) then
              defor2(i,j+1)=(utotn(i,j)*(1.-slip)*scux(i,j))**2
     .                      *scq2i(i,j+1)
            endif
          enddo
        enddo
c
c --- - vorticity, pot.vort., defor. at interior points (incl.
c --- - promontories).  defor1 = du/dx-dv/dy at mass points,
c --- - defor2 = dv/dx+du/dy at vort. points
c
        do j=-1,jj+1
          do l=1,isp(j)
          do i=max(-1,ifp(j,l)),min(ii+1,ilp(j,l))
            defor1(i,j)=((utotn(i+1,j)*scuy(i+1,j)
     .                   -utotn(i  ,j)*scuy(i  ,j))
     .                  -(vtotn(i,j+1)*scvx(i,j+1)
     .                   -vtotn(i,j  )*scvx(i,j  )))**2
     .                  *scp2i(i,j)
          enddo
          enddo
        enddo
c
        do j=1,jj+1
          do l=1,isq(j)
          do i=max(1,ifq(j,l)),min(ii+1,ilq(j,l))
            vort(i,j)=(vtotm(i,j)*scvy(i,j)-vtotm(i-1,j)*scvy(i-1,j)
     .                -utotm(i,j)*scux(i,j)+utotm(i,j-1)*scux(i,j-1))
     .                *scq2i(i,j)
            absvor(i,j,k)=vort(i,j)+corioq(i,j)
            dpvor(i,j,k)=.125*max(2.*(dp(i,j  ,km)+dp(i-1,j  ,km)+
     .                                dp(i,j-1,km)+dp(i-1,j-1,km)),
     .                            dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),
     .                            dpmx(i,j-1),dpmx(i,j+1))
            potvor(i,j)=absvor(i,j,k)/dpvor(i,j,k)
          enddo
          enddo
        enddo
c
        do j=0,jj+2
          do l=1,isq(j)
          do i=max(0,ifq(j,l)),min(ii+2,ilq(j,l))
            defor2(i,j)=(vib(i-1,j)*scvy(i,j)-via(i,j)*scvy(i-1,j)
     .                  +ujb(i,j-1)*scux(i,j)-uja(i,j)*scux(i,j-1))**2
     .                  *scq2i(i,j)
          enddo
          enddo
        enddo
c
c --- - define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to
c --- - implement biharmonic sidewall friction along near-vertical
c --- - bottom slopes.
c
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,j-1)
     .                 +wgtja(i,j)*slip*dl2u(i,j)
            dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,j+1)
     .                 +wgtjb(i,j)*slip*dl2u(i,j)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            dl2via(i,j)=(1.-wgtia(i,j))*dl2v(i-1,j)
     .                 +wgtia(i,j)*slip*dl2v(i,j)
            dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(i+1,j)
     .                 +wgtib(i,j)*slip*dl2v(i,j)
          enddo
          enddo
        enddo
c
c --- - compute area weighted kinetic energy according to Arakawa and
c --- - Lamb (1981)
c
        do j=0,jj
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ----- original
c           ke(i,j)=.25*(utotm(i,j)**2+utotm(i+1,j)**2
c    .                  +vtotm(i,j)**2+vtotm(i,j+1)**2)
c
c --- ----- no-conservation
c           ke(i,j)=.25*((utotm(i,j)+utotm(i+1,j))**2
c    .                  +(vtotm(i,j)+vtotm(i,j+1))**2)
c
c --- ----- GOLD version of ARAKAWA
            ke(i,j)=.25*(scu2(i  ,j)*utotm(i  ,j)**2
     .                  +scu2(i+1,j)*utotm(i+1,j)**2
     .                  +scv2(i,j  )*vtotm(i  ,j)**2
     .                  +scv2(i,j+1)*vtotm(i,j+1)**2)/scp2(i,j)
c
c --- ----- milicak version of ARAKAWA
c           ke(i,j)=(scu2(i  ,j)*utotm(i  ,j)**2
c    .              +scu2(i+1,j)*utotm(i+1,j)**2
c    .              +scv2(i,j  )*vtotm(i  ,j)**2
c    .              +scv2(i,j+1)*vtotm(i,j+1)**2)
c    .              /(scu2(i,j)+scu2(i+1,j)+scv2(i,j)+scv2(i,j+1))
c
c --- ----- GOLD Simple Gudonov
c           up=.5*(utotm(i  ,j)+abs(utotm(i  ,j)))
c           um=.5*(utotm(i+1,j)+abs(utotm(i+1,j)))
c           vp=.5*(vtotm(i,j  )+abs(vtotm(i,j  )))
c           vm=.5*(vtotm(i,j+1)+abs(vtotm(i,j+1)))
c           up2=up*up
c           um2=um*um
c           vp2=vp*vp
c           vm2=vm*vm
c           ke(i,j)=.5*(max(up2,um2)+max(vp2,vm2))
c
c --- ----- GOLD Area weighted Gudonov
c           up=.5*(utotm(i  ,j)+abs(utotm(i  ,j)))
c           um=.5*(utotm(i+1,j)+abs(utotm(i+1,j)))
c           vp=.5*(vtotm(i,j  )+abs(vtotm(i,j  )))
c           vm=.5*(vtotm(i,j+1)+abs(vtotm(i,j+1)))
c           up2=up*up*scu2(i,j  )
c           um2=um*um*scu2(i+1,j)
c           vp2=vp*vp*scv2(i,j  )
c           vm2=vm*vm*scv2(i,j+1)
c           ke(i,j)=.5*(max(up2,um2)+max(vp2,vm2))/scp2(i,j)
c
          enddo
          enddo
        enddo

        if (mommth.eq.'enedis') then
c
c --- -- compute min,max values for Sadourny (1975) energy conserving
c --- -- scheme with dissipation
          do j=0,jj+1
            do l=1,isu(j)
            do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
              uhc=.5*utotm(i,j)*(dp(i,j,km)+dp(i-1,j,km))
              uhm=uflux(i,j) ! mostly similar to uhc, should use uflx(i,j,k)?
              if     (abs(uhc).lt..1*abs(uhm)) then
                uhm=10.*uhc
              elseif (abs(uhc).gt.c1*abs(uhm)) then
                if     (abs(uhc).lt.c2*abs(uhm)) then
                  uhc=(3.*uhc+(1.-c2*3.)*uhm)
                elseif (abs(uhc).le.c3*abs(uhm)) then
                  uhc=uhm
                else
                  uhc=slope*uhc+(1.-c3*slope)*uhm
                endif
              endif
              if (uhc.gt.uhm) then
                uh_min(i,j)=uhm
                uh_max(i,j)=uhc
              else
                uh_max(i,j)=uhm
                uh_min(i,j)=uhc
              endif
            enddo
            enddo
            do l=1,isv(j)
            do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
              vhc=.5*vtotm(i,j)*(dp(i,j,km)+dp(i,j-1,km))
              vhm=vflux(i,j) ! mostly similar to vhc, should use vflx(i,j,k)?
              if     (abs(vhc).lt..1*abs(vhm)) then
                vhm=10.*vhc
              elseif (abs(vhc).gt.c1*abs(vhm)) then
                if     (abs(vhc).lt.c2*abs(vhm)) then 
                  vhc=(3.*vhc+(1.-c2*3.)*vhm)
                elseif (abs(vhc).le.c3*abs(vhm)) then 
                  vhc=vhm
                else 
                  vhc=slope*vhc+(1.-c3*slope)*vhm
                endif
              endif
              if (vhc.gt.vhm) then
                vh_min(i,j)=vhm
                vh_max(i,j)=vhc
              else
                vh_max(i,j)=vhm
                vh_min(i,j)=vhc
              endif
            enddo
            enddo
          enddo
c
        endif
c
c ---- compute coriolis advection terms
c
        if     (mommth.eq.'enscon') then
c
c --- --- Sadourny (1975) enstrophy conserving scheme
c
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              cau(i,j)= .125*(vflux(i  ,j)+vflux(i  ,j+1)
     .                       +vflux(i-1,j)+vflux(i-1,j+1))
     .                      *(potvor(i,j)+potvor(i,j+1))
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              cav(i,j)=-.125*(uflux(i,j  )+uflux(i+1,j  )
     .                       +uflux(i,j-1)+uflux(i+1,j-1))
     .                      *(potvor(i,j)+potvor(i+1,j))
            enddo
            enddo
          enddo
c
        elseif (mommth.eq.'enecon') then
c
c --- --- Sadourny (1975) energy conserving scheme
c
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              cau(i,j)= .25*((vflux(i ,j )+vflux(i-1,j  ))
     .                       *potvor(i,j  )
     .                      +(vflux(i,j+1)+vflux(i-1,j+1))
     .                       *potvor(i,j+1))
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              cav(i,j)=-.25*((uflux(i ,j )+uflux(i,j-1  ))
     .                       *potvor(i  ,j)
     .                      +(uflux(i+1,j)+uflux(i+1,j-1))
     .                       *potvor(i+1,j))
            enddo
            enddo
          enddo
c
        elseif (mommth.eq.'enedis') then
c
c --- --- Sadourny (1975) energy conserving scheme with a little bit
c --- --- dissipation
c
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (potvor(i,j+1)*utotm(i,j).eq.0.) then
                temp1=potvor(i,j+1)*((vh_max(i,j+1)+vh_max(i-1,j+1))
     .                              +(vh_min(i,j+1)+vh_min(i-1,j+1)))*.5
              elseif (potvor(i,j+1)*utotm(i,j).lt.0.) then !GOLD
                temp1=potvor(i,j+1)*(vh_max(i,j+1)+vh_max(i-1,j+1))
              else
                temp1=potvor(i,j+1)*(vh_min(i,j+1)+vh_min(i-1,j+1))
              endif
              if (potvor(i,j)*utotm(i,j).eq.0.) then
                temp2=potvor(i,j)*((vh_max(i,j)+vh_max(i-1,j))
     .                            +(vh_min(i,j)+vh_min(i-1,j)))*.5
              elseif (potvor(i,j)*utotm(i,j).lt.0.) then
                temp2=potvor(i,j)*(vh_max(i,j)+vh_max(i-1,j))
              else
                temp2=potvor(i,j)*(vh_min(i,j)+vh_min(i-1,j))
              endif
              cau(i,j)=.25*(temp1+temp2)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (potvor(i+1,j)*vtotm(i,j).eq.0.) then
                temp1=potvor(i+1,j)*((uh_max(i+1,j)+uh_max(i+1,j-1))
     .                              +(uh_min(i+1,j)+uh_min(i+1,j-1)))*.5
              elseif (potvor(i+1,j)*vtotm(i,j).gt.0.) then !GOLD
                temp1=potvor(i+1,j)*(uh_max(i+1,j)+uh_max(i+1,j-1))
              else
                temp1=potvor(i+1,j)*(uh_min(i+1,j)+uh_min(i+1,j-1))
              endif
              if (potvor(i,j)*vtotm(i,j).eq.0.) then
                temp2=potvor(i,j)*((uh_max(i,j)+uh_max(i,j-1))
     .                            +(uh_min(i,j)+uh_min(i,j-1)))*.5
              elseif (potvor(i,j)*vtotm(i,j).gt.0.) then
                temp2=potvor(i,j)*(uh_max(i,j)+uh_max(i,j-1))
              else
                temp2=potvor(i,j)*(uh_min(i,j)+uh_min(i,j-1))
              endif
              cav(i,j)=-.25*(temp1+temp2)
            enddo
            enddo
          enddo
c
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' mommth=',trim(mommth),
     .                        ' is unsupported!'
          endif
          call xcstop('(momtum)')
                 stop '(momtum)'
        endif
c
c --- - ----------
c --- - u equation
c --- - ----------
c
c --- - deformation-dependent eddy viscosity coefficient
c
        do j=0,jj+1
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
            q=.5*(difwgt(i-1,j)+difwgt(i,j))
            deform=sqrt(.5*(defor1(i,j)+defor1(i-1,j)
     .                     +defor2(i,j)+defor2(i,j+1)))
            vsc2(i,j)=max( q*mdv2hi+(1.-q)*mdv2lo,
     .                    (q*vsc2hi+(1.-q)*vsc2lo)*deform)
            vsc4(i,j)=max( q*mdv4hi+(1.-q)*mdv4lo,
     .                    (q*vsc4hi+(1.-q)*vsc4lo)*deform)
          enddo
          enddo
        enddo
c
        do j=1,jj
c
          do l=1,isu(j)
            i=ifu(j,l)
            if (i.gt.0   ) then
              vsc2(i-1,j)=vsc2(i,j)
              vsc4(i-1,j)=vsc4(i,j)
            endif
            i=ilu(j,l)
            if (i.lt.ii+1) then
              vsc2(i+1,j)=vsc2(i,j)
              vsc4(i+1,j)=vsc4(i,j)
            endif
          enddo
c
c --- --- longitudinal turb. momentum flux (at mass points)
c
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
            if (iu(i,j)+iu(i+1,j).gt.0) then
              dpxy=max(dpu(i  ,j,km),onemm)
              dpib=max(dpu(i+1,j,km),onemm)
              uflux1(i,j)=min(difmxp(i,j),
     .                        (vsc2(i,j)+vsc2(i+1,j))*scpy(i,j))
     .                    *hfharm(dpxy,dpib)*(utotn(i,j)-utotn(i+1,j))
     .                   +min(.125*difmxp(i,j),
     .                        (vsc4(i,j)+vsc4(i+1,j))*scpy(i,j))
     .                    *hfharm(dpxy,dpib)*(dl2u(i,j)-dl2u(i+1,j))
            endif
          enddo
          enddo
c
c --- --- lateral turb. momentum flux (at vorticity points)
c --- --- (left and right fluxes are evaluated separately because of
c --- --- sidewalls)
c
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
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
            endif
            if (iu(i,j+1).eq.0) then
              vsc2b=vsc2(i,j  )
              vsc4b=vsc4(i,j  )
            else
              vsc2b=vsc2(i,j+1)
              vsc4b=vsc4(i,j+1)
            endif
            uflux2(i,j)=min(difmxq(i,j  ),(vsc2(i,j)+vsc2a)*scqx(i,j  ))
     .                  *hfharm(dpja,dpxy)*(uja(i,j)-utotn(i,j))
     .                 +min(.125*difmxq(i,j  ),
     .                      (vsc4(i,j)+vsc4a)*scqx(i,j  ))
     .                  *hfharm(dpja,dpxy)*(dl2uja(i,j)-dl2u(i,j))
            uflux3(i,j)=min(difmxq(i,j+1),(vsc2(i,j)+vsc2b)*scqx(i,j+1))
     .                  *hfharm(dpjb,dpxy)*(utotn(i,j)-ujb(i,j))
     .                 +min(.125*difmxq(i,j+1),
     .                      (vsc4(i,j)+vsc4b)*scqx(i,j+1))
     .                  *hfharm(dpjb,dpxy)*(dl2u(i,j)-dl2ujb(i,j))
          enddo
          enddo
        enddo
c
        if (k.eq.1) then
c
c --- --- store wind forcing in -stress-
c
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              stress(i,j)=-2.*taux(i,j)*g*scux(i,j)
     .                     /(p(i,j,2)+p(i-1,j,2))
            enddo
            enddo
          enddo
        else
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              stress(i,j)=0.
            enddo
            enddo
          enddo
        endif
c
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c
            ptopl=.5*(min(pbu(i,j,m),p(i  ,j,k  ))
     .               +min(pbu(i,j,m),p(i-1,j,k  )))
            pbotl=.5*(min(pbu(i,j,m),p(i  ,j,k+1))
     .               +min(pbu(i,j,m),p(i-1,j,k+1)))
c
c --- ----- bottom boundary layer stress. stress profile is assumed
c --- ----- linear and the drag is treated implicitly
            q=.5*(drag(i,j)+drag(i-1,j))
     .        *(max(pbu(i,j,m)-thkbop,          pbotl       )
     .         -max(pbu(i,j,m)-thkbop,min(ptopl,pbotl-onemm)))
     .        /max(dpu(i,j,km),onemm)
            botstr=-utotn(i,j)*q/(1.+delt1*q)
c
c --- ----- time averaged pressure gradient term
            pgf=(1.-2.*wpgf)*pgfx(i,j,km)
     .         +wpgf*(pgfxo(i,j,k)+pgfx(i,j,kn))
c
c --- ----- time smoothing of -u- field  (part 1)
            u(i,j,km)=u(i,j,km)*(wuv1*dpu(i,j,km)+onemm)
     .               +u(i,j,kn)* wuv2*dpuold(i,j,k)
c
            u(i,j,kn)=u(i,j,kn)
     .       +delt1*(-scuxi(i,j)*(-pgf+stress(i,j)+(ke(i,j)-ke(i-1,j)))
     .               +cau(i,j)-ubrhs(i,j)+botstr
     .               -(uflux1(i,j)-uflux1(i-1,j)
     .                +uflux3(i,j)-uflux2(i  ,j))
     .                /(scu2(i,j)*max(dpu(i,j,km),onemm)))
          enddo
          enddo
        enddo
c
c --- - ----------
c --- - v equation
c --- - ----------
c
c --- - deformation-dependent eddy viscosity coefficient
c
        do j=0,jj+1
          do l=1,isv(j)
          do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
            q=.5*(difwgt(i,j-1)+difwgt(i,j))
            deform=sqrt(.5*(defor1(i,j)+defor1(i,j-1)
     .                     +defor2(i,j)+defor2(i+1,j)))
            vsc2(i,j)=max( q*mdv2hi+(1.-q)*mdv2lo,
     .                    (q*vsc2hi+(1.-q)*vsc2lo)*deform)
            vsc4(i,j)=max( q*mdv4hi+(1.-q)*mdv4lo,
     .                    (q*vsc4hi+(1.-q)*vsc4lo)*deform)
          enddo
          enddo
        enddo
c
        do i=0,ii+1
          do l=1,jsv(i)
            j=jfv(i,l)
            if (j.gt.0   ) then
              vsc2(i,j-1)=vsc2(i,j)
              vsc4(i,j-1)=vsc4(i,j)
            endif
            j=jlv(i,l)
            if (j.lt.jj+1) then
              vsc2(i,j+1)=vsc2(i,j)
              vsc4(i,j+1)=vsc4(i,j)
            endif
          enddo
        enddo
c
c --- - longitudinal turb. momentum flux (at mass points)
c
        do j=0,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (iv(i,j)+iv(i,j+1).gt.0) then
              dpxy=max(dpv(i,j  ,km),onemm)
              dpjb=max(dpv(i,j+1,km),onemm)
              vflux1(i,j)=min(difmxp(i,j),
     .                        (vsc2(i,j)+vsc2(i,j+1))*scpx(i,j))
     .                    *hfharm(dpxy,dpjb)*(vtotn(i,j)-vtotn(i,j+1))
     .                   +min(.125*difmxp(i,j),
     .                        (vsc4(i,j)+vsc4(i,j+1))*scpx(i,j))
     .                    *hfharm(dpxy,dpjb)*(dl2v(i,j)-dl2v(i,j+1))
            endif
          enddo
          enddo
        enddo
c
c --- - lateral turb. momentum flux (at vorticity points)
c --- - (left and right fluxes are evaluated separately because of
c --- - sidewalls)
c
        do j=1,jj
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
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
            endif
            if (iv(i+1,j).eq.0) then
              vsc2b=vsc2(i  ,j)
              vsc4b=vsc4(i  ,j)
            else
              vsc2b=vsc2(i+1,j)
              vsc4b=vsc4(i+1,j)
            endif
            vflux2(i,j)=min(difmxq(i  ,j),(vsc2(i,j)+vsc2a)*scqy(i  ,j))
     .                  *hfharm(dpia,dpxy)*(via(i,j)-vtotn(i,j))
     .                 +min(.125*difmxq(i  ,j),
     .                      (vsc4(i,j)+vsc4a)*scqy(i  ,j))
     .                  *hfharm(dpia,dpxy)*(dl2via(i,j)-dl2v(i,j))
            vflux3(i,j)=min(difmxq(i+1,j),(vsc2(i,j)+vsc2b)*scqy(i+1,j))
     .                  *hfharm(dpib,dpxy)*(vtotn(i,j)-vib(i,j))
     .                 +min(.125*difmxq(i+1,j),
     .                      (vsc4(i,j)+vsc4b)*scqy(i+1,j))
     .                  *hfharm(dpib,dpxy)*(dl2v(i,j)-dl2vib(i,j))
          enddo
          enddo
        enddo
c
        if (k.eq.1) then
c
c --- - store wind forcing in -stress-
c
          do j=1,jj
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              stress(i,j)=-2.*tauy(i,j)*g*scvy(i,j)
     .                     /(p(i,j,2)+p(i,j-1,2))
            enddo
            enddo
          enddo
        else
          do j=1,jj
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              stress(i,j)=0.
            enddo
            enddo
          enddo
        endif
c
        do j=1,jj
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c
            ptopl=.5*(min(pbv(i,j,m),p(i,j  ,k  ))
     .               +min(pbv(i,j,m),p(i,j-1,k  )))
            pbotl=.5*(min(pbv(i,j,m),p(i,j  ,k+1))
     .               +min(pbv(i,j,m),p(i,j-1,k+1)))
c
c --- ----- bottom boundary layer stress. stress profile is assumed
c --- ----- linear and the drag is treated implicitly
            q=.5*(drag(i,j)+drag(i,j-1))*
     .        (max(pbv(i,j,m)-thkbop,          pbotl       )
     .        -max(pbv(i,j,m)-thkbop,min(ptopl,pbotl-onemm)))
     .        /max(dpv(i,j,km),onemm)
            botstr=-vtotn(i,j)*q/(1.+delt1*q)
c
c --- ----- time averaged pressure gradient term
            pgf=(1.-2.*wpgf)*pgfy(i,j,km)
     .         +wpgf*(pgfyo(i,j,k)+pgfy(i,j,kn))
c
c --- ----- time smoothing of -v- field  (part 1)
            v(i,j,km)=v(i,j,km)*(wuv1*dpv(i,j,km)+onemm)
     .               +v(i,j,kn)* wuv2*dpvold(i,j,k)
c
            v(i,j,kn)=v(i,j,kn)
     .       +delt1*(-scvyi(i,j)*(-pgf+stress(i,j)+(ke(i,j)-ke(i,j-1)))
     .               +cav(i,j)-vbrhs(i,j)+botstr
     .               -(vflux1(i,j)-vflux1(i,j-1)
     .                +vflux3(i,j)-vflux2(i,j  ))
     .                /(scv2(i,j)*max(dpv(i,j,km),onemm)))
          enddo
          enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
c --- substitute depth-weighted averages for (u,v) at massless grid points.
c --- (scan layers in top-down direction to save time.)
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
c$OMP PARALLEL DO PRIVATE(l,i,k,km,kn,kan,q)
      do j=1,jj
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          utotn(i,j)=0.
        enddo
        do k=1,kk
          km=k+mm
          kn=k+nn
          kan=max(1,k-1)+nn
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            q=min(dpu(i,j,km),dpu(i,j,kn),onem)
            u(i,j,kn)=(u(i,j,kn)*q+u(i,j,kan)*(onem-q))/onem
            u(i,j,kn)=max(-umax(i,j),min(umax(i,j),u(i,j,kn)+ub(i,j,m)))
     .               -ub(i,j,m)
            utotn(i,j)=utotn(i,j)+u(i,j,kn)*dpu(i,j,kn)
          enddo
        enddo
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          utotn(i,j)=utotn(i,j)/pbu_p(i,j)
        enddo
        enddo
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotn(i,j)=0.
        enddo
        do k=1,kk
          km=k+mm
          kn=k+nn
          kan=max(1,k-1)+nn
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            q=min(dpv(i,j,km),dpv(i,j,kn),onem)
            v(i,j,kn)=(v(i,j,kn)*q+v(i,j,kan)*(onem-q))/onem
            v(i,j,kn)=max(-vmax(i,j),min(vmax(i,j),v(i,j,kn)+vb(i,j,m)))
     .               -vb(i,j,m)
            vtotn(i,j)=vtotn(i,j)+v(i,j,kn)*dpv(i,j,kn)
          enddo
        enddo
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotn(i,j)=vtotn(i,j)/pbv_p(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- time smoothing of -u,v- fields  (part 2)
c
      do k=1,kk
        km=k+mm
        kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            u(i,j,kn)=u(i,j,kn)-utotn(i,j)
            u(i,j,km)=(u(i,j,km)+u(i,j,kn)*wuv2*dpu(i,j,kn))
     .                /(wuv1*dpu(i,j,km)+onemm
     .                 +wuv2*(dpuold(i,j,k)+dpu(i,j,kn)))
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            v(i,j,kn)=v(i,j,kn)-vtotn(i,j)
            v(i,j,km)=(v(i,j,km)+v(i,j,kn)*wuv2*dpv(i,j,kn))
     .                /(wuv1*dpv(i,j,km)+onemm
     .                 +wuv2*(dpvold(i,j,k)+dpv(i,j,kn)))
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          utotn(i,j)=utotn(i,j)*dt1inv
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotn(i,j)=vtotn(i,j)*dt1inv
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- store 'old' interface pressures in -pu,pv- (to be used later for
c --- momentum redistribution)
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
          enddo
          enddo
        enddo
      enddo
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
