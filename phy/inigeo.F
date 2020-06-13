! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen, Alok Kumar Gupta
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

      subroutine inigeo
c
c --- ------------------------------------------------------------------
c --- Initialize the geographic environment
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_nctools
      use mod_dia, only : iotype
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
c
      real mval,fval
      parameter (mval=-1.e12,fval=-1.e13)
c
      real, dimension(itdm,jtdm) :: tmpg
      real rnwp,rmxnbp,rtnbp,rnbp,dx2,dy2,btdtmx,
     .     umaxmin,vmaxmin,umaxmax,vmaxmax
      integer i,j,k,l,kmax
c
c --- ------------------------------------------------------------------
c --- Define bathymetry, grid specification and Coriolis parameter
c --- ------------------------------------------------------------------
c
      if     (expcnf.eq.'cesm') then
        call geoenv_file
        call geoenv_cesmextra
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn'.or.
     .        expcnf.eq.'isomip1'.or.expcnf.eq.'isomip2') then
        call geoenv_file
      elseif (expcnf.eq.'test') then
        call geoenv_test
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(inigeo)')
               stop '(inigeo)'
      endif
c
c --- ------------------------------------------------------------------
c --- Compute auxilary grid parameters
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          scq2i(i,j)=1./max(1.,scq2(i,j))
          scp2i(i,j)=1./max(1.,scp2(i,j))
          scuxi(i,j)=1./max(1.,scux(i,j))
          scvyi(i,j)=1./max(1.,scvy(i,j))
          scuyi(i,j)=1./max(1.,scuy(i,j))
          scvxi(i,j)=1./max(1.,scvx(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Determine do-loop limits for u,v,p,q points
c --- ------------------------------------------------------------------
c
      call bigrid(depths)
c
c --- ------------------------------------------------------------------
c --- Update halos for parameters related to the geographic environment
c --- ------------------------------------------------------------------
c
      call xctilr(qlat, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(qlon, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(plat, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(plon, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(ulat, 1,1, nbdy,nbdy, halo_us)
      call xctilr(ulon, 1,1, nbdy,nbdy, halo_us)
      call xctilr(vlat, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(vlon, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scqx, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqy, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scpx, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scpy, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scux, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvx, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scq2, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scp2, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scu2, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scv2, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scq2i, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scp2i, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scuxi, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvyi, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scuyi, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvxi, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(corioq, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(coriop, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(betafp, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(qclat, 1,4, nbdy,nbdy, halo_qs)
      call xctilr(qclon, 1,4, nbdy,nbdy, halo_qs)
      call xctilr(pclat, 1,4, nbdy,nbdy, halo_ps)
      call xctilr(pclon, 1,4, nbdy,nbdy, halo_ps)
      call xctilr(uclat, 1,4, nbdy,nbdy, halo_us)
      call xctilr(uclon, 1,4, nbdy,nbdy, halo_us)
      call xctilr(vclat, 1,4, nbdy,nbdy, halo_vs)
      call xctilr(vclon, 1,4, nbdy,nbdy, halo_vs)
      if (expcnf.eq.'cesm') then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            util1(i,j)=cplmsk(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
        call xctilr(util1, 1,1, nbdy,nbdy, halo_ps)
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            cplmsk(i,j)=nint(util1(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- ------------------------------------------------------------------
c --- Set mask used for global sums
c --- ------------------------------------------------------------------
c
      if (nreg.eq.2.and.nproc.eq.jpr) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj-1
          do i=1-nbdy,ii+nbdy
            ips(i,j)=ip(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
        do j=jj,jj+nbdy
          do i=1-nbdy,ii+nbdy
            ips(i,j)=0
          enddo
        enddo
      else
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            ips(i,j)=ip(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          util1(i,j)=ip(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcsum(rnwp,util1,ip)
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          util1(i,j)=depths(i,j)*scp2(i,j)
          util2(i,j)=scp2(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcsum(avgbot,util1,ips)
      call xcsum(area,  util2,ips)
      avgbot=avgbot/area
      if     (mnproc.eq.1) then
        if (nwp.ne.nint(rnwp)) then
          write (lp,'(a)') ' xcsum test failed!'
          write (lp,'(a,i7)') ' number of wet points:',nwp
          write (lp,'(a,i7)') ' xcsum on ocean mask: ',nint(rnwp)
          call xchalt('(inigeo)')
                 stop '(inigeo)'
        endif
        write (lp,100) avgbot,area
        call flush(lp)
 100    format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .         -16p,f9.1)
      endif
c
c --- ------------------------------------------------------------------
c --- Set mask for grid cells connected to the world ocean
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          util1(i,j)=ips(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcsum(rnwp,util1,ips)
      if (mnproc.eq.1) then
        write (lp,*) 'Number of wet points',nint(rnwp)
        call flush(lp)
      endif
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          if (ips(i,j).eq.1) then
            util1(i,j)=fval
          else
            util1(i,j)=mval
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      k=0
      rmxnbp=0.
      rtnbp=0.
      do
        k=k+1
        call xcaget(tmpg,util1,1)
        if (mnproc.eq.1) then
          do l=1,itdm*jtdm
            j=(l-1)/itdm+1
            i=l-(j-1)*itdm
            if (tmpg(i,j).eq.fval) then
              tmpg(i,j)=k
              exit
            endif
          enddo
        endif
        call xcaput(tmpg,util1,1)
        call fill_global(mval,fval,halo_ps,util1)
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            if (util1(i,j).eq.mval.or.util1(i,j).eq.fval) then
              util2(i,j)=0.
            else
              if (nint(util1(i,j)).eq.k) then
                util2(i,j)=1.
              else
                util2(i,j)=0.
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(rnbp,util2,ips)
        if (mnproc.eq.1) then
          write (lp,*) 'Number of basin points',nint(rnbp)
          call flush(lp)
        endif
        if (rnbp.gt.rmxnbp) then
          rmxnbp=rnbp
          kmax=k
        endif
        rtnbp=rtnbp+rnbp
        if (nint(rtnbp-rnwp).eq.0) exit
      enddo
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          if (util1(i,j).eq.mval) then
            util1(i,j)=0.
          else
            if (nint(util1(i,j)).eq.kmax) then
              util1(i,j)=1.
            else
              util1(i,j)=0.
            endif
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, nbdy,nbdy, halo_ps)
c$OMP PARALLEL DO PRIVATE(i)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          ipwocn(i,j)=nint(util1(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
c      
      call ncfopn('ipwocn.nc','w','c',1,iotype)
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdefvar('ipwocn','x y',nfint,2)
      call ncedef
      call ncwrti('ipwocn','x y',ipwocn,ip,1)
      call ncfcls
c
c --- ------------------------------------------------------------------
c --- Determine upper bound of lateral diffusivity based on numerical
c --- stability concerns
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i,dx2,dy2)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          difmxp(i,j)=.9*.5*dx2*dy2/max(1.,(dx2+dy2)*(baclin+baclin))
          dx2=scqx(i,j)*scqx(i,j)
          dy2=scqy(i,j)*scqy(i,j)
          difmxq(i,j)=.9*.5*dx2*dy2/max(1.,(dx2+dy2)*(baclin+baclin))
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Estimate maximum barotropic time step
c --- ------------------------------------------------------------------
c
      btdtmx=86400.
c$OMP PARALLEL DO PRIVATE(l,i)
c$OMP+REDUCTION(min:btdtmx)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          btdtmx=min(btdtmx,
     .               scpx(i,j)*scpy(i,j)
     .               /sqrt(g*depths(i,j)*100.*(scpx(i,j)*scpx(i,j)
     .                                        +scpy(i,j)*scpy(i,j))))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(btdtmx)
      if (mnproc.eq.1) then
        write (lp,*) 'estimated max. barotropic time step:',
     .               btdtmx/sqrt(2.)
        call flush(lp)
      endif
c
c --- ------------------------------------------------------------------
c --- Set maximum velocities allowable ensuring stability of the upwind
c --- scheme
c --- ------------------------------------------------------------------
c
      umaxmin=spval
      vmaxmin=spval
      umaxmax=0.
      vmaxmax=0.
c$OMP PARALLEL DO PRIVATE(l,i)
c$OMP+REDUCTION(min:umaxmin,vmaxmin) REDUCTION(max:umaxmax,vmaxmax)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          umax(i,j)=.9*.125*min(scp2(i-1,j),scp2(i,j))
     .              /(scuy(i,j)*baclin)
          umaxmin=min(umaxmin,umax(i,j))
          umaxmax=max(umaxmax,umax(i,j))
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vmax(i,j)=.9*.125*min(scp2(i,j-1 ),scp2(i,j))
     .              /(scvx(i,j)*baclin)
          vmaxmin=min(vmaxmin,vmax(i,j))
          vmaxmax=max(vmaxmax,vmax(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(umax, 1,1, nbdy,nbdy, halo_us)
      call xctilr(vmax, 1,1, nbdy,nbdy, halo_vs)
c
      call xcmin(umaxmin)
      call xcmax(umaxmax)
      call xcmin(vmaxmin)
      call xcmax(vmaxmax)
c
      if (mnproc.eq.1) then
        write (lp,*) 'min/max umax:',umaxmin,umaxmax
        write (lp,*) 'min/max vmax:',vmaxmin,vmaxmax
        call flush(lp)
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inigeo:'
        endif
        call chksummsk(depths,ip,1,'depths')
        call chksummsk(plat,ip,1,'plat')
        call chksummsk(plon,ip,1,'plon')
        call chksummsk(pclat,ip,4,'pclat')
        call chksummsk(pclon,ip,4,'pclon')
        call chksummsk(corioq,iq,1,'corioq')
        call chksummsk(coriop,ip,1,'coriop')
        call chksummsk(betafp,ip,1,'betafp')
        call chksummsk(scqx,iq,1,'scqx')
        call chksummsk(scqy,iq,1,'scqy')
        call chksummsk(scpx,ip,1,'scpx')
        call chksummsk(scpy,ip,1,'scpy')
        call chksummsk(scux,iu,1,'scux')
        call chksummsk(scuy,iu,1,'scuy')
        call chksummsk(scvx,iv,1,'scvx')
        call chksummsk(scvy,iv,1,'scvy')
        call chksummsk(scq2,iq,1,'scq2')
        call chksummsk(scp2,ip,1,'scp2')
        call chksummsk(scu2,iu,1,'scu2')
        call chksummsk(scv2,iv,1,'scv2')
        call chksummsk(scp2i,ip,1,'scp2i')
        call chksummsk(scq2i,iq,1,'scq2i')
        call chksummsk(scuxi,iu,1,'scuxi')
        call chksummsk(scvyi,iv,1,'scvyi')
        call chksummsk(scuyi,iu,1,'scuyi')
        call chksummsk(scvxi,iv,1,'scvxi')
        call chksummsk(difmxp,ip,1,'difmxp')
        call chksummsk(difmxq,iq,1,'difmxq')
        call chksummsk(umax,iu,1,'umax')
        call chksummsk(vmax,iv,1,'vmax')
      endif
c
      return
      end
