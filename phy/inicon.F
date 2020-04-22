! ------------------------------------------------------------------------------
! Copyright (C) 2008-2018 Mats Bentsen
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

      subroutine inicon
c
c --- ------------------------------------------------------------------
c --- Define initial conditions
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_swtfrz, only: swtfrz
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_seaice.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tfrz
c
      integer i,j,k,l
      real q,dps
c
      real sig,sofsig,getlpi
      external sig,sofsig,getlpi
c
c --- ------------------------------------------------------------------
c --- Define layer interface heights and layer temperature and salinity
c --- ------------------------------------------------------------------
      if     (expcnf.eq.'cesm'.or.
     .        expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call ictsz_file
      elseif (expcnf.eq.'isomip1') then
c       call ictsz_isomip1
      elseif (expcnf.eq.'isomip2') then
c       call ictsz_isomip2
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(inicon)')
               stop '(inicon)'
      endif
c
c --- ------------------------------------------------------------------
c --- Set minimum physical temperature for each isopycnic layer
c --- ------------------------------------------------------------------
c
      call settemmin
c
c --- ------------------------------------------------------------------
c --- Initialize configuration specific variables
c --- ------------------------------------------------------------------
c
      if     (expcnf.eq.'cesm') then
        call icaux_cesm
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call icaux_ben02
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(inicon)')
               stop '(inicon)'
      endif
c
c --- ------------------------------------------------------------------
c --- Make sure layer temperature is greater than the lower physical
c --- bound and make temperature, salinity, and potential density
c --- variables consistent.
c --- ------------------------------------------------------------------
c
      do k=1,2
        tfrz(:,:)=swtfrz(p(:,:,1),saln(:,:,k))
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k)=max(tfrz(i,j),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
      do k=3,kk
        tfrz(:,:)=swtfrz(p(:,:,1),saln(:,:,k))
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k)=max(tfrz(i,j),temp(i,j,k))
            saln(i,j,k)=sofsig(sigmar(i,j,k),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar(i,j)=0.
          buoyfl(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (mnproc.eq.ptest) then
        write (lp,'('' sigmar(k)    :'',7f9.5/(15x,7f9.5))')
     .     (sigmar(itest,jtest,k),k=1,kk)
      endif
c
c --- ------------------------------------------------------------------
c --- Find layer interface pressure
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          p(i,j,1)=getlpi(temp(i,j,1),saln(i,j,1),0.,phi(i,j,1),0.)
        enddo
        enddo
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=getlpi(temp(i,j,k),saln(i,j,k),
     .                        phi(i,j,k),phi(i,j,k+1),p(i,j,k))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(p, 1,kk+1, 2,2, halo_ps)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)
      call xctilr(temp, 1,kk, 1,1, halo_ps)
      call xctilr(saln, 1,kk, 1,1, halo_ps)
c
c --- ------------------------------------------------------------------
c --- Set layer thickness and bottom pressure
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=0,jj+1
        do k=1,kk
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
            dp(i,j,k+kk)=dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=0,jj+1
        do k=1,kk
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          pb(i,j,1)=p(i,j,kk+1)
          pb(i,j,2)=pb(i,j,1)
          pb_mn(i,j,1)=pb(i,j,1)
          pb_mn(i,j,2)=pb(i,j,1)
          pb_p(i,j)=pb(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pbu(i,j,1)=min(pb(i,j,1),pb(i-1,j,1))
          pbu(i,j,2)=pbu(i,j,1)
          pbu_p(i,j)=pbu(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pbv(i,j,1)=min(pb(i,j,1),pb(i,j-1,1))
          pbv(i,j,2)=pbv(i,j,1)
          pbv_p(i,j)=pbv(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i,q)
      do j=-1,jj+2
        do k=1,kk
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,k)=
     .        .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .           +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,k)=
     .        .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .           +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            dpu(i,j,k+kk)=dpu(i,j,k)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            dpv(i,j,k+kk)=dpv(i,j,k)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Initialize potential vorticity of barotropic flow
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=0,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          q=2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,1)=corioq(i,j  )*q
          pvtrop(i,j+1,1)=corioq(i,j+1)*q
          pvtrop(i,j  ,2)=pvtrop(i,j  ,1)
          pvtrop(i,j+1,2)=pvtrop(i,j+1,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=1,jj
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
          q=2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,1)=corioq(i  ,j)*q
          pvtrop(i+1,j,1)=corioq(i+1,j)*q
          pvtrop(i  ,j,2)=pvtrop(i  ,j,1)
          pvtrop(i+1,j,2)=pvtrop(i+1,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isq(j)
        do i=max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,1)=corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  )
     .                                 +pb_p(i,j-1)+pb_p(i-1,j-1))
          pvtrop(i,j,2)=pvtrop(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Initialize fields related to the pressure gradient force
c --- ------------------------------------------------------------------
c
      call pgforc(2,1,kk,0,kk+1,1)
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            pgfx(i,j,k+kk)=pgfx(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            pgfy(i,j,k+kk)=pgfy(i,j,k)
          enddo
          enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pgfxm(i,j,2)=pgfxm(i,j,1)
          xixp(i,j,2)=xixp(i,j,1)
          xixm(i,j,2)=xixm(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pgfym(i,j,2)=pgfym(i,j,1)
          xiyp(i,j,2)=xiyp(i,j,1)
          xiym(i,j,2)=xiym(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Define first physical interior layer
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k,dps)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          k=3
          dps=0.
          do while (dp(i,j,k).lt.epsil)
            dps=dps+dp(i,j,k)
            dp(i,j,k)=0.
            dp(i,j,k+kk)=0.
            k=k+1
            if (k.gt.kk) exit
          enddo
          if (k.gt.kk) then
            dp(i,j,2)=dp(i,j,2)+dps
            dp(i,j,2+kk)=dp(i,j,2)
          else
            dp(i,j,k)=dp(i,j,k)+dps
            dp(i,j,k+kk)=dp(i,j,k)
          endif
          kfpla(i,j,1)=k
          kfpla(i,j,2)=k
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Set other time level layer thicknesses
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,k+kk)=dp(i,j,k)
            dpold(i,j,k)=dp(i,j,k)
            dpold(i,j,k+kk)=dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Initialize variables related to forcing
c --- ------------------------------------------------------------------
c
      if (sprfac) then
        prfac=1.
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            eiacc(i,j)=0.
            pracc(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          rnfres(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (mnproc.eq.ptest) then
        i=itest
        j=jtest
        write (lp,103) nstep,i0+i,j0+j,
     .    '  init.profile  temp    saln    dens   thkns    dpth',
     .    (k,temp(i,j,k),saln(i,j,k),1000.*sig(temp(i,j,k),saln(i,j,k)),
     .    dp(i,j,k)/onem,p(i,j,k+1)/onem,k=1,kk)
 103    format (i9,2i5,a/(28x,i3,3f8.2,2f8.1))
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inicon:'
        endif
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(pb,ip,3,'pb')
        call chksummsk(pbu,iu,2,'pbu')
        call chksummsk(pbv,iv,2,'pbv')
        call chksummsk(pvtrop,iq,2,'pvtrop')
        call chksummsk(pu,iu,kk+1,'pu')
        call chksummsk(pv,iv,kk+1,'pv')
      endif
c
      return
      end
