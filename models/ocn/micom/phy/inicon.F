      subroutine inicon
c
c --- version 2.8 (NERSC version)
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#ifndef CCSMCOUPLED
#  if defined(CLIM) || defined(SYN)
#    include "common_asf.h"
#  endif
#  include "common_ice.h"
#endif
#include "netcdf.inc"
c
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: start,count
      integer i,j,k,l,status,ncid,dimid,varid
      real q,dps
c
      real swtfrz,sig,sofsig,getlpi
      external swtfrz,sig,sofsig,getlpi
c
c --- ------------------------------------------------------------------
c --- mass field initialization
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading initial condition from ',
     .                    path(1:path_len)//'inicon.nc'
        call flush(lp)
c
c --- - open netcdf file
        status=nf_open(path(1:path_len)//'inicon.nc',nf_nowrite,
     .                 ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open: inicon.nc: ',
     .                     nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
c
c --- - check dimensions
        status=nf_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: x: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
        status=nf_inq_dimlen(ncid,dimid,i)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: x: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
        status=nf_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: y: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
        status=nf_inq_dimlen(ncid,dimid,j)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: y: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in inicon.nc'
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
c
      endif
c
      start(1)=1
      start(2)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
c
c --- read reference potential density
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'sigma',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: sigma: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: sigma: ',
     .                       nf_strerror(status)
            call xchalt('(inicon)')
                   stop '(inicon)'
          endif
        endif
        call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read potential temperature
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'temp',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: temp: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: temp: ',
     .                       nf_strerror(status)
            call xchalt('(inicon)')
                   stop '(inicon)'
          endif
        endif
        call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read salinity
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'saln',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: saln: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: saln: ',
     .                       nf_strerror(status)
            call xchalt('(inicon)')
                   stop '(inicon)'
          endif
        endif
        call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read layer thickness
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'dz',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: dz: ',nf_strerror(status)
          call xchalt('(inicon)')
                 stop '(inicon)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmp2d)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: dz: ',
     .                       nf_strerror(status)
            call xchalt('(inicon)')
                   stop '(inicon)'
          endif
        endif
        call xcaput(tmp2d,dp(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- convert unit of reference potential density from [kg/m^3] to
c --- [g/cm^3] and construct interface depths from layer thicknesses
c
      do k=1,kk
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=min(depths(i,j)*1.e2,p(i,j,k)+dp(i,j,k)*1.e2)
            sigmar(i,j,k)=sigmar(i,j,k)*1.e-3
          enddo
          enddo
        enddo
      enddo
c
      do j=1,jj
        do k=2,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (p(i,j,kk+1)-p(i,j,k).lt.1.e-4)
     .        p(i,j,k)=depths(i,j)*1.e2
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          p(i,j,kk+1)=depths(i,j)*1.e2
        enddo
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- set minimum physical temperature for each isopycnic layer
c --- ------------------------------------------------------------------
c
      call settemmin
c
#ifndef CCSMCOUPLED
c --- ------------------------------------------------------------------
c --- initialize variables for the thermodynamic ice model and make
c --- surface conditions consistent with the ice cover
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (ricclm(i,j,12).lt..05) then
            hicem(i,j)=0.
            ficem(i,j)=0.
            hsnwm(i,j)=0.
            tsrfm(i,j)=temp(i,j,1)+t0deg
            tsi_tda(i,j)=temp(i,j,1)+t0deg
          else
            hicem(i,j)=2.
            ficem(i,j)=min(fice_max,ricclm(i,j,12)) ! assuming start in January
            hsnwm(i,j)=0.1
            tsrfm(i,j)=sstclm(i,j,1)
            tsi_tda(i,j)=(sstclm(i,j,1)
     .                   -(1.-ficem(i,j))*(temp(i,j,1)+t0deg))
     .                   /ficem(i,j)
            iagem(i,j)=0.
          endif
c
          ticem(i,j)=tsrfm(i,j)
          tml_tda(i,j)=temp(i,j,1)+t0deg
          sml_tda(i,j)=saln(i,j,1)
          alb_tda(i,j)=albs_f*ficem(i,j)+albw_d*(1.-ficem(i,j))
          fice_tda(i,j)=ficem(i,j)
        enddo
        enddo
      enddo
c
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ssu_tda(i,j)=0.
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          ssv_tda(i,j)=0.
        enddo
        enddo
      enddo
c
      ntda=1
c
#else
c --- ------------------------------------------------------------------
c --- Initialize variables used when coupled to CCSM
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          frzpot(i,j)=0.
          mltpot(i,j)=0.
          flxco2(i,j)=0.
        enddo
        enddo
      enddo
#endif
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          do k=1,2
            temp(i,j,k)=max(swtfrz(saln(i,j,k)),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          do k=3,kk
            temp(i,j,k)=max(temmin(i,j,k),temp(i,j,k))
            saln(i,j,k   )=sofsig(sigmar(i,j,k),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          ustar(i,j)=0.
          buoyfl(i,j)=0.
        enddo
        enddo
      enddo
c
      if (mnproc.eq.ptest) then
        write (lp,'('' sigmar(k)    :'',9f7.5/(15x,9f7.5))')
     .     (sigmar(itest,jtest,k),k=1,kk)
      endif
c
c --- ------------------------------------------------------------------
c --- represent layer interfaces in units of pressure by assuming an
c --- initially flat ocean surface and set reference geopotential at the
c --- ocean floor
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          phi(i,j,1)=0.
        enddo
        enddo
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ----- geopotential at interface
            phi(i,j,k+1)=-g*p(i,j,k+1)
c
            p(i,j,k+1)=getlpi(temp(i,j,k),saln(i,j,k),
     .                        phi(i,j,k),phi(i,j,k+1),p(i,j,k))
c
          enddo
          enddo
        enddo
      enddo
c
      call xctilr(p, 1,kk+1, 2,2, halo_ps)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)
      call xctilr(temp, 1,kk, 1,1, halo_ps)
      call xctilr(saln, 1,kk, 1,1, halo_ps)
c
c --- ------------------------------------------------------------------
c --- set layer thickness and bottom pressure
c --- ------------------------------------------------------------------
c
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
c
      do j=0,jj+1
        do k=1,kk
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c
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
c
      do  j=1,jj
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
c
c$OMP PARALLEL DO PRIVATE(q)
      do j=-1,jj+2
        do k=1,kk
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,k)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i-1,j,k+1))) 
     .              -min(q,.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,k)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i,j-1,k+1)))
     .              -min(q,.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
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
c --- initialize potential vorticity of barotropic flow
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(q)
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
c$OMP PARALLEL DO PRIVATE(q)
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
c$OMP PARALLEL DO
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
c --- initialize fields related to the pressure gradient force
c --- ------------------------------------------------------------------
c
      call pgforc(2,1,kk,0,kk+1,1)
c
c$OMP PARALLEL DO
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
c --- define first physical interior layer
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,dps)
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
c --- set other time level layer thicknesses
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
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
c --- initialize variables related to forcing
c --- ------------------------------------------------------------------
c
      if (sprfac) then
        prfac=1.
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            eiacc(i,j)=0.
            pracc(i,j)=0.
          enddo
          enddo
        enddo
      endif
c
#if !defined(CCSMCOUPLED) && (defined(CLIM) || defined(SYN))
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cd_d(i,j)=cd_r
          ch_d(i,j)=ch_r
          ce_d(i,j)=ce_r
          wg2_d(i,j)=wg2_r*1.e-4
          cd_m(i,j)=cd_r
          ch_m(i,j)=ch_r
          ce_m(i,j)=ce_r
          wg2_m(i,j)=wg2_r*1.e-4
          rhoa(i,j)=rhoa_r*1.e3
          rnfres(i,j)=0.
        enddo
        enddo
      enddo
#endif
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
c
c
c> Revision history:
c>
c> May  1997 - eliminated "+ thbase" from definition of thmix (loop 10)
c> Mar. 1998 - added -th3d-
c> Aug. 2000 - added diagnostic count of static instabilities
