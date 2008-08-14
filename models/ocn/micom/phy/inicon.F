      subroutine inicon
c
c --- version 2.8 (NERSC version)
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#ifndef CCSMCOUPLED
#  include "common_forc.h"
#  if defined(CLIM) || defined(SYN)
#    include "common_asf.h"
#  endif
#  include "common_ice.h"
#endif
c
      real, dimension(itdm,jtdm) :: tmp2d
      integer i,j,k,l,kmin
      real dphi,alpl,alpu,q
c
      real swtfrz,sig,sofsig,getlpi
      external swtfrz,sig,sofsig,getlpi
c
c --- -------------------------
c --- mass field initialization
c --- -------------------------
c
      if (mnproc.eq.1) then
        write (lp,*) 'reading initial layer interfaces...'
        open (unit=nfu,file=path(1:path_len)//'sep_z.dat',status='old',
     .        action='read')
        read (nfu,'(10f8.3)') ((tmp2d(i,j),i=1,itdm),j=1,jtdm)
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          read (nfu,'(10f8.3)') ((tmp2d(i,j),i=1,itdm),j=1,jtdm)
        endif
        call xcaput(tmp2d,util1,1)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=min(depths(i,j),util1(i,j))*1.e2
          enddo
          enddo
        enddo
      enddo
      close (unit=nfu)
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
      if (mnproc.eq.1) then
        write (lp,*) 'reading initial temperature state...'
        open (unit=nfu,file=path(1:path_len)//'sep_temp.dat',
     .        status='old',action='read')
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          read (nfu,'(10f8.4)') ((tmp2d(i,j),i=1,itdm),j=1,jtdm)
        endif
        call xcaput(tmp2d,util1,1)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k)=util1(i,j)
          enddo
          enddo
        enddo
      enddo
      close (unit=nfu)
c
      if (mnproc.eq.1) then
        write (lp,*) 'reading initial salinity state...'
        open (unit=nfu,file=path(1:path_len)//'sep_saln.dat',
     .        status='old',action='read')
      endif
      if (mnproc.eq.1) then
        read (nfu,'(10f8.4)') ((tmp2d(i,j),i=1,itdm),j=1,jtdm)
      endif
      call xcaput(tmp2d,util1,1)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          saln(i,j,1)=util1(i,j)
        enddo
        enddo
      enddo
      close (unit=nfu)
c
#ifndef CCSMCOUPLED
c --- ------------------------------------------------------------------
c --- initialize variables for the thermodynamic ice model and make
c --- surface conditions consistent with the ice cover
c --- ------------------------------------------------------------------
c
      do 110 j=1,jj
      do 110 l=1,isp(j)
      do 110 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      sigma(i,j,1)=sig(temp(i,j,1),saln(i,j,1))
      if (ricclm(i,j,1).lt..5) then
      hicem(i,j)=0.
      ficem(i,j)=0.
      hsnwm(i,j)=0.
      tsrfm(i,j)=temp(i,j,1)+t0deg
      tsi_tda(i,j)=temp(i,j,1)+t0deg
      else
      hicem(i,j)=2.
      ficem(i,j)=min(fice_max,ricclm(i,j,1)) ! assuming start in January
      hsnwm(i,j)=0.1
      kmin=1
      do 120 k=kk,2,-1
        if (sigma(i,j,1)+sigjmp.lt.sigmar(k).and.
     .      p(i,j,k+1)-p(i,j,k).gt.100.) kmin=k
 120  continue
 130  temp(i,j,1)=swtfrz(saln(i,j,1))
      sigma(i,j,1)=sig(temp(i,j,1),saln(i,j,1))
      if (kmin.gt.1) then
        if (sigma(i,j,1)+sigjmp.gt.sigmar(kmin)) then
          sigma(i,j,1)=sigmar(kmin)-sigjmp
          saln(i,j,1)=sofsig(sigma(i,j,1),temp(i,j,1))
          goto 130
        endif
      endif
      tsrfm(i,j)=sstclm(i,j,1)
      tsi_tda(i,j)=(sstclm(i,j,1)-(1.-ficem(i,j))*(temp(i,j,1)+t0deg))
     .             /ficem(i,j)
      iagem(i,j)=0.
      endif
c
      ticem(i,j)=tsrfm(i,j)
      tml_tda(i,j)=temp(i,j,1)+t0deg
      sml_tda(i,j)=saln(i,j,1)
      alb_tda(i,j)=albs_f*ficem(i,j)+albw_d*(1.-ficem(i,j))
 110  fice_tda(i,j)=ficem(i,j)
c
      do 140 j=1,jj
      do 150 l=1,isu(j)
      do 150 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 150  ssu_tda(i,j)=0.
      do 140 l=1,isv(j)
      do 140 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 140  ssv_tda(i,j)=0.
c
      ntda=1
#endif
c
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 10 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      temp(i,j,1)=max(swtfrz(saln(i,j,1)),temp(i,j,1))
      tmxold(i,j,1)=temp(i,j,1)
      tmxold(i,j,2)=temp(i,j,1)
      smxold(i,j,1)=saln(i,j,1)
      smxold(i,j,2)=saln(i,j,1)
      temp(i,j,1+kk)=temp(i,j,1)
      saln(i,j,1+kk)=saln(i,j,1)
      sigma(i,j,1)=sig(temp(i,j,1),saln(i,j,1))
      sigma(i,j,1+kk)=sigma(i,j,1)
      do 10 k=2,kk
      temp(i,j,k)=max(temmin(k),temp(i,j,k))
      saln(i,j,k   )=sofsig(sigmar(k),temp(i,j,k))
      temp(i,j,k+kk)=temp(i,j,k)
      saln(i,j,k+kk)=saln(i,j,k)
      sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
      sigma(i,j,k+kk)=sigma(i,j,k)
 10   continue
c
      if (mnproc.eq.1) then
        write (lp,'('' sigmar(k)    :'',9f7.5/(15x,9f7.5))')
     .     (sigmar(k),k=1,kk)
        write (lp,'('' alphar(k)    :'',9f7.5/(15x,9f7.5))')
     .     (alphar(k),k=1,kk)
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
c$OMP PARALLEL DO
      do j=-1,jj+2
        do k=1,kk
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            dpu(i,j,k)=
     .        max(0.,min(pbu(i,j,1),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .              -min(pbu(i,j,1),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            dpv(i,j,k)=
     .        max(0.,min(pbv(i,j,1),.5*(p(i,j  ,k+1)+p(i,j-1,k+1)))
     .              -min(pbv(i,j,1),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
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
          pvtrop(i,j  ,1)=corio(i,j  )*q
          pvtrop(i,j+1,1)=corio(i,j+1)*q
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
          pvtrop(i  ,j,1)=corio(i  ,j)*q
          pvtrop(i+1,j,1)=corio(i+1,j)*q
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
          pvtrop(i,j,1)=corio(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  )
     .                                +pb_p(i,j-1)+pb_p(i-1,j-1))
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
c
      enddo
c
c --- ------------------------------------------------------------------
c --- initialize variables related to forcing
c --- ------------------------------------------------------------------
c
#if !defined(CCSMCOUPLED) && (defined(CLIM) || defined(SYN))
      do 300 j=1,jj
      do 300 l=1,isp(j)
      do 300 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      cd_d(i,j)=cd_r
      ch_d(i,j)=ch_r
      ce_d(i,j)=ce_r
      wg2_d(i,j)=wg2_r*1.e-4
      cd_m(i,j)=cd_r
      ch_m(i,j)=ch_r
      ce_m(i,j)=ce_r
      wg2_m(i,j)=wg2_r*1.e-4
      rhoa(i,j)=rhoa_r*1.e3
 300  rnfres(i,j)=0.
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
