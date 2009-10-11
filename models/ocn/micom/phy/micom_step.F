      subroutine micom_step
c
c --- ------------------------------------------------------------------
c --- integrate a model time step
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
#include "common_forc.h"
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
#endif
#include "param_bud.h"
#include "common_bud.h"
#include "timing.h"
c
      real q
      integer i,j,k,l,m,n,mm,nn,k1m,k1n
c
      real step_time,
     .     auxil_time ,
     .     getflx_time,
     .     seaice_time,
     .     tmsmt1_time,
     .     advdif_time,
     .     sfcstr_time,
     .     momtum_time,
     .     pgforc_time,
     .     barotp_time,
     .     pbtcor_time,
     .     convec_time,
     .     diapfl_time,
     .     thermf_time,
     .     mxlayr_time,
     .     tmsmt2_time,
     .     diaacc_time,
     .     io_time
c
      real gettime
      external gettime
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
c
      nstep=nstep+1
      time=nstep*baclin/86400.
c
#ifdef TRC
c --- advance tracers only every second time step
      ntrc=ntr*mod(nstep+1,2)
#  ifdef ATRC
      natrc=natr*mod(nstep+1,2)
#  endif
#endif
c
      if (cnsvdi) call stdbud(1,n,nn,k1n)
c
c --- Update calender
      call updcln(nstep)
c
c --- ------------------------------------------------------------------
c --- Set parameters for time interpolation of climatological fields
c --- ------------------------------------------------------------------
      x=(real(nday_of_year-1)+mod(time,1.))*12./real(nday_in_year)
      l3=int(x)+1
      x=x-real(l3-1)
      l1=mod(l3+ 9,12)+1
      l2=mod(l3+10,12)+1
      l4=mod(l3   ,12)+1
      l5=mod(l3+ 1,12)+1
c
      auxil_time=gettime()
c
#ifdef NEST
c --- ------------------------------------------------------------------
c --- Get boundary conditions from outer coarse model
c --- ------------------------------------------------------------------
c
      call rdotr
#endif
c
#ifndef CCSMCOUPLED
c --- ------------------------------------------------------------------
c --- Get surface fluxes
c --- ------------------------------------------------------------------
c
      if (mod(nstep-1,nfrco).eq.0) then
#  if   defined(CLIM)
      if (mnproc.eq.1) then
        write (lp,*) 'getflux...'
      endif
      call getflux_clim
#  elif defined(SYN)
      call getflux_syn
#  endif
      endif
#endif
c
      getflx_time=gettime()
c
#ifdef ICEDYN
cdiag write (lp,*) 'icedyn...'
      call icedyn(m,n,mm,nn,k1m,k1n)
#endif
      seaice_time=gettime()
c
cdiag write (lp,*) 'tmsmt1...'
      call tmsmt1(m,n,mm,nn,k1m,k1n)
      tmsmt1_time=gettime()
c
cdiag write (lp,*) 'advdif...'
      call difest(m,n,mm,nn,k1m,k1n)
      call intdif(m,n,mm,nn,k1m,k1n)
      call advect(m,n,mm,nn,k1m,k1n)
      call diffus(m,n,mm,nn,k1m,k1n)
      advdif_time=gettime()
c
      if (cnsvdi) call stdbud(2,n,nn,k1n)
      auxil_time=auxil_time+gettime()
c
cdiag write (lp,*) 'sfcstr...'
#ifdef CCSMCOUPLED
      call sfcstr_ccsm(m,n,mm,nn,k1m,k1n)
#else
      call sfcstr(m,n,mm,nn,k1m,k1n)
#endif
      sfcstr_time=gettime()
c
cdiag write (lp,*) 'pgforc...'
      call pgforc(m,n,mm,nn,k1m,k1n)
      pgforc_time=gettime()
c
cdiag write (lp,*) 'momtum...'
      call momtum(m,n,mm,nn,k1m,k1n)
      momtum_time=gettime()
c
#ifdef NEST
cdiag write (lp,*) 'bclnst...'
      call bclnst(m,n,mm,nn,k1m,k1n)
c     bclnst_time=gettime()
#endif
c
cdiag write (lp,*) 'convec...'
      call convec(m,n,mm,nn,k1m,k1n)
      convec_time=gettime()
c
      if (cnsvdi) call stdbud(3,n,nn,k1n)
      auxil_time=auxil_time+gettime()
c
cdiag write (lp,*) 'diapfl...'
      call diapfl(m,n,mm,nn,k1m,k1n)
      diapfl_time=gettime()
c
      if (cnsvdi) call stdbud(4,n,nn,k1n)
      auxil_time=auxil_time+gettime()
c
cdiag write (lp,*) 'thermf...'
#ifdef CCSMCOUPLED
      call thermf_ccsm(m,n,mm,nn,k1m,k1n)
#else
      call thermf(m,n,mm,nn,k1m,k1n)
#endif
      thermf_time=gettime()
c
cdiag write (lp,*) 'mxlayr...'
      call mxlayr(m,n,mm,nn,k1m,k1n)
      mxlayr_time=gettime()
c
#ifdef TRC
c --- update tracer due to non-passive processes
      if (mod(nstep+1,2).eq.1) call updtrc(m,n,mm,nn,k1m,k1n)
#endif
c
      if (cnsvdi) call stdbud(5,n,nn,k1n)
      auxil_time=auxil_time+gettime()
c
cdiag write (lp,*) 'barotp...'
      call barotp(m,n,mm,nn,k1m,k1n)
      barotp_time=gettime()
c
cdiag write (lp,*) 'pbtcor...'
      call pbtcor(m,n,mm,nn,k1m,k1n)
      pbtcor_time=gettime()
c
      if (cnsvdi) call stdbud(6,m,mm,k1m)
      auxil_time=auxil_time+gettime()
c
cdiag write (lp,*) 'tmsmt2...'
      call tmsmt2(m,n,mm,nn,k1m,k1n)
      tmsmt2_time=gettime()
c
      if (cnsvdi) call stdbud(7,m,mm,k1m)
c
      if (cnsvdi) then
        if (mnproc.eq.1.and.nstep.gt.nstep1+1) then
          open (unit=nfu,file='salbud',position='append')
          write (nfu,'(i8,6e12.4)') nstep-1,
     .      (sdp(2,m)-sdp(1,m))/vol0,
     .      (sdp(3,m)-sdp(2,m))/vol0,
     .      (sdp(4,m)-sdp(3,m))/vol0,
     .      (sdp(5,m)-sdp(4,m)+sf*g)/vol0,
     .      (sdp(6,m)-sdp(5,m))/vol0,
     .      (sdp(7,m)-sdp(6,m))/vol0
          close (nfu)
          open (unit=nfu,file='tembud',position='append')
          write (nfu,'(i8,6e12.4)') nstep-1,
     .      (tdp(2,m)-tdp(1,m))/vol0,
     .      (tdp(3,m)-tdp(2,m))/vol0,
     .      (tdp(4,m)-tdp(3,m))/vol0,
     .      (tdp(5,m)-tdp(4,m)+tf*g/spcifh)/vol0,
     .      (tdp(6,m)-tdp(5,m))/vol0,
     .      (tdp(7,m)-tdp(6,m))/vol0
          close (nfu)
          open (unit=nfu,file='denbud',position='append')
          write (nfu,'(i8,6e12.4)') nstep-1,
     .      (ddp(2,m)-ddp(1,m)),(ddp(3,m)-ddp(2,m)),(ddp(4,m)-ddp(3,m)),
     .      (ddp(5,m)-ddp(4,m)),(ddp(6,m)-ddp(5,m)),(ddp(7,m)-ddp(6,m))
          close (nfu)
#ifdef TRC
          if (m.eq.1) then
            open (unit=nfu,file='trcbud',position='append')
            write (nfu,'(i8,6e12.4)') nstep,
     .        (trdp(2)-trdp(1))/vol0,
     .        (trdp(3)-trdp(2))/vol0,
     .        (trdp(4)-trdp(3))/vol0,
     .        (trdp(5)-trdp(4)+trf*g)/vol0,
     .        (trdp(6)-trdp(5))/vol0,
     .        (trdp(7)-trdp(6))/vol0
            close (nfu)
            open (unit=nfu,file='trcbudtot',position='append')
            write (nfu,'(i8,7e18.10)') nstep,
     .        trdp(1)/vol0,trdp(2)/vol0,trdp(3)/vol0,
     .        trdp(4)/vol0,trdp(5)/vol0,trdp(6)/vol0,trdp(7)/vol0
            close (nfu)
          endif
#endif
        endif
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=(salflx(i,j)+salrlx(i,j))*scp2(i,j)*delt1
            util2(i,j)=(surflx(i,j)+surrlx(i,j))*scp2(i,j)*delt1
#ifdef TRC
            util3(i,j)=trflx(1,i,j)*scp2(i,j)*delt1
#endif
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
#ifdef ARCTIC
        do j=jtdm-j0,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=0.
            util2(i,j)=0.
#ifdef TRC
            util3(i,j)=0.
#endif
          enddo
          enddo
        enddo
#endif
        call xcsum(sf,util1,ip)
        call xcsum(tf,util2,ip)
#ifdef TRC
        call xcsum(trf,util3,ip)
#endif
      endif
c
      auxil_time=auxil_time+gettime()
c
c ----------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ----------------------------------------------------------------------
c
      call chkvar(m,n,mm,nn,k1m,k1n)
c
      call diaacc(m,n,mm,nn,k1m,k1n)
      diaacc_time=gettime()
c
      if (mod(nstep,nstep_in_day).eq.0.and.nday_of_year.eq.1) then
c
c --- ------------------------------------------------------------------
c --- - output diagnosed heat and salt flux
c --- ------------------------------------------------------------------
c
        call wdiflx
c
      endif
c
      if (diag2dmon.and.(mod(nstep,nstep_in_day).eq.0.and.nday.eq.1).or.
     .    .not.diag2dmon.and.mod(nstep+.5,diag2dfq).lt.1.) then
c
c --- ------------------------------------------------------------------
c --- - output 2D diagnostics
c --- ------------------------------------------------------------------
c
        call ncwrt2
c
      endif
c
      if (diag3dmon.and.(mod(nstep,nstep_in_day).eq.0.and.nday.eq.1).or.
     .    .not.diag3dmon.and.mod(nstep+.5,diag3dfq).lt.1.) then
c
c --- ------------------------------------------------------------------
c --- - output 3D diagnostics
c --- ------------------------------------------------------------------
c
        call ncwrt3
c
      endif
c
c --- update total time spent by various tasks
      auxil_total_time=auxil_total_time+auxil_time
      getflx_total_time=getflx_total_time+getflx_time
      seaice_total_time=seaice_total_time+seaice_time
      tmsmt1_total_time=tmsmt1_total_time+tmsmt1_time
      advdif_total_time=advdif_total_time+advdif_time
      sfcstr_total_time=sfcstr_total_time+sfcstr_time
      momtum_total_time=momtum_total_time+momtum_time
      pgforc_total_time=pgforc_total_time+pgforc_time
      barotp_total_time=barotp_total_time+barotp_time
      pbtcor_total_time=pbtcor_total_time+pbtcor_time
      convec_total_time=convec_total_time+convec_time
      diapfl_total_time=diapfl_total_time+diapfl_time
      thermf_total_time=thermf_total_time+thermf_time
      mxlayr_total_time=mxlayr_total_time+mxlayr_time
      tmsmt2_total_time=tmsmt2_total_time+tmsmt2_time
      diaacc_total_time=diaacc_total_time+diaacc_time
c
      if (diag3dmon.and.(mod(nstep,nstep_in_day).eq.0.and.nday.eq.1).or.
     .    .not.diag3dmon.and.
     .    mod(nstep,max(nint(diag3dfq),nstep_in_day)).eq.0) then
c
#ifndef CCSMCOUPLED
c --- ------------------------------------------------------------------
c --- - output restart files
c --- ------------------------------------------------------------------
c
        call restart_wt
c
#  ifdef ARCTIC
c --- - this is needed to reproduce results after restart with arctic
c --- - patch
        call xctilr(uflx(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_uv)
        call xctilr(vflx(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_vv)
#  endif
#endif
c
        io_time=gettime()
c
c --- ------------------------------------------------------------------
c --- - write timing diagnostics to standard out
c --- ------------------------------------------------------------------
c
        io_total_time=io_total_time+io_time
        step_time=auxil_time +getflx_time+seaice_time+tmsmt1_time
     .           +advdif_time+sfcstr_time+momtum_time+pgforc_time
     .           +barotp_time+pbtcor_time+convec_time+diapfl_time
     .           +thermf_time+mxlayr_time+tmsmt2_time+diaacc_time
     .           +io_time
        total_time=total_time+step_time
        total_xio_time=total_xio_time+step_time-io_time
c
        if (mnproc.eq.1) then
          write (lp,'(f12.4,a,i8)')
     .      step_time, '  sec for step ', nstep
          write (lp,'(f12.4,a,i8)')
     .      total_time/(nstep-nstep1),' Avg Time'
          write (lp,'(f12.4,a,i8)')
     .      total_xio_time/(nstep-nstep1),' Avg Time excluding IO'
          write (lp,'(f12.4,a,i8)')
     .      total_time,' Tot Time with contributions:'
          q=100./total_time
          write (lp,'(f12.4,a,i8)') auxil_total_time*q ,'% auxil '
          write (lp,'(f12.4,a,i8)') getflx_total_time*q,'% getflx'
          write (lp,'(f12.4,a,i8)') seaice_total_time*q,'% seaice'
          write (lp,'(f12.4,a,i8)') tmsmt1_total_time*q,'% tmsmt1'
          write (lp,'(f12.4,a,i8)') advdif_total_time*q,'% advdif'
          write (lp,'(f12.4,a,i8)') sfcstr_total_time*q,'% sfcstr'
          write (lp,'(f12.4,a,i8)') momtum_total_time*q,'% momtum'
          write (lp,'(f12.4,a,i8)') pgforc_total_time*q,'% pgforc'
          write (lp,'(f12.4,a,i8)') barotp_total_time*q,'% barotp'
          write (lp,'(f12.4,a,i8)') pbtcor_total_time*q,'% pbtcor'
          write (lp,'(f12.4,a,i8)') convec_total_time*q,'% convec'
          write (lp,'(f12.4,a,i8)') diapfl_total_time*q,'% diapfl'
          write (lp,'(f12.4,a,i8)') thermf_total_time*q,'% thermf'
          write (lp,'(f12.4,a,i8)') mxlayr_total_time*q,'% mxlayr'
          write (lp,'(f12.4,a,i8)') tmsmt2_total_time*q,'% tmsmt2'
          write (lp,'(f12.4,a,i8)') diaacc_total_time*q,'% diaacc'
          write (lp,'(f12.4,a,i8)') io_total_time*q    ,'% IO'
        endif
c
      else
c
c --- ------------------------------------------------------------------
c --- - write time spent for current time step
c --- ------------------------------------------------------------------
c
        io_time=gettime()
        io_total_time=io_total_time+io_time
        step_time=auxil_time +getflx_time+seaice_time+tmsmt1_time
     .           +advdif_time+sfcstr_time+momtum_time+pgforc_time
     .           +barotp_time+pbtcor_time+convec_time+diapfl_time
     .           +thermf_time+mxlayr_time+tmsmt2_time+diaacc_time
     .           +io_time
        total_time=total_time+step_time
        total_xio_time=total_xio_time+step_time-io_time
c
        if (mnproc.eq.1) then
          write (lp,'(f12.4,a,i8)') step_time, '  sec for step ', nstep
        endif
c
      endif
c
      delt1=baclin+baclin
c
      return
      end
