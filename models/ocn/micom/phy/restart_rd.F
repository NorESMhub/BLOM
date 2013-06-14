      subroutine restart_rd(nday1)
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia
      use mod_nctools
c
      implicit none
c
      integer nday1
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#if defined(CLIM) || defined(SYN)
#  include "common_asf.h"
#endif
#ifdef CCSMCOUPLED
#  include "common_ccsm.h"
#endif
#include "common_ice.h"
c
      real rtmp
      integer i,j,k,n,ndayr,nmonthr,nyearr,l
      character(len=80) :: rstfnm,fnm
      character(len=2) :: c2
c
      integer daydif
      external daydif
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv,iqq
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
      logical first,fexist
      real rfexist
      data first /.true./
      save first,iuu,ivv,iqq
c
#if defined TKE
#  include "common_tke.h"
#endif
#if defined(TRC) && defined(IDLAGE)
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
      if (mnproc.eq.1) then
c
#ifdef CCSMCOUPLED
c --- - get file name from rpointer.ocn
        inquire(file='rpointer.ocn',exist=fexist)
        if (fexist) then
          open (unit=nfu,file='rpointer.ocn')
          read (nfu,'(a)') rstfnm
          close (unit=nfu)
        else
          write (lp,*) 'Could not find file rpointer.ocn!'
          call xchalt('(restart_rd)')
          stop '(restart_rd)' 
        endif
#else
c --- - first try file name:
c --- -  <experiment name>_restphy_<year>.<month>.<day>_<integration day>.nc
        write (rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')
     .    runid(1:runid_len),'_restphy_',nyear,'.',nmonth,'.',nday,'_',
     .    nday1,'.nc'
#endif
        inquire(file=path2(1:path2_len)//rstfnm,exist=fexist)
        if (fexist) then
          call ncfopn(path2(1:path2_len)//rstfnm,'r',' ',1)
          call ncgeti('nday0',ndayr)
          call ncgeti('nmonth0',nmonthr)
          call ncgeti('nyear0',nyearr)
          call ncgetr('time0',time0)
          call ncgetr('time',time)
c
          i=daydif(nyear0,nmonth0,nday0,nyear,nmonth,nday)
          if (nyearr.ne.nyear0.or.nmonthr.ne.nmonth0.or.
     .        ndayr.ne.nday0) then
            time0=nday1-i
            write (lp,'(a,i4.4,a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2)')
     .        ' Warning! Initial experiment date changed from ',
     .        nyearr,'.',nmonthr,'.',ndayr,' to ',
     .        nyear0,'.',nmonth0,'.',nday0
          elseif (nint(time0).ne.nday1-i) then
            write (lp,'(2a)')
     .       ' Warning! Integration day corresponding to initial',
     .       ' experiment'
            write (lp,'(a,i6,a,i6)')
     .        ' date changed from ',nint(time0),' to ',nday1-i
            time0=nday1-i
          endif
        else
c
c --- - then try automatic selection of file with consistent integration
c --- - day and date among files named:
c --- -   <experiment name>_restphy_1.nc
c --- -   <experiment name>_restphy_2.nc
c --- -   <experiment name>_restphy_3.nc
          do i=1,4
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',i,'.nc'
            inquire(file=path2(1:path2_len)//rstfnm,exist=fexist)
            if (fexist) then
              call ncfopn(path2(1:path2_len)//rstfnm,'r',' ',1)
              call ncgeti('nday0',ndayr)
              call ncgeti('nmonth0',nmonthr)
              call ncgeti('nyear0',nyearr)
              call ncgetr('time0',time0)
              call ncgetr('time',time)
              if (nint(time).eq.nday1.and.
     .          daydif(nyearr,nmonthr,ndayr,nyear,nmonth,nday).eq.
     .          nint(time-time0)) exit 
            endif 
          enddo 
          if (i.gt.3) then 
            write (lp,*) 'Could not find proper restart file!'
            call xchalt('(restart_rd)')
            stop '(restart_rd)' 
          endif 
        endif 
c
        write (lp,'(2a)') ' reading restart file ',trim(rstfnm)
c
      endif
c
c --- Compute extended uv masks 
      if (first) then
        first=.false.
        do j=1,jj
          do i=1,ii
            if ((ip(i,j)+ip(i-1,j)).ge.1) then
              iuu(i,j)=1
            else
              iuu(i,j)=0
            endif
            if ((ip(i,j)+ip(i,j-1)).ge.1) then
              ivv(i,j)=1
            else
              ivv(i,j)=0
            endif
            if ((iu(i,j)+iv(i,j)+iu(i,j-1)+iv(i-1,j)).ge.1) then
              iqq(i,j)=1
            else
              iqq(i,j)=0
            endif
          enddo 
        enddo
      endif
c
      call ncread('u',u,iuu,1,0.)
      call ncread('v',v,ivv,1,0.)
      call ncread('dp',dp,ip,1,0.)
      call ncread('temp',temp,ip,1,0.)
      call ncread('saln',saln,ip,1,0.)
      call ncread('sigma',sigma,ip,1,0.)
      call ncread('sigmar',sigmar,ip,1,0.)
      call ncread('pgfx',pgfx,iuu,1,0.)
      call ncread('pgfy',pgfy,ivv,1,0.)
      call ncread('pb',pb,ip,1,0.)
      call ncread('pb_mn',pb_mn,ip,1,0.)
      call ncread('pb_p',pb_p,ip,1,0.)
      call ncread('pbu',pbu,iuu,1,0.)
      call ncread('pbv',pbv,ivv,1,0.)
      call ncread('pbu_p',pbu_p,iuu,1,0.)
      call ncread('pbv_p',pbv_p,ivv,1,0.)
      call ncread('ub',ub,iuu,1,0.)
      call ncread('vb',vb,ivv,1,0.)
      call ncread('uflx',uflx,iuu,1,0.)
      call ncread('utflx',utflx,iuu,1,0.)
      call ncread('usflx',usflx,iuu,1,0.)
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('umfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('umfltd',umfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''umfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (iuu(i,j).eq.1) then
                umfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('utfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('utfltd',utfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''utfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (iuu(i,j).eq.1) then
                utfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('utflld')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('utflld',utflld,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''utflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (iuu(i,j).eq.1) then
                utflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('usfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('usfltd',usfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''usfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (iuu(i,j).eq.1) then
                usfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('usflld')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('usflld',usflld,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''usflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (iuu(i,j).eq.1) then
                usflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      call ncread('vflx',vflx,ivv,1,0.)
      call ncread('vtflx',vtflx,ivv,1,0.)
      call ncread('vsflx',vsflx,ivv,1,0.)
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('vmfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('vmfltd',vmfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vmfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vmfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('vtfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('vtfltd',vtfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vtfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vtfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('vtflld')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('vtflld',vtflld,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vtflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vtflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('vsfltd')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('vsfltd',vsfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vsfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vsfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('vsflld')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('vsflld',vsflld,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vsflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,2*kk
          do j=1,jj
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vsflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      endif
c
      call ncread('ubflx',ubflx,iuu,1,0.)
      call ncread('vbflx',vbflx,ivv,1,0.)
      call ncread('ubflx_mn',ubflx_mn,iuu,1,0.)
      call ncread('vbflx_mn',vbflx_mn,ivv,1,0.)
      call ncread('ubflxs',ubflxs,iuu,1,0.)
      call ncread('vbflxs',vbflxs,ivv,1,0.)
      call ncread('ubflxs_p',ubflxs_p,iuu,1,0.)
      call ncread('vbflxs_p',vbflxs_p,ivv,1,0.)
      call ncread('ubcors_p',ubcors_p,iuu,1,0.)
      call ncread('vbcors_p',vbcors_p,ivv,1,0.)
      call ncread('pvtrop',pvtrop,iqq,1,0.)
      call ncread('pgfxm',pgfxm,iuu,1,0.)
      call ncread('pgfym',pgfym,ivv,1,0.)
      call ncread('xixp',xixp,iuu,1,0.)
      call ncread('xixm',xixm,iuu,1,0.)
      call ncread('xiyp',xiyp,ivv,1,0.)
      call ncread('xiym',xiym,ivv,1,0.)
      call ncread('phi',phi(1-nbdy,1-nbdy,kk+1),ip,1,0.)
      call ncread('sealv',sealv,ip,1,0.)
      call ncread('ustar',ustar,ip,1,0.)
      call ncread('buoyfl',buoyfl,ip,1,0.)
      call ncread('kfpla',rkfpla,ip,1,0.)
      call ncread('hicem',hicem,ip,1,0.)
      call ncread('ficem',ficem,ip,1,0.)
      call ncread('tsrfm',tsrfm,ip,1,0.)
      call ncread('hsnwm',hsnwm,ip,1,0.)
      call ncread('ticem',ticem,ip,1,0.)
      call ncread('iagem',iagem,ip,1,0.)
      call ncread('tsi_tda',tsi_tda,ip,1,0.)
      call ncread('tml_tda',tml_tda,ip,1,0.)
      call ncread('sml_tda',sml_tda,ip,1,0.)
      call ncread('alb_tda',alb_tda,ip,1,0.)
      call ncread('fice_tda',fice_tda,ip,1,0.)
c
      rtmp=0.
      if (mnproc.eq.1) then
        call ncgeti('ntda',ntda)
        rtmp=ntda
      endif
      call xcmaxr(rtmp)
      ntda=nint(rtmp)
c
      call ncread('rnfres',rnfres,ip,1,0.)
c
      if (sprfac) then
        prfac=-1.
        if (mnproc.eq.1) then
          if (ncinqa('prfac')) then
            call ncgetr('prfac',prfac)
          endif
        endif
        call xcmaxr(prfac)
        if (prfac.gt.0.) then
          call ncread('eiacc',eiacc,ip,1,0.)
          call ncread('pracc',pracc,ip,1,0.)
        else
          if (mnproc.eq.1) then
            write (lp,*)
     .        'Warning: fields needed for balancing fresh water '//
     .        'budget are not read'
            write (lp,*) 'from restart file and will be initialized.'
          endif
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
      endif
c
#if defined(CLIM) || defined(SYN)
      call ncread('cd_d',cd_d,ip,1,0.)
      call ncread('ch_d',ch_d,ip,1,0.)
      call ncread('ce_d',ce_d,ip,1,0.)
      call ncread('wg2_d',wg2_d,ip,1,0.)
      call ncread('cd_m',cd_m,ip,1,0.)
      call ncread('ch_m',ch_m,ip,1,0.)
      call ncread('ce_m',ce_m,ip,1,0.)
      call ncread('wg2_m',wg2_m,ip,1,0.)
      call ncread('rhoa',rhoa,ip,1,0.)
#endif
c
#ifdef CCSMCOUPLED
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('ustarw_da')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('ustarw_da',ustarw_da,ip,1,0.)
        call ncread('ztx_da',ztx_da,ip,1,0.)
        call ncread('mty_da',mty_da,ip,1,0.)
        call ncread('lip_da',lip_da,ip,1,0.)
        call ncread('sop_da',sop_da,ip,1,0.)
        call ncread('eva_da',eva_da,ip,1,0.)
        call ncread('rnf_da',rnf_da,ip,1,0.)
        call ncread('rfi_da',rfi_da,ip,1,0.)
        call ncread('fmltfz_da',fmltfz_da,ip,1,0.)
        call ncread('sfl_da',sfl_da,ip,1,0.)
        call ncread('swa_da',swa_da,ip,1,0.)
        call ncread('nsf_da',nsf_da,ip,1,0.)
        call ncread('hmlt_da',hmlt_da,ip,1,0.)
        call ncread('slp_da',slp_da,ip,1,0.)
        call ncread('ficem_da',ficem_da,ip,1,0.)
        call ncread('abswnd_da',abswnd_da,ip,1,0.)
        call ncread('atmco2_da',atmco2_da,ip,1,0.)
        rtmp=0.
        if (mnproc.eq.1) then
          call ncgeti('ll2',ll2)
          rtmp=ll2
        endif
        call xcmaxr(rtmp)
        ll2=nint(rtmp)
        ll1=3-ll2
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: time levels for interpolation of forcing fields is'
          write (lp,*)
     .     'not read from restart file.'
        endif
        ll1=1
        ll2=1
      endif
#endif
c
      call ncread('frzpot',frzpot,ip,1,0.)
      call ncread('mltpot',mltpot,ip,1,0.)
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('flxco2')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('flxco2',flxco2,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: Air-sea CO2 flux is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            flxco2(i,j)=0.
          enddo
          enddo
        enddo
      endif
c
#ifdef TKE
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('tke')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('tke',tke,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: TKE is not read from restart file and'
          write (lp,*)
     .     'will be initialized to tke_min.'
        endif
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('gls_psi')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('gls_psi',gls_psi,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: gls_psi is not read from restart file and'
          write (lp,*)
     .     'will be initialized to gls_psi_min.'
        endif
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('L_scale')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('L_scale',L_scale,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: L_scale is not read from restart file and'
          write (lp,*)
     .     'will be initialized to Ls_unlmt_min.'
        endif
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('difdia')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('difdia',difdia,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: difdia is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do k=1,kk
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              difdia(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('ustarb')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('ustarb',ustarb,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: ustarb is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            ustarb(i,j)=0.
          enddo
          enddo
        enddo
      endif
#endif
#if defined(TRC) && defined(IDLAGE)
      rtmp=0.
      if (mnproc.eq.1) then
        if (ncinqv('idlage')) rtmp=1.
      endif
      call xcmaxr(rtmp)
      if (rtmp.gt..5) then
        call ncread('idlage',trc(1-nbdy,1-nbdy,1,itriag),ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: Ideal age tracer is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
        call idlage_init
      endif
#endif
c
c --- read accumulated fields
      do n=1,nphy        
        rtmp=0.
        write(c2,'(i2.2)') n 
        if (mnproc.eq.1) then
          if (ncinqa('nacc_phy'//c2)) then 
            call ncgeti('nacc_phy'//c2,nacc_phy(n))
          else
            nacc_phy(n)=0
          endif 
          rtmp=nacc_phy(n)
        endif 
        call xcmaxr(rtmp)
        nacc_phy(n)=nint(rtmp)
        if (nacc_phy(n).gt.0) then
          if (ACC_UVEL(n).ne.0) call ncread('uvel_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)),iuu,1,0.)
          if (ACC_VVEL(n).ne.0) call ncread('vvel_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)),ivv,1,0.)
          if (ACC_DPU(n).ne.0) call ncread('dpu_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)),iuu,1,0.)
          if (ACC_DPV(n).ne.0) call ncread('dpv_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)),ivv,1,0.)
          if (ACC_UFLX(n).ne.0) call ncread('uflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)),iuu,1,0.)
          if (ACC_VFLX(n).ne.0) call ncread('vflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)),ivv,1,0.)
          if (ACC_UTFLX(n).ne.0) call ncread('utflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)),iuu,1,0.)
          if (ACC_VTFLX(n).ne.0) call ncread('vtflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)),ivv,1,0.)
          if (ACC_USFLX(n).ne.0) call ncread('usflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)),iuu,1,0.)
          if (ACC_VSFLX(n).ne.0) call ncread('vsflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)),ivv,1,0.)
          if (ACC_UMFLTD(n).ne.0) call ncread('umfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)),iuu,1,0.)
          if (ACC_VMFLTD(n).ne.0) call ncread('vmfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)),ivv,1,0.)
          if (ACC_UTFLTD(n).ne.0) call ncread('utfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)),iuu,1,0.)
          if (ACC_VTFLTD(n).ne.0) call ncread('vtfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)),ivv,1,0.)
          if (ACC_UTFLLD(n).ne.0) call ncread('utflld_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)),iuu,1,0.)
          if (ACC_VTFLLD(n).ne.0) call ncread('vtflld_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)),ivv,1,0.)
          if (ACC_USFLTD(n).ne.0) call ncread('usfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)),iuu,1,0.)
          if (ACC_VSFLTD(n).ne.0) call ncread('vsfltd_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)),ivv,1,0.)
          if (ACC_USFLLD(n).ne.0) call ncread('usflld_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)),iuu,1,0.)
          if (ACC_VSFLLD(n).ne.0) call ncread('vsflld_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)),ivv,1,0.)
          if (ACC_SALN(n).ne.0) call ncread('saln_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)),ip,1,0.)
          if (ACC_TEMP(n).ne.0) call ncread('temp_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)),ip,1,0.)
          if (ACC_DP(n).ne.0) call ncread('dp_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)),ip,1,0.)
          if (ACC_DZ(n).ne.0) call ncread('dz_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)),ip,1,0.)
          if (ACC_DIAFLX(n).ne.0) call ncread('diaflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DIAFLX(n)),ip,1,0.)
          if (ACC_DIFDIA(n).ne.0) call ncread('difdia_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)),ip,1,0.)
          if (ACC_DIFINT(n).ne.0) call ncread('difint_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)),ip,1,0.)
          if (ACC_DIFISO(n).ne.0) call ncread('difiso_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)),ip,1,0.)
          if (ACC_WFLX(n).ne.0) call ncread('wflx_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)),ip,1,0.)
          if (ACC_WFLX2(n).ne.0) call ncread('wflx2_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)),ip,1,0.)
#if defined TKE
          if (ACC_TKE(n).ne.0) call ncread('tke_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)),ip,1,0.)
          if (ACC_GLS_PSI(n).ne.0) call ncread('gls_psi_phy'//c2,
     .      phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)),ip,1,0.)
#endif
          if (ACC_UVELLVL(n).ne.0) call ncread('uvellvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UVELLVL(n)),iuu,1,0.)
          if (ACC_VVELLVL(n).ne.0) call ncread('vvellvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VVELLVL(n)),ivv,1,0.)
          if (ACC_UFLXLVL(n).ne.0) call ncread('uflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)),iuu,1,0.)
          if (ACC_VFLXLVL(n).ne.0) call ncread('vflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)),ivv,1,0.)
          if (ACC_UTFLXLVL(n).ne.0) call ncread('utflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UTFLXLVL(n)),iuu,1,0.)
          if (ACC_VTFLXLVL(n).ne.0) call ncread('vtflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VTFLXLVL(n)),ivv,1,0.)
          if (ACC_USFLXLVL(n).ne.0) call ncread('usflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_USFLXLVL(n)),iuu,1,0.)
          if (ACC_VSFLXLVL(n).ne.0) call ncread('vsflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VSFLXLVL(n)),ivv,1,0.)
          if (ACC_UMFLTDLVL(n).ne.0) call ncread('umfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UMFLTDLVL(n)),iuu,1,0.)
          if (ACC_VMFLTDLVL(n).ne.0) call ncread('vmfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VMFLTDLVL(n)),ivv,1,0.)
          if (ACC_UTFLTDLVL(n).ne.0) call ncread('utfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UTFLTDLVL(n)),iuu,1,0.)
          if (ACC_VTFLTDLVL(n).ne.0) call ncread('vtfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VTFLTDLVL(n)),ivv,1,0.)
          if (ACC_UTFLLDLVL(n).ne.0) call ncread('utflldlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_UTFLLDLVL(n)),iuu,1,0.)
          if (ACC_VTFLLDLVL(n).ne.0) call ncread('vtflldlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VTFLLDLVL(n)),ivv,1,0.)
          if (ACC_USFLTDLVL(n).ne.0) call ncread('usfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_USFLTDLVL(n)),iuu,1,0.)
          if (ACC_VSFLTDLVL(n).ne.0) call ncread('vsfltdlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VSFLTDLVL(n)),ivv,1,0.)
          if (ACC_USFLLDLVL(n).ne.0) call ncread('usflldlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_USFLLDLVL(n)),iuu,1,0.)
          if (ACC_VSFLLDLVL(n).ne.0) call ncread('vsflldlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_VSFLLDLVL(n)),ivv,1,0.)
          if (ACC_SALNLVL(n).ne.0) call ncread('salnlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_SALNLVL(n)),ip,1,0.)
          if (ACC_TEMPLVL(n).ne.0) call ncread('templvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_TEMPLVL(n)),ip,1,0.)
          if (ACC_DZLVL(n).ne.0) call ncread('dzlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_DZLVL(n)),ip,1,0.)
          if (ACC_DIAFLXLVL(n).ne.0) call ncread('diaflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_DIAFLXLVL(n)),ip,1,0.)
          if (ACC_DIFDIALVL(n).ne.0) call ncread('difdialvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_DIFDIALVL(n)),ip,1,0.)
          if (ACC_DIFINTLVL(n).ne.0) call ncread('difintlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_DIFINTLVL(n)),ip,1,0.)
          if (ACC_DIFISOLVL(n).ne.0) call ncread('difisolvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_DIFISOLVL(n)),ip,1,0.)
          if (ACC_WFLXLVL(n).ne.0) call ncread('wflxlvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)),ip,1,0.)
          if (ACC_WFLX2LVL(n).ne.0) call ncread('wflx2lvl_phy'//c2,
     .       phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)),ip,1,0.)
#if defined TKE
          if (ACC_TKELVL(n).ne.0) call ncread('tkelvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_TKELVL(n)),ip,1,0.)
          if (ACC_GLS_PSILVL(n).ne.0) call ncread('gls_psilvl_phy'//c2,
     .      phylvl(1-nbdy,1-nbdy,1,ACC_GLS_PSILVL(n)),ip,1,0.)
#endif
          if (ACC_UB(n).ne.0) call ncread('ub_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_UB(n)),iuu,1,0.)
          if (ACC_VB(n).ne.0) call ncread('vb_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_VB(n)),ivv,1,0.)
          if (ACC_ZTX(n).ne.0) call ncread('ztx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_ZTX(n)),iuu,1,0.)
          if (ACC_MTY(n).ne.0) call ncread('mty_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MTY(n)),ivv,1,0.)
          if (ACC_TAUX(n).ne.0) call ncread('taux_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_TAUX(n)),iuu,1,0.)
          if (ACC_TAUY(n).ne.0) call ncread('tauy_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_TAUY(n)),ivv,1,0.)
          if (ACC_MXLU(n).ne.0) call ncread('mxlu_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MXLU(n)),iuu,1,0.)
          if (ACC_MXLV(n).ne.0) call ncread('mxlv_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MXLV(n)),ivv,1,0.)
          if (ACC_MLDU(n).ne.0) call ncread('mldu_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MLDU(n)),iuu,1,0.)
          if (ACC_MLDV(n).ne.0) call ncread('mldv_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MLDV(n)),ivv,1,0.)
          if (ACC_UICE(n).ne.0) call ncread('uice_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_UICE(n)),iuu,1,0.)
          if (ACC_VICE(n).ne.0) call ncread('vice_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_VICE(n)),ivv,1,0.)
          if (ACC_IVOLU(n).ne.0) call ncread('ivolu_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_IVOLU(n)),iuu,1,0.)
          if (ACC_IVOLV(n).ne.0) call ncread('ivolv_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_IVOLV(n)),ivv,1,0.)
          if (ACC_SEALV(n).ne.0) call ncread('sealv_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SEALV(n)),ip,1,0.)
          if (ACC_SFL(n).ne.0) call ncread('sfl_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SFL(n)),ip,1,0.)
          if (ACC_SIGMX(n).ne.0) call ncread('sigmx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SIGMX(n)),ip,1,0.)
          if (ACC_HICE(n).ne.0) call ncread('hice_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_HICE(n)),ip,1,0.)
          if (ACC_HSNW(n).ne.0) call ncread('hsnw_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_HSNW(n)),ip,1,0.)
          if (ACC_FICE(n).ne.0) call ncread('fice_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_FICE(n)),ip,1,0.)
          if (ACC_TSRF(n).ne.0) call ncread('tsrf_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_TSRF(n)),ip,1,0.)
          if (ACC_TICE(n).ne.0) call ncread('tice_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_TICE(n)),ip,1,0.)
          if (ACC_IAGE(n).ne.0) call ncread('iage_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_IAGE(n)),ip,1,0.)
          if (ACC_SWA(n).ne.0) call ncread('swa_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SWA(n)),ip,1,0.)
          if (ACC_NSF(n).ne.0) call ncread('nsf_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_NSF(n)),ip,1,0.)
          if (ACC_DFL(n).ne.0) call ncread('dfl_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_DFL(n)),ip,1,0.)
          if (ACC_LIP(n).ne.0) call ncread('lip_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_LIP(n)),ip,1,0.)
          if (ACC_SOP(n).ne.0) call ncread('sop_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SOP(n)),ip,1,0.)
          if (ACC_EVA(n).ne.0) call ncread('eva_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_EVA(n)),ip,1,0.)
          if (ACC_ALB(n).ne.0) call ncread('alb_phy'//c2, 
     .      physrf(1-nbdy,1-nbdy,ACC_ALB(n)),ip,1,0.)
          if (ACC_RNFFLX(n).ne.0) call ncread('rnfflx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_RNFFLX(n)),ip,1,0.)
          if (ACC_USTAR(n).ne.0) call ncread('ustar_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_USTAR(n)),ip,1,0.)
          if (ACC_SSS(n).ne.0) call ncread('sss_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SSS(n)),ip,1,0.)
          if (ACC_SST(n).ne.0) call ncread('sst_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SST(n)),ip,1,0.)
          if (ACC_MLD(n).ne.0) call ncread('mld_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MLD(n)),ip,1,0.)
          if (ACC_MAXMLD(n).ne.0) call ncread('maxmld_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_MAXMLD(n)),ip,1,0.)
          if (ACC_SURFLX(n).ne.0) call ncread('surflx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SURFLX(n)),ip,1,0.)
          if (ACC_SURRLX(n).ne.0) call ncread('surrlx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SURRLX(n)),ip,1,0.)
          if (ACC_SALFLX(n).ne.0) call ncread('salflx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SALFLX(n)),ip,1,0.)
          if (ACC_SALRLX(n).ne.0) call ncread('salrlx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_SALRLX(n)),ip,1,0.)
          if (ACC_ABSWND(n).ne.0) call ncread('abswnd_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_ABSWND(n)),ip,1,0.)
          if (ACC_BRNFLX(n).ne.0) call ncread('brnflx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_BRNFLX(n)),ip,1,0.)
          if (ACC_BRNPD(n).ne.0) call ncread('brnpd_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_BRNPD(n)),ip,1,0.)
          if (ACC_FMLTFZ(n).ne.0) call ncread('fmltfz_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_FMLTFZ(n)),ip,1,0.)
          if (ACC_HMLTFZ(n).ne.0) call ncread('hmltfz_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_HMLTFZ(n)),ip,1,0.)
          if (ACC_RFIFLX(n).ne.0) call ncread('rfiflx_phy'//c2,
     .      physrf(1-nbdy,1-nbdy,ACC_RFIFLX(n)),ip,1,0.)
        endif
      enddo
c
      call ncfcls
c
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.1) then
            kfpla(i,j,1)=nint(rkfpla(i,j,1))
            kfpla(i,j,2)=nint(rkfpla(i,j,2))
          endif
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- set minimum physical temperature for each isopycnic layer
c --- ------------------------------------------------------------------
c
      call settemmin
c
#ifdef ICEDYN
      call restart_icerd(nday1,rstfnm,path2,path2_len)
#endif
c
#ifdef TRC
      call restart_trcrd(nday1,rstfnm)
#endif
c
      if (ditflx) then
c
c --- - read diag. heat flux restart file if available
#ifdef CCSMCOUPLED
        fnm=runid(1:runid_len)//'.micom.rtflx.'//rstfnm(runid_len+10:)
#else
        fnm=runid(1:runid_len)//'_tflx_'//rstfnm(runid_len+10:)
#endif
        inquire(file=path2(1:path2_len)//fnm,exist=fexist)
        rfexist=0.
        if (mnproc.eq.1.and.fexist) rfexist=1.
        call xcmaxr(rfexist)
        if (rfexist.gt.0.5) then
          call ncfopn(path2(1:path2_len)//fnm,'r',' ',1)
          if (mnproc.eq.1) write (lp,'(a,a)')
     .      ' reading diag. heat flux restart file ',
     .      trim(path2(1:path2_len)//fnm)
          call ncgetr('time',time)
          if (nint(time).ne.nday1.and.mnproc.eq.1) then
            write (lp,'(a,i6.6,a)')
     .        ' Integration day ',nint(time),
     .        ' in diag. heat flux restart file differs from'
            write (lp,'(a,i6.6,a)')
     .        ' start day ',nday1,' in limits file!'
            call ncfcls
            call xcstop('(restart_rd)') 
                   stop '(restart_rd)' 
          endif
          call ncread('tflxdi',tflxdi,ip,1,0.)
          call ncgeti('nflxdi',nflxdi)
        else
          if (mnproc.eq.1) write (lp,*)
     .      'Warning! No diag. heat flux restart file found'
        endif 
c
      endif
c
      if (disflx) then
c
c --- - read diag. salt flux restart file if available
#ifdef CCSMCOUPLED
        fnm=runid(1:runid_len)//'.micom.rsflx.'//rstfnm(runid_len+10:)
#else
        fnm=runid(1:runid_len)//'_sflx_'//rstfnm(runid_len+10:)
#endif
        inquire(file=path2(1:path2_len)//fnm,exist=fexist)
        rfexist=0.
        if (mnproc.eq.1.and.fexist) rfexist=1.
        call xcmaxr(rfexist)
        if (rfexist.gt.0.5) then
          call ncfopn(path2(1:path2_len)//fnm,'r',' ',1)
          if (mnproc.eq.1) write (lp,'(a,a)')
     .      ' reading diag. salt flux restart file ',
     .      trim(path2(1:path2_len)//fnm)
          call ncgetr('time',time)
          if (nint(time).ne.nday1.and.mnproc.eq.1) then
            write (lp,'(a,i6.6,a)')
     .        ' Integration day ',nint(time),
     .        ' in diag. salt flux restart file differs from'
            write (lp,'(a,i6.6,a)')
     .        ' start day ',nday1,' in limits file!'
            call ncfcls
            call xcstop('(restart_rd)') 
                   stop '(restart_rd)' 
          endif
          call ncread('sflxdi',sflxdi,ip,1,0.)
          call ncgeti('nflxdi',nflxdi)
        else
          if (mnproc.eq.1) write (lp,*)
     .      'Warning! No diag. salt flux restart file found'
        endif
c
      endif  
c
      return
      end
