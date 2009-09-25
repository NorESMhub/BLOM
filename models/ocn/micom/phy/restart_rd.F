      subroutine restart_rd(nday1)
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
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
#include "common_ice.h"
#include "common_dia.h"
c
      real rtmp
      integer i,j,ndayr,nmonthr,nyearr
      character rstfnm*80,rstfnm_ext*80,fnm*80
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
          call ncfopn(path2(1:path2_len)//rstfnm,'r')
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
            write (lp,'(a)') ' Warning! Integration day corresponding'
     .        //             ' to initial experiment'
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
              call ncfopn(path2(1:path2_len)//rstfnm,'r')
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
        rstfnm_ext=rstfnm(runid_len+9:80)
c
        write (lp,'(2a)') ' reading restart file ',rstfnm
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
      call ncread('vflx',vflx,ivv,1,0.)
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
      call ncread('frzpot',frzpot,ip,1,0.)
      call ncread('mltpot',mltpot,ip,1,0.)
c
c --- read accumulated fields
      rtmp=0.
      if (mnproc.eq.1) then
        call ncgeti('nacc_phy2d',nacc_phy2d)
        rtmp=nacc_phy2d
      endif
      call xcmaxr(rtmp)
      nacc_phy2d=nint(rtmp)
      if (nacc_phy2d.gt.0) then
        call ncread('ub_dia',ub_dia,iuu,1,0.)
        call ncread('vb_dia',vb_dia,ivv,1,0.)
        call ncread('ztx_dia',ztx_dia,iuu,1,0.)
        call ncread('mty_dia',mty_dia,ivv,1,0.)
        call ncread('taux_dia',taux_dia,iuu,1,0.)
        call ncread('tauy_dia',tauy_dia,ivv,1,0.)
        call ncread('mxlu_dia',mxlu_dia,iuu,1,0.)
        call ncread('mxlv_dia',mxlv_dia,ivv,1,0.)
        call ncread('mldu_dia',mldu_dia,iuu,1,0.)
        call ncread('mldv_dia',mldv_dia,ivv,1,0.)
        call ncread('sealv_dia',sealv_dia,ip,1,0.)
        call ncread('sigmx_dia',sigmx_dia,ip,1,0.)
        call ncread('hice_dia',hice_dia,ip,1,0.)
        call ncread('hsnw_dia',hsnw_dia,ip,1,0.)
        call ncread('fice_dia',fice_dia,ip,1,0.)
        call ncread('tsrf_dia',tsrf_dia,ip,1,0.)
        call ncread('tice_dia',tice_dia,ip,1,0.)
        call ncread('swa_dia',swa_dia,ip,1,0.)
        call ncread('nsf_dia',nsf_dia,ip,1,0.)
        call ncread('hmltfz_dia',hmltfz_dia,ip,1,0.)
        call ncread('dfl_dia',dfl_dia,ip,1,0.)
        call ncread('lip_dia',lip_dia,ip,1,0.)
        call ncread('sop_dia',sop_dia,ip,1,0.)
        call ncread('eva_dia',eva_dia,ip,1,0.)
        call ncread('fmltfz_dia',fmltfz_dia,ip,1,0.)
        call ncread('sfl_dia',sfl_dia,ip,1,0.)
        call ncread('alb_dia',alb_dia,ip,1,0.)
        call ncread('rnfflx_dia',rnfflx_dia,ip,1,0.)
        call ncread('rfiflx_dia',rfiflx_dia,ip,1,0.)
        call ncread('ustar_dia',ustar_dia,ip,1,0.)
        call ncread('abswnd_dia',abswnd_dia,ip,1,0.)
        call ncread('sss_dia',sss_dia,ip,1,0.)
        call ncread('sst_dia',sst_dia,ip,1,0.)
        call ncread('mld_dia',mld_dia,ip,1,0.)
        call ncread('maxmld_dia',maxmld_dia,ip,1,0.)
        call ncread('surflx_dia',surflx_dia,ip,1,0.)
        call ncread('salflx_dia',salflx_dia,ip,1,0.)
        call ncread('brnflx_dia',brnflx_dia,ip,1,0.)
        call ncread('uice_dia',uice_dia,iuu,1,0.)
        call ncread('vice_dia',vice_dia,ivv,1,0.)
      endif
c
      rtmp=0.
      if (mnproc.eq.1) then
        call ncgeti('nacc_phy3d',nacc_phy3d)
        rtmp=nacc_phy3d
      endif
      call xcmaxr(rtmp)
      nacc_phy3d=nint(rtmp)
      if (nacc_phy3d.gt.0) then
        call ncread('uvel_dia',uvel_dia,iuu,1,0.)
        call ncread('vvel_dia',vvel_dia,ivv,1,0.)
        call ncread('dpu_dia',dpu_dia,iuu,1,0.)
        call ncread('dpv_dia',dpv_dia,ivv,1,0.)
        call ncread('uflx_dia',uflx_dia,iuu,1,0.)
        call ncread('vflx_dia',vflx_dia,ivv,1,0.)
        call ncread('saln_dia',saln_dia,ip,1,0.)
        call ncread('temp_dia',temp_dia,ip,1,0.)
        call ncread('dp_dia',dp_dia,ip,1,0.)
        call ncread('dz_dia',dz_dia,ip,1,0.)
        call ncread('difint_dia',difint_dia,ip,1,0.)
        call ncread('difiso_dia',difiso_dia,ip,1,0.)
        call ncread('difdia_dia',difdia_dia,ip,1,0.)
        call ncread('diaflx',diaflx,ip,1,0.)
      endif
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
      call restart_icerd(runid,runid_len,nday1,rstfnm_ext,
     .                   path2,path2_len)
#endif
c
#ifdef TRC
      call restart_trcrd(nday1,rstfnm_ext)
#endif
c
      if (ditflx) then
c
c --- - read diag. heat flux restart file if available
        fnm=path2(1:path2_len)//runid(1:runid_len)//'_tflx'//rstfnm_ext
        inquire(file=fnm,exist=fexist)
        rfexist=0.
        if (mnproc.eq.1.and.fexist) rfexist=1.
        call xcmaxr(rfexist)
        if (rfexist.gt.0.5) then
          call ncfopn(fnm,'r')
          if (mnproc.eq.1) write (lp,'(a,a)')
     .      ' reading diag. heat flux restart file ',fnm
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
        fnm=path2(1:path2_len)//runid(1:runid_len)//'_sflx'//rstfnm_ext
        inquire(file=fnm,exist=fexist)
        rfexist=0.
        if (mnproc.eq.1.and.fexist) rfexist=1.
        call xcmaxr(rfexist)
        if (rfexist.gt.0.5) then
          call ncfopn(fnm,'r')
          if (mnproc.eq.1) write (lp,'(a,a)')
     .      ' reading diag. salt flux restart file ',fnm
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
