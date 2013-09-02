      subroutine restart_wt
c
c --- ------------------------------------------------------------------
c --- Write model state to restart files
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia
      use mod_nctools
c
      implicit none
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
      integer i,j,n
      character(len=80), dimension(4) :: rstdate_str
      character(len=80) :: rstfnm,fnm
      character(len=2) :: c2
      character(len=5) :: c5p,c5u,c5v,c5q
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv,iqq
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
      logical first
      data first /.true./
      save first,iuu,ivv,iqq
c
      character startdate*20,timeunits*30
      real datenum
c
#if defined TKE
#  include "common_tke.h"
#endif
#if defined(TRC) && defined(IDLAGE)
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write(startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*rstfrq/nstep_in_day-time0)*24.
c
c --- formulate restart name and open restart file
c
      if (mnproc.eq.1) then
c
#if defined(CCSMCOUPLED)
        write (rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)')
     .    runid(1:runid_len),'.micom.r.',nyear,'-',nmonth,'-',nday,'-',
     .    mod(nstep,nstep_in_day)*86400,'.nc'
c
#elif defined(CLIM) || defined(SYN)
        if (nday_of_year-max(1,nint(rstfrq/nstep_in_day)).le.0) 
     .     then
          write(rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')
     .      runid(1:runid_len),'_restphy_',nyear,'.',nmonth,'.',nday,
     .      '_',nint(time),'.nc'
        else
          if (rstmon) then
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',mod(nmonth+10,3)+1,'.nc'
          else
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',
     .        mod(nint(min(nstep/rstfrq,time))-1,3)+1,'.nc'
          endif
          open (unit=nfu,file=path2(1:path2_len)//'rstdate.txt')
          i=1
 300      read (nfu,'(a)',end=301) rstdate_str(i)
          i=i+1
          goto 300
 301      close (unit=nfu)
          write(rstdate_str(i),'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6)')
     .      rstfnm(1:runid_len+13),': date ',nyear,'.',nmonth,'.',nday,
     .      ', integration day ',nint(time)
          if (i.eq.1) then
            open (unit=nfu,file=path2(1:path2_len)//'rstdate.txt')
            write (nfu,'(a)') rstdate_str(1)(1:runid_len+54)
            close (unit=nfu)
          elseif (rstdate_str(max(1,i-2)).ne.rstdate_str(i).and.
     .            rstdate_str(i-1       ).ne.rstdate_str(i)) then
            open (unit=nfu,file=path2(1:path2_len)//'rstdate.txt')
            do j=max(1,i-2),i
              write (nfu,'(a)') rstdate_str(j)(1:runid_len+54)
            enddo
            close (unit=nfu)
          endif
        endif
#endif
c
        write (lp,'(a,a)')
     .    ' saving restart file ',trim(path2(1:path2_len)//rstfnm)
c
        if (rstfmt.eq.1) then 
          call ncfopn(path2(1:path2_len)//rstfnm,'w','6',1)
        elseif (rstfmt.eq.2) then
          call ncfopn(path2(1:path2_len)//rstfnm,'w','h',1)
        else
          call ncfopn(path2(1:path2_len)//rstfnm,'w','c',1)
        endif
        call ncputi('nday0',nday0)
        call ncputi('nmonth0',nmonth0)
        call ncputi('nyear0',nyear0)
        call ncputr('time0',time0)
        call ncputr('time',time)
c
      endif
c
c --- define spatial and time dimensions
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
      if (rstcmp.eq.1) then 
        call ncdimc('pcomp',ip,0)
        call ncdimc('qcomp',iqq,0)
        call ncdimc('ucomp',iuu,0)
        call ncdimc('vcomp',ivv,0)
        c5p='pcomp'
        c5u='ucomp'
        c5v='vcomp'
        c5q='qcomp'
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
        c5p='x y'
        c5u='x y'
        c5v='x y'
        c5q='x y'
      endif 
      call ncdims('k2',2)
      call ncdims('k3',3)
      call ncdims('kk',kk)
      call ncdims('kk2',2*kk)
      call ncdims('plev',ddm)
      call ncputr('plev',depthslev)
      call ncdims('time',1)
c
c --- output model fields to restart file
      if (sprfac) then
        call ncputr('prfac',prfac)
      endif
      call ncputi('ntda',ntda)
      do n=1,nphy 
        write(c2,'(i2.2)') n
        call ncputi('nacc_phy'//c2,nacc_phy(n))
      enddo 
#ifdef CCSMCOUPLED
      call ncputi('ll2',ll2)
#endif
c
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.1) then
            rkfpla(i,j,1)=real(kfpla(i,j,1))
            rkfpla(i,j,2)=real(kfpla(i,j,2))
          else
            rkfpla(i,j,1)=0.
            rkfpla(i,j,2)=0.
          endif
        enddo
      enddo
c
      call wrtrst('u',trim(c5u)//' kk2 time',u,iuu)
      call wrtrst('v',trim(c5v)//' kk2 time',v,ivv)
      call wrtrst('dp',trim(c5p)//' kk2 time',dp,ip)
      call wrtrst('dpold',trim(c5p)//' kk2 time',dpold,ip)
      call wrtrst('temp',trim(c5p)//' kk2 time',temp,ip)
      call wrtrst('saln',trim(c5p)//' kk2 time',saln,ip)
      call wrtrst('sigma',trim(c5p)//' kk2 time',sigma,ip)
      call wrtrst('sigmar',trim(c5p)//' kk time',sigmar,ip)
      call wrtrst('pgfx',trim(c5u)//' kk2 time',pgfx,iuu)
      call wrtrst('pgfy',trim(c5v)//' kk2 time',pgfy,ivv)
      call wrtrst('pb',trim(c5p)//' k2 time',pb,ip)
      call wrtrst('pb_mn',trim(c5p)//' k2 time',pb_mn,ip)
      call wrtrst('pb_p',trim(c5p)//' time',pb_p,ip)
      call wrtrst('pbu',trim(c5u)//' k2 time',pbu,iuu)
      call wrtrst('pbv',trim(c5v)//' k2 time',pbv,ivv)
      call wrtrst('pbu_p',trim(c5u)//' time',pbu_p,iuu)
      call wrtrst('pbv_p',trim(c5v)//' time',pbv_p,ivv)
      call wrtrst('ub',trim(c5u)//' k2 time',ub,iuu)
      call wrtrst('vb',trim(c5v)//' k2 time',vb,ivv)
      call wrtrst('uflx',trim(c5u)//' kk2 time',uflx,iuu)
      call wrtrst('utflx',trim(c5u)//' kk2 time',utflx,iuu)
      call wrtrst('usflx',trim(c5u)//' kk2 time',usflx,iuu)
      call wrtrst('umfltd',trim(c5u)//' kk2 time',umfltd,iuu)
      call wrtrst('utfltd',trim(c5u)//' kk2 time',utfltd,iuu)
      call wrtrst('utflld',trim(c5u)//' kk2 time',utflld,iuu)
      call wrtrst('usfltd',trim(c5u)//' kk2 time',usfltd,iuu)
      call wrtrst('usflld',trim(c5u)//' kk2 time',usflld,iuu)
      call wrtrst('vflx',trim(c5v)//' kk2 time',vflx,ivv)
      call wrtrst('vtflx',trim(c5v)//' kk2 time',vtflx,ivv)
      call wrtrst('vsflx',trim(c5v)//' kk2 time',vsflx,ivv)
      call wrtrst('vmfltd',trim(c5v)//' kk2 time',vmfltd,ivv)
      call wrtrst('vtfltd',trim(c5v)//' kk2 time',vtfltd,ivv)
      call wrtrst('vtflld',trim(c5v)//' kk2 time',vtflld,ivv)
      call wrtrst('vsfltd',trim(c5v)//' kk2 time',vsfltd,ivv)
      call wrtrst('vsflld',trim(c5v)//' kk2 time',vsflld,ivv)
      call wrtrst('ubflx',trim(c5u)//' k2 time',ubflx,iuu)
      call wrtrst('vbflx',trim(c5v)//' k2 time',vbflx,ivv)
      call wrtrst('ubflx_mn',trim(c5u)//' k2 time',ubflx_mn,iuu)
      call wrtrst('vbflx_mn',trim(c5v)//' k2 time',vbflx_mn,ivv)
      call wrtrst('ubflxs',trim(c5u)//' k3 time',ubflxs,iuu)
      call wrtrst('vbflxs',trim(c5v)//' k3 time',vbflxs,ivv)
      call wrtrst('ubflxs_p',trim(c5u)//' k2 time',ubflxs_p,iuu)
      call wrtrst('vbflxs_p',trim(c5v)//' k2 time',vbflxs_p,ivv)
      call wrtrst('ubcors_p',trim(c5u)//' time',ubcors_p,iuu)
      call wrtrst('vbcors_p',trim(c5v)//' time',vbcors_p,ivv)
      call wrtrst('pvtrop',trim(c5q)//' k2 time',pvtrop,iqq)
      call wrtrst('pgfxm',trim(c5u)//' k2 time',pgfxm,iuu)
      call wrtrst('pgfym',trim(c5v)//' k2 time',pgfym,ivv)
      call wrtrst('xixp',trim(c5u)//' k2 time',xixp,iuu)
      call wrtrst('xixm',trim(c5u)//' k2 time',xixm,iuu)
      call wrtrst('xiyp',trim(c5v)//' k2 time',xiyp,ivv)
      call wrtrst('xiym',trim(c5v)//' k2 time',xiym,ivv)
      call wrtrst('phi',trim(c5p)//' time',phi(1-nbdy,1-nbdy,kk+1),ip)
      call wrtrst('sealv',trim(c5p)//' time',sealv,ip)
      call wrtrst('ustar',trim(c5p)//' time',ustar,ip)
      call wrtrst('buoyfl',trim(c5p)//' time',buoyfl,ip)
      call wrtrst('kfpla',trim(c5p)//' k2 time',rkfpla,ip)
      call wrtrst('hicem',trim(c5p)//' time',hicem,ip)
      call wrtrst('ficem',trim(c5p)//' time',ficem,ip)
      call wrtrst('tsrfm',trim(c5p)//' time',tsrfm,ip)
      call wrtrst('hsnwm',trim(c5p)//' time',hsnwm,ip)
      call wrtrst('ticem',trim(c5p)//' time',ticem,ip)
      call wrtrst('iagem',trim(c5p)//' time',iagem,ip)
      call wrtrst('tsi_tda',trim(c5p)//' time',tsi_tda,ip)
      call wrtrst('tml_tda',trim(c5p)//' time',tml_tda,ip)
      call wrtrst('sml_tda',trim(c5p)//' time',sml_tda,ip)
      call wrtrst('alb_tda',trim(c5p)//' time',alb_tda,ip)
      call wrtrst('fice_tda',trim(c5p)//' time',fice_tda,ip)
      call wrtrst('rnfres',trim(c5p)//' time',rnfres,ip)
c
      if (sprfac) then
        call wrtrst('eiacc',trim(c5p)//' time',eiacc,ip)
        call wrtrst('pracc',trim(c5p)//' time',pracc,ip)
      endif
c
#if defined(CLIM) || defined(SYN)
      call wrtrst('cd_d',trim(c5p)//' time',cd_d,ip)
      call wrtrst('ch_d',trim(c5p)//' time',ch_d,ip)
      call wrtrst('ce_d',trim(c5p)//' time',ce_d,ip)
      call wrtrst('wg2_d',trim(c5p)//' time',wg2_d,ip)
      call wrtrst('cd_m',trim(c5p)//' time',cd_m,ip)
      call wrtrst('ch_m',trim(c5p)//' time',ch_m,ip)
      call wrtrst('ce_m',trim(c5p)//' time',ce_m,ip)
      call wrtrst('wg2_m',trim(c5p)//' time',wg2_m,ip)
      call wrtrst('rhoa',trim(c5p)//' time',rhoa,ip)
#endif
c
#ifdef CCSMCOUPLED
      call wrtrst('ustarw_da',trim(c5p)//' k2 time',ustarw_da,ip)
      call wrtrst('ztx_da',trim(c5p)//' k2 time',ztx_da,ip)
      call wrtrst('mty_da',trim(c5p)//' k2 time',mty_da,ip)
      call wrtrst('lip_da',trim(c5p)//' k2 time',lip_da,ip)
      call wrtrst('sop_da',trim(c5p)//' k2 time',sop_da,ip)
      call wrtrst('eva_da',trim(c5p)//' k2 time',eva_da,ip)
      call wrtrst('rnf_da',trim(c5p)//' k2 time',rnf_da,ip)
      call wrtrst('rfi_da',trim(c5p)//' k2 time',rfi_da,ip)
      call wrtrst('fmltfz_da',trim(c5p)//' k2 time',fmltfz_da,ip)
      call wrtrst('sfl_da',trim(c5p)//' k2 time',sfl_da,ip)
      call wrtrst('swa_da',trim(c5p)//' k2 time',swa_da,ip)
      call wrtrst('nsf_da',trim(c5p)//' k2 time',nsf_da,ip)
      call wrtrst('hmlt_da',trim(c5p)//' k2 time',hmlt_da,ip)
      call wrtrst('slp_da',trim(c5p)//' k2 time',slp_da,ip)
      call wrtrst('ficem_da',trim(c5p)//' k2 time',ficem_da,ip)
      call wrtrst('abswnd_da',trim(c5p)//' k2 time',abswnd_da,ip)
      call wrtrst('atmco2_da',trim(c5p)//' k2 time',atmco2_da,ip)
#endif
c
      call wrtrst('frzpot',trim(c5p)//' time',frzpot,ip)
      call wrtrst('mltpot',trim(c5p)//' time',mltpot,ip)
      call wrtrst('flxco2',trim(c5p)//' time',flxco2,ip)
c
#if defined TKE
      call wrtrst('tke',trim(c5p)//' kk2 time',tke,ip)
      call wrtrst('gls_psi',trim(c5p)//' kk2 time',gls_psi,ip)
      call wrtrst('L_scale',trim(c5p)//' kk time',L_scale,ip)
      call wrtrst('difdia',trim(c5p)//' kk time',difdia,ip)
      call wrtrst('ustarb',trim(c5p)//' time',ustarb,ip)
#endif
#if defined(TRC) && defined(IDLAGE)
      call wrtrst('idlage',trim(c5p)//' kk time',
     .            trc(1-nbdy,1-nbdy,1,itriag),ip)
#endif
c
c --- write accumulated fields
      do n=1,nphy 
        write(c2,'(i2.2)') n
        if (nacc_phy(n).ne.0) then
          if (ACC_UB(n)    .ne.0) call wrtrst('ub_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_UB(n)),iuu)
          if (ACC_VB(n)    .ne.0) call wrtrst('vb_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_VB(n)),ivv)
          if (ACC_ZTX(n)   .ne.0) call wrtrst('ztx_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_ZTX(n)),iuu)
          if (ACC_MTY(n)   .ne.0) call wrtrst('mty_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_MTY(n)),ivv)
          if (ACC_TAUX(n)  .ne.0) call wrtrst('taux_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_TAUX(n)),iuu)
          if (ACC_TAUY(n)  .ne.0) call wrtrst('tauy_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_TAUY(n)),ivv)
          if (ACC_MXLU(n)  .ne.0) call wrtrst('mxlu_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_MXLU(n)),iuu)
          if (ACC_MXLV(n)  .ne.0) call wrtrst('mxlv_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_MXLV(n)),ivv)
          if (ACC_MLDU(n)  .ne.0) call wrtrst('mldu_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_MLDU(n)),iuu)
          if (ACC_MLDV(n)  .ne.0) call wrtrst('mldv_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_MLDV(n)),ivv)
          if (ACC_SEALV(n) .ne.0) call wrtrst('sealv_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SEALV(n)),ip)
          if (ACC_SIGMX(n) .ne.0) call wrtrst('sigmx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SIGMX(n)),ip)
          if (ACC_HICE(n)  .ne.0) call wrtrst('hice_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_HICE(n)),ip)
          if (ACC_IAGE(n)  .ne.0) call wrtrst('iage_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_IAGE(n)),ip)
          if (ACC_HSNW(n)  .ne.0) call wrtrst('hsnw_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_HSNW(n)),ip)
          if (ACC_FICE(n)  .ne.0) call wrtrst('fice_phy'//c2, 
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_FICE(n)),ip)
          if (ACC_TSRF(n)  .ne.0) call wrtrst('tsrf_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_TSRF(n)),ip)
          if (ACC_TICE(n)  .ne.0) call wrtrst('tice_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_TICE(n)),ip)
          if (ACC_SWA(n)   .ne.0) call wrtrst('swa_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SWA(n)),ip)
          if (ACC_NSF(n)   .ne.0) call wrtrst('nsf_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_NSF(n)),ip)
          if (ACC_DFL(n)   .ne.0) call wrtrst('dfl_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_DFL(n)),ip)
          if (ACC_LIP(n)   .ne.0) call wrtrst('lip_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_LIP(n)),ip)
          if (ACC_SOP(n)   .ne.0) call wrtrst('sop_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SOP(n)),ip)
          if (ACC_EVA(n)   .ne.0) call wrtrst('eva_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_EVA(n)),ip)
          if (ACC_ABSWND(n)   .ne.0) call wrtrst('abswnd_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_ABSWND(n)),ip)
          if (ACC_ALB(n)   .ne.0) call wrtrst('alb_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_ALB(n)),ip)
          if (ACC_SFL(n)   .ne.0) call wrtrst('sfl_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SFL(n)),ip)
          if (ACC_BRNFLX(n)   .ne.0) call wrtrst('brnflx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_BRNFLX(n)),ip)
          if (ACC_BRNPD(n)   .ne.0) call wrtrst('brnpd_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_BRNPD(n)),ip)
          if (ACC_FMLTFZ(n)   .ne.0) call wrtrst('fmltfz_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_FMLTFZ(n)),ip)
          if (ACC_HMLTFZ(n)   .ne.0) call wrtrst('hmltfz_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_HMLTFZ(n)),ip)
          if (ACC_RFIFLX(n).ne.0) call wrtrst('rfiflx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_RFIFLX(n)),ip)
          if (ACC_RNFFLX(n).ne.0) call wrtrst('rnfflx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_RNFFLX(n)),ip)
          if (ACC_USTAR(n) .ne.0) call wrtrst('ustar_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_USTAR(n)),ip)
          if (ACC_SSS(n)   .ne.0) call wrtrst('sss_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SSS(n)),ip)
          if (ACC_SST(n)   .ne.0) call wrtrst('sst_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SST(n)),ip)
          if (ACC_MLD(n)   .ne.0) call wrtrst('mld_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_MLD(n)),ip)
          if (ACC_MAXMLD(n).ne.0) call wrtrst('maxmld_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_MAXMLD(n)),ip)
          if (ACC_SURFLX(n).ne.0) call wrtrst('surflx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SURFLX(n)),ip)
          if (ACC_SURRLX(n).ne.0) call wrtrst('surrlx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SURRLX(n)),ip)
          if (ACC_SALFLX(n).ne.0) call wrtrst('salflx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SALFLX(n)),ip)
          if (ACC_SALRLX(n).ne.0) call wrtrst('salrlx_phy'//c2,
     .      trim(c5p)//' time',physrf(1-nbdy,1-nbdy,ACC_SALRLX(n)),ip)
          if (ACC_UICE(n)  .ne.0) call wrtrst('uice_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_UICE(n)),iuu)
          if (ACC_VICE(n)  .ne.0) call wrtrst('vice_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_VICE(n)),ivv)
          if (ACC_IVOLU(n)  .ne.0) call wrtrst('ivolu_phy'//c2,
     .      trim(c5u)//' time',physrf(1-nbdy,1-nbdy,ACC_IVOLU(n)),iuu)
          if (ACC_IVOLV(n)  .ne.0) call wrtrst('ivolv_phy'//c2,
     .      trim(c5v)//' time',physrf(1-nbdy,1-nbdy,ACC_IVOLV(n)),ivv)
          if (ACC_UVEL(n)  .ne.0) call wrtrst('uvel_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)),
     .      iuu)
          if (ACC_VVEL(n)  .ne.0) call wrtrst('vvel_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)),
     .      ivv)
          if (ACC_DPU(n)   .ne.0) call wrtrst('dpu_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)),
     .      iuu)
          if (ACC_DPV(n)   .ne.0) call wrtrst('dpv_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)),
     .      ivv)
          if (ACC_UFLX(n)  .ne.0) call wrtrst('uflx_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)),
     .      iuu)
          if (ACC_VFLX(n)  .ne.0) call wrtrst('vflx_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)),
     .      ivv)
          if (ACC_UTFLX(n) .ne.0) call wrtrst('utflx_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)),
     .      iuu)
          if (ACC_VTFLX(n) .ne.0) call wrtrst('vtflx_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)),
     .      ivv)
          if (ACC_USFLX(n) .ne.0) call wrtrst('usflx_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)),
     .      iuu)
          if (ACC_VSFLX(n) .ne.0) call wrtrst('vsflx_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)),
     .      ivv)
          if (ACC_UMFLTD(n).ne.0) call wrtrst('umfltd_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)),
     .      iuu)
          if (ACC_VMFLTD(n).ne.0) call wrtrst('vmfltd_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)),
     .      ivv)
          if (ACC_UTFLTD(n).ne.0) call wrtrst('utfltd_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)),
     .      iuu)
          if (ACC_VTFLTD(n).ne.0) call wrtrst('vtfltd_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)),
     .      ivv)
          if (ACC_UTFLLD(n).ne.0) call wrtrst('utflld_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)),
     .      iuu)
          if (ACC_VTFLLD(n).ne.0) call wrtrst('vtflld_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)),
     .      ivv)
          if (ACC_USFLTD(n).ne.0) call wrtrst('usfltd_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)),
     .      iuu)
          if (ACC_VSFLTD(n).ne.0) call wrtrst('vsfltd_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)),
     .      ivv)
          if (ACC_USFLLD(n).ne.0) call wrtrst('usflld_phy'//c2,
     .      trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)),
     .      iuu)
          if (ACC_VSFLLD(n).ne.0) call wrtrst('vsflld_phy'//c2,
     .      trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)),
     .      ivv)
          if (ACC_SALN(n)  .ne.0) call wrtrst('saln_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)),
     .      ip)
          if (ACC_TEMP(n)  .ne.0) call wrtrst('temp_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)),
     .      ip)
          if (ACC_DP(n)    .ne.0) call wrtrst('dp_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)),ip)
          if (ACC_DZ(n)    .ne.0) call wrtrst('dz_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)),ip)
          if (ACC_DIAFLX(n).ne.0) call wrtrst('diaflx_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIAFLX(n)),
     .      ip)
          if (ACC_DIFDIA(n).ne.0) call wrtrst('difdia_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)),
     .      ip)
          if (ACC_DIFINT(n).ne.0) call wrtrst('difint_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)),
     .      ip)
          if (ACC_DIFISO(n).ne.0) call wrtrst('difiso_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)),
     .      ip)
          if (ACC_WFLX(n)  .ne.0) call wrtrst('wflx_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)),
     .      ip)
          if (ACC_WFLX2(n) .ne.0) call wrtrst('wflx2_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)),
     .      ip)
#if defined TKE
          if (ACC_TKE(n)  .ne.0) call wrtrst('tke_phy'//c2,
     .      trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)),
     .      ip)
          if (ACC_GLS_PSI(n).ne.0) call wrtrst('gls_psi_phy'//c2,
     .     trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)),
     .      ip)
#endif
          if (ACC_UVELLVL(n) .ne.0) call wrtrst('uvellvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UVELLVL(n)),iuu)
          if (ACC_VVELLVL(n) .ne.0) call wrtrst('vvellvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VVELLVL(n)),ivv)
          if (ACC_UFLXLVL(n) .ne.0) call wrtrst('uflxlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UFLXLVL(n)),iuu)
          if (ACC_VFLXLVL(n) .ne.0) call wrtrst('vflxlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1
     .      ,ACC_VFLXLVL(n)),ivv)
          if (ACC_UTFLXLVL(n).ne.0) call wrtrst('utflxlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UTFLXLVL(n)),iuu)
          if (ACC_VTFLXLVL(n).ne.0) call wrtrst('vtflxlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VTFLXLVL(n)),ivv)
          if (ACC_USFLXLVL(n).ne.0) call wrtrst('usflxlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_USFLXLVL(n)),iuu)
          if (ACC_VSFLXLVL(n).ne.0) call wrtrst('vsflxlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VSFLXLVL(n)),ivv)
          if (ACC_UMFLTDLVL(n).ne.0) call wrtrst('umfltdlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UMFLTDLVL(n)),iuu)
          if (ACC_VMFLTDLVL(n).ne.0) call wrtrst('vmfltdlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VMFLTDLVL(n)),ivv)
          if (ACC_UTFLTDLVL(n).ne.0) call wrtrst('utfltdlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UTFLTDLVL(n)),iuu)
          if (ACC_VTFLTDLVL(n).ne.0) call wrtrst('vtfltdlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VTFLTDLVL(n)),ivv)
          if (ACC_UTFLLDLVL(n).ne.0) call wrtrst('utflldlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_UTFLLDLVL(n)),iuu)
          if (ACC_VTFLLDLVL(n).ne.0) call wrtrst('vtflldlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VTFLLDLVL(n)),ivv)
          if (ACC_USFLTDLVL(n).ne.0) call wrtrst('usfltdlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_USFLTDLVL(n)),iuu)
          if (ACC_VSFLTDLVL(n).ne.0) call wrtrst('vsfltdlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VSFLTDLVL(n)),ivv)
          if (ACC_USFLLDLVL(n).ne.0) call wrtrst('usflldlvl_phy'//c2,
     .      trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_USFLLDLVL(n)),iuu)
          if (ACC_VSFLLDLVL(n).ne.0) call wrtrst('vsflldlvl_phy'//c2,
     .      trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_VSFLLDLVL(n)),ivv)
          if (ACC_SALNLVL(n) .ne.0) call wrtrst('salnlvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_SALNLVL(n)),ip)
          if (ACC_TEMPLVL(n) .ne.0) call wrtrst('templvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_TEMPLVL(n)),ip)
          if (ACC_DZLVL(n)    .ne.0) call wrtrst('dzlvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_DZLVL(n)),ip)
          if (ACC_DIAFLXLVL(n).ne.0) call wrtrst('diaflxlvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_DIAFLXLVL(n)),ip)
          if (ACC_DIFDIALVL(n).ne.0) call wrtrst('difdialvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_DIFDIALVL(n)),ip)
          if (ACC_DIFINTLVL(n).ne.0) call wrtrst('difintlvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_DIFINTLVL(n)),ip)
          if (ACC_DIFISOLVL(n).ne.0) call wrtrst('difisolvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_DIFISOLVL(n)),ip)
          if (ACC_WFLXLVL(n).ne.0) call wrtrst('wflxlvl_phy'//c2,
     .      trim(c5p)//' plev time',
     .      phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)),ip)
          if (ACC_WFLX2LVL(n).ne.0) call wrtrst('wflx2lvl_phy'//c2,
     .      trim(c5p)//' plev time',
     .      phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)),ip)
#if defined TKE
          if (ACC_TKELVL(n).ne.0) call wrtrst('tkelvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_TKELVL(n)),ip)
          if(ACC_GLS_PSILVL(n).ne.0) call wrtrst('gls_psilvl_phy'//c2,
     .      trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1,
     .      ACC_GLS_PSILVL(n)),ip)
#endif
        endif
      enddo
c
      call ncfcls
c
#ifdef ICEDYN
      call restart_icewt(time,rstfnm,path2,path2_len)
#endif
c
#ifdef TRC
      call restart_trcwt(rstfnm)
#endif
c
      if (ditflx) then
c
c --- - write diag. heat flux restart file
#ifdef CCSMCOUPLED
        fnm=runid(1:runid_len)//'.micom.rtflx.'//rstfnm(runid_len+10:)
#else
        fnm=runid(1:runid_len)//'_tflx_'//rstfnm(runid_len+10:)
#endif
        if (mnproc.eq.1) write (lp,'(a,a)')
     .    ' saving diag. heat flux restart file ',
     .    trim(path2(1:path2_len)//fnm)
        if (rstfmt.eq.1) then 
          call ncfopn(path2(1:path2_len)//fnm,'w','6',1)
        elseif (rstfmt.eq.2) then
          call ncfopn(path2(1:path2_len)//fnm,'w','h',1)
        else
          call ncfopn(path2(1:path2_len)//fnm,'w','c',1)
        endif
c
        if (rstcmp.eq.1) then 
          call ncdimc('pcomp',ip,0)
        else
          call ncdims('x',itdm)
          call ncdims('y',jtdm)
        endif
        call ncdims('week',48)
c
        call ncputr('time',time)
        call ncputi('nflxdi',nflxdi)
c
        call wrtrst('tflxdi',trim(c5p)//' week',tflxdi,ip)
c       
        call ncfcls   
      endif  
c
      if (disflx) then
c
c --- - write diag. salt flux restart file
#ifdef CCSMCOUPLED
        fnm=runid(1:runid_len)//'.micom.rsflx.'//rstfnm(runid_len+10:)
#else
        fnm=runid(1:runid_len)//'_sflx_'//rstfnm(runid_len+10:)
#endif
        if (mnproc.eq.1) write (lp,'(a,a)')
     .    ' saving diag. salt flux restart file ',
     .    path2(1:path2_len)//trim(fnm)
        if (rstfmt.eq.1) then 
          call ncfopn(path2(1:path2_len)//fnm,'w','6',1)
        elseif (rstfmt.eq.2) then
          call ncfopn(path2(1:path2_len)//fnm,'w','h',1)
        else
          call ncfopn(path2(1:path2_len)//fnm,'w','c',1)
        endif
c
        if (rstcmp.eq.1) then 
          call ncdimc('pcomp',ip,0)
        else
          call ncdims('x',itdm)
          call ncdims('y',jtdm)
        endif
        call ncdims('week',48)
c
        call ncputr('time',time)
        call ncputi('nflxdi',nflxdi)
c
        call wrtrst('sflxdi',trim(c5p)//' week',sflxdi,ip)
c       
        call ncfcls
      endif
c
#ifdef CCSMCOUPLED
c --- write restart filename to rpointer.ocn
      if (mnproc.eq.1) then
        open (unit=nfu,file='rpointer.ocn')
        write (nfu,'(a)') rstfnm
        close (unit=nfu)
      endif
#endif
c
      return
      end



      subroutine wrtrst(vnm,dims,fld,msk)
c
c --- ------------------------------------------------------------------
      use mod_xc
      use mod_nctools
c
      implicit none
c 
      character(len=*) :: vnm,dims
      real, dimension(*) :: fld
      integer, dimension(*) :: msk
c
c
c --- Write data in compressed or uncompressed format
      if (dims(2:5).eq.'comp') then 
        call nccomp(vnm,dims,fld,msk,1.,0.,8)
      else
        call ncwrtr(vnm,dims,fld,msk,1,1.,0.,8)
      endif
c
      end
