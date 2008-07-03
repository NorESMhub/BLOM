      subroutine restart_wt
c
c --- ------------------------------------------------------------------
c --- Write model state to restart files
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#if defined(CLIM) || defined(SYN)
#  include "common_asf.h"
#endif
#include "common_ice.h"
#ifdef CCSMCOUPLED
#  include "common_dia.h"
#endif
c
      integer i,j
      character*80, dimension(4) :: rstdate_str
      character rstfnm*80,rstfnm_ext*80,fnm*80
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv,iqq
      logical first
      data first /.true./
      save first,iuu,ivv,iqq
c
      character startdate*20,timeunits*30
      real datenum
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write(startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*diag2dfq/nstep_in_day-time0)*24.
c
c --- formulate restart name and open restart file
c
      if (mnproc.eq.1) then
c
#if defined(CCSMCOUPLED)
        write (rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')
     .    runid(1:runid_len),'_restphy_',nyear,'-',nmonth,'-',nday,'_',
     .    nint(time),'.nc'
c
#elif defined(CLIM) || defined(SYN)
        if (nday_of_year-max(1,nint(diag3dfq/nstep_in_day)).le.0) then
          write(rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')
     .      runid(1:runid_len),'_restphy_',nyear,'.',nmonth,'.',nday,
     .      '_',nint(time),'.nc'
        else
          if (diag3dmon) then
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',mod(nmonth+10,3)+1,'.nc'
          else
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',
     .        mod(nint(min(nstep/diag3dfq,time))-1,3)+1,'.nc'
          endif
          open (unit=12,file=path2(1:path2_len)//'rstdate.txt')
          i=1
 300      read (12,'(a)',end=301) rstdate_str(i)
          i=i+1
          goto 300
 301      close (unit=12)
          write(rstdate_str(i),'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6)')
     .      rstfnm(1:runid_len+13),': date ',nyear,'.',nmonth,'.',nday,
     .      ', integration day ',nint(time)
          if (i.eq.1) then
            open (unit=12,file=path2(1:path2_len)//'rstdate.txt')
            write (12,'(a)') rstdate_str(1)(1:runid_len+54)
            close (unit=12)
          elseif (rstdate_str(max(1,i-2)).ne.rstdate_str(i).and.
     .            rstdate_str(i-1       ).ne.rstdate_str(i)) then
            open (unit=12,file=path2(1:path2_len)//'rstdate.txt')
            do j=max(1,i-2),i
              write (12,'(a)') rstdate_str(j)(1:runid_len+54)
            enddo
            close (unit=12)
          endif
        endif
#endif
c
        rstfnm_ext=rstfnm(runid_len+9:80)
c
        write (lp,'(a,a)')
     .    ' saving restart file ',path2(1:path2_len)//rstfnm
c
        call ncfopn(path2(1:path2_len)//rstfnm,'w')
        call ncputi('nday0',nday0,1)
        call ncputi('nmonth0',nmonth0,1)
        call ncputi('nyear0',nyear0,1)
        call ncputr('time0',time0,1)
        call ncputr('time',time,1)
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
               end if
               if ((ip(i,j)+ip(i,j-1)).ge.1) then
                  ivv(i,j)=1
               else
                  ivv(i,j)=0
               end if
               if ((iu(i,j)+iv(i,j)+iu(i,j-1)+iv(i-1,j)).ge.1) then
                  iqq(i,j)=1
               else
                  iqq(i,j)=0
               end if
            end do 
         end do
      end if
#ifdef COMPRES     
      call ncdimc('pcomp',ip,0)
      call ncdimc('qcomp',iqq,0)
      call ncdimc('ucomp',iuu,0)
      call ncdimc('vcomp',ivv,0)
#else 
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif 
      call ncdims('k2',2)
      call ncdims('k3',3)
      call ncdims('kk',kk)
      call ncdims('kk2',2*kk)
      call nctime(datenum,calendar,timeunits,startdate)
c
c --- output model fields to restart file
      call ncputi('ntda',ntda,1)
c
#ifdef COMPRES 
      call nccomp('u','ucomp kk2 time',u,iuu)
      call nccomp('v','vcomp kk2 time',v,ivv)
      call nccomp('dp','pcomp kk2 time',dp,ip)
      call nccomp('temp','pcomp kk2 time',temp,ip)
      call nccomp('saln','pcomp kk2 time',saln,ip)
      call nccomp('sigma','pcomp kk2 time',sigma,ip)
      call nccomp('pgfx','ucomp kk2 time',pgfx,iuu)
      call nccomp('pgfxo','ucomp kk time',pgfxo,iuu)
      call nccomp('pgfy','vcomp kk2 time',pgfy,ivv)
      call nccomp('pgfyo','vcomp kk time',pgfyo,ivv)
      call nccomp('pb','pcomp k3 time',pb,ip)
      call nccomp('pb_p','pcomp time',pb_p,ip)
      call nccomp('pbu','ucomp k2 time',pbu,iuu)
      call nccomp('pbv','vcomp k2 time',pbv,ivv)
      call nccomp('pbu_p','ucomp time',pbu_p,iuu)
      call nccomp('pbv_p','vcomp time',pbv_p,ivv)
      call nccomp('ub','ucomp k2 time',ub,iuu)
      call nccomp('vb','vcomp k2 time',vb,ivv)
      call nccomp('ubflx','ucomp k3 time',ubflx,iuu)
      call nccomp('vbflx','vcomp k3 time',vbflx,ivv)
      call nccomp('ubflxs','ucomp k2 time',ubflxs,iuu)
      call nccomp('vbflxs','vcomp k2 time',vbflxs,ivv)
      call nccomp('ubflxs_p','ucomp time',ubflxs_p,iuu)
      call nccomp('vbflxs_p','vcomp time',vbflxs_p,ivv)
      call nccomp('ubcors','ucomp k2 time',ubcors,iuu)
      call nccomp('vbcors','vcomp k2 time',vbcors,ivv)
      call nccomp('ubcors_p','ucomp time',ubcors_p,iuu)
      call nccomp('vbcors_p','vcomp time',vbcors_p,ivv)
      call nccomp('pvtrop','qcomp k2 time',pvtrop,iqq)
      call nccomp('pgfxm','ucomp k2 time',pgfxm,iuu)
      call nccomp('pgfym','vcomp k2 time',pgfym,ivv)
      call nccomp('xixp','ucomp k2 time',xixp,iuu)
      call nccomp('xixm','ucomp k2 time',xixm,iuu)
      call nccomp('xiyp','vcomp k2 time',xiyp,ivv)
      call nccomp('xiym','vcomp k2 time',xiym,ivv)
      call nccomp('tmxold','pcomp k2 time',tmxold,ip)
      call nccomp('smxold','pcomp k2 time',smxold,ip)
      call nccomp('phi','pcomp time',phi(1-nbdy,1-nbdy,kk+1),ip)
      call nccomp('sealv','pcomp time',sealv,ip)
      call nccomp('hicem','pcomp time',hicem,ip)
      call nccomp('ficem','pcomp time',ficem,ip)
      call nccomp('tsrfm','pcomp time',tsrfm,ip)
      call nccomp('hsnwm','pcomp time',hsnwm,ip)
      call nccomp('ticem','pcomp time',ticem,ip)
      call nccomp('iagem','pcomp time',iagem,ip)
      call nccomp('tsi_tda','pcomp time',tsi_tda,ip)
      call nccomp('tml_tda','pcomp time',tml_tda,ip)
      call nccomp('sml_tda','pcomp time',sml_tda,ip)
      call nccomp('alb_tda','pcomp time',alb_tda,ip)
      call nccomp('fice_tda','pcomp time',fice_tda,ip)
      call nccomp('rnfres','pcomp time',rnfres,ip)
#else
      call ncwrtr('u','x y kk2 time',u,iuu,1)
      call ncwrtr('v','x y kk2 time',v,ivv,1)
      call ncwrtr('dp','x y kk2 time',dp,ip,1)
      call ncwrtr('temp','x y kk2 time',temp,ip,1)
      call ncwrtr('saln','x y kk2 time',saln,ip,1)
      call ncwrtr('sigma','x y kk2 time',sigma,ip,1)
      call ncwrtr('pgfx','x y kk time',pgfx,iuu,1)
      call ncwrtr('pgfy','x y kk time',pgfy,ivv,1)
      call ncwrtr('pb','x y k3 time',pb,ip,1)
      call ncwrtr('pb_p','x y time',pb_p,ip,1)
      call ncwrtr('pbu','x y k2 time',pbu,iuu,1)
      call ncwrtr('pbv','x y k2 time',pbv,ivv,1)
      call ncwrtr('pbu_p','x y time',pbu_p,iuu,1)
      call ncwrtr('pbv_p','x y time',pbv_p,ivv,1)
      call ncwrtr('ub','x y k2 time',ub,iuu,1)
      call ncwrtr('vb','x y k2 time',vb,ivv,1)
      call ncwrtr('ubflx','x y k3 time',ubflx,iuu,1)
      call ncwrtr('vbflx','x y k3 time',vbflx,ivv,1)
      call ncwrtr('ubflxs','x y k2 time',ubflxs,iuu,1)
      call ncwrtr('vbflxs','x y k2 time',vbflxs,ivv,1)
      call ncwrtr('ubflxs_p','x y time',ubflxs_p,iuu,1)
      call ncwrtr('vbflxs_p','x y time',vbflxs_p,ivv,1)
      call ncwrtr('ubcors','x y k2 time',ubcors,iuu,1)
      call ncwrtr('vbcors','x y k2 time',vbcors,ivv,1)
      call ncwrtr('ubcors_p','x y time',ubcors_p,iuu,1)
      call ncwrtr('vbcors_p','x y time',vbcors_p,ivv,1)
      call ncwrtr('pvtrop','x y k2 time',pvtrop,iqq,1)
      call ncwrtr('pgfxm','x y k2 time',pgfxm,iuu,1)
      call ncwrtr('pgfym','x y k2 time',pgfym,ivv,1)
      call ncwrtr('xixp','x y k2 time',xixp,iuu,1)
      call ncwrtr('xixm','x y k2 time',xixm,iuu,1)
      call ncwrtr('xiyp','x y k2 time',xiyp,ivv,1)
      call ncwrtr('xiym','x y k2 time',xiym,ivv,1)
      call ncwrtr('tmxold','x y k2 time',tmxold,ip,1)
      call ncwrtr('smxold','x y k2 time',smxold,ip,1)
      call ncwrtr('phi','x y time',phi(1-nbdy,1-nbdy,kk+1),ip,1)
      call ncwrtr('sealv','x y time',sealv,ip,1)
      call ncwrtr('hicem','x y time',hicem,ip,1)
      call ncwrtr('ficem','x y time',ficem,ip,1)
      call ncwrtr('tsrfm','x y time',tsrfm,ip,1)
      call ncwrtr('hsnwm','x y time',hsnwm,ip,1)
      call ncwrtr('ticem','x y time',ticem,ip,1)
      call ncwrtr('iagem','x y time',iagem,ip,1)
      call ncwrtr('tsi_tda','x y time',tsi_tda,ip,1)
      call ncwrtr('tml_tda','x y time',tml_tda,ip,1)
      call ncwrtr('sml_tda','x y time',sml_tda,ip,1)
      call ncwrtr('alb_tda','x y time',alb_tda,ip,1)
      call ncwrtr('fice_tda','x y time',fice_tda,ip,1)
      call ncwrtr('rnfres','x y time',rnfres,ip,1)
#endif 
c
#if defined(CLIM) || defined(SYN)
#  ifdef COMPRES 
      call nccomp('cd_d','pcomp time',cd_d,ip)
      call nccomp('ch_d','pcomp time',ch_d,ip)
      call nccomp('ce_d','pcomp time',ce_d,ip)
      call nccomp('wg2_d','pcomp time',wg2_d,ip)
      call nccomp('cd_m','pcomp time',cd_m,ip)
      call nccomp('ch_m','pcomp time',ch_m,ip)
      call nccomp('ce_m','pcomp time',ce_m,ip)
      call nccomp('wg2_m','pcomp time',wg2_m,ip)
      call nccomp('rhoa','pcomp time',rhoa,ip)
#  else
      call ncwrtr('cd_d','x y time',cd_d,ip,1)
      call ncwrtr('ch_d','x y time',ch_d,ip,1)
      call ncwrtr('ce_d','x y time',ce_d,ip,1)
      call ncwrtr('wg2_d','x y time',wg2_d,ip,1)
      call ncwrtr('cd_m','x y time',cd_m,ip,1)
      call ncwrtr('ch_m','x y time',ch_m,ip,1)
      call ncwrtr('ce_m','x y time',ce_m,ip,1)
      call ncwrtr('wg2_m','x y time',wg2_m,ip,1)
      call ncwrtr('rhoa','x y time',rhoa,ip,1)
#  endif
#endif
c
      call ncfcls
c
#ifdef ICEDYN
      call restart_icewt(runid,runid_len,time,rstfnm_ext,
     .                   path2,path2_len)
#endif
c
#ifdef TRC
      call restart_trcwt(rstfnm_ext)
#endif
c
c
#ifdef CCSMCOUPLED
      if (diag2dmon) then
c
c --- write accumulated 2d fields restart file for monthly output
      fnm=path2(1:path2_len)//runid(1:runid_len)//'_diaacc2d'//
     .    rstfnm_ext
      if (mnproc.eq.1) write (lp,'(a,a)')
     .  ' saving accumulated 2d fields restart file ',fnm
      call ncfopn(fnm,'w')
c
#  ifdef COMPRES 
      call ncdimc('pcomp',ip,0)
      call ncdimc('ucomp',iuu,0)
      call ncdimc('vcomp',ivv,0)
#  else 
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#  endif 
      call nctime(datenum,calendar,timeunits,startdate)
c
      call ncputr('time',time,1)
      call ncputi('nacc_phy2d',nacc_phy2d,1)
c
#  ifdef COMPRES 
      call nccomp('ub_dia','ucomp time',ub_dia,iuu)
      call nccomp('vb_dia','vcomp time',vb_dia,ivv)
      call nccomp('ztx_dia','ucomp time',ztx_dia,iuu)
      call nccomp('mty_dia','vcomp time',mty_dia,ivv)
      call nccomp('taux_dia','ucomp time',taux_dia,iuu)
      call nccomp('tauy_dia','vcomp time',tauy_dia,ivv)
      call nccomp('mxlu_dia','ucomp time',mxlu_dia,iuu)
      call nccomp('mxlv_dia','vcomp time',mxlv_dia,ivv)
      call nccomp('mxldu_dia','ucomp time',mxldu_dia,iuu)
      call nccomp('mxldv_dia','vcomp time',mxldv_dia,ivv)
      call nccomp('sealv_dia','pcomp time',sealv_dia,ip)
      call nccomp('sigmx_dia','pcomp time',sigmx_dia,ip)
      call nccomp('hice_dia','pcomp time',hice_dia,ip)
      call nccomp('hsnw_dia','pcomp time',hsnw_dia,ip)
      call nccomp('fice_dia','pcomp time',fice_dia,ip)
      call nccomp('tsrf_dia','pcomp time',tsrf_dia,ip)
      call nccomp('tice_dia','pcomp time',tice_dia,ip)
      call nccomp('swa_dia','pcomp time',swa_dia,ip)
      call nccomp('nsf_dia','pcomp time',nsf_dia,ip)
      call nccomp('dfl_dia','pcomp time',dfl_dia,ip)
      call nccomp('lip_dia','pcomp time',lip_dia,ip)
      call nccomp('sop_dia','pcomp time',sop_dia,ip)
      call nccomp('eva_dia','pcomp time',eva_dia,ip)
      call nccomp('alb_dia','pcomp time',alb_dia,ip)
      call nccomp('rnfflx_dia','pcomp time',rnfflx_dia,ip)
      call nccomp('ustar_dia','pcomp time',ustar_dia,ip)
      call nccomp('sss_dia','pcomp time',sss_dia,ip)
      call nccomp('sst_dia','pcomp time',sst_dia,ip)
      call nccomp('mxld_dia','pcomp time',mxld_dia,ip)
      call nccomp('surflx_dia','pcomp time',surflx_dia,ip)
      call nccomp('salflx_dia','pcomp time',salflx_dia,ip)
      call nccomp('uice_dia','ucomp time',uice_dia,iuu)
      call nccomp('vice_dia','vcomp time',vice_dia,ivv)
      call nccomp('convec_dia','pcomp time',convec_dia,ip)
#  else 
      call ncwrtr('ub_dia','x y time',ub_dia,iuu,1)
      call ncwrtr('vb_dia','x y time',vb_dia,ivv,1)
      call ncwrtr('ztx_dia','x y time',ztx_dia,iuu,1)
      call ncwrtr('mty_dia','x y time',mty_dia,ivv,1)
      call ncwrtr('taux_dia','x y time',taux_dia,iuu,1)
      call ncwrtr('tauy_dia','x y time',tauy_dia,ivv,1)
      call ncwrtr('mxlu_dia','x y time',mxlu_dia,iuu,1)
      call ncwrtr('mxlv_dia','x y time',mxlv_dia,ivv,1)
      call ncwrtr('mxldu_dia','x y time',mxldu_dia,iuu,1)
      call ncwrtr('mxldv_dia','x y time',mxldv_dia,ivv,1)
      call ncwrtr('sealv_dia','x y time',sealv_dia,ip,1)
      call ncwrtr('sigmx_dia','x y time',sigmx_dia,ip,1)
      call ncwrtr('hice_dia','x y time',hice_dia,ip,1)
      call ncwrtr('hsnw_dia','x y time',hsnw_dia,ip,1)
      call ncwrtr('fice_dia','x y time',fice_dia,ip,1)
      call ncwrtr('tsrf_dia','x y time',tsrf_dia,ip,1)
      call ncwrtr('tice_dia','x y time',tice_dia,ip,1)
      call ncwrtr('swa_dia','x y time',swa_dia,ip,1)
      call ncwrtr('nsf_dia','x y time',nsf_dia,ip,1)
      call ncwrtr('dfl_dia','x y time',dfl_dia,ip,1)
      call ncwrtr('lip_dia','x y time',lip_dia,ip,1)
      call ncwrtr('sop_dia','x y time',sop_dia,ip,1)
      call ncwrtr('eva_dia','x y time',eva_dia,ip,1)
      call ncwrtr('alb_dia','x y time',alb_dia,ip,1)
      call ncwrtr('rnfflx_dia','x y time',rnfflx_dia,ip,1)
      call ncwrtr('ustar_dia','x y time',ustar_dia,ip,1)
      call ncwrtr('sss_dia','x y time',sss_dia,ip,1)
      call ncwrtr('sst_dia','x y time',sst_dia,ip,1)
      call ncwrtr('mxld_dia','x y time',mxld_dia,ip,1)
      call ncwrtr('surflx_dia','x y time',surflx_dia,ip,1)
      call ncwrtr('salflx_dia','x y time',salflx_dia,ip,1)
      call ncwrtr('uice_dia','x y time',uice_dia,iuu,1)
      call ncwrtr('vice_dia','x y time',vice_dia,ivv,1)
      call ncwrtr('convec_dia','x y time',convec_dia,ip,1)
#  endif
c
      call ncfcls       
      endif
c
      if (diag3dmon) then
c
c --- write accumulated 3d fields restart file for monthly output
      fnm=path2(1:path2_len)//runid(1:runid_len)//'_diaacc3d'//
     .    rstfnm_ext
      if (mnproc.eq.1) write (lp,'(a,a)')
     .  ' saving accumulated 3d fields restart file ',fnm
      call ncfopn(fnm,'w')
c
#  ifdef COMPRES
      call ncdimc('pcomp',ip,0)
      call ncdimc('ucomp',iuu,0)
      call ncdimc('vcomp',ivv,0)
#  else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#  endif 
      call ncdims('z',kk)
      call nctime(datenum,calendar,timeunits,startdate)
c
      call ncputr('time',time,1)
      call ncputi('nacc_phy3d',nacc_phy3d,1)
c
#  ifdef COMPRES
      call nccomp('uvel_dia','ucomp z time',uvel_dia,iuu)
      call nccomp('vvel_dia','vcomp z time',vvel_dia,ivv)
      call nccomp('dpu_dia','ucomp z time',dpu_dia,iuu)
      call nccomp('dpv_dia','vcomp z time',dpv_dia,ivv)
      call nccomp('uflx_dia','ucomp z time',uflx_dia,iuu)
      call nccomp('vflx_dia','vcomp z time',vflx_dia,ivv)
      call nccomp('saln_dia','pcomp z time',saln_dia,ip)
      call nccomp('temp_dia','pcomp z time',temp_dia,ip)
      call nccomp('dp_dia','pcomp z time',dp_dia,ip)
      call nccomp('dz_dia','pcomp z time',dz_dia,ip)
      call nccomp('diaflx','pcomp z time',diaflx,ip)
#  else
      call ncwrtr('uvel_dia','x y z time',uvel_dia,iuu,1)
      call ncwrtr('vvel_dia','x y z time',vvel_dia,ivv,1)
      call ncwrtr('dpu_dia','x y z time',dpu_dia,iuu,1)
      call ncwrtr('dpv_dia','x y z time',dpv_dia,ivv,1)
      call ncwrtr('uflx_dia','x y z time',uflx_dia,iuu,1)
      call ncwrtr('vflx_dia','x y z time',vflx_dia,ivv,1)
      call ncwrtr('saln_dia','x y z time',saln_dia,ip,1)
      call ncwrtr('temp_dia','x y z time',temp_dia,ip,1)
      call ncwrtr('dp_dia','x y z time',dp_dia,ip,1)
      call ncwrtr('dz_dia','x y z time',dz_dia,ip,1)
      call ncwrtr('diaflx','x y z time',diaflx,ip,1)
#  endif 
c
      call ncfcls       
      endif
#endif
c
      if (ditflx) then
c
c --- write diag. heat flux restart file
      fnm=path2(1:path2_len)//runid(1:runid_len)//'_tflx'//rstfnm_ext
      if (mnproc.eq.1) write (lp,'(a,a)')
     .  ' saving diag. heat flux restart file ',fnm
      call ncfopn(fnm,'w')
c
#ifdef COMPRES
      call ncdimc('pcomp',ip,0)
#else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif 
      call ncdims('week',48)
c
      call ncputr('time',time,1)
      call ncputi('nflxdi',nflxdi,48)
c
#ifdef COMPRES
      call nccomp('tflxdi','pcomp week',tflxdi,ip)
#else
      call ncwrtr('tflxdi','x y week',tflxdi,ip,1)
#endif 
c       
      call ncfcls   
      endif  
c
      if (disflx) then
c
c --- write diag. salt flux restart file
      fnm=path2(1:path2_len)//runid(1:runid_len)//'_sflx'//rstfnm_ext
      if (mnproc.eq.1) write (lp,'(a,a)')
     .  ' saving diag. salt flux restart file ',fnm
      call ncfopn(fnm,'w')
c
#ifdef COMPRES
      call ncdimc('pcomp',ip,0)
#else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif 
      call ncdims('week',48)
c
      call ncputr('time',time,1)
      call ncputi('nflxdi',nflxdi,48)
c
#ifdef COMPRES
      call nccomp('sflxdi','pcomp week',sflxdi,ip)
#else
      call ncwrtr('sflxdi','x y week',sflxdi,ip,1)
#endif 
c       
      call ncfcls
      endif
c
#ifdef CCSMCOUPLED
c --- write restart filename to rpointer.ocn
      if (mnproc.eq.1) then
        open (unit=12,file='rpointer.ocn')
        write (12,'(a)') rstfnm
        close (unit=12)
      endif
#endif
c
      return
      end
