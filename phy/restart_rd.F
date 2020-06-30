! ------------------------------------------------------------------------------
! Copyright (C) 2006-2020 Mats Bentsen, Mehmet Ilicak, Alok Kumar Gupta,
!                         Jerry Tjiputra
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

      subroutine restart_rd
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia
      use mod_ben02, only: cd_d,ch_d,ce_d,wg2_d,
     .                     cd_m,ch_m,ce_m,wg2_m,rhoa,
     .                     tsi_tda,tml_tda,sml_tda,alb_tda,fice_tda,ntda
      use mod_thdysi, only: tsrfm,ticem
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_cesm.h"
#include "common_seaice.h"
c
      integer i,j,k,n,ndayr,nmonthr,nyearr,l
      character(len=256) :: rstfnm,fnm
      character(len=2) :: c2
c
      integer daydif
      external daydif
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv,iqq
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
      logical first,fexist,vexist
      data first /.true./
      save first,iuu,ivv,iqq
c
#if defined(TRC) && (defined(TKE) || defined(IDLAGE))
#  include "param_trc.h"
#  include "common_trc.h"
#  ifdef TKE
#    include "common_tke.h"
#  endif
#endif
c
c --- open restart file and adjust time information if needed
c
      if     (nday1+nint(time0).eq.0) then
c
c --- - open restart file for initial conditions and adjust integration
c --- - time corresponding to start date
        rstfnm=icfile
        if (mnproc.eq.1) inquire(file=rstfnm,exist=fexist)
        call xcbcst(fexist)
        if (fexist) then
          call ncfopn(rstfnm,'r',' ',1,iotype)
          call ncgetr('time',time)
          time0=time
        else
          if (mnproc.eq.1) then
            write (lp,*)
     .'restart_rd: could not find restart file for initial conditions!'
          endif
          call xcstop('(restart_rd)')
                 stop '(restart_rd)'
        endif
c
      elseif (expcnf.eq.'cesm') then
c
c --- - get restart file name from rpointer.ocn
        if (mnproc.eq.1) inquire(file='rpointer.ocn',exist=fexist)
        call xcbcst(fexist)
        if (fexist) then
          if (mnproc.eq.1) then
            open (unit=nfu,file='rpointer.ocn')
            read (nfu,'(a)') rstfnm
            close (unit=nfu)
          endif
          call xcbcst(rstfnm)
        else
          if (mnproc.eq.1) then
            write (lp,*) 'restart_rd: could not find file rpointer.ocn!'
          endif
          call xcstop('(restart_rd)')
                 stop '(restart_rd)'
        endif
c
c --- - open restart file
        if (mnproc.eq.1) inquire(file=rstfnm,exist=fexist)
        call xcbcst(fexist)
        if (.not.fexist) then
          if (mnproc.eq.1) then
            write (lp,*) 'restart_rd: could not find restart file!'
          endif
          call xcstop('(restart_rd)')
                 stop '(restart_rd)'
        endif
        call ncfopn(rstfnm,'r',' ',1,iotype)
c
      else
c
c --- - first try file name:
c --- -  <experiment name>_restphy_<year>.<month>.<day>_<integration day>.nc
        write (rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')
     .    runid(1:runid_len),'_restphy_',nyear,'.',nmonth,'.',nday,
     .    '_',nday1,'.nc'
c
        if (mnproc.eq.1) inquire(file=rstfnm,exist=fexist)
        call xcbcst(fexist)
        if (fexist) then
          call ncfopn(rstfnm,'r',' ',1,iotype)
          call ncgeti('nday0',ndayr)
          call ncgeti('nmonth0',nmonthr)
          call ncgeti('nyear0',nyearr)
          call ncgetr('time0',time0)
          call ncgetr('time',time)
          if (nyearr.ne.nyear0.or.nmonthr.ne.nmonth0.or.
     .        ndayr.ne.nday0) then
            if (mnproc.eq.1) then
              write (lp,'(2a)')
     .' restart_rd: expected identical initial experiment date in',
     .' namelist and restart but found:'
              write (lp,'(a,i4.4,a,i2.2,a,i2.2)')
     .          ' restart_rd: initial date namelist: ',
     .          nyear0,'.',nmonth0,'.',nday0
              write (lp,'(a,i4.4,a,i2.2,a,i2.2)')
     .          ' restart_rd: initial date restart: ',
     .          nyearr,'.',nmonthr,'.',ndayr
            endif
            call xcstop('(restart_rd)')
                   stop '(restart_rd)'
          endif
c
        else
c
c --- --- then try automatic selection of file with consistent
c --- --- integration day and date among files named:
c --- ---   <experiment name>_restphy_1.nc
c --- ---   <experiment name>_restphy_2.nc
c --- ---   <experiment name>_restphy_3.nc
          do i=1,4
            write (rstfnm,'(2a,i1,a)')
     .        runid(1:runid_len),'_restphy_',i,'.nc'
            inquire(file=rstfnm,exist=fexist)
            if (fexist) then
              call ncfopn(rstfnm,'r',' ',1,iotype)
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
            if (mnproc.eq.1) then
              write (lp,*)
     .          'restart_rd: Could not find proper restart file!'
            endif
            call xcstop('(restart_rd)')
                   stop '(restart_rd)' 
          endif 
        endif 
c
      endif 
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' restart_rd: reading restart file ',
     .                    trim(rstfnm)
      endif
c
c --- Compute extended uv masks 
      if (first) then
        first=.false.
c$OMP PARALLEL DO PRIVATE(i)
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
c$OMP END PARALLEL DO
      endif
c
      call ncread('u',u,iuu,1,0.)
      call ncread('v',v,ivv,1,0.)
      call ncread('dp',dp,ip,1,0.)
      call ncread('dpold',dpold,ip,1,0.)
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
      vexist=ncinqv('umfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('umfltd',umfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''umfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (iuu(i,j).eq.1) then
                umfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('utfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('utfltd',utfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''utfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (iuu(i,j).eq.1) then
                utfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('utflld')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('utflld',utflld,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''utflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (iuu(i,j).eq.1) then
                utflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('usfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('usfltd',usfltd,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''usfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (iuu(i,j).eq.1) then
                usfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('usflld')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('usflld',usflld,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''usflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (iuu(i,j).eq.1) then
                usflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      call ncread('vflx',vflx,ivv,1,0.)
      call ncread('vtflx',vtflx,ivv,1,0.)
      call ncread('vsflx',vsflx,ivv,1,0.)
c
      vexist=ncinqv('vmfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vmfltd',vmfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vmfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vmfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vtfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vtfltd',vtfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vtfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vtfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vtflld')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vtflld',vtflld,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vtflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vtflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vsfltd')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vsfltd',vsfltd,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vsfltd'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vsfltd(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vsflld')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vsflld',vsflld,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vsflld'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1,jj
          do k=1,2*kk
            do i=1,ii
              if (ivv(i,j).eq.1) then
                vsflld(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
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
      call ncread('ficem',ficem,ip,1,0.)
c
      vexist=ncinqv('uml')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('uml',uml,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''uml'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,4
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              uml(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vml')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vml',vml,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vml'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj+1
          do k=1,4
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vml(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('umlres')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('umlres',umlres,iuu,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''umlres'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,2
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              umlres(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('vmlres')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('vmlres',vmlres,ivv,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: field ''vmlres'' is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj+1
          do k=1,2
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vmlres(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sprfac) then
        vexist=ncinqa('prfac')
!        call xcbcst(vexist)
        if (vexist) then
          call ncgetr('prfac',prfac)
          call xcbcst(prfac)
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
      endif
c
      if (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call ncread('cd_d',cd_d,ip,1,0.)
        call ncread('ch_d',ch_d,ip,1,0.)
        call ncread('ce_d',ce_d,ip,1,0.)
        call ncread('wg2_d',wg2_d,ip,1,0.)
        call ncread('cd_m',cd_m,ip,1,0.)
        call ncread('ch_m',ch_m,ip,1,0.)
        call ncread('ce_m',ce_m,ip,1,0.)
        call ncread('wg2_m',wg2_m,ip,1,0.)
        call ncread('rhoa',rhoa,ip,1,0.)
        call ncread('tsi_tda',tsi_tda,ip,1,0.)
        call ncread('tml_tda',tml_tda,ip,1,0.)
        call ncread('sml_tda',sml_tda,ip,1,0.)
        call ncread('alb_tda',alb_tda,ip,1,0.)
        call ncread('fice_tda',fice_tda,ip,1,0.)
c
        call ncgeti('ntda',ntda)
        call xcbcst(ntda)
c
        call ncread('hicem',hicem,ip,1,0.)
        call ncread('tsrfm',tsrfm,ip,1,0.)
        call ncread('hsnwm',hsnwm,ip,1,0.)
        call ncread('ticem',ticem,ip,1,0.)
        call ncread('iagem',iagem,ip,1,0.)
        call ncread('rnfres',rnfres,ip,1,0.)
      endif
c
      if (expcnf.eq.'cesm') then
        vexist=ncinqv('ustarw_da')
!        call xcbcst(vexist)
        if (vexist) then
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
          call ncgeti('l2ci',l2ci)
          call xcbcst(l2ci)
          l1ci=3-l2ci
        else
          if (mnproc.eq.1) then
            write (lp,*)
     .'Warning: time levels for interpolation of forcing fields is not'
            write (lp,*)
     .       'read from restart file.'
          endif
          l1ci=1
          l2ci=1
        endif
        call ncread('frzpot',frzpot,ip,1,0.)
        call ncread('mltpot',mltpot,ip,1,0.)
        vexist=ncinqv('flxco2')
!        call xcbcst(vexist)
        if (vexist) then
          call ncread('flxco2',flxco2,ip,1,0.)
        else
          if (mnproc.eq.1) then
            write (lp,*)
     .'Warning: Air-sea CO2 flux is not read from restart file and'
            write (lp,*)
     .       'will be initialized to zero.'
          endif
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              flxco2(i,j)=0.
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        vexist=ncinqv('flxdms')
!        call xcbcst(vexist)
        if (vexist) then
          call ncread('flxdms',flxdms,ip,1,0.)
        else
          if (mnproc.eq.1) then
            write (lp,*)
     .'Warning: DMS flux is not read from restart file and'
            write (lp,*)
     .       'will be initialized to zero.'
          endif
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              flxdms(i,j)=0.
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
      endif
c
#ifdef TRC
#  ifdef TKE
      vexist=ncinqv('tke')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('tke',trc(1-nbdy,1-nbdy,1,itrtke),ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: TKE is not read from restart file and'
          write (lp,*)
     .     'will be initialized to tke_min.'
        endif
      endif
c
      vexist=ncinqv('gls_psi')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('gls_psi',trc(1-nbdy,1-nbdy,1,itrgls),ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: gls_psi is not read from restart file and'
          write (lp,*)
     .     'will be initialized to gls_psi_min.'
        endif
      endif
c
      vexist=ncinqv('L_scale')
!      call xcbcst(vexist)
      if (vexist) then
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
      vexist=ncinqv('difdia')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('difdia',difdia,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: difdia is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              difdia(i,j,k)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      vexist=ncinqv('ustarb')
!      call xcbcst(vexist)
      if (vexist) then
        call ncread('ustarb',ustarb,ip,1,0.)
      else
        if (mnproc.eq.1) then
          write (lp,*)
     .     'Warning: ustarb is not read from restart file and'
          write (lp,*)
     .     'will be initialized to zero.'
        endif
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            ustarb(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
#  endif
#  ifdef IDLAGE
      vexist=ncinqv('idlage')
!      call xcbcst(vexist)
      if (vexist) then
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
#  endif
#endif
c
c --- read accumulated fields
      if (nstep1.gt.nstep0) then
        do n=1,nphy        
          write(c2,'(i2.2)') n 
            if (ncinqa('nacc_phy'//c2)) then 
              call ncgeti('nacc_phy'//c2,nacc_phy(n))
            else
              nacc_phy(n)=0
            endif 
          call xcbcst(nacc_phy(n))
          if (nacc_phy(n).gt.0) then
            if (ACC_UB(n).ne.0) call ncread('ub_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_UB(n)),iuu,1,0.)
            if (ACC_VB(n).ne.0) call ncread('vb_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_VB(n)),ivv,1,0.)
            if (ACC_UBFLXS(n).ne.0) call ncread('ubflxs_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(n)),iuu,1,0.)
            if (ACC_VBFLXS(n).ne.0) call ncread('vbflxs_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(n)),ivv,1,0.)
            if (ACC_ZTX(n).ne.0) call ncread('ztx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_ZTX(n)),iuu,1,0.)
            if (ACC_MTY(n).ne.0) call ncread('mty_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTY(n)),ivv,1,0.)
            if (ACC_TAUX(n).ne.0) call ncread('taux_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_TAUX(n)),iuu,1,0.)
            if (ACC_TAUY(n).ne.0) call ncread('tauy_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_TAUY(n)),ivv,1,0.)
            if (ACC_MXLU(n).ne.0) call ncread('mxlu_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MXLU(n)),iuu,1,0.)
            if (ACC_MXLV(n).ne.0) call ncread('mxlv_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MXLV(n)),ivv,1,0.)
            if (ACC_MLDU(n).ne.0) call ncread('mldu_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLDU(n)),iuu,1,0.)
            if (ACC_MLDV(n).ne.0) call ncread('mldv_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLDV(n)),ivv,1,0.)
            if (ACC_UICE(n).ne.0) call ncread('uice_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_UICE(n)),iuu,1,0.)
            if (ACC_VICE(n).ne.0) call ncread('vice_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_VICE(n)),ivv,1,0.)
            if (ACC_IVOLU(n).ne.0) call ncread('ivolu_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_IVOLU(n)),iuu,1,0.)
            if (ACC_IVOLV(n).ne.0) call ncread('ivolv_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_IVOLV(n)),ivv,1,0.)
            if (ACC_PSRF(n).ne.0) call ncread('psrf_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_PSRF(n)),ip,1,0.)
            if (ACC_PBOT(n).ne.0) call ncread('pbot_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_PBOT(n)),ip,1,0.)
            if (ACC_SEALV(n).ne.0) call ncread('sealv_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SEALV(n)),ip,1,0.)
            if (ACC_SLVSQ(n).ne.0) call ncread('slvsq_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SLVSQ(n)),ip,1,0.)
            if (ACC_SSS(n).ne.0) call ncread('sss_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SSS(n)),ip,1,0.)
            if (ACC_SSSSQ(n).ne.0) call ncread('ssssq_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SSSSQ(n)),ip,1,0.)
            if (ACC_SBOT(n).ne.0) call ncread('sbot_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SBOT(n)),ip,1,0.)
            if (ACC_SST(n).ne.0) call ncread('sst_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SST(n)),ip,1,0.)
            if (ACC_SSTSQ(n).ne.0) call ncread('sstsq_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SSTSQ(n)),ip,1,0.)
            if (ACC_TBOT(n).ne.0) call ncread('tbot_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_TBOT(n)),ip,1,0.)
            if (ACC_SIGMX(n).ne.0) call ncread('sigmx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SIGMX(n)),ip,1,0.)
            if (ACC_MLD(n).ne.0) call ncread('mld_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLD(n)),ip,1,0.)
            if (ACC_MAXMLD(n).ne.0) call ncread('maxmld_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MAXMLD(n)),ip,1,0.)
            if (ACC_MLTS(n).ne.0) call ncread('mlts_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLTS(n)),ip,1,0.)
            if (ACC_MLTSMN(n).ne.0) call ncread('mltsmn_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLTSMN(n)),ip,1,0.)
            if (ACC_MLTSMX(n).ne.0) call ncread('mltsmx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLTSMX(n)),ip,1,0.)
            if (ACC_MLTSSQ(n).ne.0) call ncread('mltssq_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MLTSSQ(n)),ip,1,0.)
            if (ACC_T20D(n).ne.0) call ncread('t20d_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_T20D(n)),ip,1,0.)
            if (ACC_ALB(n).ne.0) call ncread('alb_phy'//c2, 
     .        phyh2d(1-nbdy,1-nbdy,ACC_ALB(n)),ip,1,0.)
            if (ACC_SWA(n).ne.0) call ncread('swa_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SWA(n)),ip,1,0.)
            if (ACC_NSF(n).ne.0) call ncread('nsf_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_NSF(n)),ip,1,0.)
            if (ACC_DFL(n).ne.0) call ncread('dfl_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_DFL(n)),ip,1,0.)
            if (ACC_LIP(n).ne.0) call ncread('lip_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_LIP(n)),ip,1,0.)
            if (ACC_SOP(n).ne.0) call ncread('sop_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SOP(n)),ip,1,0.)
            if (ACC_EVA(n).ne.0) call ncread('eva_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_EVA(n)),ip,1,0.)
            if (ACC_RNFFLX(n).ne.0) call ncread('rnfflx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_RNFFLX(n)),ip,1,0.)
            if (ACC_RFIFLX(n).ne.0) call ncread('rfiflx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_RFIFLX(n)),ip,1,0.)
            if (ACC_SFL(n).ne.0) call ncread('sfl_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SFL(n)),ip,1,0.)
            if (ACC_BRNFLX(n).ne.0) call ncread('brnflx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_BRNFLX(n)),ip,1,0.)
            if (ACC_BRNPD(n).ne.0) call ncread('brnpd_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_BRNPD(n)),ip,1,0.)
            if (ACC_SURFLX(n).ne.0) call ncread('surflx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SURFLX(n)),ip,1,0.)
            if (ACC_SURRLX(n).ne.0) call ncread('surrlx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SURRLX(n)),ip,1,0.)
            if (ACC_SALFLX(n).ne.0) call ncread('salflx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SALFLX(n)),ip,1,0.)
            if (ACC_SALRLX(n).ne.0) call ncread('salrlx_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_SALRLX(n)),ip,1,0.)
            if (ACC_ABSWND(n).ne.0) call ncread('abswnd_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_ABSWND(n)),ip,1,0.)
            if (ACC_USTAR(n).ne.0) call ncread('ustar_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_USTAR(n)),ip,1,0.)
            if (ACC_USTAR3(n).ne.0) call ncread('ustar3_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_USTAR3(n)),ip,1,0.)
            if (ACC_IDKEDT(n).ne.0) call ncread('idkedt_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_IDKEDT(n)),ip,1,0.)
            if (ACC_MTKEUS(n).ne.0) call ncread('mtkeus_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKEUS(n)),ip,1,0.)
            if (ACC_MTKENI(n).ne.0) call ncread('mtkeni_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKENI(n)),ip,1,0.)
            if (ACC_MTKEBF(n).ne.0) call ncread('mtkebf_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKEBF(n)),ip,1,0.)
            if (ACC_MTKERS(n).ne.0) call ncread('mtkers_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKERS(n)),ip,1,0.)
            if (ACC_MTKEPE(n).ne.0) call ncread('mtkepe_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKEPE(n)),ip,1,0.)
            if (ACC_MTKEKE(n).ne.0) call ncread('mtkeke_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_MTKEKE(n)),ip,1,0.)
            if (ACC_FMLTFZ(n).ne.0) call ncread('fmltfz_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_FMLTFZ(n)),ip,1,0.)
            if (ACC_HMLTFZ(n).ne.0) call ncread('hmltfz_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_HMLTFZ(n)),ip,1,0.)
            if (ACC_HICE(n).ne.0) call ncread('hice_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_HICE(n)),ip,1,0.)
            if (ACC_HSNW(n).ne.0) call ncread('hsnw_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_HSNW(n)),ip,1,0.)
            if (ACC_FICE(n).ne.0) call ncread('fice_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_FICE(n)),ip,1,0.)
            if (ACC_TSRF(n).ne.0) call ncread('tsrf_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_TSRF(n)),ip,1,0.)
            if (ACC_TICE(n).ne.0) call ncread('tice_phy'//c2,
     .        phyh2d(1-nbdy,1-nbdy,ACC_TICE(n)),ip,1,0.)
            if (ACC_UVEL(n).ne.0) call ncread('uvel_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)),iuu,1,0.)
            if (ACC_VVEL(n).ne.0) call ncread('vvel_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)),ivv,1,0.)
            if (ACC_DPU(n).ne.0) call ncread('dpu_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)),iuu,1,0.)
            if (ACC_DPV(n).ne.0) call ncread('dpv_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)),ivv,1,0.)
            if (ACC_UFLX(n).ne.0) call ncread('uflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)),iuu,1,0.)
            if (ACC_VFLX(n).ne.0) call ncread('vflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)),ivv,1,0.)
            if (ACC_UTFLX(n).ne.0) call ncread('utflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)),iuu,1,0.)
            if (ACC_VTFLX(n).ne.0) call ncread('vtflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)),ivv,1,0.)
            if (ACC_USFLX(n).ne.0) call ncread('usflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)),iuu,1,0.)
            if (ACC_VSFLX(n).ne.0) call ncread('vsflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)),ivv,1,0.)
            if (ACC_UMFLTD(n).ne.0) call ncread('umfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)),iuu,1,0.)
            if (ACC_VMFLTD(n).ne.0) call ncread('vmfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)),ivv,1,0.)
            if (ACC_UTFLTD(n).ne.0) call ncread('utfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)),iuu,1,0.)
            if (ACC_VTFLTD(n).ne.0) call ncread('vtfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)),ivv,1,0.)
            if (ACC_UTFLLD(n).ne.0) call ncread('utflld_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)),iuu,1,0.)
            if (ACC_VTFLLD(n).ne.0) call ncread('vtflld_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)),ivv,1,0.)
            if (ACC_USFLTD(n).ne.0) call ncread('usfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)),iuu,1,0.)
            if (ACC_VSFLTD(n).ne.0) call ncread('vsfltd_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)),ivv,1,0.)
            if (ACC_USFLLD(n).ne.0) call ncread('usflld_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)),iuu,1,0.)
            if (ACC_VSFLLD(n).ne.0) call ncread('vsflld_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)),ivv,1,0.)
            if (ACC_SALN(n).ne.0) call ncread('saln_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)),ip,1,0.)
            if (ACC_TEMP(n).ne.0) call ncread('temp_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)),ip,1,0.)
            if (ACC_DP(n).ne.0) call ncread('dp_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)),ip,1,0.)
            if (ACC_DZ(n).ne.0) call ncread('dz_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)),ip,1,0.)
            if (ACC_BFSQ(n).ne.0) call ncread('bfsq_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_BFSQ(n)),ip,1,0.)
            if (ACC_DIFDIA(n).ne.0) call ncread('difdia_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)),ip,1,0.)
            if (ACC_DIFINT(n).ne.0) call ncread('difint_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)),ip,1,0.)
            if (ACC_DIFISO(n).ne.0) call ncread('difiso_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)),ip,1,0.)
            if (ACC_WFLX(n).ne.0) call ncread('wflx_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)),ip,1,0.)
            if (ACC_WFLX2(n).ne.0) call ncread('wflx2_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)),ip,1,0.)
            if (ACC_AVDSG(n).ne.0) call ncread('avdsg_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_AVDSG(n)),ip,1,0.)
            if (ACC_DPVOR(n).ne.0) call ncread('dpvor_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_DPVOR(n)),ip,1,0.)
#if defined(TRC) && defined(TKE)
            if (ACC_TKE(n).ne.0) call ncread('tke_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)),ip,1,0.)
            if (ACC_GLS_PSI(n).ne.0) call ncread('gls_psi_phy'//c2,
     .        phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)),ip,1,0.)
#endif
            if (ACC_UVELLVL(n).ne.0) call ncread('uvellvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_UVELLVL(n)),iuu,1,0.)
            if (ACC_VVELLVL(n).ne.0) call ncread('vvellvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VVELLVL(n)),ivv,1,0.)
            if (ACC_UFLXLVL(n).ne.0) then
              call ncread('uflxlvl_phy'//c2,
     .                    phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)),
     .                    iuu,1,0.)
              call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)),
     .                    1,ddm, 1,1, halo_uv)
            endif
            if (ACC_VFLXLVL(n).ne.0) then
              call ncread('vflxlvl_phy'//c2,
     .                    phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)),
     .                    ivv,1,0.)
              call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)),
     .                    1,ddm, 1,1, halo_vv)
            endif
            if (ACC_UTFLXLVL(n).ne.0) call ncread('utflxlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_UTFLXLVL(n)),iuu,1,0.)
            if (ACC_VTFLXLVL(n).ne.0) call ncread('vtflxlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VTFLXLVL(n)),ivv,1,0.)
            if (ACC_USFLXLVL(n).ne.0) call ncread('usflxlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_USFLXLVL(n)),iuu,1,0.)
            if (ACC_VSFLXLVL(n).ne.0) call ncread('vsflxlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VSFLXLVL(n)),ivv,1,0.)
            if (ACC_UMFLTDLVL(n).ne.0) call ncread('umfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_UMFLTDLVL(n)),iuu,1,0.)
            if (ACC_VMFLTDLVL(n).ne.0) call ncread('vmfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VMFLTDLVL(n)),ivv,1,0.)
            if (ACC_UTFLTDLVL(n).ne.0) call ncread('utfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_UTFLTDLVL(n)),iuu,1,0.)
            if (ACC_VTFLTDLVL(n).ne.0) call ncread('vtfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VTFLTDLVL(n)),ivv,1,0.)
            if (ACC_UTFLLDLVL(n).ne.0) call ncread('utflldlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_UTFLLDLVL(n)),iuu,1,0.)
            if (ACC_VTFLLDLVL(n).ne.0) call ncread('vtflldlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VTFLLDLVL(n)),ivv,1,0.)
            if (ACC_USFLTDLVL(n).ne.0) call ncread('usfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_USFLTDLVL(n)),iuu,1,0.)
            if (ACC_VSFLTDLVL(n).ne.0) call ncread('vsfltdlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VSFLTDLVL(n)),ivv,1,0.)
            if (ACC_USFLLDLVL(n).ne.0) call ncread('usflldlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_USFLLDLVL(n)),iuu,1,0.)
            if (ACC_VSFLLDLVL(n).ne.0) call ncread('vsflldlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_VSFLLDLVL(n)),ivv,1,0.)
            if (ACC_SALNLVL(n).ne.0) call ncread('salnlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_SALNLVL(n)),ip,1,0.)
            if (ACC_TEMPLVL(n).ne.0) call ncread('templvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_TEMPLVL(n)),ip,1,0.)
            if (ACC_DZLVL(n).ne.0) call ncread('dzlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_DZLVL(n)),ip,1,0.)
            if (ACC_BFSQLVL(n).ne.0) call ncread('bfsqlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_BFSQLVL(n)),ip,1,0.)
            if (ACC_DIFDIALVL(n).ne.0) call ncread('difdialvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_DIFDIALVL(n)),ip,1,0.)
            if (ACC_DIFINTLVL(n).ne.0) call ncread('difintlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_DIFINTLVL(n)),ip,1,0.)
            if (ACC_DIFISOLVL(n).ne.0) call ncread('difisolvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_DIFISOLVL(n)),ip,1,0.)
            if (ACC_WFLXLVL(n).ne.0) call ncread('wflxlvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)),ip,1,0.)
            if (ACC_WFLX2LVL(n).ne.0) call ncread('wflx2lvl_phy'//c2,
     .         phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)),ip,1,0.)
            if (ACC_PVLVL(n).ne.0) call ncread('pvlvl_phy'//c2,
     .         phylvl(1-nbdy,1-nbdy,1,ACC_PVLVL(n)),ip,1,0.)
#if defined(TRC) && defined(TKE)
            if (ACC_TKELVL(n).ne.0) call ncread('tkelvl_phy'//c2,
     .        phylvl(1-nbdy,1-nbdy,1,ACC_TKELVL(n)),ip,1,0.)
            if (ACC_GLS_PSILVL(n).ne.0)
     .        call ncread('gls_psilvl_phy'//c2,
     .                    phylvl(1-nbdy,1-nbdy,1,ACC_GLS_PSILVL(n)),
     .                    ip,1,0.)
#endif
          endif
        enddo
      endif
c
      call ncfcls
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.1) then
            kfpla(i,j,1)=nint(rkfpla(i,j,1))
            kfpla(i,j,2)=nint(rkfpla(i,j,2))
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- set minimum physical temperature for each isopycnic layer
c --- ------------------------------------------------------------------
c
      call settemmin
c
#ifdef TRC
      call restart_trcrd(rstfnm)
#endif
c
      if (ditflx) then
c
c --- - read diag. heat flux restart file if available
        if (expcnf.eq.'cesm') then
          fnm=runid(1:runid_len)//'.blom.rtflx.'//rstfnm(runid_len+10:)
        else
          fnm=runid(1:runid_len)//'_resttflx_'//rstfnm(runid_len+10:)
        endif
        if(mnproc .eq.1) inquire(file=fnm,exist=fexist)
        call xcbcst(fexist)
        if (fexist) then
          call ncfopn(fnm,'r',' ',1,iotype)
          if (mnproc.eq.1)
     .      write (lp,'(a,a)') ' reading diag. heat flux restart file ',
     .                         trim(fnm)
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
        if (expcnf.eq.'cesm') then
          fnm=runid(1:runid_len)//'.blom.rsflx.'//rstfnm(runid_len+10:)
        else
          fnm=runid(1:runid_len)//'_restsflx_'//rstfnm(runid_len+10:)
        endif
        if (mnproc.eq.1) inquire(file=fnm,exist=fexist)
        call xcbcst(fexist)
        if (fexist) then
          call ncfopn(fnm,'r',' ',1,iotype)
          if (mnproc.eq.1)
     .      write (lp,'(a,a)') ' reading diag. salt flux restart file ',
     .                         trim(fnm)
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
          call ncfcls
        else
          if (mnproc.eq.1) write (lp,*)
     .      'Warning! No diag. salt flux restart file found'
        endif
c
      endif  
c
      return
      end
