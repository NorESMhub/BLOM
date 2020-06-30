! ------------------------------------------------------------------------------
! Copyright (C) 2015-2018 Mats Bentsen, Alok Kumar Gupta
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

      subroutine getfrc_cesm
c
c --- ------------------------------------------------------------------
c --- Interpolate CESM forcing fields
c --- ------------------------------------------------------------------
#define DIAG
#undef DIAG
c
      use mod_xc
#ifdef DIAG
      use mod_nctools
#endif
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_seaice.h"
#include "common_cesm.h"
#include "common_clndr.h"
c
      integer i,j,l
      real w1,w2
c
      if (smtfrc) then
        w1=.5*(1.+cos((mod(nstep-1,nstep_in_day)+1)*pi/nstep_in_day))
      else
        w1=0.
      endif
      w2=1.-w1
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
           ustarw(i,j)=w1*ustarw_da(i,j,l1ci)+w2*ustarw_da(i,j,l2ci)
           lip(i,j)=w1*lip_da(i,j,l1ci)+w2*lip_da(i,j,l2ci)
           sop(i,j)=w1*sop_da(i,j,l1ci)+w2*sop_da(i,j,l2ci)
           eva(i,j)=w1*eva_da(i,j,l1ci)+w2*eva_da(i,j,l2ci)
           rnf(i,j)=w1*rnf_da(i,j,l1ci)+w2*rnf_da(i,j,l2ci)
           rfi(i,j)=w1*rfi_da(i,j,l1ci)+w2*rfi_da(i,j,l2ci)
           fmltfz(i,j)=w1*fmltfz_da(i,j,l1ci)+w2*fmltfz_da(i,j,l2ci)
           sfl(i,j)=w1*sfl_da(i,j,l1ci)+w2*sfl_da(i,j,l2ci)
           swa(i,j)=w1*swa_da(i,j,l1ci)+w2*swa_da(i,j,l2ci)
           nsf(i,j)=w1*nsf_da(i,j,l1ci)+w2*nsf_da(i,j,l2ci)
           hmlt(i,j)=w1*hmlt_da(i,j,l1ci)+w2*hmlt_da(i,j,l2ci)
           slp(i,j)=w1*slp_da(i,j,l1ci)+w2*slp_da(i,j,l2ci)
           ficem(i,j)=w1*ficem_da(i,j,l1ci)+w2*ficem_da(i,j,l2ci)
           abswnd(i,j)=w1*abswnd_da(i,j,l1ci)+w2*abswnd_da(i,j,l2ci)
           atmco2(i,j)=w1*atmco2_da(i,j,l1ci)+w2*atmco2_da(i,j,l2ci)
        enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
           ztx(i,j)=w1*ztx_da(i,j,l1ci)+w2*ztx_da(i,j,l2ci)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
           mty(i,j)=w1*mty_da(i,j,l1ci)+w2*mty_da(i,j,l2ci)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
#ifdef DIAG
      call ncfopn('getfrc_cesm.nc','w','c',1,iotype)
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdefvar('ustarw_da','x y',ndouble,8) 
      call ncdefvar('lip_da','x y',ndouble,8)
      call ncdefvar('sop_da','x y',ndouble,8)
      call ncdefvar('eva_da','x y',ndouble,8)
      call ncdefvar('rnf_da','x y',ndouble,8)
      call ncdefvar('rfi_da','x y',ndouble,8)
      call ncdefvar('fmltfz_da','x y',ndouble,8)
      call ncdefvar('sfl_da','x y',ndouble,8)
      call ncdefvar('swa_da','x y',ndouble,8)
      call ncdefvar('nsf_da','x y',ndouble,8)
      call ncdefvar('hmlt_da','x y',ndouble,8)
      call ncdefvar('slp_da','x y',ndouble,8)
      call ncdefvar('ficem_da','x y',ndouble,8)
      call ncdefvar('abswnd_da','x y',ndouble,8)
      call ncdefvar('atmco2_da','x y',ndouble,8)
      call ncdefvar('ztx_da','x y',ndouble,8)
      call ncdefvar('mty_da','x y',ndouble,8)
      call ncedef

      call ncwrtr('ustarw_da','x y',ustarw_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('lip_da','x y',lip_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('sop_da','x y',sop_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('eva_da','x y',eva_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('rnf_da','x y',rnf_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('rfi_da','x y',rfi_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('fmltfz_da','x y',fmltfz_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('sfl_da','x y',sfl_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('swa_da','x y',swa_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('nsf_da','x y',nsf_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('hmlt_da','x y',hmlt_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('slp_da','x y',slp_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('ficem_da','x y',ficem_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('abswnd_da','x y',abswnd_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('atmco2_da','x y',atmco2_da(1-nbdy,1-nbdy,l2ci),
     .            ip,1,1.,0.,8)
      call ncwrtr('ztx_da','x y',ztx_da(1-nbdy,1-nbdy,l2ci),
     .            iu,1,1.,0.,8)
      call ncwrtr('mty_da','x y',mty_da(1-nbdy,1-nbdy,l2ci),
     .            iv,1,1.,0.,8)
      call ncfcls
      call xcstop('(getfrc_cesm)')
             stop '(getfrc_cesm)'
#endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'getfrc_cesm:'
        endif
        call chksummsk(ustarw,ip,1,'ustarw')
        call chksummsk(ztx,iu,1,'ztx')
        call chksummsk(mty,iv,1,'mty')
        call chksummsk(lip,ip,1,'lip')
        call chksummsk(sop,ip,1,'sop')
        call chksummsk(eva,ip,1,'eva')
        call chksummsk(rnf,ip,1,'rnf')
        call chksummsk(rfi,ip,1,'rfi')
        call chksummsk(fmltfz,ip,1,'fmltfz')
        call chksummsk(sfl,ip,1,'sfl')
        call chksummsk(swa,ip,1,'swa')
        call chksummsk(nsf,ip,1,'nsf')
        call chksummsk(hmlt,ip,1,'hmlt')
        call chksummsk(slp,ip,1,'slp')
        call chksummsk(ficem,ip,1,'ficem')
        call chksummsk(abswnd,ip,1,'abswnd')
        call chksummsk(atmco2,ip,1,'abswnd')
      endif
c
      return
      end
