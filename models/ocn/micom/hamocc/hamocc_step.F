      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
cka initialisation of HAMOCC      
c
      use mod_xc
      use mo_bgcmean
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_geo.h"
#include "common_ice.h"
#include "common_trc.h"
c
      integer i,j,k,l,kn,nt
      integer m,n,mm,nn,k1m,k1n
      real q

      call micom2hamocc(nn,n)
c
      do l=1,nbgc
        bgcwrt(l)=0. 
        if (diagann_bgc(l).and.(mod(nstep,nstep_in_day).le.1.and.
     .    nday_of_year.eq.1).or.diagmon_bgc(l).and.
     .    (mod(nstep,nstep_in_day).le.1.and.nday.eq.1).or.
     .    .not.(diagann_bgc(l).or.diagmon_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.2. ) then
            bgcwrt(l)=1.
        endif
      enddo 
c
      call hamocc4bcm(idm,jdm,kdm,nbdy,
     .            bgc_swr,bgc_fice,bgc_t,bgc_s,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,trc,ntr,ntrbgc,itrbgc,nday_in_year)
c
c$OMP PARALLEL DO 
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          flxco2(i,j)=bgc_flxco2(i,j)
        enddo
        enddo
      enddo
cc$OMP END PARALLEL DO
c
      return
      end
