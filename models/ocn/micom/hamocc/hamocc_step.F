      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- perform one HAMOCC step     
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mo_bgcmean
c
      implicit none
c
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
c
      integer :: m,n,mm,nn,k1m,k1n
      integer :: l
c
      call trc_limitc
c
      call micom2hamocc(m,n,mm,nn)
c
      do l=1,nbgc 
        bgcwrt(l)=0. 
        if (((diagann_bgc(l).and.nday_of_year.eq.1.or.diagmon_bgc(l)
     .    .and.nday.eq.1).and.mod(nstep,nstep_in_day).eq.0).or..not.
     .    (diagann_bgc(l).or.diagmon_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.1.) 
     .    bgcwrt(l)=1.
      enddo
c
      call hamocc4bcm(idm,jdm,kdm,pglat,bgc_swr,bgc_fice,
     .            bgc_t,bgc_s,bgc_slp,bgc_rho,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,nday_in_year)
c
      call hamocc2micom(m,n,mm,nn)
c
      return
      end
