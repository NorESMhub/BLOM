      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
cka initialisation of HAMOCC      
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
#include "common_forc.h"
#include "common_geo.h"
#include "common_ice.h"
#include "common_trc.h"
c
      integer i,j,k,l,kn,nt
      integer m,n,mm,nn,k1m,k1n
      real q

      do k=1,kk
         kn=k+nn
cc$OMP PARALLEL DO 
         do j=1,jj
         do l=1,isp(j)
         do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c      write(lp,*) i,j,k,kn
cka layer thicknesses from last BGC timestep 
         dp_dia_bgc(i,j,k)=dp_dia_bgc(i,j,k)+dp(i,j,kn)
         enddo
         enddo
         enddo
cc$OMP END PARALLEL DO
      enddo
c
      call micom2hamocc(nn,n)
c
      bgc3dwrt=0. 
      if (diag3dmon.and.(mod(nstep,nstep_in_day).eq.0.and.nday.eq.1).or.
     .    diag3dmon.and.(mod(nstep,nstep_in_day).eq.1.and.nday.eq.1).or.
     .    .not.diag3dmon.and.mod(nstep+.5,diag3dfq).lt.1.  .or.
     .    .not.diag3dmon.and.mod(nstep-.5,diag3dfq).lt.1.) then
          bgc3dwrt=1.
      endif
c
      bgc2dwrt=0. 
      if (diag2dmon.and.(mod(nstep,nstep_in_day).eq.0.and.nday.eq.1).or.
     .    diag2dmon.and.(mod(nstep,nstep_in_day).eq.1.and.nday.eq.1).or.
     .    .not.diag2dmon.and.mod(nstep+.5,diag2dfq).lt.1.  .or.
     .    .not.diag2dmon.and.mod(nstep-.5,diag2dfq).lt.1.) then
          bgc2dwrt=1.
      endif
c
      call hamocc4bcm(idm,jdm,kdm,nbdy,
     .            bgc_swr,bgc_fice,bgc_t,bgc_s,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,trc,ntr,bgc3dwrt,bgc2dwrt,nday_in_year)
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
