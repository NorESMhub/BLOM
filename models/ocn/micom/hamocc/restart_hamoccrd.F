      subroutine restart_hamoccrd(nday1,rstfnm_ocn)
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
#include "common_geo.h"
#include "common_trc.h"
c
      integer i,j,k,l,nt,nday1
      character rstfnm_ocn*(*)
      real q

#ifdef RESTART_BGC
cka initialise arrays for dynamic part of BGC
      call inibgc
c
cka initialise HAMOCC
      do j=1,jj
        do i=1,ii
          pgila(2*i-1,2*j-1)=plon(i,j)
          pgila(2*i  ,2*j  )=plon(i,j)
          pgiph(2*i-1,2*j-1)=plat(i,j)
          pgiph(2*i  ,2*j  )=plat(i,j)
        enddo
      enddo

      pmonts=(nstep2-nstep1)/nstep_in_day/30.

      call micom2hamocc(0,1)
c
      call INI_HAMOCC(1,1,bgcdt,kpndtrun,idm,jdm,kdm,nbdy
     .               ,bgc_dp,bgc_t,bgc_s,bgc_dx,bgc_dy,bgc_pu,bgc_pw
     .               ,nyear,nmonth,nday,nstep1/nphys,pmonts
c     .               ,plon,plat,omask,ocetra,nocetra,th3d
     .               ,pgila,pgiph,omask,trc,ntr,ntrbgc,itrbgc
     .               ,rstfnm_ocn,path,path_len
     .               ,path2,path2_len)
c
#else
cka initialise arrays for dynamic part of BGC
      call inibgc
cka initialise HAMOCC
      do j=1,jj
        do i=1,ii
          pgila(2*i-1,2*j-1)=plon(i,j)
          pgila(2*i  ,2*j  )=plon(i,j)
          pgiph(2*i-1,2*j-1)=plat(i,j)
          pgiph(2*i  ,2*j  )=plat(i,j)
        enddo
      enddo

      pmonts=(nstep2-nstep1)/nstep_in_day/30.

      call micom2hamocc(0,1)
c
      call INI_HAMOCC(0,1,bgcdt,kpndtrun,idm,jdm,kdm,nbdy
     .               ,bgc_dp,bgc_t,bgc_s,bgc_dx,bgc_dy,bgc_pu,bgc_pw
     .               ,nyear,nmonth,nday,nstep1/nphys,pmonts
c     .               ,plon,plat,omask,ocetra,nocetra,th3d
     .               ,pgila,pgiph,omask,trc,ntr,ntrbgc,itrbgc
     .               ,rstfnm_ocn,path,path_len
     .               ,path2,path2_len)
c
#endif
c
      return
      end
