c--------------------------------------------------------------------
c arrays to feed BCM variables to HAMOCC
c
      c o m m o n  /bgc_hamocc/
     . bgc_dx  (idm,jdm)      ,bgc_dy  (idm,jdm)
     .,bgc_dp  (idm,jdm,kdm)  ,bgc_dpio(idm,jdm,kdm)
     .,bgc_pu  (idm,jdm,kdm+1),bgc_pw  (idm,jdm,kdm+1)
     .,bgc_t   (idm,jdm,kdm)  ,bgc_s   (idm,jdm,kdm)
     .,omask   (idm,jdm)
     .,bgc_swr (idm,jdm)      ,bgc_fice(idm,jdm)
     .,bgc_awnd(idm,jdm) 
     .,bgc_atmco2(idm,jdm)    ,bgc_flxco2(idm,jdm)
c
      real bgc_dx,bgc_dy,bgc_dp,bgc_dpio
      real bgc_pu,bgc_pw,bgc_t,bgc_s,omask
      real bgc_swr,bgc_fice,bgc_awnd
      real bgc_atmco2,bgc_flxco2
c
      c o m m o n /bgc_hamocc_b/
     . ldtday,ldtmonth,kpndtrun,pmonts
c
      integer ldtday,ldtmonth,kpndtrun,pmonts
c
c----------------------------------------------------------------------
      c o m m o n/bgcc/
     . bgcdt
c
      real bgcdt
c
c----------------------------------------------------------------------
c number of physical timesteps in every bgc timestep
c
      integer nphys
      parameter (nphys=2)
c
c----------------------------------------------------------------------

      c o m m o n /bgc_hamocc2/
     . pgila,pgiph,bgc3dwrt,bgc2dwrt

      real pgila(2*idm,2*jdm),pgiph(2*idm,2*jdm)
      real bgc3dwrt,bgc2dwrt

