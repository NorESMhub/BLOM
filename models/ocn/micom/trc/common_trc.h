c
c --- common blocks related to tracers
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ntr) :: trc
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  uflxtr,vflxtr,trflx
      integer ntrc
c
      common /trc1/ trc,uflxtr,vflxtr,trflx,ntrc
#ifdef ATRC
c
c --- common blocks related to age tracers
c
      integer natrc
c
      common /atrc1/ natrc
c
#endif
