c
c --- common blocks related to tracers
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm,ntr) :: trc
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ntr) :: trcold
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  uflxtr,vflxtr,trflx
c
      common /trc1/ trc,trcold,uflxtr,vflxtr,trflx
