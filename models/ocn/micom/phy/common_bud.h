c
c --- common blocks related to budget computations
c
      c o m m o n  /bud1/
     . sdp(ncalls,2),tdp(ncalls,2),ddp(ncalls,2)
#ifdef TRC
     .,trdp(ncalls)
#endif
     .,vol0,sf,tf
c
      real sdp,tdp,ddp,vol0,sf,tf
#ifdef TRC
      real trdp,trf
#endif
