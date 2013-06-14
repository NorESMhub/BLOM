c
c --- common blocks related to budget computations
c
      c o m m o n  /bud1/
     . sdp(ncalls,2),tdp(ncalls,2)
#ifdef TKE
     .,tkedp(ncalls,2)
#  ifdef GLS
     .,glsdp(ncalls,2)
#  endif
#endif
#ifdef TRC
     .,trdp(ncalls)
#endif
     .,mass0,sf,tf
#ifdef TRC
     .,trf
#endif
c
      real sdp,tdp
#ifdef TKE
     .    ,tkedp
#  ifdef GLS
     .    ,glsdp
#  endif
#endif
#ifdef TRC
     .    ,trdp
#endif
     .    ,mass0,sf,tf
#ifdef TRC
     .    ,trf
#endif
