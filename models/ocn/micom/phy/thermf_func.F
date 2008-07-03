      real function qsatw(t,p)
c
c --- computes saturation specific humidity [kg/kg] over open water,
c --- Buck (1981) JAM 20, 1527-1532
c
      implicit none
c
c --- input variables:
c ---   t      - air temperature [K]
c ---   p      - atmospheric pressure [Pa]
c
      real t,p
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c
      real eps,c1,c2,c3,c4,t0,t1
      parameter(eps=0.62197,c1=611.21,c2=1.0007,c3=3.46e-8,c4=17.502,
     .          t0=273.15,t1=32.19)
c
      real tl,e
c
      tl=max(150.,t)
      e=c1*(c2+c3*p)*exp(c4*(tl-t0)/(tl-t1))
      qsatw=eps*e/(p-(1.-eps)*e)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dqsatw(t,p)
c
c --- computes the derivative of saturation specific humidity [kg/kg/K]
c --- over open water with respect to temperature,
c --- Buck (1981) JAM 20, 1527-1532
c
      implicit none
c
c --- input variables:
c ---   t      - air temperature [K]
c ---   p      - atmospheric pressure [Pa]
c
      real t,p
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c
      real eps,c1,c2,c3,c4,t0,t1
      parameter(eps=0.62197,c1=611.21,c2=1.0007,c3=3.46e-8,c4=17.502,
     .          t0=273.15,t1=32.19)
c
      real tl,e,dedt
c
      tl=max(150.,t)
      e=c1*(c2+c3*p)*exp(c4*(tl-t0)/(tl-t1))
      dedt=e*c4*(t0-t1)/(tl-t1)**2
      dqsatw=dedt*eps*p/(p-(1.-eps)*e)**2
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function qsati(t,p)
c
c --- computes saturation specific humidity [kg/kg] over sea ice,
c --- Parkinson and Washington (1979) JGR 84, 311-337
c
      implicit none
c
c --- input variables:
c ---   t      - air temperature [K]
c ---   p      - atmospheric pressure [Pa]
c
      real t,p
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c
      real eps,c1,c2,t0,t1
      parameter(eps=0.62197,c1=611.,c2=9.5,t0=273.15,t1=7.66)
c
      real tl,e
c
      tl=max(150.,t)
      e=c1*10.**(c2*(tl-t0)/(tl-t1))
      qsati=eps*e/(p-(1.-eps)*e)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function dqsati(t,p)
c
c --- computes the derivative of saturation specific humidity [kg/kg/k]
c --- over sea ice with respect to temperature,
c --- Parkinson and Washington (1979) JGR 84, 311-337
c
      implicit none
c
c --- input variables:
c ---   t      - air temperature [K]
c ---   p      - atmospheric pressure [Pa]
c
      real t,p
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c
      real eps,c1,c2,t0,t1
      parameter(eps=0.62197,c1=611.,c2=9.5,t0=273.15,t1=7.66)
c
      real tl,e,dedt
c
      tl=max(150.,t)
      e=c1*10.**(c2*(tl-t0)/(tl-t1))
      dedt=e*c2*(t0-t1)*log(10.)/(tl-t1)**2
      dqsati=dedt*eps*p/(p-(1.-eps)*e)**2
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function rhoair(t,q,p)
c
c --- computes air density [kg/m^3]
c
      implicit none
c
c --- input variables:
c ---   t      - air temperature [K]
c ---   q      - specific humidity of air [kg/kg]
c ---   p      - atmospheric pressure [Pa]
c
      real t,q,p
c
c --- parameters:
c ---   eps    - molecular weight ratio of dry air and water vapour
c ---   rgas   - gas constant for dry air [J/(kg*K)]
c
      real eps,cv,rgas
      parameter(eps=.62197,cv=1./eps-1.,rgas=287.04)
c
      real tv
c
c --- virtual temperature
      tv=t*(1.+cv*q)
c
c --- moist air density [kg/m^3]
      rhoair=p/(rgas*tv)
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      real function swtfrz(s)
c
c --- computes sea water freezing point temperature [deg. celsius],
c
      implicit none
c
c --- input variables:
c ---   s      - salinity [psu]
c
      real s
c
#ifdef CCSMCOUPLED
      swtfrz=-1.8
#else
      swtfrz=-0.0547*s
#endif
c
      return
      end
