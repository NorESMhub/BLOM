      real function hfharm(a,b)
c
c --- harmonic average divided by 2
c
      implicit none
      real a,b
c
      hfharm=a*b/(a+b)
      return
      end
