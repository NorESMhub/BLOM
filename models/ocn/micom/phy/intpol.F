      real function intpol(d1,d2,d3,d4,d5,x)
c
c --- interpolate data
c
      real d1,d2,d3,d4,d5,x
c
      real a1,a2,a3,b1,b2,b3,b4,b5,c1,c2
      parameter(a1=-3./7.,a2=-15./7.,a3= 3./2.,
     .          b1= 4./7.,b2=-16./7.,b3=15./7.,b4=-5./7.,b5=2./7.,
     .          c1=-1./7.,c2=  9./14)
      real a,b,c
c
      a=a1*(d1+d5)+a2*d3+a3*(d2+d4)
      b=b1*d1+b2*d2+b3*d3+b4*d4+b5*d5
      c=c1*(d1+d4)+c2*(d2+d3)
c
      intpol=(a*x+b)*x+c
c
      return
      end
