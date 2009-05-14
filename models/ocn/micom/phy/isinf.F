      logical function isinf(a) 
c
c --- check if a real number is Inf
c
      implicit none
c
      real a 
c
      if ((a*0).ne.0) then 
        isinf=.true. 
      else 
        isinf=.false. 
      endif 
c
      return 
      end
