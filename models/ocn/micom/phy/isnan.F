      logical function isnan(a) 
c
c --- check if a real number is NaN
c
      implicit none
c
      real a 
c
      if (a.ne.a) then 
        isnan=.true. 
      else 
        isnan=.false. 
      endif 
c
      return 
      end 
