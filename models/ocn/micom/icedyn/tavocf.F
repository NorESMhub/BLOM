      subroutine tavocf
c
c --- time average of ocean forcing
c
      implicit none
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      integer i,j
c
c$OMP PARALLEL DO
      do j=0,M
        do i=1,L
          uwat(i,j)=.5*(uwat2l(i,j,1)+uwat2l(i,j,2))
          vwat(i,j)=.5*(vwat2l(i,j,1)+vwat2l(i,j,2))
          gradhx(i,j)=.5*(gradhx2l(i,j,1)+gradhx2l(i,j,2))
          gradhy(i,j)=.5*(gradhy2l(i,j,1)+gradhy2l(i,j,2))
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
