      subroutine iniconicedyn

      use mod_xc, only: mnproc

      implicit none

#include "param_icedyn.h"
#include "common_icedyn.h"

      integer i,j

      if (mnproc.eq.1) then

c-----------------------------------------------------------------------
c Initial values of running indices
c-----------------------------------------------------------------------
        lold=1
        lnew=2
C----------------------------------------------------------------------
c Initial conditions for vector points
C----------------------------------------------------------------------
c$OMP PARALLEL DO
        do 105 j=0,M
          do 105 i=1,L
c----------------------------------------------------------------------
c Velocities
c----------------------------------------------------------------------
            u(i,j,1) = 0.001*vm(i,j)    ! small initial velocities
            v(i,j,1) = 0.001*vm(i,j)    ! small initial velocities
            u(i,j,2) = 0.0
            v(i,j,2) = 0.0
            u(i,j,3) = 0.0
            v(i,j,3) = 0.0
            uwat2l(i,j,1) = 0.0
            vwat2l(i,j,1) = 0.0
            uwat2l(i,j,2) = 0.0
            vwat2l(i,j,2) = 0.0
            gradhx2l(i,j,1) = 0.0
            gradhy2l(i,j,1) = 0.0
            gradhx2l(i,j,2) = 0.0
            gradhy2l(i,j,2) = 0.0
  105   continue
c$OMP END PARALLEL DO

      endif

      return
      end
