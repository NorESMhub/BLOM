c***********************************************************************
      subroutine drift
c***********************************************************************
c Solve momentum equation
c
c Code rewritten and documented by Markus Harder, 1994
c-----------------------------------------------------------------------
c The momentum equation is solved in 2 parts:
c      1)      Predictive time step dt/2 :
c            Calculate centered velocities U(i,j,3) using the viscosities
c            derived from the old velocities U(i,j,lold)
c            Then the viscosities are recalculated with the velocities
c            from the predictive time step, and finally:
c      2)      Centered time step dt :
c            Calculate new velocities U(i,j,lnew) using the viscosities
c            derived from the centered velocities U(i,j,3)
c
c llhs (left hand side)  denotes variables that are recalculated
c lrhs (right hand side) denotes variables that are given for the actual step
c
c Modification by Markus Harder, 14-FEB-95:
c      4 iteration steps
c-----------------------------------------------------------------------
      implicit none
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer nloops,ntotal,jjc,lrhs,llhs

c-----------------------------------------------------------------------
c  Time average of ocean forcing
c-----------------------------------------------------------------------
      call tavocf

c-----------------------------------------------------------------------
c  Half number of iterations in outer loop
c-----------------------------------------------------------------------
      nloops=10


c-----------------------------------------------------------------------
c Multiply by 2 to get an even number of iterations
c-----------------------------------------------------------------------
      nloops=2*nloops
c-----------------------------------------------------------------------
c Set the indices for the actual step
c
c In centered scheme:
c      jjc=1      predictive step
c      jjc=2      centered step
c-----------------------------------------------------------------------
      ntotal=0
      do 60 jjc=1,nloops
        if (mod(jjc,2).eq.0) then
          lrhs=3
          llhs=lnew
        else
          lrhs=lnew
          llhs=3
        endif

        if (jjc.eq.1) then
          lrhs=lold
          llhs=3
        endif

c-----------------------------------------------------------------------
c Define error margin for the solution
c-----------------------------------------------------------------------
      vrmax=((nloops+1)-jjc)*1.E-4

c-----------------------------------------------------------------------
c 1.) Calculate ice pressure as function of ice thickness and concentration.
c Thicknesses and concentrations are taken from the previous time step.
c
c Hibler (1979):   P = Pstar * h * exp(C * (1-A))
c
c
c 2.) Calculate the viscosities as a function of the ice strength
c and of the strain rates, derived from the rhs velocities
c
c Hibler (1979):
c      zeta = 0.5 * P / delta(strain rates)
c      eta  = zeta / e**2
c
c where e is the eccentricity of the elliptical yield curve.
c (Other concepts may be applied to describe the viscosities.)
c-----------------------------------------------------------------------
        call plast (lrhs)
c-----------------------------------------------------------------------
c Calculate those parts of the momentum equation
c that do not depend on the new velocity,
c i. e. set up the coefficients for the equation
c to be solved by subroutine  relax.
c
c Attention is required for the oceanic drag:
c The oceanic drag DOES depend on the new velocities,
c but it is calculated here using the old velocities.
c-----------------------------------------------------------------------
        call relcon (lrhs)
c-----------------------------------------------------------------------
c Solve the momentum equation in an iterative scheme.
c To obtain a faster approach to the solution, over-relaxation is used.
c (If numerical instabilities occur, over-relaxation should
c be reduced or turned off by setting the over-relaxation factor WTA
c to a value closer to 1 or identical to 1.)
c-----------------------------------------------------------------------
        call relax (lrhs, llhs)
c-----------------------------------------------------------------------
c Store the characteristics of the numerical solution
c For diagnostics only
c-----------------------------------------------------------------------
        ntotal=ntotal+mrelax
c-----------------------------------------------------------------------
   60      continue
c
      write(*,*) 'average no. of iterations:',real(ntotal)/nloops
c=======================================================================
c Solution of momentum equation completed !
c=======================================================================

      return
      end            ! subroutine drift
