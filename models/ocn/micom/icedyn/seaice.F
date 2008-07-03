c***********************************************************************
      subroutine seaice
c
c This main subroutine of the sea ice model includes all sea ice physics.
c
c Markus Harder, March 1995
c***********************************************************************
      implicit none
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
c-----------------------------------------------------------------------
c Calculate ice drift
c-----------------------------------------------------------------------
      call drift
c-----------------------------------------------------------------------
      lold=lnew
      lnew=3-lold
c
c-----------------------------------------------------------------------
c Finished
c-----------------------------------------------------------------------
c
      return
      end             ! seaice
