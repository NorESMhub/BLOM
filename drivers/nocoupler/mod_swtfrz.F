      module mod_swtfrz
c
c --- ------------------------------------------------------------------
c --- This module contains routines for computing the freezing point of
c --- sea water.
c --- ------------------------------------------------------------------
c
      implicit none
c
      private
c
      public :: swtfrz
c
      interface swtfrz
         module procedure swtfrz_0d
         module procedure swtfrz_1d
         module procedure swtfrz_2d
      end interface swtfrz
c
      contains
c
c --- ------------------------------------------------------------------
c
      function swtfrz_0d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p      ! Pressure [g cm-1 s-2]
      real, intent(in) :: s      ! Salinity [g kg-1]
      real :: swtfrz
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_0d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_1d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p(:)   ! Pressure [g cm-1 s-2]
      real, intent(in) :: s(:)   ! Salinity [g kg-1]
      real :: swtfrz(size(s))
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_1d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_2d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p(:,:) ! Pressure [g cm-1 s-2]
      real, intent(in) :: s(:,:) ! Salinity [g kg-1]
      real :: swtfrz(size(s,1),size(s,2))
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_2d
c
c --- ------------------------------------------------------------------
c
      end module mod_swtfrz
