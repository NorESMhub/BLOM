      module mod_swtfrz
c
c --- ------------------------------------------------------------------
c --- This module contains routines for computing the freezing point of
c --- sea water.
c --- ------------------------------------------------------------------
c
      use types, only : r8
      use shr_frz_mod, only: shr_frz_freezetemp
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
c --- Retrieve freezing temperature from shared CESM function.
c --- ------------------------------------------------------------------
c
      implicit none
c
      real (r8), intent(in) :: p      ! Pressure [g cm-1 s-2]
      real (r8), intent(in) :: s      ! Salinity [g kg-1]
      real (r8) :: swtfrz
c
      swtfrz=shr_frz_freezetemp(s)
c
      end function swtfrz_0d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_1d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Retrieve freezing temperature from shared CESM function.
c --- ------------------------------------------------------------------
c
      implicit none
c
      real (r8), intent(in) :: p(:)   ! Pressure [g cm-1 s-2]
      real (r8), intent(in) :: s(:)   ! Salinity [g kg-1]
      real (r8) :: swtfrz(size(s))
c
      swtfrz(:)=shr_frz_freezetemp(s(:))
c
      end function swtfrz_1d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_2d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Retrieve freezing temperature from shared CESM function.
c --- ------------------------------------------------------------------
c
      implicit none
c
      real (r8), intent(in) :: p(:,:) ! Pressure [g cm-1 s-2]
      real (r8), intent(in) :: s(:,:) ! Salinity [g kg-1]
      real (r8) :: swtfrz(size(s,1),size(s,2))
c
      swtfrz(:,:)=shr_frz_freezetemp(s(:,:))
c
      end function swtfrz_2d
c
c --- ------------------------------------------------------------------
c
      end module mod_swtfrz
