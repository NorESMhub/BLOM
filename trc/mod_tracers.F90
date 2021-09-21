! ------------------------------------------------------------------------------
! Copyright (C) 2007-2020 Mats Bentsen, Jörg Schwinger, Jerry Tjiputra,
!                         Alok Kumar Gupta
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_tracers
! ------------------------------------------------------------------------------
! This module declares variables related to tracers.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_xc

   implicit none

   private

   ! Number of ocean tracers, not including turbulent kinetic energy, generic
   ! length scale, ideal age and HAMOCC tracers.
   integer, parameter :: ntrocn = 0

   ! Number of age tracers.
   integer, parameter :: natr = 0

   ! Number of turbulent kinetic energy and generic length scale tracers.
#ifdef TKE
   integer, parameter :: ntrtke = 1, ntrgls = 1
   ! Indices of first turbulent kinetic energy and generic length scale tracers.
   integer, parameter :: itrtke = ntrocn - natr + 1
   integer, parameter :: itrgls = ntrocn - natr + ntrtke + 1
#else
   integer, parameter :: ntrtke = 0,  ntrgls = 0
   integer, parameter :: itrtke = -1, itrgls = -1
#endif

   ! Number of ideal age tracer.
#ifdef IDLAGE
   integer, parameter :: ntriag = 1
   ! Index of first ideal age tracer.
   integer, parameter :: itriag = ntrocn - natr + ntrtke + ntrgls + 1
#else
   integer, parameter :: ntriag = 0
   integer, parameter :: itriag = -1
#endif

   ! HAMOCC tracers.
   integer, protected :: ntrbgc = -999   ! Number of HAMOCC tracers.
   integer, protected :: itrbgc = -999   ! Index of first HAMOCC tracer.

   ! Total number of tracers. (Return error if default value is used)
   integer, protected :: ntr = -999

   real(r8), allocatable, dimension(:,:,:,:) ::  &
      trc,      & ! Tracer array.
      trcold      ! Tracer array at old time level.

   real(r8), allocatable, dimension(:,:,:) ::  &
      uflxtr,   & ! u-component of tracer flux.
      vflxtr,   & ! v-component of tracer flux.
      trflx       ! Surface flux of tracer.

   public :: ntrocn, ntrtke, ntrgls, ntriag, ntrbgc, ntr, natr,                 &
             itrtke, itrgls, itriag, itrbgc,                                    &
             trc, trcold, uflxtr, vflxtr, trflx,                                &
             inivar_tracers

contains

   subroutine allocate_tracers
   ! ---------------------------------------------------------------------------
   ! If using HAMOCC, set number of bgc tracers : ntrbgc, itrbgc
   ! Set total number of tracers : ntr
   ! Allocate tracer arrays.
   ! ---------------------------------------------------------------------------

#ifdef HAMOCC
     use mo_param1_bgc, only: nocetra
#endif

     implicit none

     integer :: errstat
     integer :: num_bgc_tracers = 0

     ! Number of HAMOCC tracers.
#ifdef HAMOCC
     ntrbgc = nocetra
     itrbgc = ntrocn - natr + ntrtke + ntrgls + ntriag + 1
#else
     ntrbgc = 0
     itrbgc = -1
#endif

     ! Total number of tracers.
     if (mnproc.eq.1) then
        ntr = ntrocn + ntrtke + ntrgls + ntriag + ntrbgc

        write(lp,'(A,1X,I8)') 'Total tracer count: ', ntr
        if (ntr < 0) then
           write(lp,'(a)') 'Number of tracers must be non-negative.'
           call xchalt('(allocate_tracers)')
           stop '(allocate_tracers)'
        endif
     endif
     call xcbcst(ntr)

     ! Tracer array.
     allocate (trc(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 2*kdm, ntr),       &
          stat = errstat)
     if(errstat.ne.0) then
        write(lp,'(a)') 'Not enough memory for variable trc.'
        call xchalt('(allocate_tracers)')
        stop '(allocate_tracers)'
     endif

     ! Tracer array at old time level.
     allocate (trcold(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm, ntr),      &
          stat = errstat)
     if(errstat.ne.0) then
        write(lp,'(a)') 'Not enough memory for variable trcold.'
        call xchalt('(allocate_tracers)')
        stop '(allocate_tracers)'
     endif

     ! u-component of tracer flux.
     allocate (uflxtr(ntr, 1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy),           &
          stat = errstat)
     if(errstat.ne.0) then
        write(lp,'(a)') 'Not enough memory for variable uflxtr.'
        call xchalt('(allocate_tracers)')
        stop '(allocate_tracers)'
     endif

     ! v-component of tracer flux.
     allocate (vflxtr(ntr, 1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy),           &
          stat = errstat)
     if(errstat.ne.0) then
        write(lp,'(a)') 'Not enough memory for variable vflxtr.'
        call xchalt('(allocate_tracers)')
        stop '(allocate_tracers)'
     endif

     ! Surface flux of tracer.
     allocate (trflx(ntr, 1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy),            &
          stat = errstat)
     if(errstat.ne.0) then
        write(lp,'(a)') 'allocate_tracers: Not enough memory for variable trflx.'
        call xchalt('(allocate_tracers)')
        stop '(allocate_tracers)'
     endif
   end subroutine allocate_tracers

   subroutine inivar_tracers
   ! ---------------------------------------------------------------------------
   ! Initialize variables related to tracers.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
         uflux, vflux

      integer :: i, j, k, l, nt

      ! Allocate tracer arrays
      call allocate_tracers

   !$omp parallel do private(i, k, nt)
      do j = 1 - nbdy, jj + nbdy
         do i = 1 - nbdy, ii + nbdy
            uflux(i, j) = spval
            vflux(i, j) = spval
            do nt = 1, ntr
               trflx(nt, i, j) = spval
            enddo
            do k = 1, kk
               do nt = 1, ntr
                  trc(i ,j ,k     , nt) = spval
                  trc(i ,j ,k + kk, nt) = spval
               enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! Initialize uflxtr at points located upstream and downstream (in i
      ! direction) of p-points.
   !$omp parallel do private(l, i, k)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l) + 1)
            uflux(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(uflux, 1, 1, nbdy, nbdy, halo_us)
   !$omp parallel do private(i, nt)
      do j = 1 - nbdy, jdm + nbdy
         do i = 1 - nbdy, idm + nbdy
            do nt = 1, ntr
               uflxtr(nt, i, j) = uflux(i, j)
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! Initialize vflxtr at points located upstream and downstream (in j
      ! direction) of p-points.
   !$omp parallel do private(l, j, k)
      do i = 1, ii
         do l = 1, jsp(i)
         do j = max(1, jfp(i, l)), min(jj, jlp(i, l) + 1)
            vflux(i, j)=0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(vflux, 1, 1, nbdy, nbdy, halo_vs)
   !$omp parallel do private(i, nt)
      do j = 1 - nbdy, jdm + nbdy
         do i = 1 - nbdy, idm + nbdy
            do nt = 1, ntr
               vflxtr(nt, i, j) = vflux(i, j)
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! Initialize surface flux of tracer.
   !$omp parallel do private(l, i, nt)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            do nt = 1, ntr
               trflx(nt, i, j) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_tracers

end module mod_tracers
