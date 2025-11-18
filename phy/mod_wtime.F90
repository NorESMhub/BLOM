! ------------------------------------------------------------------------------
! Copyright (C) 2005 HYCOM Consortium and contributors
! Copyright (C) 2006 Mats Bentsen
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

module mod_wtime

  use mod_types,  only: r8

  implicit none
  public

contains

#if defined(AIX)

  real(r8) function wtime()
    ! use the function  rtc  to return wall time.
    real(r8) :: rtc
    wtime = rtc()
  end function wtime

#elif defined(MPI)

  real(r8) function wtime()
    ! use the mpi function  mpi_wtime  to return wall time.
    real(r8) :: mpi_wtime
    wtime = mpi_wtime()
  end function wtime

#else

  real(r8) function wtime()
    ! use the f90 intrinsic  system_clock  to return wall time.
    ! will fail if the count is ever negative, but the standard
    ! says that it is aways non-negative if a clock exists.
    ! not thread-safe, unless lcount and iover are threadprivate.
    real(r8), parameter :: zero= 0.0
    real(r8), parameter :: one = 1.0
    integer :: count, mcount, rate
    real(r8) , save :: offsec, offset, persec
    integer, save :: icount, iover,  lcount, ncount
    data iover, lcount / -1, -1 /
    !
    call system_clock(count)
    !
    if (count.lt.lcount) then
      ! count is supposed to be non-decreasing except when it wraps,
      ! but some implementations don''t do this.  so ignore any
      ! decrease of less than one percent of the range.
      if (lcount-count.lt.ncount) then
        count  = lcount
      else
        iover  = iover + 1
        offset = offset + offsec
      endif
    endif
    lcount = count

    if     (iover.eq.0) then
      ! first cycle, for accuracy with 64-bit counts.
      wtime = (count - icount) * persec
    elseif (iover.gt.0) then
      ! all other cycles.
      wtime = count * persec + offset
    else
      ! initialization.
      call system_clock(icount, rate, mcount)
      ncount =  mcount/100
      persec =  one/rate
      offsec =  mcount * persec
      offset = -icount * persec
      iover  =  0
      wtime  =  zero
    endif
  end function wtime

#endif  ! MPI:else

end module mod_wtime
