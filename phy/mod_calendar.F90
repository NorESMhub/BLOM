! ------------------------------------------------------------------------------
! Copyright (C) 2020 Mats Bentsen
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

module mod_calendar
  ! ------------------------------------------------------------------------------
  ! This module contains calendar routines.
  !
  ! The supported calendars are the same as in the NetCDF Climate and Forecast
  ! (CF) Metadata Conventions:
  !
  !    - 'gregorian' or 'standard': Mixed Gregorian/Julian calendar as defined by
  !                                 UDUNITS.
  !    - 'proleptic_gregorian':     A Gregorian calendar extended to dates before
  !                                 15 Oct 1582. That is, a year is a leap year if
  !                                 either (i) it is divisible by 4 but not by 100
  !                                 or (ii) it is divisible by 400.
  !    - 'julian':                  Julian calendar.
  !    - 'noleap' or '365_day':     Gregorian calendar without leap years, i.e.,
  !                                 all years are 365 days long.
  !    - 'all_leap' or '366_day':   Gregorian calendar with every year being a
  !                                 leap year, i.e., all years are 366 days long.
  !    - '360_day':                 All years are 360 days divided into 30 day
  !                                 months.
  !
  ! For Julian and Gregorian calendars, the Chronological Julian Day Number (CJDN)
  ! is used (algorithms at https://www.aa.quae.nl/en/reken/juliaansedag.html).
  ! Zero CJDN corresponds to 1 Jan -4712 and in the Julian calendar and 24 Nov
  ! -4713 in the proleptic Gregorian calendar. The other calendars uses a day
  ! number where zero day number corresponds to 1 Jan 1.
  !
  ! The available functions are invoked by:
  !
  !    errstat = date_to_daynum(calendar, date, daynum)
  !    errstat = daynum_to_date(calendar, daynum, date)
  !    errstat = daynum_diff(calendar, date1, date2, dndiff)
  !    errstat = date_offset(calendar, date, dnoffset)
  !    errstat = date_check(calendar, date)
  !    errstr = calendar_errstr(errstat)
  !
  ! The 'calendar' argument (character(len = *)) gives the calendar type. The
  ! arguments 'date', 'date1', 'date2' are of type 'date_type' defined by the
  ! module. All other arguments are of default integer type. The first 5 functions
  ! returns an integer error status value which is equal to 'calendar_noerr' in
  ! case of no error. The function 'calendar_errstr' returns a static reference to
  ! an error message string (character(len = 80)) corresponding to an integer
  ! error status.
  !
  ! With default 4-byte integer type, the valid range of the various calendar
  ! types are:
  !
  !    - 'gregorian' or 'standard':
  !       - Day number range: [     -535149630,      538592031]
  !       - Date range:       [ 1 Mar -1469872, 18 Oct 1469902]
  !    - 'proleptic_gregorian':
  !       - Day number range: [     -535148831,      538592031]
  !       - Date range:       [ 1 Mar -1469900, 18 Oct 1469902]
  !    - 'julian':
  !       - Day number range: [     -535149630,      538592029]
  !       - Date range:       [ 1 Mar -1469872,  8 Nov 1469872]
  !    - 'noleap' or '365_day':
  !       - Day number range: [    -2147483648,     2147483341]
  !       - Date range:       [27 Feb -5883516,  2 Jan 5883517]
  !    - 'all_leap' or '366_day':
  !       - Day number range: [    -2147483648,     2147483341]
  !       - Date range:       [ 4 May -5867441, 28 Oct 5867441]
  !    - '360_day':
  !       - Day number range: [    -2147483648,     2147483647]
  !       - Date range:       [23 Aug -5965232,  8 May 5965233]
  ! ------------------------------------------------------------------------------

  implicit none
  private

  integer, parameter :: &
       calendar_noerr                  = 0, &
       calendar_unsupported            = 1, &
       calendar_invalid_date           = 2, &
       calendar_invalid_gregorian_date = 3, &
       calendar_daynum_overflow        = 4, &
       calendar_daynum_offset_overflow = 5, &
       errmsg_num                      = 5, &
       last_julian_daynum              = 2299160

  character(len = 80), dimension(errmsg_num), parameter :: errmsg = &
       (/"Calendar type is not supported!                              ", &
         "Date argument is incorrect or causes overflow!               ", &
         "Date argument is invalid for mixed Julian/Gregorian calendar!", &
         "Day number causes overflow!                                  ", &
         "Day number offset causes overflow!                           "/)

  type date_type
    integer :: year, month, day
  end type date_type

  interface operator(==)
    module procedure dates_equal
  end interface operator(==)

  interface operator(<)
    module procedure date1_lt_date2
  end interface operator(<)

  interface operator(>)
    module procedure date1_gt_date2
  end interface operator(>)

  interface operator(/=)
    module procedure dates_not_equal
  end interface operator(/=)

  interface operator(<=)
    module procedure date1_le_date2
  end interface operator(<=)

  interface operator(>=)
    module procedure date1_ge_date2
  end interface operator(>=)

  public :: date_type, date_to_daynum, daynum_to_date, daynum_diff, &
       date_offset, date_check, calendar_noerr, calendar_errstr, &
       operator(==), operator(<), operator(>), operator(/=), &
       operator(<=), operator(>=)

contains

  ! ---------------------------------------------------------------------------
  ! Private procedures.
  ! ---------------------------------------------------------------------------

  pure function intdivfloor(a, b)
    ! ---------------------------------------------------------------------------
    ! Returns integer quotient of integer division, rounded towards minus
    ! infinity. This function assumes a positive denominator.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: a, b

    integer :: intdivfloor

    intdivfloor = a/b
    if (mod(a, b) < 0) intdivfloor = intdivfloor - 1

  end function intdivfloor

  subroutine date_to_daynum_julian(date, daynum)
    ! ---------------------------------------------------------------------------
    ! Convert from Julian calendar date to Chronological Julian Day Number
    ! (CJDN).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date ! Julian calendar date.
    integer, intent(out) :: daynum      ! CJDN.

    integer :: c0

    c0 = intdivfloor(date%month - 3, 12)

    daynum = intdivfloor(1461*(date%year + c0), 4) &
         + (153*date%month - 1836*c0 - 457)/5 &
         + date%day + 1721117

  end subroutine date_to_daynum_julian

  subroutine date_to_daynum_gregorian(date, daynum)
    ! ---------------------------------------------------------------------------
    ! Convert from proleptic Gregorian calendar date to Chronological Julian Day
    ! Number (CJDN).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date ! Proleptic Gregorian calendar date.
    integer, intent(out) :: daynum      ! CJDN.

    integer :: c0, k1, q1

    c0 = intdivfloor(date%month - 3, 12)
    k1 = date%year + c0
    q1 = intdivfloor(k1, 100)

    daynum = intdivfloor(146097*q1, 4) + 36525*(k1 - q1*100)/100 &
         + (153*date%month - 1836*c0 - 457)/5 &
         + date%day + 1721119

  end subroutine date_to_daynum_gregorian

  subroutine date_to_daynum_365_day(date, daynum)
    ! ---------------------------------------------------------------------------
    ! Convert from 'noleap' or '365_day' calendar date to day number.
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date ! 'noleap' or '365_day' calendar date.
    integer, intent(out) :: daynum      ! Day number.

    integer :: c0

    c0 = intdivfloor(date%month - 3, 12)

    daynum = 365*(date%year + c0) &
         + (153*date%month - 1836*c0 - 457)/5 &
         + date%day - 307

  end subroutine date_to_daynum_365_day

  subroutine date_to_daynum_366_day(date, daynum)
    ! ---------------------------------------------------------------------------
    ! Convert from 'all_leap' or '366_day' calendar date to day number.
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date ! 'all_leap' or '366_day' calendar
    ! date.
    integer, intent(out) :: daynum      ! Day number.

    integer :: c0

    c0 = intdivfloor(date%month - 3, 12)

    daynum = 366*(date%year + c0) &
         + (153*date%month - 1836*c0 - 457)/5 &
         + date%day - 307

  end subroutine date_to_daynum_366_day

  subroutine date_to_daynum_360_day(date, daynum)
    ! ---------------------------------------------------------------------------
    ! Convert from '360_day' calendar date to day number.
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date ! '360_day' calendar date.
    integer, intent(out) :: daynum      ! Day number.

    daynum = 360*(date%year - 1) + 30*(date%month - 1) + date%day - 1

  end subroutine date_to_daynum_360_day

  subroutine daynum_to_date_julian(daynum, date)
    ! ---------------------------------------------------------------------------
    ! Convert from Chronological Julian Day Number (CJDN) to Julian calendar
    ! date.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: daynum        ! CJDN.
    type(date_type), intent(out) :: date ! Julian calendar date.

    integer :: k1, k2, q1, q2, c0

    k2 = 4*daynum - 6884469
    q2 = intdivfloor(k2, 1461)
    k1 = 5*((k2 - q2*1461)/4) + 2
    q1 = k1/153
    c0 = (q1 + 2)/12

    date = date_type(q2 + c0, q1 - 12*c0 + 3, (k1 - q1*153)/5 + 1)

  end subroutine daynum_to_date_julian

  subroutine daynum_to_date_gregorian(daynum, date)
    ! ---------------------------------------------------------------------------
    ! Convert from Chronological Julian Day Number (CJDN) to proleptic Gregorian
    ! calendar date.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: daynum        ! CJDN.
    type(date_type), intent(out) :: date ! Proleptic Gregorian calendar date.

    integer :: k1, k2, k3, q1, q2, q3, c0

    k3 = 4*daynum - 6884477
    q3 = intdivfloor(k3, 146097)
    k2 = 100*((k3 - q3*146097)/4) + 99
    q2 = k2/36525
    k1 = 5*((k2 - q2*36525)/100) + 2
    q1 = k1/153
    c0 = (q1 + 2)/12

    date = date_type(100*q3 + q2 + c0, q1 - 12*c0 + 3, (k1 - q1*153)/5 + 1)

  end subroutine daynum_to_date_gregorian

  subroutine daynum_to_date_365_day(daynum, date)
    ! ---------------------------------------------------------------------------
    ! Convert from day number to 'noleap' or '365_day' calendar date.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: daynum        ! Day number.
    type(date_type), intent(out) :: date ! 'noleap' or '365_day' calendar.
    ! date.

    integer :: k1, k2, q1, q2, c0

    k2 = daynum + 306
    q2 = intdivfloor(k2, 365)
    k1 = 5*(k2 - q2*365) + 2
    q1 = k1/153
    c0 = (q1 + 2)/12

    date = date_type(q2 + c0, q1 - 12*c0 + 3, (k1 - q1*153)/5 + 1)

  end subroutine daynum_to_date_365_day

  subroutine daynum_to_date_366_day(daynum, date)
    ! ---------------------------------------------------------------------------
    ! Convert from day number to 'all_leap' or '366_day' calendar date.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: daynum        ! Day number.
    type(date_type), intent(out) :: date ! 'all_leap' or '366_day' calendar.
    ! date.

    integer :: k1, k2, q1, q2, c0

    k2 = daynum + 306
    q2 = intdivfloor(k2, 366)
    k1 = 5*(k2 - q2*366) + 2
    q1 = k1/153
    c0 = (q1 + 2)/12

    date = date_type(q2 + c0, q1 - 12*c0 + 3, (k1 - q1*153)/5 + 1)

  end subroutine daynum_to_date_366_day

  subroutine daynum_to_date_360_day(daynum, date)
    ! ---------------------------------------------------------------------------
    ! Convert from day number to '360_day' calendar date.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: daynum        ! Day number.
    type(date_type), intent(out) :: date ! '360_day' calendar date.

    integer :: r

    date%year = intdivfloor(daynum, 360)
    r = daynum - date%year*360
    date%month = r/30
    date%day = r - date%month*30 + 1
    date%month = date%month + 1
    date%year = date%year + 1

  end subroutine daynum_to_date_360_day

  ! ---------------------------------------------------------------------------
  ! Public procedures.
  ! ---------------------------------------------------------------------------

  function date_to_daynum(calendar, date, daynum) result(errstat)
    ! ---------------------------------------------------------------------------
    ! Convert from calendar date to day number.
    ! ---------------------------------------------------------------------------

    character(len = *), intent(in) :: calendar ! Calendar type.
    type(date_type), intent(in) :: date        ! Calendar date.
    integer, intent(out) :: daynum             ! Day number.

    integer :: errstat                         ! Error status.

    type(date_type) :: date_tmp

    errstat = calendar_noerr

    select case (trim(calendar))
    case ('gregorian', 'standard')
      call date_to_daynum_gregorian(date, daynum)
      if (daynum > last_julian_daynum) then
        call daynum_to_date_gregorian(daynum, date_tmp)
        if (date_tmp /= date) then
          errstat = calendar_invalid_date
        endif
      else
        call date_to_daynum_julian(date, daynum)
        call daynum_to_date_julian(daynum, date_tmp)
        if (date_tmp /= date) then
          errstat = calendar_invalid_date
        else
          if (daynum > last_julian_daynum) then
            errstat = calendar_invalid_gregorian_date
          endif
        endif
      endif
    case ('proleptic_gregorian')
      call date_to_daynum_gregorian(date, daynum)
      call daynum_to_date_gregorian(daynum, date_tmp)
      if (date_tmp /= date) then
        errstat = calendar_invalid_date
      endif
    case ('julian')
      call date_to_daynum_julian(date, daynum)
      call daynum_to_date_julian(daynum, date_tmp)
      if (date_tmp /= date) then
        errstat = calendar_invalid_date
      endif
    case ('noleap', '365_day')
      call date_to_daynum_365_day(date, daynum)
      call daynum_to_date_365_day(daynum, date_tmp)
      if (date_tmp /= date) then
        errstat = calendar_invalid_date
      endif
    case ('all_leap', '366_day')
      call date_to_daynum_366_day(date, daynum)
      call daynum_to_date_366_day(daynum, date_tmp)
      if (date_tmp /= date) then
        errstat = calendar_invalid_date
      endif
    case ('360_day')
      call date_to_daynum_360_day(date, daynum)
      call daynum_to_date_360_day(daynum, date_tmp)
      if (date_tmp /= date) then
        errstat = calendar_invalid_date
      endif
    case default
      daynum = 0
      errstat = calendar_unsupported
    end select

  end function date_to_daynum

  function daynum_to_date(calendar, daynum, date) result(errstat)
    ! ---------------------------------------------------------------------------
    ! Convert from day number to calendar date.
    ! ---------------------------------------------------------------------------

    character(len = *), intent(in) :: calendar ! Calendar type.
    integer, intent(in) :: daynum              ! Day number.
    type(date_type), intent(out) :: date       ! Calendar date.

    integer :: errstat                         ! Error status.

    integer :: daynum_tmp

    errstat = calendar_noerr

    select case (trim(calendar))
    case ('gregorian', 'standard')
      if (daynum > last_julian_daynum) then
        call daynum_to_date_gregorian(daynum, date)
        call date_to_daynum_gregorian(date, daynum_tmp)
      else
        call daynum_to_date_julian(daynum, date)
        call date_to_daynum_julian(date, daynum_tmp)
      endif
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case ('proleptic_gregorian')
      call daynum_to_date_gregorian(daynum, date)
      call date_to_daynum_gregorian(date, daynum_tmp)
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case ('julian')
      call daynum_to_date_julian(daynum, date)
      call date_to_daynum_julian(date, daynum_tmp)
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case ('noleap', '365_day')
      call daynum_to_date_365_day(daynum, date)
      call date_to_daynum_365_day(date, daynum_tmp)
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case ('all_leap', '366_day')
      call daynum_to_date_366_day(daynum, date)
      call date_to_daynum_366_day(date, daynum_tmp)
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case ('360_day')
      call daynum_to_date_360_day(daynum, date)
      call date_to_daynum_360_day(date, daynum_tmp)
      if (daynum_tmp /= daynum) then
        errstat = calendar_daynum_overflow
      endif
    case default
      date = date_type(0, 0, 0)
      errstat = calendar_unsupported
    end select

  end function daynum_to_date

  function daynum_diff(calendar, date1, date2, dndiff) result(errstat)
    ! ---------------------------------------------------------------------------
    ! Calculate number of days between two calendar dates.
    ! ---------------------------------------------------------------------------

    character(len = *), intent(in) :: calendar ! Calendar type.
    type(date_type), intent(in) :: date1       ! First date.
    type(date_type), intent(in) :: date2       ! Second date.
    integer, intent(out) :: dndiff             ! Day number difference.

    integer :: errstat                         ! Error status.

    integer :: errstat1, errstat2, dn1, dn2

    errstat1 = date_to_daynum(calendar, date1, dn1)
    errstat2 = date_to_daynum(calendar, date2, dn2)

    if (errstat1 == calendar_noerr .and. errstat2 == calendar_noerr) then
      dndiff = dn2 - dn1
      errstat = calendar_noerr
    else
      dndiff = 0
      errstat = min(errstat1, errstat2)
    endif

  end function daynum_diff

  function date_offset(calendar, date, dnoffset) result(errstat)
    ! ---------------------------------------------------------------------------
    ! Offset a date by a specified number of days.
    ! ---------------------------------------------------------------------------

    character(len = *), intent(in) :: calendar ! Calendar type.
    type(date_type), intent(inout) :: date     ! Calendar date.
    integer, intent(in) :: dnoffset            ! Day number offset.

    integer :: errstat                         ! Error status.

    integer :: dn

    errstat = date_to_daynum(calendar, date, dn)

    if (errstat == calendar_noerr) then
      dn = dn + dnoffset
      errstat = daynum_to_date(calendar, dn, date)
      if (errstat == calendar_daynum_overflow) then
        errstat = calendar_daynum_offset_overflow
      endif
    endif

  end function date_offset

  function date_check(calendar, date) result(errstat)
    ! ---------------------------------------------------------------------------
    ! Check validity of date.
    ! ---------------------------------------------------------------------------

    character(len = *), intent(in) :: calendar ! Calendar type.
    type(date_type), intent(inout) :: date     ! Calendar date.

    integer :: errstat                         ! Error status.

    integer :: daynum

    errstat = date_to_daynum(calendar, date, daynum)

  end function date_check

  pure function calendar_errstr(errstat) result(errstr)
    ! ---------------------------------------------------------------------------
    ! Returns static reference to an error message string corresponding to a
    ! calendar error status.
    ! ---------------------------------------------------------------------------

    integer, intent(in) :: errstat ! Error status.

    character(len = 80) :: errstr  ! Error message string.

    if (errstat > 0 .and. errstat <= errmsg_num) then
      errstr = errmsg(errstat)
    else
      errstr = 'Unknown error status!'
    endif

  end function calendar_errstr

  pure logical function dates_equal(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if dates are equal. For overloading operator (==).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    dates_equal = date1%year  == date2%year  .and. &
         date1%month == date2%month .and. &
         date1%day   == date2%day

  end function dates_equal

  pure logical function date1_lt_date2(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if date1 is less than date2. For overloading operator (<).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    date1_lt_date2 =  date1%year  <  date2%year   .or.  &
         (date1%year  == date2%year   .and. &
         date1%month <  date2%month) .or.  &
         (date1%year  == date2%year   .and. &
         date1%month == date2%month  .and. &
         date1%day   <  date2%day)

  end function date1_lt_date2

  pure logical function date1_gt_date2(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if date1 is greater than date2. For overloading operator (>).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    date1_gt_date2 =  date1%year  >  date2%year   .or.  &
         (date1%year  == date2%year   .and. &
         date1%month >  date2%month) .or.  &
         (date1%year  == date2%year   .and. &
         date1%month == date2%month  .and. &
         date1%day   >  date2%day)

  end function date1_gt_date2

  pure logical function dates_not_equal(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if dates are not equal. For overloading operator (/=).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    dates_not_equal = .not. dates_equal(date1, date2)

  end function dates_not_equal

  pure logical function date1_le_date2(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if date1 is less or equal date2. For overloading operator
    ! (<=).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    date1_le_date2 = .not. date1_gt_date2(date1, date2)

  end function date1_le_date2

  pure logical function date1_ge_date2(date1, date2)
    ! ---------------------------------------------------------------------------
    ! Returns true if date1 is greater or equal date2. For overloading operator
    ! (>=).
    ! ---------------------------------------------------------------------------

    type(date_type), intent(in) :: date1, date2

    date1_ge_date2 = .not. date1_lt_date2(date1, date2)

  end function date1_ge_date2

end module mod_calendar
