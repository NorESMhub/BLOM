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

module mod_time
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to time.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_config, only: expcnf
   use mod_constants, only: epsil
   use mod_calendar, only: date_type, daynum_diff, date_offset, &
                           calendar_noerr, calendar_errstr
   use mod_xc, only: lp, mnproc, xcstop

   implicit none

   private

   type(date_type) :: &
      date0, &        ! Experiment start date.
      date            ! Current experiment date.

   character(len = 19) :: &
      calendar        ! The calendar used (see mod_calendar.F90 for supported
                      ! calendars)

   integer :: &
      nday1, &        ! First day number of current integration, counting from
                      ! experiment start.
      nday2, &        ! Last day number of current integration, counting from
                      ! experiment start.
      nday_in_year, & ! Number of days in current year.
      nday_of_year, & ! Current number of days since start of year.
      nstep0, &       ! Time step number at experiment start. Always zero in
                      ! standalone configurations and non-zero when coupled to
                      ! CESM.
      nstep1, &       ! First time step number of current integration, counting
                      ! from experiment start.
      nstep2, &       ! Last time step number of current integration, counting
                      ! from experiment start (only used in standalone
                      ! configurations)
      nstep, &        ! Current time step number, counting from experiment
                      ! start.
      lstep, &        ! Number of barotropic time steps per baroclinic time
                      ! step.
      nstep_in_day    ! Number of time steps in a day.

   real(r8) :: &
      time0, &        ! Integration time in days at experiment start ('date0').
      time, &         ! Current integration time in days, including integration
                      ! time before experiment start.
      baclin, &       ! Time step between baroclinic time levels.
      batrop, &       ! Requested barotropic time step.
      delt1, &        ! Baroclinic time step equal to 'baclin' the first time
                      ! step (forward) and equal to '2*baclin' otherwise
                      ! (leap-frog).
      dlt             ! Resolved barotropic time step.

   ! Interpolation parameters for monthly climatological fields.
   real(r8) :: xmi
   integer :: l1mi, l2mi, l3mi, l4mi, l5mi

   public :: date0, date, calendar, nday1, nday2, nday_in_year, nday_of_year, &
             nstep0, nstep1, nstep2, nstep, lstep, nstep_in_day, time0, time, &
             baclin, batrop, delt1, dlt, xmi, l1mi, l2mi, l3mi, l4mi, l5mi, &
             init_timevars, set_day_of_year, step_time, blom_time

contains

   subroutine init_timevars
   ! ---------------------------------------------------------------------------
   ! Initialize various time variables. Assumes that 'baclin' and 'date' are
   ! already defined.
   ! ---------------------------------------------------------------------------

      ! Set calendar type to be used.
      select case (trim(expcnf))
         case ('cesm')               
            calendar = 'noleap'
         case ('ben02clim')          
            calendar = '360_day'
         case ('ben02syn')           
            calendar = 'standard'
         case ('fuk95')           
            calendar = '360_day'
         case ('channel')
            calendar = '360_day'
         case ('single_column')           
            calendar = '360_day'
         case ('isomip1', 'isomip2') 
            calendar = '360_day'
         case default
            if (mnproc == 1) then
               write (lp,'(3a)') ' init_timevars: expcnf = ', trim(expcnf), &
                                 ' is unsupported!'
            endif
            call xcstop('(init_timevars)')
                   stop '(init_timevars)'
      end select

      ! Get number of baroclinic time steps per day and verify that an integer
      ! number of steps fits in a day.
      nstep_in_day = nint(86400._r8/baclin)
      if (abs(86400._r8/baclin - nstep_in_day) > epsil) then
         if (mnproc == 1) then
            write (lp, *) &
               'init_timevars: '// &
               'must have an integer number of baroclinic time steps pr. day!'
         endif
         call xcstop('(init_timevars)')
                stop '(init_timevars)'
      endif
      if (mnproc == 1) then
         write (lp, *) &
            'init_timevars: number of baroclinic time steps per day: ', &
            nstep_in_day
      endif

      ! Number of barotropic time steps per baroclinic time step, 'lstep', must
      ! be even.
      lstep = 2*ceiling(.5_r8*baclin/batrop)

      ! Set barotropic time step.
      dlt = baclin/lstep

      ! Set number of days in year and number of days since start of year.
      call set_day_of_year

   end subroutine init_timevars

   subroutine set_day_of_year
   ! ---------------------------------------------------------------------------
   ! Set number of days in current year ('nday_in_year') and current number of
   ! days since start of year ('nday_of_year').
   ! ---------------------------------------------------------------------------

      integer :: errstat

      ! Set number of days in current year.
      errstat = daynum_diff(calendar, &
                            date_type(date%year    , 1, 1), &
                            date_type(date%year + 1, 1, 1), &
                            nday_in_year)
      if (errstat /= calendar_noerr) then
         write (lp, '(2a)') ' set_day_of_year: daynum_diff error: ', &
                            trim(calendar_errstr(errstat))
         call xcstop('(set_day_of_year)')
                stop '(set_day_of_year)'
      endif
 
      ! Set current number of days since start of year.
      errstat = daynum_diff(calendar, date_type(date%year, 1, 1), date, &
                            nday_of_year)
      if (errstat /= calendar_noerr) then
         write (lp, '(2a)') ' set_day_of_year: daynum_diff error: ', &
                            trim(calendar_errstr(errstat))
         call xcstop('(set_day_of_year)')
                stop '(set_day_of_year)'
      endif
      nday_of_year = nday_of_year + 1

   end subroutine set_day_of_year

   subroutine step_time
   ! ---------------------------------------------------------------------------
   ! Increment time step and ensure all time variables are updated accordingly.
   ! ---------------------------------------------------------------------------

      integer :: errstat

      ! Increment time step counter and set integration time.
      nstep = nstep + 1
      time = time0 + nstep*baclin/86400._r8

      if (mod(nstep, nstep_in_day) == 0) then

         ! Increment date.
         errstat = date_offset(calendar, date, 1)
         if (errstat /= calendar_noerr) then
            write (lp, '(2a)') ' step_time: date_offset error: ', &
                               trim(calendar_errstr(errstat))
            call xcstop('(step_time)')
                   stop '(step_time)'
         endif

         call set_day_of_year

      endif

      ! Set parameters for time interpolation of climatological fields.
      xmi = ( nday_of_year - 1 &
            + mod(nstep, nstep_in_day)/real(nstep_in_day, r8)) &
            *12._r8/real(nday_in_year, r8)
      l3mi = int(xmi) + 1
      xmi = xmi - real(l3mi-1, r8)
      l1mi = mod(l3mi +  9,12) + 1
      l2mi = mod(l3mi + 10,12) + 1
      l4mi = mod(l3mi     ,12) + 1
      l5mi = mod(l3mi +  1,12) + 1

   end subroutine step_time

   subroutine blom_time(ymd, tod)
   ! ---------------------------------------------------------------------------
   ! Through arguments, provide current date in integer format ('ymd') and time
   ! of day in seconds ('tod').
   ! ---------------------------------------------------------------------------

      integer, intent(out) :: ymd, tod

      ymd = date%year*10000 + date%month*100 + date%day
      tod = nint(mod(nstep, nstep_in_day)*baclin)

   end subroutine blom_time

end module mod_time
