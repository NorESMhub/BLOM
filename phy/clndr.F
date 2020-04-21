! ------------------------------------------------------------------------------
! Copyright (C) 2004-2017 Mats Bentsen
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

      subroutine inicln
c
c --- ------------------------------------------------------------------
c --- Initialize the calendar
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
c
      integer i
c
c --- Set calendar type to be used
      if     (expcnf.eq.'cesm') then
        calendar='noleap'
      elseif (expcnf.eq.'ben02clim') then
        calendar='360_day'
      elseif (expcnf.eq.'ben02syn') then
        calendar='standard'
      elseif (expcnf.eq.'isomip1'.or.expcnf.eq.'isomip2') then
        calendar='360_day'
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(clndr)')
               stop '(clndr)'
      endif
c
c --- Set days in months and days in the year
      if     (calendar(1:3).eq.'360') then
        nd_in_m( 1)=30
        nd_in_m( 2)=30
        nd_in_m( 3)=30
        nd_in_m( 4)=30
        nd_in_m( 5)=30
        nd_in_m( 6)=30
        nd_in_m( 7)=30
        nd_in_m( 8)=30
        nd_in_m( 9)=30
        nd_in_m(10)=30
        nd_in_m(11)=30
        nd_in_m(12)=30
        nday_in_year=360
      elseif (calendar(1:3).eq.'sta') then
        nd_in_m( 1)=31
        if (mod(nyear,4).eq.0.and.
     .      (mod(nyear,100).ne.0.or.mod(nyear,400).eq.0)) then
          nd_in_m( 2)=29
          nday_in_year=366
        else
          nd_in_m( 2)=28
          nday_in_year=365
        endif
        nd_in_m( 3)=31
        nd_in_m( 4)=30
        nd_in_m( 5)=31
        nd_in_m( 6)=30
        nd_in_m( 7)=31
        nd_in_m( 8)=31
        nd_in_m( 9)=30
        nd_in_m(10)=31
        nd_in_m(11)=30
        nd_in_m(12)=31
      elseif (calendar(1:3).eq.'nol') then
        nd_in_m( 1)=31
        nd_in_m( 2)=28
        nd_in_m( 3)=31
        nd_in_m( 4)=30
        nd_in_m( 5)=31
        nd_in_m( 6)=30
        nd_in_m( 7)=31
        nd_in_m( 8)=31
        nd_in_m( 9)=30
        nd_in_m(10)=31
        nd_in_m(11)=30
        nd_in_m(12)=31
        nday_in_year=365
      endif
c
c --- Get day of year
      nday_of_year=0
      do i=1,nmonth-1
        nday_of_year=nday_of_year+nd_in_m(i)
      enddo
      nday_of_year=nday_of_year+nday
c
      return
      end
c
      subroutine updcln(nstep)
c
c --- ------------------------------------------------------------------
c --- Update the calendar
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_clndr.h"
c
      integer nstep
c
      if (mod(nstep,nstep_in_day).eq.0) then
c
c --- - get new date
        nday=nday+1
        nday_of_year=nday_of_year+1
        if (nday.gt.nd_in_m(nmonth)) then
          nday=1
          nmonth=nmonth+1
          if (nmonth.gt.12) then
            nday_of_year=1
            nmonth=1
            nyear=nyear+1
            if (calendar(1:3).eq.'sta') then
              if (mod(nyear,4).eq.0.and.
     .            (mod(nyear,100).ne.0.or.mod(nyear,400).eq.0)) then
                nd_in_m(2)=29
                nday_in_year=366
              else
                nd_in_m(2)=28
                nday_in_year=365
              endif
            endif
          endif
        endif
c
      endif
c
      if ((calendar(1:3).eq.'sta'.or.calendar(1:3).eq.'mix'.or.
     .     calendar(1:3).eq.'gre').and.nyear.le.1582) then
        if (mnproc.eq.1) then
          write (lp,*)
     .'Do not use mixed Julian/Gregorian calendar before Oct 10th 1582!'
        endif
        call xcstop('(updcln)')
      endif
c
      return
      end
c
      integer function daydif(nyears,nmonths,ndays,nyeare,nmonthe,ndaye)
c
c --- ------------------------------------------------------------------
c --- return number of days between dates
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_clndr.h"
c
      integer nyears,nmonths,ndays,nyeare,nmonthe,ndaye
c
      integer nd_in_m2,ny1,nm1,nd1,ny2,nm2,nd2
c
      if (nyears.lt.nyeare.or.
     .    (nyears.eq.nyeare.and.nmonths.lt.nmonthe).or.
     .    (nyears.eq.nyeare.and.nmonths.eq.nmonthe.and.ndays.lt.ndaye))
     .  then
        ny1=nyears
        nm1=nmonths
        nd1=ndays
        ny2=nyeare
        nm2=nmonthe
        nd2=ndaye
      else
        ny1=nyeare
        nm1=nmonthe
        nd1=ndaye
        ny2=nyears
        nm2=nmonths
        nd2=ndays
      endif
c
      nd_in_m2=nd_in_m(2)
      if (calendar(1:3).eq.'sta') then
        if (mod(ny1,4).eq.0.and.
     .      (mod(ny1,100).ne.0.or.mod(ny1,400).eq.0)) then
          nd_in_m( 2)=29
        endif
      endif
c
      daydif=0
c
      do while (ny1.lt.ny2.or.(ny1.eq.ny2.and.nm1.lt.nm2).or.
     .          (ny1.eq.ny2.and.nm1.eq.nm2.and.nd1.lt.nd2))
        daydif=daydif+1
        nd1=nd1+1
        if (nd1.gt.nd_in_m(nm1)) then
          nd1=1
          nm1=nm1+1
          if (nm1.gt.12) then
            nm1=1
            ny1=ny1+1
            if (calendar(1:3).eq.'sta') then
              if (mod(ny1,4).eq.0.and.
     .            (mod(ny1,100).ne.0.or.mod(ny1,400).eq.0)) then
                nd_in_m(2)=29
              else
                nd_in_m(2)=28
              endif
            endif
          endif
        endif
      enddo
c
      nd_in_m(2)=nd_in_m2
c
      return
      end
c
      subroutine pstdat(ny,nm,nd,ndays)
c
c --- ------------------------------------------------------------------
c --- return date  ndays  in the past
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_clndr.h"
c
      integer ny,nm,nd,ndays
c
      integer nd_in_m2,i
c
      nd_in_m2=nd_in_m(2)
      if (calendar(1:3).eq.'sta') then
        if (mod(ny,4).eq.0.and.
     .      (mod(ny,100).ne.0.or.mod(ny,400).eq.0)) then
          nd_in_m( 2)=29
        else
          nd_in_m( 2)=28
        endif
      endif
c
c --- Go back ndays
      do i=1,ndays
c
        nd=nd-1
        if (nd.eq.0) then
          nm=nm-1
          if (nm.eq.0) then
            ny=ny-1
            if (calendar(1:3).eq.'sta') then
              if (mod(ny,4).eq.0.and.
     .            (mod(ny,100).ne.0.or.mod(ny,400).eq.0)) then
                nd_in_m(2)=29
              else
                nd_in_m(2)=28
              endif
            endif
            nm=12
          endif
          nd=nd_in_m(nm)
        endif
c
      enddo
c
      nd_in_m(2)=nd_in_m2
c
      return
      end
