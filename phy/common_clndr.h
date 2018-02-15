c --- ------------------------------------------------------------------
c --- common block related to model calendar and time management
c --- ------------------------------------------------------------------
c
      real time0
      integer, dimension(12) :: nd_in_m
      integer nday1,nday2,nday,nmonth,nyear,nday0,nmonth0,nyear0,
     .        nday_in_year,nday_of_year,nstep_in_day
      character*19 calendar
c
      common /clndr/ time0,nd_in_m,nday1,nday2,nday,nmonth,nyear,nday0,
     .               nmonth0,nyear0,nday_in_year,nday_of_year,
     .               nstep_in_day,calendar
