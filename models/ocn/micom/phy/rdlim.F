      subroutine rdlim(nday1,nday2)
c
c --- ------------------------------------------------------------------
c --- Read limits file
c --- ------------------------------------------------------------------
c
      use mod_xc
#ifdef CCSMCOUPLED
      use data_mct, only: runid_mct, runtype_mct
#endif
c
      implicit none
c
      integer nday1,nday2
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
#include "common_forc.h"
#ifdef SYN
#  include "param_atm.h"
#  include "common_atmsyn.h"
#endif
c
      integer i,j,m,n
#ifdef CCSMCOUPLED
      character*80 rstfnm
      logical fexist
#endif
c
#ifdef CCSMCOUPLED
      open (unit=nfu,file='ocn_in',status='old',action='read')
#else
      open (unit=nfu,file='limits',status='old',action='read')
#endif
      read (nfu,'(i6)') nday1       ! First day of integration
      read (nfu,'(i6)') nday2       ! Last day of integration
      read (nfu,'(i4,x,i2,x,i2)') nyear ,nmonth ,nday  ! Model date
      read (nfu,'(i4,x,i2,x,i2)') nyear0,nmonth0,nday0 ! Initial experiment date
      read (nfu,'(a)') runid        ! Experiment name
      read (nfu,'(f8.2)') baclin    ! Baroclinic time step
      read (nfu,'(f8.2)') batrop    ! Barotropic time step
      read (nfu,'(i3)' ) i          ! time between 2D diagnostics
      read (nfu,'(i3)' ) j          ! time between 3D diagnostics
      read (nfu,'(f8.2)') trxday    ! Relax. time scale for SST, if 0 no relax.
      read (nfu,'(f8.2)') srxday    ! Relax. time scale for SSS, if 0 no relax.
      read (nfu,'(l5)') aptflx      ! Apply diagnosed heat flux flag
      read (nfu,'(l5)') apsflx      ! Apply diagnosed freshwater flux flag
      read (nfu,'(l5)') ditflx      ! Diagnose heat flux flag
      read (nfu,'(l5)') disflx      ! Diagnose freshwater flux flag
      read (nfu,'(a)') path         ! Path to input files
      read (nfu,'(a)') path1        ! Path to diagnostic files
      read (nfu,'(a)') path2        ! Path to restart files
#ifdef SYN
      read (nfu,'(a)') atm_path     ! Path to synoptic atmospheric forcing files
#endif
      close (unit=nfu)
c
#ifdef CCSMCOUPLED
c --- override some information read from the micom namelist with
c --- information received from the coupler and from the restart file
c --- pointer if it exists
      runid(1:80)=runid_mct(1:80)
      inquire(file='rpointer.ocn',exist=fexist)
      if (fexist) then
        open (unit=nfu,file='rpointer.ocn')
        read (nfu,'(a)') rstfnm
        close (unit=nfu)
        n=1
        do while (n.lt.80.and.rstfnm(n:n).ne.' ')
          n=n+1
        enddo
        n=n-1
        if (runtype_mct.eq.'initial') then
          read (rstfnm(n-8 :n-3 ),'(i6)') nday1
        else
          read (rstfnm(n-19:n-16),'(i4)') nyear
          read (rstfnm(n-14:n-13),'(i2)') nmonth
          read (rstfnm(n-11:n-10),'(i2)') nday
          read (rstfnm(n-8 :n-3 ),'(i6)') nday1
        endif
      endif
#endif
c
      if (nday1.lt.0.or.nday2.lt.0) then
        if (mnproc.eq.1) then
          write (lp,*) 'Integrations days must be positive!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nday2.lt.nday1) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' First day of integration must be less ',
     .      'than or equal to last day!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nday1.eq.0.and.
     .    (nyear.ne.nyear0.or.nmonth.ne.nmonth0.or.nday.ne.nday0)) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' When first integration day is zero, ',
     .      'model date and initial experiment date '
          write (lp,'(a)') ' must be equal!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (nyear.lt.nyear0.or.(nyear.eq.nyear0.and.nmonth.lt.nmonth0).or.
     .  (nyear.eq.nyear0.and.nmonth.eq.nmonth0.and.nday.lt.nday0)) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' Model date must be greater or equal to ',
     .      'initial experiment date!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (mod(86400./baclin+epsil,1.).gt.2.*epsil) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' Must have an integer number of ',
     .     'baroclinic time steps pr. day!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (trxday.eq.0.and.ditflx) then
        if (mnproc.eq.1) then
          write (lp,*) 'trxday=0. and ditflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (srxday.eq.0.and.disflx) then
        if (mnproc.eq.1) then
          write (lp,*) 'srxday=0. and disflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
c --- trim the experiment id and paths
#ifndef CCSMCOUPLED
      n=1
      do while (runid(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          if (mnproc.eq.1) then
            write (lp,'(2a)') 'Place a colon after the end of useable ',
     .        'characters in runid!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
      enddo
      do m=n,80
        runid(m:m)=' '
      enddo
#endif
c
      n=1
      do while (path(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          if (mnproc.eq.1) then
            write (lp,'(2a)') 'Place a colon after the end of useable ',
     .        'characters in path!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
      enddo
      do m=n,80
        path(m:m)=' '
      enddo
c
      n=1
      do while (path1(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          if (mnproc.eq.1) then
            write (lp,'(2a)') 'Place a colon after the end of useable ',
     .        'characters in path1!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
      enddo
      do m=n,80
        path1(m:m)=' '
      enddo
c
      n=1
      do while (path2(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          if (mnproc.eq.1) then
            write (lp,'(2a)') 'Place a colon after the end of useable ',
     .        'characters in path2!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
      enddo
      do m=n,80
        path2(m:m)=' '
      enddo
#ifdef SYN
c
      n=1
      do while (atm_path(n:n).ne.':')
        n=n+1
        if (n.gt.80) then
          if (mnproc.eq.1) then
            write (lp,'(2a)') 'Place a colon after the end of useable ',
     .        'characters in atm_path!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
      enddo
      do m=n,80
        atm_path(m:m)=' '
      enddo
#endif
c
c --- get time steps in a day
      nstep_in_day=nint(86400./baclin)
c
c --- represent time between diagnostics in time steps
      if (i.lt.0) then
        diag2dfq=-real(nstep_in_day)/i
      else
        diag2dfq=nstep_in_day*max(1,i)
      endif
      if (i.eq.30) then
        diag2dmon=.true.
      else
        diag2dmon=.false.
      endif
c
      if (j.lt.0) then
        diag3dfq=-real(nstep_in_day)/j
      else
        diag3dfq=nstep_in_day*max(1,j)
      endif
      if (j.eq.30) then
        diag3dmon=.true.
      else
        diag3dmon=.false.
      endif
c
      if (mnproc.eq.1) then
        write (lp,'(a,i10)')   'nday1:             ',nday1
        write (lp,'(a,i10)')   'nday2:             ',nday2
        write (lp,'(a,i4.4,a,i2.2,a,i2.2)') 'nyear ,nmonth ,nday : ',
     .                             nyear ,'.',nmonth ,'.',nday
        write (lp,'(a,i4.4,a,i2.2,a,i2.2)') 'nyear0,nmonth0,nday0: ',
     .                             nyear0,'.',nmonth0,'.',nday0
        write (lp,'(a,a60)')   'runid:             ',runid
        write (lp,'(a,f10.4)') 'baclin:            ',baclin
        write (lp,'(a,f10.4)') 'batrop:            ',batrop
        if (diag2dmon) then
          write (lp,'(a,f10.4)') '2D diag. interval:    monthly'
        else
          write (lp,'(a,f10.4)') '2D diag. interval: ',
     .                             diag2dfq/nstep_in_day
        endif
        if (diag3dmon) then
          write (lp,'(a,f10.4)') '3D diag. interval:    monthly'
        else
          write (lp,'(a,f10.4)') '3D diag. interval: ',
     .                             diag3dfq/nstep_in_day
        endif
        write (lp,'(a,f10.4)') 'trxday:            ',trxday
        write (lp,'(a,f10.4)') 'srxday:            ',srxday
        write (lp,'(a,l10)')   'aptflx:            ',aptflx
        write (lp,'(a,l10)')   'apsflx:            ',apsflx
        write (lp,'(a,l10)')   'ditflx:            ',ditflx
        write (lp,'(a,l10)')   'disflx:            ',disflx
        write (lp,'(a,a60)')   'path:      ',path
        write (lp,'(a,a60)')   'path1:     ',path1
        write (lp,'(a,a60)')   'path2:     ',path2
#ifdef SYN
        write (lp,'(a,a60)')   'atm_path:  ',atm_path
#endif
c
        write (lp,101) thkdff,temdff,veldff,viscos,diapyc,vertmx
 101    format (' turb. flux parameters:',1p/
     .    ' thkdff,temdff,veldff =',3e9.2/
     .    ' viscos,diapyc,vertmx =',3e9.2)
        call flush(lp)
      endif
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
      lstep=2*ceiling(.5*baclin/batrop)
      dlt=baclin/lstep
      if (mnproc.eq.1) then
        write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .    lstep
        call flush(lp)
      endif
c
#ifndef CCSMCOUPLED
c --- 'nfrco' = number of baroclinic time steps between surface forcing
c --- updates
#  if   defined(CLIM)
      nfrco=nstep_in_day  ! once pr. day with climatological forcing
#  elif defined(SYN)
      nfrco=1             ! every time step with synoptic forcing
#  endif
#endif
c
c --- set experiment id and path lengths
      runid_len=1
      do while (runid_len.lt.80.and.runid(runid_len:runid_len).ne.' ')
        runid_len=runid_len+1
      enddo
      runid_len=runid_len-1
      path_len=1
      do while (path_len.lt.80.and.path(path_len:path_len).ne.' ')
        path_len=path_len+1
      enddo
      path_len=path_len-1
      path1_len=1
      do while (path1_len.lt.80.and.path1(path1_len:path1_len).ne.' ')
        path1_len=path1_len+1
      enddo
      path1_len=path1_len-1
      path2_len=1
      do while (path2_len.lt.80.and.path2(path2_len:path2_len).ne.' ')
        path2_len=path2_len+1
      enddo
      path2_len=path2_len-1
#ifdef SYN
      atm_path_len=1
      do while (atm_path_len.lt.80.and.
     .          atm_path(atm_path_len:atm_path_len).ne.' ')
        atm_path_len=atm_path_len+1
      enddo
      atm_path_len=atm_path_len-1
#endif
c
c --- model is to be integrated from time step 'nstep1' to 'nstep2'
      nstep1=nday1*nstep_in_day
      nstep2=nday2*nstep_in_day
c
      if (csdiag) then
        nstep2=nstep1+2
      endif
c
      return
      end
