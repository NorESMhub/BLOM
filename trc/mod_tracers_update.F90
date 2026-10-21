! ------------------------------------------------------------------------------
! Copyright (C) 2007-2024 Mats Bentsen, Jörg Schwinger, Jerry Tjiputra,
!                         Alok Kumar Gupta, Mariana Vertenstein
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

module mod_tracers_update

  use mod_types,           only: r8
  use mod_grid,            only: plon, plat
  use mod_tracers,         only: ntrocn, natr, ntr, trc
  use mod_idlage,          only: idlage_init
  use mod_idlage,          only: idlage_step
#ifdef HAMOCC
  use mo_param1_bgc,       only: init_indices, nocetra
  use mo_hamocc_init,      only: hamocc_init
  use mo_hamocc_step,      only: hamocc_step
  use mo_restart_hamoccwt, only: restart_hamoccwt
#endif
  use mod_constants,       only: spval
  use mod_calendar,        only: date_type, operator(/=)
  use mod_time,            only: date0, time0, time
  use mod_dia,             only: iotype, rstfmt, rstcmp
  use mod_config,          only: expcnf
  use mod_ifdefs,          only: use_ATRC, use_IDLAGE
  use mod_utility,         only: fnmlen
  use mod_nctools
  use mod_xc

  implicit none
  private

  ! Public routines
  public :: initrc
  public :: updtrc
  public :: restart_trcwt
  public :: restart_trcrd

  ! Private routines
  private :: restart_ocntrcwt
  private :: restart_ocntrcrd
  private :: restart_getfile

contains

  subroutine initrc()
    ! ------------------------------------------------------------------
    ! initialization of ocean tracers
    ! ------------------------------------------------------------------

    integer :: i,j,k,l,nt,nat

    ! ------------------------------------------------------------------
    ! if no ocean tracers are defined
    ! ------------------------------------------------------------------

    if (ntrocn /= 0) then

      ! ------------------------------------------------------------------
      ! if number of age tracers is greater than zero, CPP flag ATRC must
      ! be defined
      ! ------------------------------------------------------------------
      if (.not. use_ATRC) then
        if (natr > 0) then
          if (mnproc == 1) then
            write (lp,'(2a)') ' Since number of age tracers is greater ', &
                 'than zero, ATRC must be defined!'
          end if
          call xcstop('(ocntrc_init)')
          stop '(ocntrc_init)'
        end if
      end if

      ! ------------------------------------------------------------------
      ! check number of age tracers
      ! ------------------------------------------------------------------

      if (natr < 0.or.2*natr > ntrocn) then
        if (mnproc == 1) then
          write (lp,'(3a)') ' Number of age tracers must be greater ', &
               'than zero and less or equal half the total number of ', &
               'ocean tracers!'
        end if
        call xcstop('(ocntrc_init)')
        stop '(ocntrc_init)'
      end if

      ! ------------------------------------------------------------------
      ! initialization of tracers
      ! ------------------------------------------------------------------

      do nt = 1,ntrocn-natr
        !$omp parallel do private(k,l,i)
        do j = 1,jj
          do k = 1,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                trc(i,j,k,nt)= &
                     (mod(k,5)+1)*(plat(i,j)+90._r8)/(5._r8*180._r8)+nt
                trc(i,j,k+kk,nt) = trc(i,j,k,nt)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      end do

      ! ------------------------------------------------------------------
      ! initialization of age tracers
      ! ------------------------------------------------------------------

      do nt = 1,natr
        nat = ntr-natr+nt
        !$omp parallel do private(k,l,i)
        do j = 1,jj
          do k = 1,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                trc(i,j,k,nat)= &
                     (mod(k,5)*(plon(i,j)+180._r8)/(4._r8*360._r8)+nat) &
                     *trc(i,j,k,nt)
                trc(i,j,k+kk,nat) = trc(i,j,k,nat)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      end do

    end if ! if (ntrocn != 0)

#ifdef HAMOCC
    call hamocc_init(0,'c')
#endif
    if (use_IDLAGE) then
      call idlage_init
    end if

  end subroutine initrc

  ! ============================================================================

  subroutine updtrc(m,n,mm,nn,k1m,k1n)
    ! ------------------------------------------------------------------
    ! update tracers due to non-passive processes
    ! ------------------------------------------------------------------
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

#ifdef HAMOCC
    call hamocc_step(m,n,mm,nn,k1m,k1n)
#endif
#ifndef OFFLINE_SEDIMENT_SPINUP
    if (use_IDLAGE) then
      call idlage_step(m,n,mm,nn,k1m,k1n)
    end if
#endif
  end subroutine updtrc

  ! ============================================================================

  subroutine restart_trcwt(rstfnm_ocn)
    ! ------------------------------------------------------------------
    ! Write tracer state to restart files
    ! ------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: rstfnm_ocn

    ! Local variables
    logical :: error
    character(len=fnmlen) :: rstfnm_ocntrc
    character(len=fnmlen) :: rstfnm_hamocc

    ! ------------------------------------------------------------------
    ! Generate file name
    ! ------------------------------------------------------------------

    if (mnproc == 1) then
      if (expcnf == 'cesm') then
        call restart_getfile(rstfnm_ocn, 'rtrc', rstfnm_ocntrc, error)
        if (error) then
          write(lp,*) 'restart_trcwt: could not generate rstfnm_ocntrc file!'
          call xcstop('(restat_trcwt)')
          stop '(restart_trcwt)'
        endif
        call restart_getfile(rstfnm_ocn, 'rbgc', rstfnm_hamocc, error)
        if (error) then
          write(lp,*) 'restart_trcwt: could not generate rstfnm_hamocc file!'
          call xcstop('(restat_trcwt)')
          stop '(restart_trcwt)'
        endif
      else
        call restart_getfile(rstfnm_ocn, 'resttrc', rstfnm_ocntrc, error)
        if (error) then
          write(lp,*) 'restart_trcwt: could not generate rstfnm_ocntrc file!'
          call xcstop('(restat_trcwt)')
          stop '(restart_trcwt)'
        endif
        call restart_getfile(rstfnm_ocn, 'restbgc', rstfnm_hamocc, error)
        if (error) then
          write(lp,*) 'restart_trcwt: could not generate rstfnm_hamocc file!'
          call xcstop('(restat_trcwt)')
          stop '(restart_trcwt)'
        endif
      endif
    endif

    call xcbcst(rstfnm_ocntrc)
    call xcbcst(rstfnm_hamocc)

    ! ------------------------------------------------------------------
    ! Write restart files
    ! ------------------------------------------------------------------
    !
    call restart_ocntrcwt(rstfnm_ocntrc)
#ifdef HAMOCC
    call restart_hamoccwt(rstfnm_hamocc)
#endif

  end subroutine restart_trcwt

  ! ============================================================================

  subroutine restart_trcrd(rstfnm_ocn)
    ! ------------------------------------------------------------------
    ! Read tracer state from restart files
    ! ------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: rstfnm_ocn

    ! Local variables
    logical :: error
    character(len=fnmlen) :: rstfnm_ocntrc
    character(len=fnmlen) :: rstfnm_hamocc

    ! ------------------------------------------------------------------
    ! Generate file name
    ! ------------------------------------------------------------------
    if (mnproc == 1) then
      if (expcnf == 'cesm') then
        call restart_getfile(rstfnm_ocn, 'rtrc', rstfnm_ocntrc, error)
        if (error) then
          write(lp,*) 'restart_trcrd: could not generate rstfnm_ocntrc file!'
          call xcstop('(restat_trcrd)')
          stop '(restart_trcrd)'
        endif
        call restart_getfile(rstfnm_ocn, 'rbgc', rstfnm_hamocc, error)
        if (error) then
          write(lp,*) 'restart_trcrd: could not generate rstfnm_hamocc file!'
          call xcstop('(restat_trcrd)')
          stop '(restart_trcrd)'
        endif
      else
        call restart_getfile(rstfnm_ocn, 'resttrc', rstfnm_ocntrc, error)
        if (error) then
          write(lp,*) 'restart_trcrd: could not generate rstfnm_ocntrc file!'
          call xcstop('(restat_trcrd)')
          stop '(restart_trcrd)'
        endif
        call restart_getfile(rstfnm_ocn, 'restbgc', rstfnm_hamocc, error)
        if (error) then
          write(lp,*) 'restart_trcrd: could not generate rstfnm_hamocc file!'
          call xcstop('(restat_trcrd)')
          stop '(restart_trcrd)'
        endif
      endif
    endif

    call xcbcst(rstfnm_ocntrc)
    call xcbcst(rstfnm_hamocc)

    ! ------------------------------------------------------------------
    ! Read restart files
    ! ------------------------------------------------------------------
    call restart_ocntrcrd(rstfnm_ocntrc)
#ifdef HAMOCC
    call hamocc_init(1,rstfnm_hamocc)
#endif

  end subroutine restart_trcrd

  ! ============================================================================

  subroutine restart_ocntrcwt(rstfnm)
    ! ------------------------------------------------------------------
    ! Write ocean tracer state to restart file
    ! ------------------------------------------------------------------

    ! Arguments
    character :: rstfnm*(*)

    ! Local variables
    integer :: nt,nat
    character(len=fnmlen) :: trcnm

    ! ------------------------------------------------------------------
    ! if no ocean tracers are defined, return
    ! ------------------------------------------------------------------
    if (ntrocn == 0) return

    ! ------------------------------------------------------------------
    ! Create file
    ! ------------------------------------------------------------------
    if (mnproc == 1) then
      write (lp,'(2a)') ' saving ocean tracer restart file ', &
           trim(rstfnm)
    end if
    if (rstfmt == 1) then
      call ncfopn(rstfnm,'w','6',1,iotype)
    else if (rstfmt == 2) then
      call ncfopn(rstfnm,'w','h',1,iotype)
    else
      call ncfopn(rstfnm,'w','c',1,iotype)
    end if

    ! ------------------------------------------------------------------
    ! Create attributes and dimensions
    ! ------------------------------------------------------------------

    call ncputi('nday0',date0%day)
    call ncputi('nmonth0',date0%month)
    call ncputi('nyear0',date0%year)
    call ncputr('time0',time0)
    call ncputr('time',time)
    if (rstcmp == 1) then
      call ncdimc('pcomp',ip,0)
    else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
    end if
    call ncdims('kk2',2*kk)
    call ncdims('time',1)

    ! ------------------------------------------------------------------
    ! Write tracer data to file
    ! ------------------------------------------------------------------

    do nt = 1,ntrocn-natr
      write (trcnm,'(a,i3.3)') 'trc',nt
      if (rstcmp == 1) then
        call ncdefvar(trim(trcnm),'pcomp kk2 time', &
             ndouble,8)
      else
        call ncdefvar(trim(trcnm),'x y kk2 time', &
             ndouble,8)
      end if
    end do
    do nt = 1,natr
      nat = ntr-natr+nt
      write (trcnm,'(a,i3.3)') 'atrc',nt
      if (rstcmp == 1) then
        call ncdefvar(trim(trcnm),'pcomp kk2 time', &
             ndouble,8)
      else
        call ncdefvar(trim(trcnm),'x y kk2 time', &
             ndouble,8)
      end if
    end do

    call ncedef

    do nt = 1,ntrocn-natr
      write (trcnm,'(a,i3.3)') 'trc',nt
      if (rstcmp == 1) then
        call nccomp(trim(trcnm),'pcomp kk2 time', &
             trc(1-nbdy,1-nbdy,1,nt),ip,1.,0.,8)
      else
        call ncwrtr(trim(trcnm),'x y kk2 time', &
             trc(1-nbdy,1-nbdy,1,nt),ip,1,1.,0.,8)
      end if
    end do
    do nt = 1,natr
      nat = ntr-natr+nt
      write (trcnm,'(a,i3.3)') 'atrc',nt
      if (rstcmp == 1) then
        call nccomp(trim(trcnm),'pcomp kk2 time', &
             trc(1-nbdy,1-nbdy,1,nat),ip,1.,0.,8)
      else
        call ncwrtr(trim(trcnm),'x y kk2 time', &
             trc(1-nbdy,1-nbdy,1,nat),ip,1,1.,0.,8)
      end if
    end do

    call ncfcls

  end subroutine restart_ocntrcwt

  ! ============================================================================

  subroutine restart_ocntrcrd(rstfnm)
    ! ------------------------------------------------------------------
    ! Read ocean tracer state from restart file
    ! ------------------------------------------------------------------

    ! Arguments
    character :: rstfnm*(*)

    ! Local variables
    type(date_type) :: date_rest
    integer :: nt,nat
    real :: time0r,timer
    logical :: fexist
    character(len=fnmlen) :: trcnm

    ! ------------------------------------------------------------------
    ! If no ocean tracers are defined, return
    ! ------------------------------------------------------------------

    if (ntrocn == 0) return

    ! ------------------------------------------------------------------
    ! Check for file existence
    ! ------------------------------------------------------------------

    inquire(file=rstfnm,exist = fexist)

    call xcbcst(fexist)

    ! ------------------------------------------------------------------
    ! If file exists, read tracer data from file
    ! ------------------------------------------------------------------

    if (.not.fexist) then
      if (mnproc == 1) then
        write (lp,*) &
             'Warning! No tracer restart file found. Calling ocntrc_init...'
      end if
      call initrc
    else
      if (mnproc == 1) then
        write (lp,'(2a)') ' reading ocean tracer restart file ', &
             trim(rstfnm)
      end if
      call ncfopn(rstfnm,'r',' ',1,iotype)
      call ncgeti('nday0',date_rest%day)
      call ncgeti('nmonth0',date_rest%month)
      call ncgeti('nyear0',date_rest%year)
      call ncgetr('time0',time0r)
      call ncgetr('time',timer)
      if (mnproc == 1) then
        if (date_rest /= date0 .or. &
             time0r  /=  time0 .or. timer /= time) then
          write (lp,'(2a)') &
               ' Warning! The time information of the model and', &
               ' restart file is inconsistent'
          write (lp,'(a)') '                 model          file'
          write (lp,'(a,i4.4,2(i2.2),a,i4.4,2(i2.2))') &
               ' date0:      ',date0, '      ',date_rest
          write (lp,'(a,2f14.4)') &
               ' time0: ',time0,time0r
          write (lp,'(a,2f14.4)') &
               ' time:  ',time,timer
        end if
      end if
      do nt = 1,ntrocn-natr
        write (trcnm,'(a,i3.3)') 'trc',nt
        call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nt),ip,1,0.)
      end do
      do nt = 1,natr
        nat = ntr-natr+nt
        write (trcnm,'(a,i3.3)') 'atrc',nt
        call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nat),ip,1,0.)
      end do
      call ncfcls
    end if

  end subroutine restart_ocntrcrd

  ! ============================================================================

  subroutine restart_getfile(rstfnm_in, rstlabel, rstfnm_out, rstfnm_err)
    ! ------------------------------------------------------------------
    ! Generate filename for restart files to read or write tracer fields
    ! ------------------------------------------------------------------

    ! Argument
    character(len=*), intent(in)  :: rstfnm_in     ! Original restart file name
    character(len=*), intent(in)  :: rstlabel      ! Label to insert in new file
    character(len=*), intent(out) :: rstfnm_out    ! New restart file name
    logical,          intent(out) :: rstfnm_err    ! Error flag

    ! Local variables
    integer :: i_suffix, i_time, i_restart
    character(len=:), allocatable  :: str_suffix, str_time, str_restart

    rstfnm_err = .false.

    if (expcnf.eq.'cesm') then
      ! Assume file format: <str_restart.>'r'<.str_timestamp><.str_suffix>
      ! Search for '.' starting from end of "rstfnm_in" filename
      ! File suffix
      i_suffix = index(rstfnm_in, '.', back=.true.)
      str_suffix = trim(rstfnm_in(i_suffix:))
      ! File timestamp
      i_time = index(rstfnm_in(:(i_suffix-1)), '.', back=.true.)
      str_time = rstfnm_in(i_time:(i_suffix-1))
      ! File without original restart label
      i_restart = index(rstfnm_in(:(i_time-1)), '.', back=.true.)
      str_restart = rstfnm_in(:i_restart)

      if (i_suffix == 0 .or. i_time == 0 .or. i_restart == 0) then
        rstfnm_err = .true.
      else
        rstfnm_out = str_restart // trim(rstlabel) // str_time // str_suffix
      end if
    else
      ! Assume file format: <str_restart_>'restphy'<_str_suffix>
      ! Search for '_' starting from end of "rstfnm_in" filename
      ! File suffix
      i_suffix = index(rstfnm_in, '_', back=.true.)
      str_suffix = trim(rstfnm_in(i_suffix:))
      ! File without original restart label
      i_restart = index(rstfnm_in(:(i_suffix-1)), '_', back=.true.)
      str_restart = rstfnm_in(:i_restart)

      if (i_suffix == 0 .or. i_restart == 0) then
        rstfnm_err = .true.
      else
        rstfnm_out = str_restart // trim(rstlabel) // str_suffix
      end if
    end if

  end subroutine restart_getfile

end module mod_tracers_update
