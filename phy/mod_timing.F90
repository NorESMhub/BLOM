! ------------------------------------------------------------------------------
! Copyright (C) 2025 Mats Bentsen
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

module mod_timing
! ------------------------------------------------------------------------------
! This module contains routines and data structures for timers.
! ------------------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env, only: int32
   use mod_types, only: r8
   use mod_crc32, only: crc32
   use mod_wtime, only: wtime
   use mod_xc   , only: xchalt, xcstop, xcmax, xcsync, mnproc, lp

   implicit none

   private

   integer, parameter :: &
      maxlen_name = 52

   type :: group_struct
      logical :: used
      integer(int32) :: checksum
      integer :: mode, num_timers, first_tmridx, last_tmridx
      real(r8) :: last_time
   end type group_struct

   type :: timer_struct
      logical :: used
      integer(int32) :: checksum
      character(len = maxlen_name) :: name
      integer :: acc_num, next_tmridx
      real(r8) :: acc_dtime
   end type timer_struct

   type(group_struct), allocatable, dimension(:) :: group
   type(timer_struct), allocatable, dimension(:) :: timer
   integer :: maxnum_timers, maxnum_groups, hash_mask_group, hash_mask_timer
   logical :: initialized = .false.

   public :: timer_init, timer_start, timer_stop, timer_reset, &
             timer_group_time, timer_statistics

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   function group_index(group_name) result(grpidx)
   ! ---------------------------------------------------------------------------
   ! Return group array index associated with argument group_name.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name

      integer :: grpidx

      integer(int32) :: checksum
      integer :: grpidx_start

      checksum = crc32(trim(group_name))
      grpidx_start = iand(checksum, hash_mask_group) + 1

      grpidx = grpidx_start
      do
         if (group(grpidx)%used) then
            if (group(grpidx)%checksum == checksum) return
         else
            group(grpidx)%checksum = checksum
            return
         endif
         grpidx = iand(grpidx, hash_mask_group) + 1
         if (grpidx == grpidx_start) then
            if (mnproc == 1) &
               write(lp,*) &
                  'Too many groups! Increase maxnum_groups_log2 in call '// &
                  'to timers_init'
            call xcstop('(group_index)')
                   stop '(group_index)'
         endif
      enddo

   end function group_index

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine timer_init(maxnum_groups_log2, maxnum_timers_log2)
   ! ---------------------------------------------------------------------------
   ! Initialize timers. Arguments:
   !    maxnum_groups_log2: Number of groups will be 2**maxnum_groups_log2.
   !    maxnum_timers_log2: Number of timers will be 2**maxnum_timers_log2.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: maxnum_groups_log2, maxnum_timers_log2

      integer :: errstat

      if (initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer already initialized!'
         call xcstop('(timer_init)')
                stop '(timer_init)'
      endif

      maxnum_groups = 2**maxnum_groups_log2
      maxnum_timers = 2**maxnum_timers_log2
      hash_mask_group = maskr(maxnum_groups_log2)
      hash_mask_timer = maskr(maxnum_timers_log2)

      allocate(group(maxnum_groups), &
               timer(maxnum_timers), &
               stat = errstat)
      if (errstat /= 0) then
         write(lp,*) 'Failed to allocate timing arrays!'
         call xchalt('(timer_init)')
                stop '(timer_init)'
      endif

      group(:)%used = .false.
      timer(:)%used = .false.

      initialized = .true.

   end subroutine timer_init

   subroutine timer_start(group_name)
   ! ---------------------------------------------------------------------------
   ! Record the current time for the timer group associated with argument
   ! group_name.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name

      integer :: grpidx

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_start)')
                stop '(timer_start)'
      endif

      grpidx = group_index(group_name)

      if (.not.group(grpidx)%used) then
         group(grpidx)%used = .true.
         group(grpidx)%mode = 0
         group(grpidx)%num_timers = 0
      endif

      call xcsync(.false.)
      group(grpidx)%last_time = wtime()

   end subroutine timer_start

   subroutine timer_stop(group_name, timer_name)
   ! ---------------------------------------------------------------------------
   ! Record the elapsed time since previous call to either timer_start or
   ! timer_stop for the timer group associated with the argument group_name. If
   ! the optional timer_name argument is not provided, the timer group is in
   ! single-timer mode. In multi-timer mode, each unique timer_name is
   ! associated with its own timer within the group. For additional calls to
   ! timer_stop with previously used group_name and timer_name arguments,
   ! elapsed time is accumulated in the existing timer. For multiple timers
   ! where a combined timer statistic output is preferred, it is possible to use
   ! timer_name of the form "<prefix>_sg<suffix>", which defines a subgroup
   ! within the timer group. The <prefix> must be the same for all timers in the
   ! subgroup, while the <suffix> must be unique.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name
      character(len = *), optional, intent(in) :: timer_name

      real(r8) :: curr_time
      integer(int32) :: checksum
      integer :: grpidx, tmridx_start, tmridx

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_stop)')
                stop '(timer_stop)'
      endif

      grpidx = group_index(group_name)

      if (.not.group(grpidx)%used) then
         if (mnproc == 1) &
            write(lp,*) 'Timer group '''//trim(group_name)// &
                        ''' not started, call timer_start first!'
         call xcstop('(timer_stop)')
                stop '(timer_stop)'
      endif

      if (present(timer_name)) then
         if (group(grpidx)%mode == 1) then
            if (mnproc == 1) &
               write(lp,*) 'Timer group '''//trim(group_name)// &
                           '''is already in single-timer mode!'
            call xcstop('(timer_stop)')
                   stop '(timer_stop)'
         endif
         group(grpidx)%mode = 2
         checksum = crc32(trim(group_name)//trim(timer_name))
      else
         if (group(grpidx)%mode == 2) then
            if (mnproc == 1) &
               write(lp,*) 'Timer group '''//trim(group_name)// &
                           '''is already in multi-timer mode!'
            call xcstop('(timer_stop)')
                   stop '(timer_stop)'
         endif
         group(grpidx)%mode = 1
         checksum = crc32(trim(group_name))
      endif

      tmridx_start = iand(checksum, hash_mask_timer) + 1
      tmridx = tmridx_start
      do
         if (timer(tmridx)%used) then
            if (timer(tmridx)%checksum == checksum) then
               timer(tmridx)%acc_num = timer(tmridx)%acc_num + 1
               curr_time = wtime()
               timer(tmridx)%acc_dtime = timer(tmridx)%acc_dtime &
                                       + curr_time - group(grpidx)%last_time
               group(grpidx)%last_time = curr_time
               return
            endif
         else
            timer(tmridx)%used = .true.
            timer(tmridx)%checksum = checksum
            timer(tmridx)%name = timer_name
            timer(tmridx)%acc_num = 1
            timer(tmridx)%next_tmridx = 0
            group(grpidx)%num_timers = group(grpidx)%num_timers + 1
            if (group(grpidx)%num_timers == 1) then
               group(grpidx)%first_tmridx = tmridx
            else
               timer(group(grpidx)%last_tmridx)%next_tmridx = tmridx
            endif
            group(grpidx)%last_tmridx = tmridx
            curr_time = wtime()
            timer(tmridx)%acc_dtime = curr_time - group(grpidx)%last_time
            group(grpidx)%last_time = curr_time
            return
         endif
         tmridx = iand(tmridx, hash_mask_timer) + 1
         if (tmridx == tmridx_start) then
            if (mnproc == 1) &
               write(lp,*) &
                  'Too many timers! Increase maxnum_timers_log2 in call '// &
                  'to timers_init'
            call xcstop('(timer_stop)')
                   stop '(timer_stop)'
         endif
      enddo

   end subroutine timer_stop

   subroutine timer_reset(group_name)
   ! ---------------------------------------------------------------------------
   ! Reset the timer group associated with argument group_name.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name

      integer :: grpidx, tmridx, grpidx_start, i, j

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_reset)')
                stop '(timer_reset)'
      endif

      grpidx = group_index(group_name)

      if (.not.group(grpidx)%used) return

      if (group(grpidx)%num_timers > 0) then
         tmridx = group(grpidx)%first_tmridx
         do while (tmridx /= 0)
            timer(tmridx)%used = .false.
            tmridx = timer(tmridx)%next_tmridx
         enddo
      endif

      group(grpidx)%used = .false.
      i = iand(grpidx, hash_mask_group) + 1
      do while (i /= grpidx)
         grpidx_start = iand(group(i)%checksum, hash_mask_group) + 1
         j = grpidx_start
         do
            if (group(j)%used) then
               if (group(j)%checksum == group(i)%checksum) exit
            else
               group(j) = group(i)
               group(i)%used = .false.
               exit
            endif
            j = iand(j, hash_mask_group) + 1
         enddo
         i = iand(i, hash_mask_group) + 1
      enddo

   end subroutine timer_reset

   subroutine timer_statistics(group_name)
   ! ---------------------------------------------------------------------------
   ! Output statistics of the timer group, associated with argument group_name,
   ! to stdout.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name

      real(r8), allocatable, dimension(:) :: cumsum_time
      logical, allocatable, dimension(:) :: contributed
      real(r8) :: group_time, sg_time
      integer :: grpidx, tmridx, acc_num, errstat, i, sg_index, tmrjdx, j
      logical :: equal_acc_nums

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_statistics)')
                stop '(timer_statistics)'
      endif

      grpidx = group_index(group_name)

      if     (.not.group(grpidx)%used) then
         if (mnproc == 1) &
            write(lp,*) 'No timer group associated with name '''// &
                        trim(group_name)//''''

         return
      elseif (group(grpidx)%num_timers == 0) then
         if (mnproc == 1) &
            write(lp,*) 'No timers of group '''//trim(group_name)//''''
         return
      endif

      if (mnproc == 1) then
         write(lp,*) ''
         write(lp,'(a)') 'Timer statistics of group '''//trim(group_name)//''':'
      endif

      if (group(grpidx)%mode == 1) then

         tmridx = group(grpidx)%first_tmridx
         group_time = timer(tmridx)%acc_dtime
         acc_num = timer(tmridx)%acc_num
         call xcmax(group_time)

         if (mnproc == 1) then
            if (acc_num > 1) &
               write(lp,'(f12.4,a,i6,a)') &
                  group_time/acc_num, '  sec average time for ', &
                  acc_num, ' accumulations'
            write(lp,'(f12.4,a)') group_time, '  sec total time'
         endif

      else

          allocate(cumsum_time(group(grpidx)%num_timers+1), &
                   contributed(group(grpidx)%num_timers), &
                   stat = errstat)
          if (errstat /= 0) then
             write(lp,*) 'Failed to allocate timing arrays!'
             call xchalt('(timer_statistics)')
                    stop '(timer_statistics)'
          endif

         tmridx = group(grpidx)%first_tmridx
         acc_num = timer(tmridx)%acc_num
         i = 1
         cumsum_time(1) = 0._r8
         equal_acc_nums = .true.
         do while (tmridx /= 0)
            cumsum_time(i+1) = cumsum_time(i) + timer(tmridx)%acc_dtime
            if (timer(tmridx)%acc_num /= acc_num) equal_acc_nums = .false.
            tmridx = timer(tmridx)%next_tmridx
            i = i + 1
         enddo
         call xcmax(cumsum_time)

         if (mnproc == 1) then
            group_time = cumsum_time(group(grpidx)%num_timers+1)
            if (acc_num > 1 .and. equal_acc_nums) &
               write(lp,'(f12.4,a,i6,a)') &
                  group_time/acc_num, ' sec average time for ', &
                  acc_num, ' accumulations'
            write(lp,'(f12.4,a)') group_time, &
                                  ' sec total time with contributions:'
            contributed(:) = .false.
            tmridx = group(grpidx)%first_tmridx
            i = 1
            do while (tmridx /= 0)
               if (.not.contributed(i)) then
                  sg_index = index(timer(tmridx)%name, '_sg')
                  if (sg_index <= 1) then
                     write(lp,'(f12.4,a,f8.4,a)') &
                         cumsum_time(i+1) - cumsum_time(i),' sec ', &
                        (cumsum_time(i+1) - cumsum_time(i)) &
                        /group_time*100._r8, &
                        '%  '//trim(timer(tmridx)%name)
                     contributed(i) = .true.
                  else
                     sg_time = cumsum_time(i+1) - cumsum_time(i)
                     contributed(i) = .true.
                     tmrjdx = timer(tmridx)%next_tmridx
                     j = i + 1
                     do while (tmrjdx /= 0)
                        if (timer(tmrjdx)%name(1:sg_index+2) == &
                            timer(tmridx)%name(1:sg_index+2)) then
                           sg_time = sg_time + cumsum_time(j+1) - cumsum_time(j)
                           contributed(j) = .true.
                        endif
                        tmrjdx = timer(tmrjdx)%next_tmridx
                        j = j + 1
                     enddo
                     write(lp,'(f12.4,a,f8.4,a)') &
                        sg_time, ' sec ', &
                        sg_time/group_time*100._r8, &
                        '%  '//timer(tmridx)%name(1:sg_index-1)
                  endif
               endif
               tmridx = timer(tmridx)%next_tmridx
               i = i + 1
            enddo
         endif

      endif

   end subroutine timer_statistics

   subroutine timer_group_time(group_name, group_time)
   ! ---------------------------------------------------------------------------
   ! Provide the total accumulated time for the timer group, associated with
   ! argument group_name, in group_time.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: group_name
      real(r8), intent(out) :: group_time

      integer :: grpidx, tmridx

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_group_time)')
                stop '(timer_group_time)'
      endif

      grpidx = group_index(group_name)
      
      group_time = 0._r8

      if     (.not.group(grpidx)%used) then
         return
      elseif (group(grpidx)%num_timers == 0) then
         return
      endif

      tmridx = group(grpidx)%first_tmridx
      do while (tmridx /= 0)
         group_time = group_time + timer(tmridx)%acc_dtime
         tmridx = timer(tmridx)%next_tmridx
      enddo
      call xcmax(group_time)

   end subroutine timer_group_time

end module mod_timing
