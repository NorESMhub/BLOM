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
   use mod_xc   , only: xchalt, xcstop, xcmax, mnproc, lp

   implicit none

   private

   integer, parameter :: &
      maxlen_name = 52

   type :: timer_struct
      logical :: used
      integer(int32) :: checksum
      character(len = maxlen_name) :: name
      integer :: acc_num
      real(r8) :: acc_time
   end type timer_struct

   type(timer_struct), allocatable, dimension(:,:) :: timer
   integer, allocatable, dimension(:,:) :: timer_idx_list
   integer, allocatable, dimension(:) :: timer_num, timer_mode
   real(r8), allocatable, dimension(:) :: last_time
   integer :: maxnum_timers, maxnum_groups, hash_mask
   logical :: initialized = .false.

   public :: timer_init, timer_start, timer_stop, timer_reset, &
             timer_group_time, timer_statistics

contains

   subroutine timer_init(maxnum_timers_log2, maxnum_groups_arg)
   ! ---------------------------------------------------------------------------
   ! Initialize timers. Arguments:
   ! - maxnum_timers_log2: Number of timers will be 2**maxnum_timers_log2.
   ! - maxnum_groups_arg : Maximum number of timer groups.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: maxnum_timers_log2, maxnum_groups_arg

      integer :: errstat
      integer :: checksum, hash

      maxnum_groups = maxnum_groups_arg
      maxnum_timers = 2**maxnum_timers_log2
      hash_mask = maskr(maxnum_timers_log2)

      allocate(timer(maxnum_timers,maxnum_groups), &
               timer_idx_list(maxnum_timers,maxnum_groups), &
               timer_num(maxnum_groups), &
               timer_mode(maxnum_groups), &
               last_time(maxnum_groups), &
               stat = errstat)
      if (errstat /= 0) then
         write(lp,*) 'Failed to allocate timing arrays!'
         call xchalt('(timer_init)')
                stop '(timer_init)'
      endif

      timer(:,:)%used = .false.
      timer_num(:) = 0
      timer_mode(:) = 0

      initialized = .true.

   end subroutine timer_init

   subroutine timer_start(group)
   ! ---------------------------------------------------------------------------
   ! Record the current time for the timer group, provided as argument.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: group

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_start)')
                stop '(timer_start)'
      endif
      if (group < 1 .or. group > maxnum_groups) then
         if (mnproc == 1) write(lp,*) 'Group number out of bounds!'
         call xcstop('(timer_start)')
                stop '(timer_start)'
      endif

      last_time(group) = wtime()

   end subroutine timer_start

   subroutine timer_stop(group, name)
   ! ---------------------------------------------------------------------------
   ! Record the elapsed time since previous call to either timer_start or
   ! timer_stop for the timer group, provided as argument. If the optional name
   ! argument is not provided, the timer group is in single-timer mode. In
   ! multi-timer mode, each unique name is associated with its own timer within
   ! the group. For additional calls to timer_stop with previously used group
   ! and name arguments, elapsed time is accumulated in the existing timer. For
   ! multiple timers where a combined timer statistic output is preferred, it is
   ! possible to use names of the form "<prefix>_sg<suffix>", which defines a
   ! subgroup within the timer group. The <prefix> must be the same for all
   ! timers in the subgroup, while the <suffix> must be unique.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: group
      character(len = *), optional, intent(in) :: name

      real(r8) :: curr_time
      integer :: checksum, idx_start, idx

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_stop)')
                stop '(timer_stop)'
      endif
      if (group < 1 .or. group > maxnum_groups) then
         if (mnproc == 1) write(lp,*) 'Group number out of bounds!'
         call xcstop('(timer_stop)')
                stop '(timer_stop)'
      endif

      if (present(name)) then

         if (timer_mode(group) == 1) then
            if (mnproc == 1) &
               write(lp,*) 'Timer is already in single-timer mode!'
            call xcstop('(timer_stop)')
                   stop '(timer_stop)'
            
         endif

         checksum = crc32(name)
         idx_start = iand(checksum, hash_mask) + 1

         idx = idx_start
         do
            if (timer(idx,group)%used) then
               if (timer(idx,group)%checksum == checksum) then
                  timer(idx,group)%acc_num = timer(idx,group)%acc_num + 1
                  curr_time = wtime()
                  timer(idx,group)%acc_time = timer(idx,group)%acc_time &
                                            + curr_time - last_time(group)
                  last_time(group) = curr_time
                  return
               endif
            else
               timer(idx,group)%used = .true.
               timer(idx,group)%checksum = checksum
               timer(idx,group)%name = name
               timer(idx,group)%acc_num = 1
               timer_num(group) = timer_num(group) + 1
               timer_mode(group) = 2
               timer_idx_list(timer_num(group),group) = idx
               curr_time = wtime()
               timer(idx,group)%acc_time = curr_time - last_time(group)
               last_time(group) = curr_time
               return
            endif
            idx = iand(idx, hash_mask) + 1
            if (idx == idx_start) then
               if (mnproc == 1) &
                  write(lp,*) &
                     'Too many timers! Increase maxnum_timers_log2 in call '// &
                     'to timers_init'
               call xcstop('(timer_stop)')
                      stop '(timer_stop)'
            endif
         enddo

      else

         if (timer_mode(group) == 2) then
            if (mnproc == 1) &
               write(lp,*) 'Timer is already in multi-timer mode!'
            call xcstop('(timer_stop)')
                   stop '(timer_stop)'
            
         endif

         idx = 1
         if (timer(idx,group)%used) then
            timer(idx,group)%acc_num = timer(idx,group)%acc_num + 1
            curr_time = wtime()
            timer(idx,group)%acc_time = timer(idx,group)%acc_time &
                                      + curr_time - last_time(group)
            last_time(group) = curr_time
         else
            timer(idx,group)%used = .true.
            timer(idx,group)%acc_num = 1
            timer_num(group) = 1
            timer_mode(group) = 1
            timer_idx_list(timer_num(group),group) = idx
            curr_time = wtime()
            timer(idx,group)%acc_time = curr_time - last_time(group)
            last_time(group) = curr_time
            return
         endif

      endif

   end subroutine timer_stop

   subroutine timer_reset(group)
   ! ---------------------------------------------------------------------------
   ! Reset the timer group, provided as argument.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: group

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_reset)')
                stop '(timer_reset)'
      endif
      if (group < 1 .or. group > maxnum_groups) then
         if (mnproc == 1) write(lp,*) 'Group number out of bounds!'
         call xcstop('(timer_reset)')
                stop '(timer_reset)'
      endif

      timer(:,group)%used = .false.
      timer_num(group) = 0
      timer_mode(group) = 0

   end subroutine timer_reset

   subroutine timer_statistics(group)
   ! ---------------------------------------------------------------------------
   ! Output statistics of the timer group, provided as argument, to stdout.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: group

      real(r8), dimension(timer_num(group)+1) :: cumsum_time
      logical, dimension(timer_num(group)+1) :: contributed
      real(r8) :: group_time, sg_time
      integer :: i, idx, acc_num, sg_index, j, jdx
      logical :: equal_acc_nums

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_statistics)')
                stop '(timer_statistics)'
      endif
      if (group < 1 .or. group > maxnum_groups) then
         if (mnproc == 1) write(lp,*) 'Group number out of bounds!'
         call xcstop('(timer_statistics)')
                stop '(timer_statistics)'
      endif

      if (timer_num(group) == 0) then
         if (mnproc == 1) write(lp,*) 'No timer statistics of group:', group
         return
      endif

      if (mnproc == 1) then
         write(lp,*) ''
         write(lp,'(a,i3,a)') 'Timer statistics of group ', group, ':'
      endif

      if (timer_mode(group) == 1) then

         idx = 1
         group_time = timer(idx,group)%acc_time
         acc_num = timer(idx,group)%acc_num
         call xcmax(group_time)

         if (mnproc == 1) then
            if (acc_num > 1) &
               write(lp,'(f12.4,a,i6,a)') &
                  group_time/acc_num, '  sec average time for ', &
                  acc_num, ' accumulations'
            write(lp,'(f12.4,a)') group_time, '  sec total time'
            write(lp,*) ''
         endif

      else

         cumsum_time(1) = 0._r8
         acc_num = 0
         equal_acc_nums = .true.
         do i = 1, timer_num(group)
            idx = timer_idx_list(i,group)
            cumsum_time(i+1) = cumsum_time(i) + timer(idx,group)%acc_time
            if (i > 1 .and. timer(idx,group)%acc_num /= acc_num) &
               equal_acc_nums = .false.
            acc_num = timer(idx,group)%acc_num
         enddo
         call xcmax(cumsum_time)

         if (mnproc == 1) then
            group_time = cumsum_time(timer_num(group)+1)
            if (acc_num > 1 .and. equal_acc_nums) &
               write(lp,'(f12.4,a,i6,a)') &
                  group_time/acc_num, ' sec average time for ', &
                  acc_num, ' accumulations'
            write(lp,'(f12.4,a)') group_time, &
                                  ' sec total time with contributions:'
            contributed(:) = .false.
            do i = 1, timer_num(group)
               if (.not.contributed(i)) then
                  idx = timer_idx_list(i,group)
                  sg_index = index(timer(idx,group)%name, '_sg')
                  if (sg_index <= 1) then
                     write(lp,'(f12.4,a,f8.4,a)') &
                         cumsum_time(i+1) - cumsum_time(i),' sec ', &
                        (cumsum_time(i+1) - cumsum_time(i)) &
                        /group_time*100._r8, &
                        '%  '//trim(timer(idx,group)%name)
                     contributed(i) = .true.
                  else
                     sg_time = cumsum_time(i+1) - cumsum_time(i)
                     contributed(i) = .true.
                     do j = i + 1, timer_num(group)
                        jdx = timer_idx_list(j,group)
                        if (timer(jdx,group)%name(1:sg_index+2) == &
                            timer(idx,group)%name(1:sg_index+2)) then
                           sg_time = sg_time + cumsum_time(j+1) - cumsum_time(j)
                           contributed(j) = .true.
                        endif
                     enddo
                     write(lp,'(f12.4,a,f8.4,a)') &
                        sg_time, ' sec ', &
                        sg_time/group_time*100._r8, &
                        '%  '//timer(idx,group)%name(1:sg_index-1)
                  endif
               endif
            enddo
            write(lp,*) ''
         endif

      endif

   end subroutine timer_statistics

   subroutine timer_group_time(group, group_time)
   ! ---------------------------------------------------------------------------
   ! Provide the total accumulated time for the timer group, provided as
   ! argument, in group_time.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: group
      real(r8), intent(out) :: group_time

      real(r8), dimension(timer_num(group)+1) :: cumsum_time
      logical, dimension(timer_num(group)+1) :: contributed
      real(r8) :: sg_time
      integer :: i, idx, acc_num, sg_index, j, jdx
      logical :: equal_acc_nums

      if (.not.initialized) then
         if (mnproc == 1) &
            write(lp,*) 'Timer not initialized, call timer_init first!'
         call xcstop('(timer_group_time)')
                stop '(timer_group_time)'
      endif
      if (group < 1 .or. group > maxnum_groups) then
         if (mnproc == 1) write(lp,*) 'Group number out of bounds!'
         call xcstop('(timer_group_time)')
                stop '(timer_group_time)'
      endif

      group_time = 0._r8
      do i = 1, timer_num(group)
         idx = timer_idx_list(i,group)
         group_time = group_time + timer(idx,group)%acc_time
      enddo
      call xcmax(group_time)

   end subroutine timer_group_time

end module mod_timing
