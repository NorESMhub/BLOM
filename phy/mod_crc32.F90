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

module mod_crc32
! ------------------------------------------------------------------------------
! This module contains functions that return a CRC32 (Cyclic Redundancy Check
! 32-bit) checksum/hash for character, scalar or array arguments, where arrays
! can be up to rank 3. For scalar and array arguments, data types int8, int16,
! int32, int64, real32, real64 are supported as defined by the intrinsic
! iso_fortran_env module. The CRC32 implementation is based on
! https://rosettacode.org/wiki/CRC-32#Fortran.
! ------------------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, &
                                            real32, real64

   implicit none

   private

   integer(int32) :: crc_table(0:255)
   logical :: table_initialized = .false.

   ! ---------------------------------------------------------------------------
   ! Polymorphic interface of the crc32 function.
   ! ---------------------------------------------------------------------------

   interface crc32
      module procedure &
         crc32_char, &
         crc32_int8_rank0, crc32_int8_rank1, &
         crc32_int8_rank2, crc32_int8_rank3, &
         crc32_int16_rank0, crc32_int16_rank1, &
         crc32_int16_rank2, crc32_int16_rank3, &
         crc32_int32_rank0, crc32_int32_rank1, &
         crc32_int32_rank2, crc32_int32_rank3, &
         crc32_int64_rank0, crc32_int64_rank1, &
         crc32_int64_rank2, crc32_int64_rank3, &
         crc32_real32_rank0, crc32_real32_rank1, &
         crc32_real32_rank2, crc32_real32_rank3, &
         crc32_real64_rank0, crc32_real64_rank1, &
         crc32_real64_rank2, crc32_real64_rank3
   end interface crc32

   public :: crc32

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   subroutine init_table

      integer :: i, j
      integer(int32) :: k

      do i = 0, 255
         k = i
         do j = 1, 8
            if (btest(k, 0)) then
               k = ieor(shiftr(k, 1), - 306674912)
            else
               k = shiftr(k, 1)
            endif
         enddo
         crc_table(i) = k
      enddo

      table_initialized = .true.

   end subroutine init_table

   ! ---------------------------------------------------------------------------
   ! Type specific crc32 functions.
   ! ---------------------------------------------------------------------------

   function crc32_char(a, crc_init) result(crc)

      character(len = *), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer :: n, i

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = len(a)
      do i = 1, n
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, iachar(a(i:i), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_char

   function crc32_int8_rank0(a, crc_init) result(crc)

      integer(int8), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      crc = ieor(shiftr(crc, 8), &
                 crc_table(iand(ieor(crc, int(a, int32)), 255)))
      crc = not(crc)

   end function crc32_int8_rank0

   function crc32_int8_rank1(a, crc_init) result(crc)

      integer(int8), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer :: n, i

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(a(i), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_int8_rank1

   function crc32_int8_rank2(a, crc_init) result(crc)

      integer(int8), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int8_rank1(a_ptr, crc_init)
      else
         crc = crc32_int8_rank1(a_ptr)
      endif

   end function crc32_int8_rank2

   function crc32_int8_rank3(a, crc_init) result(crc)

      integer(int8), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int8_rank1(a_ptr, crc_init)
      else
         crc = crc32_int8_rank1(a_ptr)
      endif

   end function crc32_int8_rank3

   function crc32_int16_rank0(a, crc_init) result(crc)

      integer(int16), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(2) :: bytes
      integer :: j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      bytes = transfer(a, bytes)
      do j = 1, 2
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_int16_rank0

   function crc32_int16_rank1(a, crc_init) result(crc)

      integer(int16), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(2) :: bytes
      integer :: n, i, j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         bytes = transfer(a(i), bytes)
         do j = 1, 2
            crc = ieor(shiftr(crc, 8), &
                       crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
         enddo
      enddo
      crc = not(crc)

   end function crc32_int16_rank1

   function crc32_int16_rank2(a, crc_init) result(crc)

      integer(int16), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int16), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int16_rank1(a_ptr, crc_init)
      else
         crc = crc32_int16_rank1(a_ptr)
      endif

   end function crc32_int16_rank2

   function crc32_int16_rank3(a, crc_init) result(crc)

      integer(int16), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int16), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int16_rank1(a_ptr, crc_init)
      else
         crc = crc32_int16_rank1(a_ptr)
      endif

   end function crc32_int16_rank3

   function crc32_int32_rank0(a, crc_init) result(crc)

      integer(int32), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(4) :: bytes
      integer :: j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      bytes = transfer(a, bytes)
      do j = 1, 4
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_int32_rank0

   function crc32_int32_rank1(a, crc_init) result(crc)

      integer(int32), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(4) :: bytes
      integer :: n, i, j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         bytes = transfer(a(i), bytes)
         do j = 1, 4
            crc = ieor(shiftr(crc, 8), &
                       crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
         enddo
      enddo
      crc = not(crc)

   end function crc32_int32_rank1

   function crc32_int32_rank2(a, crc_init) result(crc)

      integer(int32), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int32), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int32_rank1(a_ptr, crc_init)
      else
         crc = crc32_int32_rank1(a_ptr)
      endif

   end function crc32_int32_rank2

   function crc32_int32_rank3(a, crc_init) result(crc)

      integer(int32), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int32), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int32_rank1(a_ptr, crc_init)
      else
         crc = crc32_int32_rank1(a_ptr)
      endif

   end function crc32_int32_rank3

   function crc32_int64_rank0(a, crc_init) result(crc)

      integer(int64), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(8) :: bytes
      integer :: j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      bytes = transfer(a, bytes)
      do j = 1, 8
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_int64_rank0

   function crc32_int64_rank1(a, crc_init) result(crc)

      integer(int64), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(8) :: bytes
      integer :: n, i, j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         bytes = transfer(a(i), bytes)
         do j = 1, 8
            crc = ieor(shiftr(crc, 8), &
                       crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
         enddo
      enddo
      crc = not(crc)

   end function crc32_int64_rank1

   function crc32_int64_rank2(a, crc_init) result(crc)

      integer(int64), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int64), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int64_rank1(a_ptr, crc_init)
      else
         crc = crc32_int64_rank1(a_ptr)
      endif

   end function crc32_int64_rank2

   function crc32_int64_rank3(a, crc_init) result(crc)

      integer(int64), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int64), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_int64_rank1(a_ptr, crc_init)
      else
         crc = crc32_int64_rank1(a_ptr)
      endif

   end function crc32_int64_rank3

   function crc32_real32_rank0(a, crc_init) result(crc)

      real(real32), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(4) :: bytes
      integer :: j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      bytes = transfer(a, bytes)
      do j = 1, 4
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_real32_rank0

   function crc32_real32_rank1(a, crc_init) result(crc)

      real(real32), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(4) :: bytes
      integer :: n, i, j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         bytes = transfer(a(i), bytes)
         do j = 1, 4
            crc = ieor(shiftr(crc, 8), &
                       crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
         enddo
      enddo
      crc = not(crc)

   end function crc32_real32_rank1

   function crc32_real32_rank2(a, crc_init) result(crc)

      real(real32), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      real(real32), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_real32_rank1(a_ptr, crc_init)
      else
         crc = crc32_real32_rank1(a_ptr)
      endif

   end function crc32_real32_rank2

   function crc32_real32_rank3(a, crc_init) result(crc)

      real(real32), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      real(real32), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_real32_rank1(a_ptr, crc_init)
      else
         crc = crc32_real32_rank1(a_ptr)
      endif

   end function crc32_real32_rank3

   function crc32_real64_rank0(a, crc_init) result(crc)

      real(real64), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(8) :: bytes
      integer :: j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      bytes = transfer(a, bytes)
      do j = 1, 8
         crc = ieor(shiftr(crc, 8), &
                    crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
      enddo
      crc = not(crc)

   end function crc32_real64_rank0

   function crc32_real64_rank1(a, crc_init) result(crc)

      real(real64), dimension(:), intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      integer(int8), dimension(8) :: bytes
      integer :: n, i, j

      if (.not. table_initialized) call init_table

      if (present(crc_init)) then
         crc = crc_init
      else
         crc = 0
      endif
      
      crc = not(crc)
      n = size(a)
      do i = 1, n
         bytes = transfer(a(i), bytes)
         do j = 1, 8
            crc = ieor(shiftr(crc, 8), &
                       crc_table(iand(ieor(crc, int(bytes(j), int32)), 255)))
         enddo
      enddo
      crc = not(crc)

   end function crc32_real64_rank1

   function crc32_real64_rank2(a, crc_init) result(crc)

      real(real64), dimension(:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      real(real64), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_real64_rank1(a_ptr, crc_init)
      else
         crc = crc32_real64_rank1(a_ptr)
      endif

   end function crc32_real64_rank2

   function crc32_real64_rank3(a, crc_init) result(crc)

      real(real64), dimension(:,:,:), contiguous, target, intent(in) :: a
      integer(int32), optional, intent(in) :: crc_init

      integer(int32) :: crc

      real(real64), dimension(:), pointer :: a_ptr

      a_ptr(1:size(a)) => a

      if (present(crc_init)) then
         crc = crc32_real64_rank1(a_ptr, crc_init)
      else
         crc = crc32_real64_rank1(a_ptr)
      endif

   end function crc32_real64_rank3

end module mod_crc32
