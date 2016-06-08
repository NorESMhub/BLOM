module types

   ! -------------------------------------------------------------------
   ! This module defines numerical data types
   ! -------------------------------------------------------------------

   implicit none

   integer, parameter, public :: &
      i1 = selected_int_kind(2), &
      i2 = selected_int_kind(4), &
      i4 = selected_int_kind(9), &
      i8 = selected_int_kind(18), &
      r4 = selected_real_kind(6,37), &
      r8 = selected_real_kind(13,307)

end module types
