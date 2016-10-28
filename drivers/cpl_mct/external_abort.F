      subroutine external_abort(message)

      ! Uses modules

      use shr_sys_mod

      implicit none

      ! Input/output arguments

      character*(*) message

      ! ----------------------------------------------------------------
      ! Call CESM shared abort routine
      ! ----------------------------------------------------------------
      call shr_sys_abort(message)

      end subroutine external_abort
