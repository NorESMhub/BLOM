module ocn_pio_share
   ! Initialize PIO

   use pio, only : iosystem_desc_t

   implicit none
   public

   type(iosystem_desc_t), pointer :: pio_subsystem => null()     ! pio info
   integer                        :: io_type                     ! pio info
   integer                        :: io_format                   ! pio info

! contains

!    subroutine init_pio
!       use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat

!       pio_subsystem => shr_pio_getiosys('OCN')
!       io_type       =  shr_pio_getiotype('OCN')
!       io_format     =  shr_pio_getioformat('OCN')
!    end subroutine init_pio

end module ocn_pio_share
