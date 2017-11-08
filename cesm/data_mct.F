      module data_mct

      ! Data set by mct to be used by micom

      implicit none

      private

      character (len=256) :: runid_mct  ! case name
      character (len=256) :: runtyp_mct ! run type
      integer :: ocn_cpl_dt_mct         ! Coupling frequency

      public :: runid_mct, runtyp_mct, ocn_cpl_dt_mct

      end module data_mct
