      module data_mct

      ! Data set by mct to be used by micom

      integer, public :: mpicom_mct              ! mpi communicator
      character (len=256), public :: runid_mct   ! case name
      character (len=256), public :: runtype_mct ! run type

      end module data_mct
