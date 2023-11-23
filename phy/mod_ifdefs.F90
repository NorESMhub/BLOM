module mod_ifdefs

  implicit none
  public
 
#ifdef TRC
  logical, parameter :: use_TRC = .true.
#else
  logical, parameter :: use_TRC = .false.
#endif
#ifdef TKE
  logical, parameter :: use_TKE = .true.
#else
  logical, parameter :: use_TKE = .false.
#endif
#ifdef IDLAGE
  logical, parameter :: use_IDLAGE = .true.
#else
  logical, parameter :: use_IDLAGE = .false.
#endif

end module mod_ifdefs
