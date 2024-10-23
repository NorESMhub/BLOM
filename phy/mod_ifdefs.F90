module mod_ifdefs

  implicit none
  public

#ifdef TRC
  logical :: use_TRC = .true.
#else
  logical :: use_TRC = .false.
#endif
#ifdef ATRC
  logical :: use_ATRC = .true.
#else
  logical :: use_ATRC = .false.
#endif
#ifdef TKE
  logical :: use_TKE = .true.
#else
  logical :: use_TKE = .false.
#endif
#ifdef TKEADV
  logical :: use_TKEADV = .true.
#else
  logical :: use_TKEADV = .false.
#endif
#ifdef TKEIDF
  logical :: use_TKEIDF = .true.
#else
  logical :: use_TKEIDF = .false.
#endif
#ifdef GLS
  logical :: use_GLS = .true.
#else
  logical :: use_GLS = .false.
#endif
#ifdef IDLAGE
  logical :: use_IDLAGE = .true.
#else
  logical :: use_IDLAGE = .false.
#endif
#ifdef MKS
  logical :: use_MKS = .true.
#else
  logical :: use_MKS = .false.
#endif

  ! Namelist input
  logical :: use_diag = .false.

end module mod_ifdefs
