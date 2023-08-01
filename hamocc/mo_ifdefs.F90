module mo_ifdefs

  implicit none
  public
 
#ifdef BROMO
  logical, parameter :: use_BROMO = .true.
#else
  logical, parameter :: use_BROMO = .false.
#endif
#ifdef AGG
  logical, parameter :: use_AGG = .true.
#else
  logical, parameter :: use_AGG = .false.
#endif
#ifdef WLIN
  logical, parameter :: use_WLIN = .true.
#else
  logical, parameter :: use_WLIN = .false.
#endif
#ifdef natDIC
  logical, parameter :: use_natDIC = .true.
#else
  logical, parameter :: use_natDIC = .false.
#endif
#ifdef CFC
  logical, parameter :: use_CFC = .true.
#else
  logical, parameter :: use_CFC = .false.
#endif
#ifdef cisonew
  logical, parameter :: use_cisonew = .true.
#else
  logical, parameter :: use_cisonew = .false.
#endif
#ifdef PBGC_OCNP_TIMESTEP
  logical, parameter :: use_PBGC_OCNP_TIMESTEP = .true.
#else
  logical, parameter :: use_PBGC_OCNP_TIMESTEP = .false.
#endif
#ifdef PBGC_CK_TIMESTEP
  logical, parameter :: use_PBGC_CK_TIMESTEP = .true.
#else
  logical, parameter :: use_PBGC_CK_TIMESTEP = .false.
#endif
#ifdef FB_BGC_OCE
  logical, parameter :: use_FB_BGC_OCE = .true.
#else
  logical, parameter :: use_FB_BGC_OCE = .false.
#endif
#ifdef BOXATM
  logical, parameter :: use_BOXATM = .true.
#else
  logical, parameter :: use_BOXATM = .false.
#endif
#ifdef sedbypass
  logical, parameter :: use_sedbypass = .true.
#else
  logical, parameter :: use_sedbypass = .false.
#endif
#ifdef PROGCO2
  logical, parameter :: use_PROGCO2 = .true.
#else
  logical, parameter :: use_PROGCO2 = .false.
#endif
#ifdef DIAGCO2
  logical, parameter :: use_DIAGCO2 = .true.
#else
  logical, parameter :: use_DIAGCO2 = .false.
#endif

end module mo_ifdefs
