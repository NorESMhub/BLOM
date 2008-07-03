c
c --- ------------------------------------------------------------------
c --- Parameters related to tracers:
c ---   ntr  - number of tracers
c --- ------------------------------------------------------------------
c
      integer ntr
      parameter (ntr=2)
#ifdef ATRC
c
c --- ------------------------------------------------------------------
c --- Parameters related to age tracers:
c ---   natr  - number of age tracers
c --- ------------------------------------------------------------------
c
      integer natr
      parameter (natr=1)
#endif
