c
c --- ------------------------------------------------------------------
c --- Parameters related to tracers:
c ---   ntrhamocc  - number of tracers for HAMOCC
c ---   ntriage    - number of tracers for ideal age tracer
c ---   ntr        - total number of tracers
c --- ------------------------------------------------------------------
c
      integer ntrhamocc,ntriage,ntr
c
c --- ------------------------------------------------------------------
c --- HAMOCC tracers
c --- ------------------------------------------------------------------
c
#ifdef HAMOCC
      integer i_base_adv,i_iso_adv,i_cfc_adv,i_agg_adv,ntraad,i_base,
     .        i_iso
c
c --- Advected HAMOCC tracers
      parameter (i_base_adv=14)      
#  ifdef __c_isotopes
      parameter (i_iso_adv=2)      
#  else 
      parameter (i_iso_adv=0)      
#  endif
#  ifdef PCFC  
      parameter (i_cfc_adv=3)
#  else 
      parameter (i_cfc_adv=0)
#  endif
#  ifdef AGG
      parameter (i_agg_adv=2)
#  else 
      parameter (i_agg_adv=0)
#  endif
      parameter (ntraad=i_base_adv+i_iso_adv+i_cfc_adv+i_agg_adv)
c
c --- Non-advected (fast sinking) HAMOCC tracers
      parameter (i_base=3) 
#  ifdef __c_isotopes
      parameter (i_iso=4)  
#  else
      parameter (i_iso=0)  
#  endif
c
c --- Total number of HAMOCC tracers
      parameter (ntrhamocc=ntraad+i_base+i_iso)
#else
      parameter (ntrhamocc=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Ideal age tracers
c --- ------------------------------------------------------------------
c
#ifdef IAGE
      parameter (ntriage=1)
#else
      parameter (ntriage=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Total number of tracers
c --- ------------------------------------------------------------------
c
      parameter (ntr=ntrhamocc+ntriage)
c
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
