c
c --- ------------------------------------------------------------------
c --- Parameters related to tracers:
c ---   ntrocn  - number of ocean tracers
c ---   ntrbgc  - number of HAMOCC tracers
c ---   ntriag  - number of ideal age tracers
c ---   ntr     - total number of tracers
c --- ------------------------------------------------------------------
c
      integer ntrocn,ntrbgc,ntriag,ntr
c
c --- ------------------------------------------------------------------
c --- Ocean tracers, not including HAMOCC and ideal age tracers
c --- ------------------------------------------------------------------
c
      parameter (ntrocn=0)
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
#  ifdef CFC  
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
      parameter (ntrbgc=ntraad+i_base+i_iso)
#else
      parameter (ntrbgc=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Ideal age tracer
c --- ------------------------------------------------------------------
c
#ifdef IDLAGE
      parameter (ntriag=1)
#else
      parameter (ntriag=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Total number of tracers
c --- ------------------------------------------------------------------
c
      parameter (ntr=ntrocn+ntrbgc+ntriag)
c
#ifdef ATRC
c
c --- ------------------------------------------------------------------
c --- Parameters related to age tracers:
c ---   natr  - number of age tracers
c --- ------------------------------------------------------------------
c
      integer natr
      parameter (natr=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Set tracer indexes of first HAMOCC tracer and ideal age tracer
c --- ------------------------------------------------------------------
c
#ifdef HAMOCC
      integer itrbgc
#  ifdef ATRC
      parameter (itrbgc=ntrocn-natr+1)
#  else
      parameter (itrbgc=ntrocn+1)
#  endif
#endif
c
#ifdef IDLAGE
      integer itriag
#  ifdef ATRC
      parameter (itriag=ntrocn-natr+ntrbgc+1)
#  else
      parameter (itriag=ntrocn+ntrbgc+1)
#  endif
#endif
