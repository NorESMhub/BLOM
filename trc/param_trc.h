! ------------------------------------------------------------------------------
! Copyright (C) 2007-2015 Mats Bentsen, Jörg Schwinger, Jerry Tjiputra,
!                         Alok Kumar Gupta
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

c
c --- ------------------------------------------------------------------
c --- Parameters related to tracers:
c ---   ntrocn  - number of ocean tracers
c ---   ntrtke  - number of turbulent kinetic energy tracers
c ---   ntrgls  - number of generic length scale tracers
c ---   ntriag  - number of ideal age tracers
c ---   ntrbgc  - number of HAMOCC tracers
c ---   ntr     - total number of tracers
c ---   itrtke  - index of first turbulent kinetic energy tracer
c ---   itrgls  - index of first generic length scale tracer
c ---   itriag  - index of first ideal age tracer
c ---   itrbgc  - index of first HAMOCC tracer
c --- ------------------------------------------------------------------
c
      integer ntrocn,ntrtke,ntrgls,ntriag,ntrbgc,ntr,
     .        itrtke,itrgls,itriag,itrbgc
c
c --- ------------------------------------------------------------------
c --- Ocean tracers, not including turbulent kinetic energy, generic
c --- length scale, ideal age and HAMOCC tracers
c --- ------------------------------------------------------------------
c
      parameter (ntrocn=0)
c
c --- ------------------------------------------------------------------
c --- Turbulent kinetic energy and generic length scale tracer
c --- ------------------------------------------------------------------
c
#ifdef TKE
      parameter (ntrtke=1,ntrgls=1)
#else
      parameter (ntrtke=0,ntrgls=0)
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
c --- HAMOCC tracers
c --- ------------------------------------------------------------------
c
#ifdef HAMOCC
      integer i_base,i_iso,i_cfc,i_agg,i_nat_dic
c
c --- Advected HAMOCC tracers
      parameter (i_base=22)      
#  ifdef cisonew
      parameter (i_iso=12)      
#  else 
      parameter (i_iso=0)      
#  endif
#  ifdef CFC  
      parameter (i_cfc=3)
#  else 
      parameter (i_cfc=0)
#  endif
#  ifdef AGG
      parameter (i_agg=2)
#  else 
      parameter (i_agg=0)
#  endif
#  ifdef natDIC
      parameter (i_nat_dic=3)
#  else
      parameter (i_nat_dic=0)
#  endif
c --- Total number of HAMOCC tracers
      parameter (ntrbgc=i_base+i_iso+i_cfc+i_agg+i_nat_dic)
#else
      parameter (ntrbgc=0)
#endif
c
c --- ------------------------------------------------------------------
c --- Total number of tracers
c --- ------------------------------------------------------------------
c
      parameter (ntr=ntrocn+ntrtke+ntrgls+ntriag+ntrbgc)
c
#ifdef ATRC
c --- ------------------------------------------------------------------
c --- Parameters related to age tracers:
c ---   natr  - number of age tracers
c --- ------------------------------------------------------------------
c
      integer natr
      parameter (natr=0)
c
#endif
c --- ------------------------------------------------------------------
c --- Set tracer indexes of first turbulent kinetic energy, generic
c --- length scale, ideal age and HAMOCC tracer
c --- ------------------------------------------------------------------
c
#ifdef TKE
#  ifdef ATRC
      parameter (itrtke=ntrocn-natr+1,
     .           itrgls=ntrocn-natr+ntrtke+1)
#  else
      parameter (itrtke=ntrocn+1,
     .           itrgls=ntrocn+ntrtke+1)
#  endif
#else
      parameter (itrtke=-1,itrgls=-1)
#endif
c
#ifdef IDLAGE
#  ifdef ATRC
      parameter (itriag=ntrocn-natr+ntrtke+ntrgls+1)
#  else
      parameter (itriag=ntrocn+ntrtke+ntrgls+1)
#  endif
#else
      parameter (itriag=-1)
#endif
c
#ifdef HAMOCC
#  ifdef ATRC
      parameter (itrbgc=ntrocn-natr+ntrtke+ntrgls+ntriag+1)
#  else
      parameter (itrbgc=ntrocn+ntrtke+ntrgls+ntriag+1)
#  endif
#else
      parameter (itrbgc=-1)
#endif
