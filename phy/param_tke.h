! ------------------------------------------------------------------------------
! Copyright (C) 2013 Mehmet Ilicak, Mats Bentsen
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
c --- Parameters related to tke equation computation:
c ---   glc_cmu0 ==> cmu0
c ---   Pr_t ==> Turbulent Prandtl number
c ---   tke_min ==> minimum TKE
c ---   gls_psi_min ==> min psi value
c --- ------------------------------------------------------------------
c
      real gls_cmu0,Pr_t,tke_min,zos,gls_psi_min,gls_p,gls_n,gls_m,
     .     gls_c1,gls_c2,gls_c3plus,gls_c3minus,gls_Gh0,gls_Ghmin,
     .     gls_Ghcri,gls_L1,gls_L2,gls_L3,gls_L4,gls_L5,gls_L6,gls_L7,
     .     gls_L8,vonKar,Ls_unlmt_min

      parameter (gls_cmu0=.527,Pr_t=1.,tke_min=7.6e-4,zos=.0002,
     .           gls_psi_min=1.e-10)
      parameter (gls_p=3.,gls_m=1.5,gls_n=-1.)
      parameter (gls_c1=1.44,gls_c2=1.92,gls_c3plus=1.,gls_c3minus=-.63)
      parameter (gls_L1=.107,gls_L2=.0032,gls_L3=.0864,gls_L4=.12,
     .           gls_L5=11.9,gls_L6=.4,gls_L7=.0,gls_L8=.48)
      parameter (gls_Gh0=.0329,gls_Ghmin=-.28,gls_Ghcri=.03)
      parameter (vonKar=.4,Ls_unlmt_min=1.e-6)
