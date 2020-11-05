! ------------------------------------------------------------------------------
! Copyright (C) 2013-2018 Mehmet Ilicak, Mats Bentsen
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

      subroutine initke
c
c --- ------------------------------------------------------------------
c --- Initialization of second order turbulence closure.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_trc.h"
#include "common_trc.h"
#include "param_tke.h"
#include "common_tke.h"
c
      integer i,j,k,l
c
#if defined(TRC) && defined(TKE)
c
c --- ------------------------------------------------------------------
c --- Initialize fields holding turbulent kinetic energy, generic length
c --- scale, and other fields used in the turbulence closure.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1-nbdy,jj+nbdy
        do k=1,2*kdm
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            trc(i,j,k,itrtke)=tke_min
            trc(i,j,k,itrgls)=gls_psi_min
          enddo
          enddo
        enddo
        do k=1,kk
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            difdia(i,j,k)=0.
            L_scale(i,j,k)=Ls_unlmt_min
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          ustarb(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Precompute various coefficients
c --- ------------------------------------------------------------------
c
      sqrt2=sqrt(2.)
      cmu_fac1=gls_cmu0**(-gls_p/gls_n)
      cmu_fac2=gls_cmu0**(3.+gls_p/gls_n)
      cmu_fac3=sqrt2
      tke_exp1=gls_m/gls_n
      gls_exp1=1./gls_n
      gls_fac6=8./gls_cmu0**6
      gls_s0=1.5*gls_L1*gls_L5**2
      gls_s1=-gls_L4*(gls_L6+gls_L7)
     .       +2.*gls_L4*gls_L5*(gls_L1-1./3.*gls_L2-gls_L3)
     .       +1.5*gls_L1*gls_L5*gls_L8
      gls_s2=-3./8.*gls_L1*(gls_L6**2-gls_L7**2)
      gls_s4=2.*gls_L5
      gls_s5=2.*gls_L4
      gls_s6=2./3.*gls_L5*(3.*gls_L3**2-gls_L2**2)
     .      -.5*gls_L5*gls_L1*(3.*gls_L3-gls_L2)
     .      +.75*gls_L1*(gls_L6-gls_L7)
      gls_b0=3.*gls_L5**2
      gls_b1=gls_L5*(7.*gls_L4+3.*gls_L8)
      gls_b2=gls_L5**2*(3.*gls_L3**2-gls_L2**2)
     .      -.75*(gls_L6**2-gls_L7**2)
      gls_b3=gls_L4*(4.*gls_L4+3.*gls_L8)
      gls_b4=gls_L4*(gls_L2*gls_L6-3.*gls_L3*gls_L7
     .      -gls_L5*(gls_L2**2-gls_L3**2))
     .      +gls_L5*gls_L8*(3.*gls_L3**2-gls_L2**2)
      gls_b5=.25*(gls_L2**2-3.*gls_L3**2)*(gls_L6**2-gls_L7**2)
c
#endif
c
      return
      end
