! ------------------------------------------------------------------------------
! Copyright (C) 2020-2022 Mats Bentsen
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

module mod_diffusion
! ------------------------------------------------------------------------------
! This module contains variables related to diffusion.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval, epsil
   use mod_xc

   implicit none

   private

   ! Variables to be set by namelist.
   real(r8) :: &
      egc, &    ! The parameter 'c' in the Eden and Greatbatch (2008)
                ! parameterization [].
      eggam, &  ! The parameter 'gamma' in the Eden and Greatbatch (2008)
                ! parameterization [].
      eglsmn, & ! Minimum eddy length scale in the Eden and Greatbatch (2008)
                ! parameterization [cm].
      egmndf, & ! Minimum diffusivity in the Eden and Greatbatch (2008)
                ! parameterization [cm2 s-1].
      egmxdf, & ! Maximum diffusivity in the Eden and Greatbatch (2008)
                ! parameterization [cm2 s-1].
      egidfq, & ! Factor relating the isopycnal diffusivity to the layer
                ! interface diffusivity in the Eden and Greatbatch (2008)
                ! parameterization (egidfq = difint/difiso) [].
      ri0, &    ! Critical gradient richardson number for shear driven vertical
                ! mixing [].
      bdmc1, &  ! Background diapycnal diffusivity times buoyancy frequency
                ! [cm2 s-2].
      bdmc2, &  ! Background diapycnal diffusivity [cm2 s-1].
      tkepf     ! Fraction of surface TKE that penetrates beneath mixed layer
                ! [].
   integer :: &
      bdmtyp    ! Type of background diapycnal mixing. If bdmtyp = 1 the
                ! background diffusivity is a constant divided by the
                ! Brunt-Vaisala frequency, if bdmtyp = 2 the background
                ! diffusivity is constant [].
   logical :: &
      edsprs    ! If true, apply eddy mixing suppression away from steering
                ! level.
   character(len = 80) :: &
      eitmth, & ! Eddy-induced transport parameterization method. Valid
                ! methods: 'intdif', 'gm'.
      edritp, & ! Type of Richardson number used in eddy diffusivity
                ! computation. Valid types: 'shear', 'large scale'.
      edwmth    ! Method to estimate eddy diffusivity weight as a function of
                ! the ration of Rossby radius of deformation to the horizontal
                ! grid spacing. Valid methods: 'smooth', 'step'.
   logical :: &
      ntrdif = .false.

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kdm) :: &
      difint, & ! Layer interface diffusivity [cm2 s-1].
      difiso, & ! Isopycnal diffusivity [cm2 s-1].
      difdia    ! Diapycnal diffusivity [cm2 s-1].

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kdm+1) :: &
      Kvisc_m, &      ! momentum eddy viscosity [cm2 s-1].
      Kdiff_t, &      ! temperature eddy diffusivity [cm2 s-1].
      Kdiff_s, &      ! salinity eddy  diffusivity [cm2 s-1].
      t_ns_nonloc, &  ! Non-local transport term that is the fraction of
                      ! non-shortwave flux passing a layer interface [].
      s_nonloc        ! Non-local transport term that is the fraction of
                      ! material tracer flux passing a layer interface [].

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      difmxp, & ! Maximum lateral diffusivity at p-points [cm2 s-1].
      difmxq, & ! Maximum lateral diffusivity at q-points [cm2 s-1].
      difwgt    ! Eddy diffusivity weight [].


   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, 2*kdm) :: &
      umfltd, & ! u-component of horizontal mass flux due to thickness diffusion
                ! [g cm s-2].
      vmfltd, & ! v-component of horizontal mass flux due to thickness diffusion
                ! [g cm s-2].
      utfltd, & ! u-component of horizontal heat flux due to thickness diffusion
                ! [K g cm s-2].
      vtfltd, & ! v-component of horizontal heat flux due to thickness diffusion
                ! [K g cm s-2].
      utflld, & ! u-component of horizontal heat flux due to lateral diffusion
                ! [K g cm s-2].
      vtflld, & ! v-component of horizontal heat flux due to lateral diffusion
                ! [K g cm s-2].
      usfltd, & ! u-component of horizontal salt flux due to thickness diffusion
                ! [g2 cm kg-1 s-2].
      vsfltd, & ! v-component of horizontal salt flux due to thickness diffusion
                ! [g2 cm kg-1 s-2].
      usflld, & ! u-component of horizontal salt flux due to lateral diffusion
                ! [g2 cm kg-1 s-2].
      vsflld    ! v-component of horizontal salt flux due to lateral diffusion
                ! [g2 cm kg-1 s-2].

   public :: egc, eggam, eglsmn, egmndf, egmxdf, egidfq, ri0, bdmc1, bdmc2, &
             tkepf, bdmtyp, edsprs, eitmth, edritp, edwmth, ntrdif, &
             difint, difiso, difdia, difmxp, difmxq, difwgt, &
             umfltd, vmfltd, utfltd, vtfltd, utflld, vtflld, &
             usfltd, vsfltd, usflld, vsflld, &
             Kvisc_m, Kdiff_t, Kdiff_s, t_ns_nonloc, s_nonloc, &
             inivar_diffusion

contains

   subroutine inivar_diffusion
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l

   !$omp parallel do private(i, k)
      do j = 1 - nbdy, jj + nbdy
         do i = 1 - nbdy, ii + nbdy
            difmxp(i, j) = spval
            difmxq(i, j) = spval
            difwgt(i, j) = spval
         enddo
         do k = 1, kk
            do i = 1 - nbdy, ii + nbdy
               difint(i, j, k) = spval
               difiso(i, j, k) = spval
               difdia(i, j, k) = spval
            enddo
         enddo
         do k = 1, 2*kk
            do i = 1 - nbdy, ii + nbdy
               umfltd(i, j, k) = spval
               vmfltd(i, j, k) = spval
               utfltd(i, j, k) = spval
               vtfltd(i, j, k) = spval
               utflld(i, j, k) = spval
               vtflld(i, j, k) = spval
               usfltd(i, j, k) = spval
               vsfltd(i, j, k) = spval
               usflld(i, j, k) = spval
               vsflld(i, j, k) = spval
            enddo
         enddo
         do k = 1, kk+1
            do i = 1 - nbdy, ii + nbdy
               Kvisc_m(i, j, k) = epsil
               Kdiff_t(i, j, k) = epsil
               Kdiff_s(i, j, k) = epsil
            enddo
         enddo
      enddo
   !$omp end parallel do

   ! Initialize isopycnal diffusivity .
      do j = 1, jj
         do k = 1, kk
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               difiso(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(difiso, 1, kk, nbdy, nbdy, halo_ps)

   ! Initialize diffusive fluxes at points located upstream and downstream (in
   ! i-direction) of p-points.
   !$omp parallel do private(k, l, i)
      do j = 1, jj
         do k = 1, 2*kk
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l) + 1)
               umfltd(i, j, k) = 0._r8
               utfltd(i, j, k) = 0._r8
               utflld(i, j, k) = 0._r8
               usfltd(i, j, k) = 0._r8
               usflld(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(umfltd, 1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(utflld, 1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(usflld, 1, 2*kk, nbdy, nbdy, halo_us)

   ! Initialize diffusive fluxes at points located upstream and downstream (in
   ! j-direction) of p-points.
   !$omp parallel do private(k, l, j)
      do i = 1, ii
         do k = 1, 2*kk
            do l = 1, jsp(i)
            do j = max(1, jfp(i, l)), min(jj, jlp(i, l) + 1)
               vmfltd(i, j, k) = 0._r8
               vtfltd(i, j, k) = 0._r8
               vtflld(i, j, k) = 0._r8
               vsfltd(i, j, k) = 0._r8
               vsflld(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(vmfltd, 1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vtflld, 1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vsflld, 1, 2*kk, nbdy, nbdy, halo_vs)

   end subroutine inivar_diffusion

end module mod_diffusion
