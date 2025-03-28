! ------------------------------------------------------------------------------
! Copyright (C) 2021-2025 Mats Bentsen, Mehmet Ilicak
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

module mod_vcoord
! ------------------------------------------------------------------------------
! This module contains parameter, variables and procedures related to the
! vertical coordinate.
! ------------------------------------------------------------------------------

   use mod_types,     only: r8
   use mod_config,    only: inst_suffix
   use mod_constants, only: spval, onem
   use mod_xc

   implicit none
   private

   ! Parameters:
   integer, parameter :: &
      vcoord_isopyc_bulkml = 1, & ! Vertical coordinate type: bulk surface
                                  ! mixed layer with isopycnic layers below.
      vcoord_cntiso_hybrid = 2, & ! Vertical coordinate type: Hybrid
                                  ! coordinate with pressure coordinates towards
                                  ! the surface and continuous isopycnal below.
      vcoord_plevel        = 3, & ! Vertical coordinate type: pressure
                                  ! coordinate.
      kdm_max              = 1000 ! Maximum anticipated vertical dimension.

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      vcoord_type            = 'isopyc_bulkml', &
      sigref_spec            = 'inicon', &
      plevel_spec            = 'inflation'
   real(r8) :: &
      dpmin_surface          = 1.5_r8, &
      dpmin_inflation_factor = 1._r8
   real(r8), dimension(kdm_max) :: &
      sigref                 = spval, &
      plevel                 = spval

   ! Options derived from string options.
   integer :: &
      vcoord_tag

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
      sigmar                      ! Reference potential density [kg m-3].

   public :: vcoord_tag, vcoord_isopyc_bulkml, vcoord_cntiso_hybrid, &
             vcoord_plevel, sigref_spec, sigmar, sigref, plevel, &
             readnml_vcoord, inivar_vcoord

contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine readnml_vcoord
   ! ---------------------------------------------------------------------------
   ! Read variables in the namelist group 'vcoord' and resolve options.
   ! ---------------------------------------------------------------------------

      character(len = 80) :: nml_fname
      real(r8) :: dpmin
      integer :: nfu, ios, k
      logical :: fexist

      namelist /vcoord/ &
         vcoord_type, dpmin_surface, dpmin_inflation_factor, &
         sigref_spec, plevel_spec, sigref, plevel

      ! Read variables in the namelist group 'vcoord'.
      if (mnproc == 1) then
         nml_fname = 'ocn_in'//trim(inst_suffix)
         inquire(file = nml_fname, exist = fexist)
         if (fexist) then
            open (newunit = nfu, file = nml_fname, status = 'old', &
                 action = 'read')
         else
            nml_fname = 'limits'//trim(inst_suffix)
            inquire(file = nml_fname, exist = fexist)
            if (fexist) then
               open (newunit = nfu, file = nml_fname, status = 'old', &
                    action = 'read')
            else
               write (lp,*) 'readnml_vcoord: could not find namelist file!'
               call xchalt('(readnml_vcoord)')
               stop '(readnml_vcoord)'
            endif
         endif
         read (unit = nfu, nml = vcoord, iostat = ios)
         close (unit = nfu)
      endif
      call xcbcst(ios)
      if (ios /= 0) then
         if (mnproc == 1) &
            write (lp,*) 'readnml_vcoord: No vertical coordinate variable '//  &
                         'group found in namelist. Using defaults.'
      else
         call xcbcst(vcoord_type)
         call xcbcst(dpmin_surface)
         call xcbcst(dpmin_inflation_factor)
         call xcbcst(sigref_spec)
         call xcbcst(plevel_spec)
         call xcbcst(sigref)
         call xcbcst(plevel)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_vcoord: vertical coordinate variables:'
         write (lp,*) '  vcoord_type =            ', trim(vcoord_type)
         write (lp,*) '  dpmin_surface =          ', dpmin_surface
         write (lp,*) '  dpmin_inflation_factor = ', dpmin_inflation_factor
         write (lp,*) '  sigref_spec =            ', trim(sigref_spec)
         write (lp,*) '  plevel_spec =            ', trim(plevel_spec)
      endif

      ! Change units from [m] to [kg m-1 s-2] of depth interval variables.
      dpmin_surface = dpmin_surface*onem

      ! Resolve options.
      select case (trim(vcoord_type))
         case ('isopyc_bulkml')
            vcoord_tag = vcoord_isopyc_bulkml
         case ('cntiso_hybrid')
            vcoord_tag = vcoord_cntiso_hybrid
         case ('plevel')
            vcoord_tag = vcoord_plevel
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') ' readnml_vcoord: vcoord_type = ', &
                                 trim(vcoord_type), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
            stop '(readnml_vcoord)'
      end select
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
         select case (trim(sigref_spec))
            case ('inicon')
            case ('namelist')
               k = 1
               do while (sigref(k) /= spval)
                  k = k + 1
                  if (k > kdm_max) exit
               enddo
               if (k /= kdm + 1) then
                  if (mnproc == 1) &
                     write (lp,'(3a)') &
                        ' readnml_vcoord: number of sigref values does not match vertical dimension!'
                  call xcstop('(readnml_vcoord)')
                  stop '(readnml_vcoord)'
               endif
            case default
               if (mnproc == 1) &
                  write (lp,'(3a)') ' readnml_vcoord: sigref_spec = ', &
                                    trim(sigref_spec), ' is unsupported!'
               call xcstop('(readnml_vcoord)')
               stop '(readnml_vcoord)'
         end select
         select case (trim(plevel_spec))
            case ('inflation')
               dpmin = dpmin_surface
               plevel(1) = 0._r8
               do k = 1, kk - 1
                  plevel(k+1) = plevel(k) + dpmin
                  dpmin = dpmin*dpmin_inflation_factor
               enddo
            case ('namelist')
               k = 1
               do while (plevel(k) /= spval)
                  k = k + 1
                  if (k > kdm_max) exit
               enddo
               if (k /= kdm + 1) then
                  if (mnproc == 1) &
                     write (lp,'(3a)') &
                        ' readnml_vcoord: number of plevel values does not match vertical dimension!'
                  call xcstop('(readnml_vcoord)')
                  stop '(readnml_vcoord)'
               endif
               ! Change units from [m] to [kg m-1 s-2].
               plevel(:) = plevel(:)*onem
            case default
               if (mnproc == 1) &
                  write (lp,'(3a)') ' readnml_vcoord: plevel_spec = ', &
                                    trim(plevel_spec), ' is unsupported!'
               call xcstop('(readnml_vcoord)')
               stop '(readnml_vcoord)'
         end select
      endif

   end subroutine readnml_vcoord

   subroutine inivar_vcoord
   ! ---------------------------------------------------------------------------
   ! Initialize arrays and data structures.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k

      if (vcoord_tag == vcoord_isopyc_bulkml .or. &
          trim(sigref_spec) == 'inicon') then
         !$omp parallel do private(i, k)
         do j = 1-nbdy, jj+nbdy
            do k = 1, kk
               do i = 1-nbdy, ii+nbdy
                  sigmar(i,j,k) = spval
               enddo
            enddo
         enddo
         !$omp end parallel do
      else
         !$omp parallel do private(i, k)
         do j = 1-nbdy, jj+nbdy
            do k = 1, kk
               do i = 1-nbdy, ii+nbdy
                  sigmar(i,j,k) = sigref(k)
               enddo
            enddo
         enddo
         !$omp end parallel do
      endif

   end subroutine inivar_vcoord

end module mod_vcoord
