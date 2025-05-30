! ------------------------------------------------------------------------------
! Copyright (C) 2008-2025 Mats Bentsen
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

      program blom
c
c --- ------------------------------------------------------------------
c --- Bergen Layered Ocean Model (BLOM)
c --- ------------------------------------------------------------------
c
      use mod_time, only: nstep2, nstep
      use mod_xc
      use mod_state, only: dp
      use mod_checksum, only: chksummsk
      use mod_blom_init, only: blom_init_phase1, blom_init_phase2
      use mod_blom_step, only: blom_step
c
      implicit none

      integer :: nfu
c
c --- initialize the model
      call blom_init_phase1
      call blom_init_phase2
c
c --- advance the model from time step nstep1 to nstep2
      blom_loop: do
        call blom_step
        if (nstep.eq.nstep2) exit blom_loop
      enddo blom_loop
c
c --- write check sum of layer thickness
      call chksummsk(dp(1-nbdy,1-nbdy,1+mod(nstep2,2)*kk),ip,kk,'dp')
c
      if (mnproc.eq.1) then
        open (newunit=nfu,file='run.status',status='unknown')
        write (nfu,*) 'success'
        close (unit=nfu)
      endif
      call xcstop('(normal)')
             stop '(normal)'
c
      end
