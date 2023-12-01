! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen
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
      use mod_time, only: nstep1, nstep2, nstep
      use mod_xc
      use mod_state, only: dp
      use mod_cesm
      use mod_checksum, only: chksummsk
      use mod_blom_init, only: blom_init
      use mod_blom_step, only: blom_step
c
      implicit none
c
#ifdef MPI
#  include "mpif.h"
#endif
c
      integer i
c
#ifdef MPI
      call mpi_init
      mpicom_external=mpi_comm_world
#endif
c
      open (unit=nfu,file='ocn_in',status='old',action='read')
      read (nfu,'(i6)') i
      read (nfu,*)
      read (nfu,*)
      read (nfu,*)
      read (nfu,'(a)') runid_cesm   ! Experiment name
      close (unit=nfu)
      if (i.eq.0) then
        runtyp_cesm='initial'
      else
        runtyp_cesm='continue'
      endif
c
c --- initialize the model
      call blom_init
c
c --- advance the model from time step nstep1 to nstep2
      blom_loop: do
        call blom_step
c       if (nstep.eq.nstep2) exit blom_loop
        if (nstep.eq.nstep1+2) exit blom_loop
      enddo blom_loop
c
c --- write check sum of layer thickness
      call chksummsk(dp(1-nbdy,1-nbdy,1+mod(nstep2,2)*kk),ip,1,'dp')
c
      if (mnproc.eq.1) then
        open (unit=nfu,file='run.status',status='unknown')
        write (nfu,*) 'success'
        close (unit=nfu)
      endif
      call xcstop('(normal)')
             stop '(normal)'
c
      end
