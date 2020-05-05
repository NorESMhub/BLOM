! ------------------------------------------------------------------------------
! Copyright (C) 2008 Mats Bentsen
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

      program micom
c
c --- ---------------------------------------------
c --- miami isopycnic coordinate ocean model
c
c --- NERSC version 2.0, based on MICOM version 2.8
c --- ---------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
c
c --- initialize the model
      call micom_init
c
c --- advance the model from time step nstep1 to nstep2
      micom_loop: do
        call micom_step
        if (nstep.eq.nstep2) exit micom_loop
      enddo micom_loop
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
