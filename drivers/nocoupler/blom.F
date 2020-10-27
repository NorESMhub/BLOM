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
      use mod_xc
      use cudafor
      use mytests
c
      implicit none
c
#include "common_blocks.h"

c --- testing CUDA Fortran
c --- sample code from https://docs.nvidia.com/hpc-sdk/compilers/
c ---     cuda-fortran-prog-guide/index.html#cfpg-cuda-fort-host-dev-code
      integer, parameter :: n = 100
      integer, allocatable, device :: iarr(:)
      integer h(n)
      integer istat
      integer i

      if (0) then
        istat = cudaSetDevice(0)
        allocate(iarr(n))
        h = 0; iarr = h
        call test1<<<1,n>>> (iarr)
        h = iarr
        print *, "Errors: ", count(h.ne.(/ (i,i=1,n) /))
        deallocate(iarr)

        stop '(normal)'
      endif

c
c --- initialize the model
      call blom_init
c
c --- advance the model from time step nstep1 to nstep2
      blom_loop: do
        call blom_step
        if (nstep.eq.nstep2) exit blom_loop
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
