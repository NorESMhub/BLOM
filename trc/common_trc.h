! ------------------------------------------------------------------------------
! Copyright (C) 2007-2018 Mats Bentsen
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
c --- common blocks related to tracers
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm,ntr) :: trc
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ntr) :: trcold
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  uflxtr,vflxtr,trflx
c
      common /trc1/ trc,trcold,uflxtr,vflxtr,trflx
