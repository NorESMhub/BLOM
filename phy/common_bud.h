! ------------------------------------------------------------------------------
! Copyright (C) 2007-2015 Mats Bentsen
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
c --- common blocks related to budget computations
c
      c o m m o n  /bud1/
     . sdp(ncalls,2),tdp(ncalls,2)
#ifdef TKE
     .,tkedp(ncalls,2)
#  ifdef GLS
     .,glsdp(ncalls,2)
#  endif
#endif
#ifdef TRC
     .,trdp(ncalls,2)
#endif
     .,mass0,sf,tf
#ifdef TRC
     .,trf
#endif
c
      real sdp,tdp
#ifdef TKE
     .    ,tkedp
#  ifdef GLS
     .    ,glsdp
#  endif
#endif
#ifdef TRC
     .    ,trdp
#endif
     .    ,mass0,sf,tf
#ifdef TRC
     .    ,trf
#endif
