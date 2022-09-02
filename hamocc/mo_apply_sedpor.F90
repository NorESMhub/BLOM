! Copyright (C) 2020  S. Gao, I. Bethke, J. Tjiputra, J. Schwinger
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


module mo_apply_sedpor
!*****************************************************************************
! Purpose
! -------
!   - Routine for applying sediment porosity
!
! Description
! -----------
!
! Subroutine apply_sedpor applies formerly read lon-lat-sediment depth variable
! sediment porosity  
!*****************************************************************************

implicit none

private

public :: apply_sedpor

contains 

subroutine apply_sedpor(kpie,kpje,ks,omask,sed_por)
 use mo_sedmnt,      only: porwat

 implicit none

 integer,intent(in) :: kpie,kpje,ks
 real,   intent(in) :: omask(kpie,kpje)
 real,   intent(in) :: sed_por(kpie,kpje,ks)

 ! local variables
 integer :: i,j,k

 do j=1,kpje
 do i=1,kpie
 do k=1,ks
   if(omask(i,j).gt. 0.5)then
     porwat(i,j,k) = sed_por(i,j,k)
   endif
 enddo
 enddo
 enddo

end subroutine apply_sedpor
end module mo_apply_sedpor


