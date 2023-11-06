! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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


subroutine powadi(j,kpie,kpje,solrat,sedb1,sediso,omask)
  !**********************************************************************
  !
  !**** *POWADI* - vertical diffusion with simultaneous dissolution.
  !
  !     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
  !
  !     Modified
  !     --------
  !     S.Legutke,        *MPI-MaD, HH*    10.04.01
  !
  !     Purpose
  !     -------
  !     .
  !
  !     Method
  !     -------
  !     implicit discretisation.
  !
  !**   Interface.
  !     ----------
  !
  !     *CALL*       *POWADI(j,solrat,sedb1,sediso)*
  !
  !     Input  solrat : dissolution rate
  !     =====       j : zonal grid index
  !             sedb1 : tracer at entry
  !
  !     Output: sediso: diffused tracer at exit
  !     ======
  !
  !     *PARAMETER*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
  !
  !     Externals
  !     ---------
  !     none.
  !
  !**********************************************************************

  use mo_sedmnt,     only: porwah,porwat,seddw,seddzi
  use mo_param_bgc,  only: sedict
  use mo_param1_bgc, only: ks
  use mo_vgrid,      only: bolay

  implicit none

  integer, intent(in) :: j, kpie, kpje
  real, dimension(kpie,ks),   intent(in)    :: solrat
  real, dimension(kpie,0:ks), intent(inout) :: sedb1, sediso
  real, dimension(kpie,kpje), intent(in)    :: omask

  ! Local variables
  integer :: i,k,l
  real :: asu, alo
  real, dimension(kpie,0:ks,3) :: tredsy

  !**********************************************************************

  do k = 1, ks
    do i = 1, kpie
      asu = sedict * seddzi(k) * porwah(i,j,k)
      alo = 0.
      if(k < ks) alo = sedict * seddzi(k+1) * porwah(i,j,k+1)
      tredsy(i,k,1) = -asu
      tredsy(i,k,3) = -alo
      tredsy(i,k,2) = seddw(k) * porwat(i,j,k) - tredsy(i,k,1)                   &
           &          - tredsy(i,k,3) + solrat(i,k) * porwat(i,j,k) * seddw(k)
    enddo
  enddo

  k = 0
  asu = 0.
  do i = 1, kpie
    alo = sedict * seddzi(1) * porwah(i,j,1)
    if(omask(i,j) > 0.5) then
      tredsy(i,k,1) = -asu
      tredsy(i,k,3) = -alo
      tredsy(i,k,2) = bolay(i,j) - tredsy(i,k,1) - tredsy(i,k,3)
    else
      tredsy(i,k,1) = 0
      tredsy(i,k,3) = 0
      tredsy(i,k,2) = 0
    endif
  enddo

  do k = 1, ks
    do i = 1, kpie
      if(omask(i,j) > 0.5) then
        tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)
        tredsy(i,k,2)   = tredsy(i,k,2)                                     &
             &            - tredsy(i,k-1,3) * tredsy(i,k,1) / tredsy(i,k-1,2)
      endif
    enddo
  enddo

  do k = 1, ks
    do i = 1, kpie
      sedb1(i,k) = sedb1(i,k) - tredsy(i,k-1,1) * sedb1(i,k-1)
    enddo
  enddo

  k = ks
  do i = 1, kpie
    if(omask(i,j) > 0.5) sediso(i,k) = sedb1(i,k) / tredsy(i,k,2)
  enddo

  do k = 1, ks
    l = ks - k
    do i = 1, kpie
      if(omask(i,j) > 0.5) then
        sediso(i,l) = ( sedb1(i,l) - tredsy(i,l,3) * sediso(i,l+1) )        &
             &        / tredsy(i,l,2)
      endif
    enddo
  enddo

end subroutine powadi
