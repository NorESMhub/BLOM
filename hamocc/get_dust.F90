! Copyright (C) 2002  I. Kriest, Patrick Wetzel
! Copyright (C) 2020  J. Tjiputra, J. Schwinger
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


      SUBROUTINE GET_DUST(kpie,kpje,kpke,omask,path)
!****************************************************************
!
!**** *GET_DUST* - 
!
!     Iris Kriest,    *MPI-Met, HH*    18.10.02
!
!     Modified
!     --------
!     Patrick Wetzel : read data in NetCDF format.
!     
!     Purpose
!     -------
!     - read monthly dust fluxes from nudged data set by C. Timmreck.
!
!     Method
!     -------
!     - read monthly dust fluxes [kg/m2/month] data into array dustin(i,j,k). 
!     - if there is a wet cell write value into array dusty.
!     - if there is a dry cell assign 0.
!
!**   Interface.
!     ----------
!     *CALL*       *GET_DUST*
!     *MODULE*     *MO_CARBCH* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER*   *kpie*  - 1st dimension of model grid.
!     *INTEGER*   *kpje*  - 2nd dimension of model grid.
!     *INTEGER*   *kpke*  - 3rd (vertical) dimension of model grid.
!     *REAL*      *omask* - ocean mask.
!     *CHARACTER* *path*  - path to dust data file
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
      USE netcdf
      USE mo_carbch
      USE mod_xc 

      implicit none
      INTEGER         :: kpie,kpje,kpke
      REAL            :: omask(kpie,kpje)
      character(len=*):: path

! Local variables
      INTEGER         :: i,j,k,l,ncid,ncstat,ncvarid
      REAL            :: dustin(kpie,kpje,12)

!
! Open netCDF data file
!      
       IF(mnproc==1) THEN
       ncstat = NF90_OPEN(trim(path)//'INPDUST_mhw.nc', &
     &           NF90_NOWRITE, ncid)
        IF (ncstat.NE.NF90_NOERR ) THEN
         CALL xchalt('(get_dust: Problem with netCDF1)')
                stop '(get_dust: Problem with netCDF1)'
        END IF
       END IF
!
! Read  data
       call read_netcdf_var(ncid,'DUST',dustin(1,1,1),12,0,0)

!
! Close file
       IF(mnproc==1) THEN
        ncstat = NF90_CLOSE(ncid)
        IF ( ncstat .NE. NF90_NOERR ) THEN
         CALL xchalt('(get_dust: Problem with netCDF200)')
                stop '(get_dust: Problem with netCDF200)'
        END IF
       END IF


! set missings over land
      do l=1,12
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j).gt.0.5) then 
              dusty(i,j,l) = dustin(i,j,l)
            else 
              dusty(i,j,l) = 0.
            endif    
          enddo
         enddo
      enddo

      RETURN
      END
