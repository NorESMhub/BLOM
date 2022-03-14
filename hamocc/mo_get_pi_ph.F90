! Copyright (C) 2020  J. Tjiputra
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

      MODULE mo_get_pi_ph

      implicit none
      private
      public :: get_pi_ph,pi_ph_path

      ! Path to input data, set through namelist 
      ! in hamocc_init.F 
      character(len=256),save    :: pi_ph_path = ''

      CONTAINS     
 
      SUBROUTINE GET_PI_PH(kpie,kpje,kpke,omask)
!**********************************************************************

      USE mo_carbch, only: pi_ph 
      USE mo_control_bgc, only: io_stdo_bgc 
      USE netcdf, only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open 
      USE mod_xc, only: mnproc,xchalt

      implicit none
      INTEGER, INTENT(in) :: kpie,kpje,kpke
      INTEGER ::i,j,l
  
      REAL,intent(in) ::omask(kpie,kpje)

! define the fields

      REAL :: pi_ph_in(kpie,kpje,12)

      INTEGER ncid,ncstat
!
! Open netCDF data file
!      
       IF(mnproc==1) THEN
        ncstat = NF90_OPEN(trim(pi_ph_path)//'MONTHLY_PI_PH.nc',   &
     &                   NF90_NOWRITE, ncid)
        write(io_stdo_bgc,*) 'HAMOCC: opening MONTHLY_PI_PH file'
        IF (ncstat.NE.NF90_NOERR ) THEN
         CALL xchalt('(get_pi_ph: Problem with netCDF1)')
                stop '(get_pi_ph: Problem with netCDF1)'
        END IF
       END IF
!
! Read  data
       call read_netcdf_var(ncid,'pH',pi_ph_in(1,1,1),12,0,0)

!
! Close file
       IF(mnproc==1) THEN
        ncstat = NF90_CLOSE(ncid)
        IF ( ncstat .NE. NF90_NOERR ) THEN
         CALL xchalt('(get_pi_ph: Problem with netCDF200)')
                stop '(get_pi_ph: Problem with netCDF200)'
        END IF
       END IF

! set missings over land
      do l=1,12
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j).gt.0.5) then 
              pi_ph(i,j,l) = pi_ph_in(i,j,l)
            else 
              pi_ph(i,j,l) = 0.
            endif    
          enddo
         enddo
      enddo

      RETURN
      END

      END MODULE mo_get_pi_ph
