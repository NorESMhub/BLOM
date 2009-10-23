      SUBROUTINE GET_DUST(kpie,kpje,kpke,omask,path,path_len)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/get_dust.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\
!$Name:  $\\
!
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
!
!     *CALL*       *GET_DUST*
!
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      use mo_param1_bgc 

      USE mod_xc 

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l
  
      REAL ::omask(kpie,kpje)
      character*(*) path
      integer path_len

! define the fields

      REAL :: dustin(kpie,kpje,12)

#ifdef PNETCDF
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid
!
! Open netCDF data file
!      
       IF(mnproc==1) THEN
        ncstat = NF_OPEN(path(1:path_len)//'INPDUST_mhw.nc',   &
     &                   NF_NOWRITE, ncid)
        IF (ncstat.NE.NF_NOERR ) THEN
         CALL xchalt('(get_dust: Problem with netCDF1)')
                stop '(get_dust: Problem with netCDF1)'
        END IF
       END IF
!
! Read  data
       call read_netcdf_var(ncid,'DUST',dustin(1,1,1),12,0)

!
! Close file
       IF(mnproc==1) THEN
        ncstat = NF_CLOSE(ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
         CALL xchalt('(get_dust: Problem with netCDF200)')
                stop '(get_dust: Problem with netCDF200)'
        END IF
       END IF
#else 
     
! define logical io unit (is this save, so that it doesn't override the
! io definitions in other parts of the program?)      
      INTEGER :: io_inbgc5

      io_inbgc5=175

      IF(p_pe == p_io) THEN
       open(io_inbgc5,file='INPDUST',status='unknown',    &
     &        access='sequential',form='unformatted')
      ENDIF

      do l=1,12
       call read_slice(io_inbgc5,dustin(:,:,l))
      enddo
      
      IF(p_pe == p_io) close(io_inbgc5)

#endif /*PNETCDF*/

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

!      IF(65-p_ioff>=1 .AND. 65-p_ioff<=kpie .AND. &
!     &   35-p_joff>=1 .AND. 35-p_joff<=kpje) THEN
!        write(io_stdo_bgc,*) 'DUST input at NABE'
!        write(io_stdo_bgc,*) 'i= 65 (47.34N) j=35 (19.29N)'
!        DO l=1,12
!          write(io_stdo_bgc,*) dustin(65-p_ioff,35-p_joff,l)
!        ENDDO
!      ENDIF

      RETURN
      END
