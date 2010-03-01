      SUBROUTINE NETCDF_DEF_VARDB                                       &
     & (kcid,kshort,yshort,kdims,kcdims,kcvarid,                        &
     &  kunitl,yunit,klong,ylong,pmissing,klabel,kunit)
! ****************************************************************
! 
! **** *NETCDF_DEF_VAR* - define NetCDF variable.
! 
!     S.Legutke,        *MPI-MaD, HH*    10.10.01
! 
!     Modified
!     --------
! 
!     Purpose
!     -------
!     Interface to NETCDF routines.
! 
!     Method
!     -------
!     
! 
!**   Interface.
!     ----------
! 
!     *CALL*       *NETCDF_DEF_VARDB(kcid,kshort,yshort,kdims,kcdims,kcvarid,
!                              kunitl,yunit,klong,ylong,pmissing,klabel,kunit)*
! 
! 
! **   Interface to calling routine (parameter list):
!     ----------------------------------------------
! 
!     *INTEGER*   *kcid*       - file ID.
!     *INTEGER*   *kshort*     - length of short name.
!     *INTEGER*   *kdims*      - number of dimensions.
!     *INTEGER*   *kcdims*     - dimensions.
!     *INTEGER*   *kcvarid*    - variable ID.
!     *INTEGER*   *kunitl*     - length of unit string.
!     *INTEGER*   *klong*      - length of long name.
!     *INTEGER*   *klabel*     - label for abort identification.
!     *INTEGER*   *kunit*      - stdout unit.
!     *REAL*      *pmissing*   - missing value.
!     *CHARACTER* *yshort*     - short name.
!     *CHARACTER* *yunit*      - unit string.
!     *CHARACTER* *ylong*      - long name.
! 
! 
!     Externals
!     ---------
!     none.
! 
! **************************************************************************

      USE mod_xc 

      INCLUDE 'netcdf.inc'

      INTEGER ncstat

      INTEGER kcid,kcvarid,kdims,kcdims(kdims)                          &
     &       ,kunitl,klong,kshort,klabel,kunit

      REAL pmissing

      CHARACTER*(*) yshort, yunit, ylong

      CHARACTER*24 ystring

      ystring(1:21)='NETCDF stop at label '
! 
!  Define variable
! 
      ncstat =                                                          &
     &NF_DEF_VAR(kcid,yshort(1:kshort),NF_DOUBLE,kdims,kcdims,kcvarid)  
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of NetCDF variable:'
         WRITE(kunit,*) 'kcid           : ',kcid
         WRITE(kunit,*) 'kshort         : ',kshort
         WRITE(kunit,*) 'yshort(kshort) : ',yshort(1:kshort),'---'
         WRITE(kunit,*) 'kdims          : ',kdims
         WRITE(kunit,*) 'kcdims         : ',(kcdims(k),k=1,kdims)
         WRITE(kunit,*) 'kcvarid        : ',kcvarid     
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      ENDIF
! 
!  Set unit
! 
      ncstat =                                                          &
     &NF_PUT_ATT_TEXT(kcid,kcvarid,'units',kunitl,yunit(1:kunitl))
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of unit:'
         WRITE(kunit,*) 'kcid          : ',kcid
         WRITE(kunit,*) 'kcvarid       : ',kcvarid     
         WRITE(kunit,*) 'kunitl        : ',kunitl
         WRITE(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      ENDIF

! 
!  Set long name
! 
      ncstat =                                                          &
     &NF_PUT_ATT_TEXT(kcid,kcvarid,'long_name',klong,ylong(1:klong))
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of long name:'
         WRITE(kunit,*) 'kcid         : ',kcid
         WRITE(kunit,*) 'kcvarid      : ',kcvarid     
         WRITE(kunit,*) 'klong        : ',klong
         WRITE(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      ENDIF

! 
!  Set missing value
! 

      ncstat = NF_PUT_ATT_DOUBLE                                        &
     &(kcid,kcvarid,'missing_value',NF_DOUBLE,1,pmissing)
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of missing value:'
         WRITE(kunit,*) 'kcid     : ',kcid
         WRITE(kunit,*) 'kcvarid  : ',kcvarid     
         WRITE(kunit,*) 'pmissing : ',pmissing
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      ENDIF



      RETURN
      END
