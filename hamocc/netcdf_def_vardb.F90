! Copyright (C) 2001  S. Legutke
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
      USE netcdf
      USE mod_xc 
      USE mod_dia, only:iotype
      implicit none
#ifdef PNETCDF
#include <pnetcdf.inc>
#endif
#include <mpif.h>      
       
      INTEGER ncstat

      INTEGER kcid,kcvarid,kdims,kcdims(kdims)                          &
     &       ,kunitl,klong,kshort,klabel,kunit,k

      REAL pmissing

      CHARACTER*(*) yshort, yunit, ylong

      CHARACTER*24 ystring
      integer(kind=MPI_OFFSET_KIND) clen
      ystring(1:21)='NETCDF stop at label '


! 
!  Define variable
! 
      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat =                                                          &
     &NF90_DEF_VAR(kcid,yshort(1:kshort),NF90_DOUBLE,kcdims,kcvarid)  
      IF ( ncstat .NE. NF90_NOERR ) THEN
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
     &NF90_PUT_ATT(kcid,kcvarid,'units',yunit(1:kunitl))
      IF ( ncstat .NE. NF90_NOERR ) THEN
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
     &NF90_PUT_ATT(kcid,kcvarid,'long_name',ylong(1:klong))
      IF ( ncstat .NE. NF90_NOERR ) THEN
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

      ncstat = NF90_PUT_ATT                                        &
     &(kcid,kcvarid,'missing_value',pmissing)
      IF ( ncstat .NE. NF90_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of missing value:'
         WRITE(kunit,*) 'kcid     : ',kcid
         WRITE(kunit,*) 'kcvarid  : ',kcvarid     
         WRITE(kunit,*) 'pmissing : ',pmissing
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      ENDIF
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
! 
!  Define variable
! 
      ncstat = nfmpi_def_var(kcid,yshort(1:kshort),nf_double,kdims,  &
     & kcdims,kcvarid)

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
         CALL xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      ENDIF
! 
!  Set unit
!    
      clen=len(trim(yunit(1:kunitl)))
      ncstat =                                                          &
     &NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'units',clen,yunit(1:kunitl))
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of unit:'
         WRITE(kunit,*) 'kcid          : ',kcid
         WRITE(kunit,*) 'kcvarid       : ',kcvarid
         WRITE(kunit,*) 'kunitl        : ',kunitl
         WRITE(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      ENDIF

! 
!  Set long name
!  
      clen=len(trim(ylong(1:klong)))
      ncstat =                                                          &
     &NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'long_name',clen,ylong(1:klong))
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of long name:'
         WRITE(kunit,*) 'kcid         : ',kcid
         WRITE(kunit,*) 'kcvarid      : ',kcvarid
         WRITE(kunit,*) 'klong        : ',klong
         WRITE(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      ENDIF

! 
!  Set missing value
! 
      clen=1
      ncstat = NFMPI_PUT_ATT_DOUBLE                                        &
     &(kcid,kcvarid,'missing_value',NF_DOUBLE,clen,pmissing)
      IF ( ncstat .NE. NF_NOERR ) THEN
         WRITE(kunit,*) 'Problems with definition of missing value:'
         WRITE(kunit,*) 'kcid     : ',kcid
         WRITE(kunit,*) 'kcvarid  : ',kcvarid
         WRITE(kunit,*) 'pmissing : ',pmissing
         WRITE(ystring(22:24),'(I3)') klabel
         WRITE(kunit,*) ystring
         CALL xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      ENDIF

#endif      
      ENDIF
      RETURN
      END
