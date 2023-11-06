! Copyright (C) 2001  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
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
!     *integer*   *kcid*       - file ID.
!     *integer*   *kshort*     - length of short name.
!     *integer*   *kdims*      - number of dimensions.
!     *integer*   *kcdims*     - dimensions.
!     *integer*   *kcvarid*    - variable ID.
!     *integer*   *kunitl*     - length of unit string.
!     *integer*   *klong*      - length of long name.
!     *integer*   *klabel*     - label for abort identification.
!     *integer*   *kunit*      - stdout unit.
!     *real*      *pmissing*   - missing value.
!     *CHARACTER* *yshort*     - short name.
!     *CHARACTER* *yunit*      - unit string.
!     *CHARACTER* *ylong*      - long name.
! 
! 
!     Externals
!     ---------
!     NONE.
! 
! **************************************************************************
      use netcdf,  only: nf90_double,nf90_noerr,nf90_put_att,nf90_def_var 
      use mod_xc,  only: mnproc,xchalt 
      use mod_dia, only:iotype
      IMPLICIT NONE
#ifdef PNETCDF
#include <pnetcdf.inc>
#include <mpif.h>      
#endif
       
      integer ncstat

      integer kcid,kcvarid,kdims,kcdims(kdims)                          &
     &       ,kunitl,klong,kshort,klabel,kunit,k

      real pmissing

      CHARACTER*(*) yshort, yunit, ylong

      CHARACTER*24 ystring
#ifdef PNETCDF
      integer(kind=MPI_OFFSET_KIND) clen
#endif
      ystring(1:21)='NETCDF stop at label '


! 
!  Define variable
! 
      if (mnproc==1 .AND. IOTYPE==0) then
      ncstat =                                                          &
     &NF90_DEF_VAR(kcid,yshort(1:kshort),NF90_DOUBLE,kcdims,kcvarid)  
      if ( ncstat .NE. NF90_NOERR ) then
         write(kunit,*) 'Problems with definition of NetCDF variable:'
         write(kunit,*) 'kcid           : ',kcid
         write(kunit,*) 'kshort         : ',kshort
         write(kunit,*) 'yshort(kshort) : ',yshort(1:kshort),'---'
         write(kunit,*) 'kdims          : ',kdims
         write(kunit,*) 'kcdims         : ',(kcdims(k),k=1,kdims)
         write(kunit,*) 'kcvarid        : ',kcvarid     
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      endif
! 
!  Set unit
! 
      ncstat =                                                          &
     &NF90_PUT_ATT(kcid,kcvarid,'units',yunit(1:kunitl))
      if ( ncstat .NE. NF90_NOERR ) then
         write(kunit,*) 'Problems with definition of unit:'
         write(kunit,*) 'kcid          : ',kcid
         write(kunit,*) 'kcvarid       : ',kcvarid     
         write(kunit,*) 'kunitl        : ',kunitl
         write(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      endif

! 
!  Set long name
! 
      ncstat =                                                          &
     &NF90_PUT_ATT(kcid,kcvarid,'long_name',ylong(1:klong))
      if ( ncstat .NE. NF90_NOERR ) then
         write(kunit,*) 'Problems with definition of long name:'
         write(kunit,*) 'kcid         : ',kcid
         write(kunit,*) 'kcvarid      : ',kcvarid     
         write(kunit,*) 'klong        : ',klong
         write(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      endif

! 
!  Set missing value
! 

      ncstat = NF90_PUT_ATT                                        &
     &(kcid,kcvarid,'missing_value',pmissing)
      if ( ncstat .NE. NF90_NOERR ) then
         write(kunit,*) 'Problems with definition of missing value:'
         write(kunit,*) 'kcid     : ',kcid
         write(kunit,*) 'kcvarid  : ',kcvarid     
         write(kunit,*) 'pmissing : ',pmissing
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(netcdf_def_vardb)')
                stop '(netcdf_def_vardb)'
      endif
      else if (IOTYPE==1) then
#ifdef PNETCDF
! 
!  Define variable
! 
      ncstat = nfmpi_def_var(kcid,yshort(1:kshort),nf_double,kdims,  &
     & kcdims,kcvarid)

      if ( ncstat .NE. NF_NOERR ) then
         write(kunit,*) 'Problems with definition of NetCDF variable:'
         write(kunit,*) 'kcid           : ',kcid
         write(kunit,*) 'kshort         : ',kshort
         write(kunit,*) 'yshort(kshort) : ',yshort(1:kshort),'---'
         write(kunit,*) 'kdims          : ',kdims
         write(kunit,*) 'kcdims         : ',(kcdims(k),k=1,kdims)
         write(kunit,*) 'kcvarid        : ',kcvarid
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      endif
! 
!  Set unit
!    
      clen=len(trim(yunit(1:kunitl)))
      ncstat =                                                          &
     &NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'units',clen,yunit(1:kunitl))
      if ( ncstat .NE. NF_NOERR ) then
         write(kunit,*) 'Problems with definition of unit:'
         write(kunit,*) 'kcid          : ',kcid
         write(kunit,*) 'kcvarid       : ',kcvarid
         write(kunit,*) 'kunitl        : ',kunitl
         write(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      endif

! 
!  Set long name
!  
      clen=len(trim(ylong(1:klong)))
      ncstat =                                                          &
     &NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'long_name',clen,ylong(1:klong))
      if ( ncstat .NE. NF_NOERR ) then
         write(kunit,*) 'Problems with definition of long name:'
         write(kunit,*) 'kcid         : ',kcid
         write(kunit,*) 'kcvarid      : ',kcvarid
         write(kunit,*) 'klong        : ',klong
         write(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      endif

! 
!  Set missing value
! 
      clen=1
      ncstat = NFMPI_PUT_ATT_DOUBLE                                        &
     &(kcid,kcvarid,'missing_value',NF_DOUBLE,clen,pmissing)
      if ( ncstat .NE. NF_NOERR ) then
         write(kunit,*) 'Problems with definition of missing value:'
         write(kunit,*) 'kcid     : ',kcid
         write(kunit,*) 'kcvarid  : ',kcvarid
         write(kunit,*) 'pmissing : ',pmissing
         write(ystring(22:24),'(I3)') klabel
         write(kunit,*) ystring
         call xchalt('(pnetcdf_def_vardb)')
                stop '(pnetcdf_def_vardb)'
      endif

#endif      
      endif
      return
      END
