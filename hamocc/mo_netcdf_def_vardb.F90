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

module mo_netcdf_def_vardb

  implicit none
  private

  public :: netcdf_def_vardb

contains

  subroutine netcdf_def_vardb (kcid,kshort,yshort,kdims,kcdims,kcvarid, &
                               kunitl,yunit,klong,ylong,pmissing,klabel,kunit)

    ! ****************************************************************
    !  Interface to NETCDF routines - define NetCDF variable.
    !
    !  S.Legutke,        *MPI-MaD, HH*    10.10.01
    ! **************************************************************************

    use netcdf,  only: nf90_double,nf90_noerr,nf90_put_att,nf90_def_var
    use mod_xc,  only: mnproc,xchalt
    use mod_dia, only: iotype
#ifdef PNETCDF
#include <pnetcdf.inc>
#include <mpif.h>
#endif

    ! Arguments
    integer,          intent(in)  :: kcid          ! file ID.
    integer,          intent(in)  :: kshort        ! length of short name.
    integer,          intent(in)  :: kdims         ! number of dimensions.
    integer,          intent(in)  :: kcdims(kdims) ! dimensions.
    integer,          intent(out) :: kcvarid       ! variable ID.
    integer,          intent(in)  :: kunitl        ! length of unit string.
    integer,          intent(in)  :: klong         ! length of long name.
    integer,          intent(in)  :: klabel        ! label for abort identification.
    integer,          intent(in)  :: kunit         ! stdout unit.
    character(len=*), intent(in)  :: yshort        ! short name.
    character(len=*), intent(in)  :: yunit         ! unit string.
    character(len=*), intent(in)  :: ylong         ! long name.

    ! Local variables
    integer           :: k
    real              :: pmissing
    character(len=24) :: ystring
    integer           ::ncstat
#ifdef PNETCDF
    integer(kind=MPI_OFFSET_KIND) :: clen
#endif

    ystring(1:21)='NETCDF stop at label '

    !
    !  Define variable
    !
    if (mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_DEF_VAR(kcid,yshort(1:kshort),NF90_DOUBLE,kcdims,kcvarid)
      if ( ncstat .NE. NF90_NOERR ) THEN
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
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'units',yunit(1:kunitl))
      if ( ncstat .NE. NF90_NOERR ) THEN
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
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'long_name',ylong(1:klong))
      if ( ncstat .NE. NF90_NOERR ) THEN
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
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'missing_value',pmissing)
      if ( ncstat .NE. NF90_NOERR ) THEN
        write(kunit,*) 'Problems with definition of missing value:'
        write(kunit,*) 'kcid     : ',kcid
        write(kunit,*) 'kcvarid  : ',kcvarid
        write(kunit,*) 'pmissing : ',pmissing
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(netcdf_def_vardb)')
        stop '(netcdf_def_vardb)'
      endif
    else if (IOTYPE==1) THEN
#ifdef PNETCDF
      !
      !  Define variable
      !
      ncstat = nfmpi_def_var(kcid,yshort(1:kshort),nf_double,kdims,kcdims,kcvarid)
      if ( ncstat .NE. NF_NOERR ) THEN
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
      ncstat = NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'units',clen,yunit(1:kunitl))
      if ( ncstat .NE. NF_NOERR ) THEN
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
      ncstat = NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'long_name',clen,ylong(1:klong))
      if ( ncstat .NE. NF_NOERR ) THEN
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
      ncstat = NFMPI_PUT_ATT_DOUBLE(kcid,kcvarid,'missing_value',NF_DOUBLE,clen,pmissing)
      if ( ncstat .NE. NF_NOERR ) THEN
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

  end subroutine netcdf_def_vardb

end module mo_netcdf_def_vardb
