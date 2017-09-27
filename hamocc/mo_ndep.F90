module mo_ndep
!********************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  - Routines for reading and applying nitrogen deposition fluxes 
!
!
! Description:
! ------------
!  
!  The routine n_deposition reads nitrogen deposition from file and applies it 
!  to the top-most model layer. 
!
!  The file name is optionally read from HAMOCC's bgcnml namelist into the
!  namelist variable ndepfile. The default name ndep_CMIP6_1850.nc is assumed if 
!  ndepfile is not specified in bgcnml. 
! 
!  The deposition flux is set to zero if the file is not found or if the namelist 
!  variable ndepfile is set to empty (i.e. ' ').  
! 
!  The input data must be already pre-interpolated to the ocean grid and stored 
!  in the same folder with MICOM's grid information.   
!
!
! Changes: 
! --------
!     Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to alkalinity (minus 1 mol)
!  
!
!********************************************************************************

character(len=512), save :: ndepfile='ndep_CMIP6_1850.nc',ndepdir=' '  

!********************************************************************************
contains

subroutine n_deposition(kpie,kpje,kpke,kplyear,kplmon,pddpo,omask)

use netcdf,         only: nf90_open,nf90_close,nf90_nowrite
use dimensions,     only: idm,jdm 
use mod_dia,        only: iotype
use mod_xc,         only: mnproc
use mod_nctools,    only: ncfopn,ncgeti,ncfcls
use mo_control_bgc, only: io_stdo_bgc,dtb
use mo_param1_bgc,  only: iano3,ialkali
use mo_carbch,      only: ocetra

implicit none 

integer :: kpie,kpje,kpke,kplyear,kplmon,kplday 
real, dimension(kpie,kpje,kpke) :: pddpo
real, dimension(kpie,kpje) :: omask

! local variables 

integer :: i,j,month_in_file,startyear,endyear,dummymask(2),ncstat,ncid
integer, save :: oldmonth=0 
logical, save :: first=.true.,file_exists=.false.  
real, save, dimension(idm,jdm) :: ndep 

! check if nitrogen deposition specified. if not, exit routine 
if (first) then 
  if (len_trim(ndepfile).eq.0) then  
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'No N deposition file specified. Setting ndep to zero.' 
    end if 
  else      
    inquire(file=trim(ndepdir)//trim(ndepfile),exist=file_exists)
    if (.not.file_exists.and.mnproc.eq.1) then 
      write(io_stdo_bgc,*) 'Cannot find N deposition file ' & 
      & //trim(ndepdir)//trim(ndepfile)//'. Setting ndep to zero.' 
    end if 
  end if 
  first=.false.
end if  
if (.not.file_exists) return 
  
! read ndep data from file 
if (kplmon.ne.oldmonth) then 
  if (oldmonth.eq.0) then 
    call ncfopn(trim(ndepdir)//trim(ndepfile),'r',' ',1,iotype)
    call ncgeti('startyear',startyear)
    call ncgeti('endyear',endyear)
    call ncfcls
  end if 
  month_in_file=(max(startyear,min(endyear,kplyear))-startyear)*12+kplmon 
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) 'Read N deposition month ',month_in_file, & 
      & ' from file ',trim(ndepdir)//trim(ndepfile) 
  end if 
  ncstat=nf90_open(trim(ndepdir)//trim(ndepfile),nf90_nowrite,ncid)
  call read_netcdf_var(ncid,'ndep',ndep,1,month_in_file,0) 
  ncstat=nf90_close(ncid)   
  oldmonth=kplmon 
end if 

! deposite N in topmost layer 
do j=1,kpje
  do i=1,kpie
    if (omask(i,j).gt.0.5) then
      ocetra(i,j,1,iano3)=ocetra(i,j,1,iano3)+ndep(i,j)*dtb/365./pddpo(i,j,1)
      ocetra(i,j,1,ialkali)=ocetra(i,j,1,ialkali)-ndep(i,j)*dtb/365./pddpo(i,j,1)
   end if
  enddo
enddo

end subroutine n_deposition 

!********************************************************************************
end module mo_ndep 
