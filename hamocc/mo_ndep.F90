module mo_ndep
!********************************************************************************
!
!     S.Gao             *Gfi, Bergen*    2017-08-19
!
! Modified
! --------
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - re-organised this module into an initialisation routine and
!       a routine that does the deposition, otherwise unchanged functionality
!
! Purpose
! -------
!  -Routines for reading and applying nitrogen deposition fluxes 
!
!
! Description:
! ------------
!  
!  The routine n_deposition reads nitrogen deposition from file and applies it 
!  to the top-most model layer. 
!
!  The file name is optionally read from HAMOCC's bgcnml namelist into the
!  namelist variable ndepfname. The default name ndep_CMIP6_1850.nc is assumed if 
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
!  Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to 
!    alkalinity (minus 1 mol)
!  
!
!********************************************************************************

character(len=512), save :: ndepfname='ndep_CMIP6_1850.nc'
character(len=512), save :: ndepfile
integer,            save :: startyear,endyear  
real, allocatable,  save :: ndep(:,:)
logical,            save :: lini = .false.
logical,            save :: file_exists=.false.  

!********************************************************************************
contains



subroutine ini_ndep(kpie,kpje,path)
!********************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Set filename and check existence of n-deposition input file. This is a 
!   preliminary version, functionality of reading filename from namelist needs
!   to be added.
!
! Changes: 
! --------
!
!********************************************************************************
use mod_xc,         only: mnproc
use mo_control_bgc, only: io_stdo_bgc
use mod_dia,        only: iotype
use mod_nctools,    only: ncfopn,ncgeti,ncfcls

implicit none 

integer,          intent(in) :: kpie,kpje
character(len=*), intent(in) :: path


! check if nitrogen deposition specified. if not, exit routine 
if (.not. lini) then 
  if (len_trim(ndepfname).eq.0) then  
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'No N deposition file specified. Setting ndep to zero.' 
    end if 
  else
    ndepfile=trim(path)//trim(ndepfname)
    inquire(file=ndepfile,exist=file_exists)
    if (.not.file_exists.and.mnproc.eq.1) then 
      write(io_stdo_bgc,*) 'Cannot find N deposition file ' & 
      & //trim(ndepfile)//'. Setting ndep to zero.' 
    end if 
  end if 

  if (file_exists) then
    ! read start and end year of n-deposition file
    call ncfopn(trim(ndepfile),'r',' ',1,iotype)
    call ncgeti('startyear',startyear)
    call ncgeti('endyear',endyear)
    call ncfcls
    ! allocate the field to hold one month of n-deposition data
    allocate(ndep(kpie,kpje))
    ndep(:,:) = 0.0
    if (mnproc.eq.1) then 
      write(io_stdo_bgc,*) 'Using N deposition file '//trim(ndepfile) 
    end if
  endif

  lini=.true.
end if  


end



subroutine n_deposition(kpie,kpje,kpke,kplyear,kplmon,pddpo,omask)
!********************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Read n-deposition data and apply to the top-most model layer.
!
! Changes: 
! --------
!  Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to 
!    alkalinity (minus 1 mol)
!
!********************************************************************************
use mod_xc,         only: mnproc
use netcdf,         only: nf90_open,nf90_close,nf90_nowrite
use mo_control_bgc, only: io_stdo_bgc,dtb
use mo_param1_bgc,  only: iano3,ialkali
use mo_carbch,      only: ocetra

integer :: kpie,kpje,kpke,kplyear,kplmon,kplday 
real, dimension(kpie,kpje,kpke) :: pddpo
real, dimension(kpie,kpje) :: omask

! local variables 
integer :: i,j,month_in_file,ncstat,ncid
integer, save :: oldmonth=0 

if (.not.file_exists) return 

  
! read ndep data from file 
if (kplmon.ne.oldmonth) then 
  month_in_file=(max(startyear,min(endyear,kplyear))-startyear)*12+kplmon 
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) 'Read N deposition month ',month_in_file, & 
      & ' from file ',trim(ndepfile) 
  end if 
  ncstat=nf90_open(trim(ndepfile),nf90_nowrite,ncid)
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
