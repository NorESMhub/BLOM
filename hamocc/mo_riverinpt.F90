module mo_riverinpt
!********************************************************************************
!     C. Bernard
!     J. Schwinger,        *Gfi, Bergen*    28.09.2012
!
! Modified
! --------
!     
! Purpose
! -------
!  - Routines for reading and applying riverine nutrient and carbon input data
!
!
! Description:
! ------------
!  Public routines and variable of this module:
!
!  -read_gnews
!    read gnews riverine nutrient and cabon data 
!
!  -subroutine riverinpt
!    apply riverine input to the ocean tracer fields
!
!********************************************************************************
use dimensions,     only: itdm,jtdm,idm,jdm
use mod_xc ,        only: mnproc,nbdy,xcaput,xchalt
use mo_control_bgc, only: io_stdo_bgc
implicit none

private
public :: riverinpt,ini_riverinpt

logical, save              :: lini=.false. 
integer,         parameter :: nrivmax=10000
character(len=*),parameter :: infile = 'riv_in.txt'


! Global "Tuning factor" for riverine input (GNEWS input data
! is multiplied by this factor). Test with the low resolution
! NorESM version showed that an increase of riverine nutrient
! input by ~10% reduces model drift.
real,            parameter :: tfac = 1.1


! arrays for riverine inputs on the model grid
real,save,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIN2d, riv_DIP2d,   &
                                                        riv_DSI2d, riv_DIC2d,   &
                                                        riv_idet2d,riv_idoc2d,  &
                                                        riv_DFe2d


!********************************************************************************
contains


subroutine ini_riverinpt(path)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!   Read river input and compute update
!   all Global_NEWS nutrients given in Mg/yr,
!   Fe in MG/yr(1.45 Tg from Chester 1990)
!
! Description:
! ------------
! 
!
! Arguments:
! ----------
!
!--------------------------------------------------------------------------------
use mo_biomod, only: rcar,rnit

implicit none

character(len=*)            :: path

integer                     :: stat,nriv,riv_indice,i
integer, dimension(nrivmax) :: riv_i,riv_j
real                        :: distmin,lon_gn,lat_gn,lon_noresm,lat_noresm
real                        :: res_DOP, res_DON, res_DOC
real                        :: res_PP,  res_PN,  res_POC
real,    dimension(nrivmax) :: riv_DIN, riv_DIP, riv_DON, riv_DOP, riv_DOC,     &
                               riv_PN,  riv_PP,  riv_POC, riv_DIC, riv_DSi,     &
                               riv_idet,riv_idoc,riv_DFe
       

! Temporary global fields 
real,allocatable,dimension(:,:) :: rivt_DIN2d,rivt_DIP2d, rivt_DSI2d,           &
                                   rivt_DIC2d,rivt_idet2d,rivt_idoc2d,          &
                                   rivt_DFe2d
   

if( lini ) then 
  write(io_stdo_bgc,*) 'Warning: Riverine nutrients data already initialised'
  return
end if

! allocate temporary global fields
allocate(  rivt_DIN2d(itdm,jtdm),  &
           rivt_DIP2d(itdm,jtdm),  &
           rivt_DSI2d(itdm,jtdm),  &
           rivt_DIC2d(itdm,jtdm),  &
           rivt_idet2d(itdm,jtdm), &
           rivt_idoc2d(itdm,jtdm), &
           rivt_DFe2d(itdm,jtdm) )

rivt_DIN2d(:,:)  = 0.0
rivt_DIP2d(:,:)  = 0.0
rivt_DSI2d(:,:)  = 0.0
rivt_DIC2d(:,:)  = 0.0
rivt_idet2d(:,:) = 0.0
rivt_idoc2d(:,:) = 0.0
rivt_DFe2d(:,:)  = 0.0


! Read and process the data on proc 1
if (mnproc.eq.1) then

  write(io_stdo_bgc,*) 'Read GNEWS riverine nutrients inputs data'
  open(99,file=trim(path)//infile)
  read(99,*)
  read(99,*)

  nriv=0
  read_data: do
    
    read(99,*,iostat=stat) riv_indice,      &
                           riv_j(nriv+1),   &
                           riv_i(nriv+1),   &
                           distmin,         &
                           Lon_GN,          &
                           Lat_GN,          &
                           Lon_NorESM,      &
                           Lat_NorESM,      &
                           riv_DIN(nriv+1), &
                           riv_DIP(nriv+1), &
                           riv_DON(nriv+1), &
                           riv_DOP(nriv+1), &
                           riv_DOC(nriv+1), &
                           riv_PN(nriv+1),  &
                           riv_PP(nriv+1),  &
                           riv_POC(nriv+1), &
                           riv_DIC(nriv+1), &
                           riv_DSi(nriv+1), &
                           riv_DFe(nriv+1)

    if (stat.ge.0) then
      nriv=nriv+1
      riv_DIN(nriv) = riv_DIN(nriv)*1E6/14. /1E3*tfac 
      riv_DIP(nriv) = riv_DIP(nriv)*1E6/31. /1E3*tfac
      riv_DON(nriv) = riv_DON(nriv)*1E6/14. /1E3*tfac
      riv_DOP(nriv) = riv_DOP(nriv)*1E6/31. /1E3*tfac
      riv_DOC(nriv) = riv_DOC(nriv)*1E6/12. /1E3*tfac
      riv_PN(nriv)  = riv_PN(nriv) *1E6/14. /1E3*tfac
      riv_PP(nriv)  = riv_PP(nriv) *1E6/31. /1E3*tfac
      riv_POC(nriv) = riv_POC(nriv)*1E6/12. /1E3*tfac
      riv_DIC(nriv) = riv_DIC(nriv)*1E6/12. /1E3*tfac  
      riv_DSi(nriv) = riv_DSi(nriv)*1E6/28. /1E3*tfac
      riv_DFe(nriv) = riv_DFe(nriv)*1E6/55.8/1E3*tfac
    else
      exit read_data
    endif
    
  end do read_data
  
  close(99)


  do i=1,nriv
    
    ! Inorganic tracers
    rivt_DIN2d(riv_i(i),riv_j(i))=                         &
    rivt_DIN2d(riv_i(i),riv_j(i))+riv_DIN(i)
    rivt_DIP2d(riv_i(i),riv_j(i))=                         &
    rivt_DIP2d(riv_i(i),riv_j(i))+riv_DIP(i)
    rivt_DSI2d(riv_i(i),riv_j(i))=                         &
    rivt_DSI2d(riv_i(i),riv_j(i))+riv_DSI(i)
    rivt_DIC2d(riv_i(i),riv_j(i))=                         &
    rivt_DIC2d(riv_i(i),riv_j(i))+riv_DIC(i)
    rivt_DFe2d(riv_i(i),riv_j(i))=                         &
    rivt_DFe2d(riv_i(i),riv_j(i))+riv_DFe(i)

    ! Dissolved organic matter is assumed to be in Redfield ratio in the
    ! model. For this reason DOC is assembled according to Redfield here,
    ! and excess in any of the components P, N or C is put into the 
    ! inorganic part
    riv_idoc(i) = min( riv_DOP(i), riv_DON(i)/rnit, riv_DOC(i)/rcar )   
    res_DOP = max( riv_DOP(i)      - riv_idoc(i) , 0.0 ) 
    res_DON = max( riv_DON(i)/rnit - riv_idoc(i) , 0.0 ) 
    res_DOC = max( riv_DOC(i)/rcar - riv_idoc(i) , 0.0 ) 

    rivt_idoc2d(riv_i(i),riv_j(i))=                          &
    rivt_idoc2d(riv_i(i),riv_j(i))+riv_idoc(i)

    ! Add excess dissolved organic matter to the inorganic part
    rivt_DIP2d(riv_i(i),riv_j(i))=                             &
    rivt_DIP2d(riv_i(i),riv_j(i))+res_DOP
    rivt_DIN2d(riv_i(i),riv_j(i))=                             &
    rivt_DIN2d(riv_i(i),riv_j(i))+res_DON*rnit
    rivt_DIC2d(riv_i(i),riv_j(i))=                             &
    rivt_DIC2d(riv_i(i),riv_j(i))+res_DOC*rcar


    ! Particulate organic matter is assumed to be in Redfield ratio in the
    ! model. For this reason POC is assembled according to Redfield here,
    ! and excess in any of the components P, N or C is put into the 
    ! inorganic part
    riv_idet(i) = min(riv_PP(i), riv_PN(i)/rnit, riv_POC(i)/rcar)
    res_PP  = max( riv_PP(i)       - riv_idet(i) , 0.0 ) 
    res_PN  = max( riv_PN(i)/rnit  - riv_idet(i) , 0.0 )  
    res_POC = max( riv_POC(i)/rcar - riv_idet(i) , 0.0 )  


    rivt_idet2d(riv_i(i),riv_j(i))=                          &
    rivt_idet2d(riv_i(i),riv_j(i))+riv_idet(i)

    ! Add excess particulate organic matter to the inorganic part
    rivt_DIP2d(riv_i(i),riv_j(i))=                             &
    rivt_DIP2d(riv_i(i),riv_j(i))+res_PP
    rivt_DIN2d(riv_i(i),riv_j(i))=                             &
    rivt_DIN2d(riv_i(i),riv_j(i))+res_PN*rnit
    rivt_DIC2d(riv_i(i),riv_j(i))=                             &
    rivt_DIC2d(riv_i(i),riv_j(i))+res_POC*rcar

  end do


endif ! mnproc == 1


! Send data to all procs
call xcaput(rivt_DIN2d, riv_DIN2d,1)
call xcaput(rivt_DIP2d, riv_DIP2d,1)
call xcaput(rivt_DSI2d, riv_DSI2d,1)
call xcaput(rivt_DIC2d, riv_DIC2d,1)
call xcaput(rivt_idoc2d,riv_idoc2d,1)
call xcaput(rivt_idet2d,riv_idet2d,1)
call xcaput(rivt_DFe2d, riv_DFe2d,1)



deallocate( rivt_DIN2d, rivt_DIP2d, rivt_DSI2d,  &
            rivt_DIC2d, rivt_idet2d,rivt_idoc2d, &
            rivt_DFe2d )

lini = .true.

!--------------------------------------------------------------------------------
end subroutine ini_riverinpt


subroutine riverinpt(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,omask)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Apply riverine input to oceanic tracer fields
!
! Description:
! ------------
! 
!
! Arguments:
! ----------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *omask*   - ocean mask
!
!--------------------------------------------------------------------------------
use mo_param1_bgc,  only: iano3,iphosph,isilica,isco212,iiron,idoc,idet
use mo_control_bgc, only: dtb
use mo_carbch,      only: ocetra

implicit none


integer,intent(in) :: kpie,kpje,kpke
real,   intent(in) :: pddpo(kpie,kpje,kpke)
real,   intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
real,   intent(in) :: omask(kpie,kpje)

! local variables
integer            :: i,j
real               :: fdt,volij

! initialisation is now called from ini_hamocc (changed by Ingo Bethke on 2013.06.07)
if( .not. lini ) then
  write(io_stdo_bgc,*) 'Error: Riverine nutrients data not initialised'
  call xchalt('(mo_riverinpt)')
  stop '(mo_riverinpt)'
end if


!$OMP PARALLEL DO PRIVATE(fdt,volij)
DO j=1,kpje
DO i=1,kpie
  IF(omask(i,j).GT.0.5) THEN

    fdt   = dtb/365.

    ! Distribute riverine inputs over the top two model layers
    ! THIS MIGHT BE DEPENDENT ON FUTURE CHANGES IN MICOM!
    volij = (pddpo(i,j,1)+pddpo(i,j,2))*pdlxp(i,j)*pdlyp(i,j)

    ! Inorganic elements
    ocetra(i,j,1:2,iano3)   = ocetra(i,j,1:2,iano3)   + riv_DIN2d(i,j)*fdt/volij
    ocetra(i,j,1:2,iphosph) = ocetra(i,j,1:2,iphosph) + riv_DIP2d(i,j)*fdt/volij
    ocetra(i,j,1:2,isilica) = ocetra(i,j,1:2,isilica) + riv_DSI2d(i,j)*fdt/volij
    ocetra(i,j,1:2,isco212) = ocetra(i,j,1:2,isco212) + riv_DIC2d(i,j)*fdt/volij
    ocetra(i,j,1:2,iiron)   = ocetra(i,j,1:2,iiron)   + riv_DFe2d(i,j)*fdt/volij

    ! Dissolved organic matter
    ocetra(i,j,1:2,idoc)    = ocetra(i,j,1:2,idoc)    + riv_idoc2d(i,j)*fdt/volij

    ! Particulate organic matter
    ocetra(i,j,1:2,idet)    = ocetra(i,j,1:2,idet)    + riv_idet2d(i,j)*fdt/volij

  ENDIF
ENDDO
ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------------------
end subroutine riverinpt



!********************************************************************************
end module mo_riverinpt
