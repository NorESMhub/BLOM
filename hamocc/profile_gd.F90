subroutine profile_gd(kpie,kpje,kpke,pglon,pglat,ptiestw,omask,path)

!********************************************************************************
!     J.Schwinger,        *Gfi, Bergen*           2011-05-19
!
!     Modified
!     --------
!     J.Schwinger,        *Uni Climate, BCCR*     2017-07-07
!     - moved conversion from mumol to mol to mod_gdata_read
!     - changed linear interpolation from data-levels to model levels to propper
!       mapping of data profile to model-levels
!
!     Purpose
!     -------
!     - initialise HAMOCC fields with gridded (1x1 deg) WOA and GLODAP 
!       data using the module mo_Gdata_read. Note that the routine get_profile
!       returns the mean of all data profiles within a rectangular region 
!       ("smoothing region") of dxy x dxy degrees extent, where dxy is an 
!       adjustable parameter.
!
!
!********************************************************************************
   
use mod_xc,         only: mnproc,xchalt
use mo_carbch,      only: ocetra
use mo_Gdata_read,  only: set_Gdata,clean_Gdata,get_profile,nzmax,nz,zlev_bnds
use mo_control_bgc, only: io_stdo_bgc
use mo_param1_bgc


implicit none

integer,         intent(in) :: kpie,kpje,kpke
real,            intent(in) :: ptiestw(kpie,kpje,kpke+1)
real,            intent(in) :: omask(kpie,kpje)
real,            intent(in) :: pglon(kpie,kpje)
real,            intent(in) :: pglat(kpie,kpje)
character(len=*),intent(in) :: path

! Local variables
integer         :: i,j,k,l,ll,n
integer         :: idx,izmax
real            :: prf(nzmax),wgt(nzmax),zbnds(2,nzmax),zmax,clon,clat
real, parameter :: fillval = -1.0e34

! Extent of "smoothing region"
real,             parameter :: dxy = 5.0

#ifdef natDIC
! Number of data fields to read
integer,          parameter :: nflds = 8
! Names of data fields to read
character(len=3), parameter :: vname(nflds) = (/ 'dic', 'alk', 'pho', 'nit', 'sil', 'oxy', 'dic', 'alk'/)
! Index of fields in ocetra
integer,          parameter :: ifld(nflds)   = (/ isco212,ialkali,iphosph,iano3,isilica,ioxygen,inatsco212,inatalkali/)
#else
! Number of data fields to read
integer,          parameter :: nflds = 6
! Names of data fields to read
character(len=3), parameter :: vname(nflds) = (/ 'dic', 'alk', 'pho', 'nit', 'sil', 'oxy'/)
! Index of fields in ocetra
integer,          parameter :: ifld(nflds)   = (/ isco212,ialkali,iphosph,iano3,isilica,ioxygen/)
#endif


do n = 1, nflds  ! Loop over tracer

   call set_Gdata(path,vname(n),dxy)

   do j=1,kpje
      do i=1,kpie

         If(omask(i,j) > 0.5) THEN
            
            clon = pglon(i,j)
            clat = pglat(i,j)
            idx  = ifld(n)
            call get_profile(clon,clat,prf)

            ! Find depest z-level with valid data
            izmax=nz
            do l=2,nz
               if( prf(l) < 0 ) then
                  izmax = l-1
                  exit
               endif
            enddo
            ! Set data level-boundaries for this profile
            zbnds         = fillval
            zbnds(:,1:nz) = zlev_bnds
            zbnds(1,1)    = 0.0                          ! make sure that upper data bnd is 0
            if(zbnds(2,izmax) < ptiestw(i,j,kpke+1)) then
               zbnds(2,izmax) = ptiestw(i,j,kpke+1)+10.0 ! extend lower bound of bottom layer
            endif
            
            Do k=1,kpke

               wgt(:)=0.0
                  
               loop_obs: do l=1,izmax

                  ! 1st case: Model layer completely within data-layer
                  if(zbnds(1,l) <= ptiestw(i,j,k) .and. zbnds(2,l) >=  ptiestw(i,j,k+1)) then  
                     ocetra(i,j,k,idx)=prf(l)
                     exit loop_obs
                  endif

                  ! 2nd case: one (or both) data-layer boundary are within model layer
                  
                  ! a) The lower data level-boundary is lower than the upper model level-interface.
                  !    and the upper data level-boundary is higher than the lower model 
                  !    level-interface => some overlap between data and model level exists.
                  !    Calculate the corresponding weight.
                  if(zbnds(2,l) > ptiestw(i,j,k) .and. zbnds(1,l) <= ptiestw(i,j,k+1))   &
                     wgt(l) =     zbnds(2,l)-ptiestw(i,j,k)                              &
                            - max(zbnds(1,l)-ptiestw(i,j,k),  0.0)                       &
                            - max(zbnds(2,l)-ptiestw(i,j,k+1),0.0)

                  ! a) The upper data level-boundary is lower than the lower model level-interface 
                  !    => all weights have been calculated, calculate concentration and exit
                  if(zbnds(1,l) > ptiestw(i,j,k+1) .or. l==izmax) then 
                     wgt(:) = wgt(:)/(ptiestw(i,j,k+1)-ptiestw(i,j,k))
                     if( abs(sum(wgt(:))-1.0) > 1.0e-6 ) then
                        write(io_stdo_bgc,*) 'profile_gd error: inconsisten weihts'
                        write(io_stdo_bgc,*) 'profile_gd error: ', k,l,abs(sum(wgt(:))-1.0)
                        write(io_stdo_bgc,*) 'profile_gd error: ', wgt(1:izmax)
                        write(io_stdo_bgc,*) 'profile_gd error: ', ptiestw(i,j,k),ptiestw(i,j,k+1)
                        call flush(io_stdo_bgc)
                        call xchalt('(profile_gd)')
                     endif
                     do ll=1,l
                        ocetra(i,j,k,idx) =  ocetra(i,j,k,idx) + prf(ll)*wgt(ll)
                     enddo
                     exit loop_obs
                  endif


               enddo loop_obs
               
            ENDDO ! k=1,kpke
            
         ENDIF ! omask > 0.5

      ENDDO 
   ENDDO

   call clean_Gdata()

end do  ! Loop over fields

RETURN

!********************************************************************************
END subroutine profile_gd
