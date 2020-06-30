! Copyright (C) 2020  J. Schwinger
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


subroutine profile_gd(kpie,kpje,kpke,kbnd,pglon,pglat,omask)
!*******************************************************************************
!     J.Schwinger,      *Gfi, Bergen*            2011-05-19
!
!     Modified
!     --------
!     J.Schwinger,      *Uni Climate, BCCR*      2017-07-07
!     - moved conversion from mumol to mol to mod_gdata_read
!     - changed linear interpolation from data-levels to model levels to propper
!       mapping of data profile to model-levels
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - adaptions for reading c-isotope initial values as d13C and d14C
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
!*******************************************************************************
   
use mod_xc,         only: xchalt
use mo_carbch,      only: ocetra
use mo_Gdata_read,  only: set_Gdata,clean_Gdata,get_profile,nzmax,nz,zlev_bnds,fillval
use mo_control_bgc, only: io_stdo_bgc
use mo_vgrid,       only: ptiestw
use mo_param1_bgc


implicit none

integer, intent(in) :: kpie,kpje,kpke,kbnd
real,    intent(in) :: omask(kpie,kpje)
real,    intent(in) :: pglon(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
real,    intent(in) :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)

! Local variables
integer             :: i,j,k,l,ll,n
integer             :: idx,izmax
real                :: prf(nzmax),wgt(nzmax),zbnds(2,nzmax),clon,clat

! Extent of "smoothing region"
real,     parameter :: dxy = 5.0

! Number of fields to read
integer,  parameter :: nread_base = 6
integer,  parameter :: nread_ndic = 2
integer,  parameter :: nread_ciso = 2
integer,  parameter :: maxflds    = nread_base+nread_ndic+nread_ciso

integer             :: nflds, no
integer             :: ifld(maxflds) 
character(len=3)    :: vname(maxflds)

nflds = nread_base
vname( 1:nflds) = (/ 'dic',  'alk',  'pho',  'nit','sil',  'oxy'  /)
 ifld( 1:nflds) = (/ isco212,ialkali,iphosph,iano3,isilica,ioxygen/)

#ifdef natDIC
no    = nflds+1
nflds = nflds+nread_ndic
vname(no:nflds) = (/'dic',     'alk'/)
 ifld(no:nflds) = (/inatsco212,inatalkali/)
#endif

#ifdef cisonew
no    = nflds+1
nflds = nflds+nread_ciso
vname(no:nflds) = (/'d13', 'd14'/)
 ifld(no:nflds) = (/isco213,isco214/)
#endif


do n = 1, nflds  ! Loop over tracer

   call set_Gdata(vname(n),dxy)

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
               if( prf(l) < fillval*0.1 ) then
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

                  ! b) The upper data level-boundary is lower than the lower model level-interface 
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
