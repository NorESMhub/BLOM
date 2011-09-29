subroutine profile_gd(kpie,kpje,kpke,pgila,pgiph,ptiestu,omask,path)

!********************************************************************************
!     J.Schwinger,        *Gfi, Bergen*    19.05.2011
!
!     Modified
!     --------
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
   
use mo_carbch,     only: ocetra
use mo_Gdata_read, only: set_Gdata,clean_Gdata,get_profile,nz,zlev
use mo_param1_bgc


implicit none

integer,         intent(in) :: kpie,kpje,kpke
real,            intent(in) :: ptiestu(kpie,kpje,kpke+1)
real,            intent(in) :: omask(kpie,kpje)
real,            intent(in) :: pgila(kpie*2,kpje*2)
real,            intent(in) :: pgiph(kpie*2,kpje*2)
character(len=*),intent(in) :: path

! Local variables
integer :: i,j,k,l,n
integer :: idx,izmax
real    :: prf(nz),zmax,clon,clat
real    :: d,d1,d2


! Extent of "smoothing region"
real,             parameter :: dxy = 10.0
! Number of data fields to read
integer,          parameter :: nflds = 6
! Names of data fields to read
character(len=3), parameter :: vname(nflds) = (/ 'dic', 'alk', 'pho', 'nit', 'sil', 'oxy'/)
! Index of fields in ocetra
integer,          parameter :: ifld(nflds)   = (/ isco212,ialkali,iphosph,iano3,isilica,ioxygen/)


do n = 1, nflds  ! Loop over tracer

   call set_Gdata(path,vname(n),dxy)

   do j=1,kpje
      do i=1,kpie
         
         If(omask(i,j) > 0.5) THEN
            
            clon = pgila(i*2,j*2)
            clat = pgiph(i*2,j*2)
            idx  = ifld(n)
            call get_profile(clon,clat,prf)

            ! Find depest z-level with valid data
            izmax=nz
            do l=2,nz
               if( prf(l) < 0 ) then
                  izmax = l-1
                  exit
               end if
            end do
            zmax = zlev(izmax)
            
            ! Surface layer
            ocetra(i,j,1,idx)=prf( 1)*1.e-6
            
            Do k=2,kpke
               
               IF(ptiestu(i,j,k).ge. zmax) THEN
                  
                  ocetra(i,j,k,idx)=prf(izmax)*1.e-6
                  
               ELSE 
                  
                  do l=1,izmax-1
                     
                     if(zlev(l  ).le.ptiestu(i,j,k).and.                 &
                        zlev(l+1).ge.ptiestu(i,j,k)) then
                        
                        d =zlev(l+1)-zlev(l)
                        d1=ptiestu(i,j,k)-zlev(l)
                        d2=zlev(l+1)-ptiestu(i,j,k)
                        
                        ocetra(i,j,k,idx)=((d2*prf(l)+d1*prf(l+1))/d)*1.e-6
                        
                     endif

                  enddo
                  
               ENDIF
               
            ENDDO
            
         ENDIF ! omask > 0.5

      ENDDO 
   ENDDO

   call clean_Gdata()

end do  ! Loop over fields

RETURN

!********************************************************************************
END subroutine profile_gd
