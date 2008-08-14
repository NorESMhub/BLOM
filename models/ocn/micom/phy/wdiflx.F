      subroutine wdiflx
c
c --- Write accumulated diagnosed heat and salt fluxes
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
c
      integer i2max
      parameter(i2max=32767)
c
      real*8 fldmax,fldmin,scale_factor,offset
      character fname*18
      integer i,j,k
      integer*2 fldout(idm,jdm)
c
      if (ditflx) then
        write (fname,'(2a,i4.4,a)')
     .    runid(1:runid_len),'_tflxdi_',nyear-1,'.uf'
        write (lp,'(2a)') 'Writing diagnosed heat flux to ',
     .    path1(1:path1_len)//fname
        open (unit=nfu,file=path1(1:path1_len)//fname,
     .        form='unformatted')
        do k=1,48
          if (nflxdi(k).gt.0) then
            fldmin=tflxdi(1,1,k)
            fldmax=tflxdi(1,1,k)
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
            do j=1,jj
              do i=1,ii
                fldmin=min(fldmin,tflxdi(i,j,k))
                fldmax=max(fldmax,tflxdi(i,j,k))
              enddo
            enddo
c$OMP END PARALLEL DO
            scale_factor=(fldmax-fldmin)/real(2*i2max)
            offset=.5*(fldmax+fldmin)
c$OMP PARALLEL DO
            do j=1,jj
              do i=1,ii
                fldout(i,j)=nint((tflxdi(i,j,k)-offset)/scale_factor)
              enddo
            enddo
c$OMP END PARALLEL DO
            i=nflxdi(k)
            write (nfu) k,i,scale_factor,offset,fldout
          endif
        enddo
        close (unit=nfu)
c$OMP PARALLEL DO
        do k=1,48
          nflxdi(k)=0
          do j=1,jj
            do i=1,ii
              tflxdi(i,j,k)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (disflx) then
        write (fname,'(2a,i4.4,a)')
     .    runid(1:runid_len),'_sflxdi_',nyear-1,'.uf'
        write (lp,'(2a)') 'Writing diagnosed salt flux to ',
     .    path1(1:path1_len)//fname
        open (unit=nfu,file=path1(1:path1_len)//fname,
     .        form='unformatted')
        do k=1,48
          if (nflxdi(k).gt.0) then
            fldmin=sflxdi(1,1,k)
            fldmax=sflxdi(1,1,k)
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
            do j=1,jj
              do i=1,ii
                fldmin=min(fldmin,sflxdi(i,j,k))
                fldmax=max(fldmax,sflxdi(i,j,k))
              enddo
            enddo
c$OMP END PARALLEL DO
            scale_factor=(fldmax-fldmin)/real(2*i2max)
            offset=.5*(fldmax+fldmin)
c$OMP PARALLEL DO
            do j=1,jj
              do i=1,ii
                fldout(i,j)=nint((sflxdi(i,j,k)-offset)/scale_factor)
              enddo
            enddo
c$OMP END PARALLEL DO
            i=nflxdi(k)
            write (nfu) k,i,scale_factor,offset,fldout
          endif
        enddo
        close (unit=nfu)
c$OMP PARALLEL DO
        do k=1,48
          nflxdi(k)=0
          do j=1,jj
            do i=1,ii
              sflxdi(i,j,k)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      return
      end
