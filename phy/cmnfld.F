      subroutine cmnfld(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute fields that are used by several subsequent routines
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      integer i,j,l
c
c --- ------------------------------------------------------------------
c --- Update halos of various fields.
c --- ------------------------------------------------------------------
c
      call xctilr(temp, 1,2*kk, 3,3, halo_ps)
      call xctilr(saln, 1,2*kk, 3,3, halo_ps)
c     call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
c     call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kfpla(i,j,n)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 2,2, halo_ps)
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=-1,jj+2
        do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
          kfpla(i,j,n)=nint(util1(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute fields depending on selection of physics.
c --- ------------------------------------------------------------------
c
      if (edritp.eq.'large scale'.or.eitmth.eq.'gm') then
c
c --- - Compute filtered buoyancy frequency squared.
        call cmnfld_bfsqf(m,n,mm,nn,k1m,k1n)
      endif
c
      if (edritp.eq.'large scale'.or.eitmth.eq.'gm') then
c
c --- - Estimate slope of local neutral surface.
        call cmnfld_nslope(m,n,mm,nn,k1m,k1n)
      endif
c
      return
      end
