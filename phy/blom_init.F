! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen, Mehmet Ilicak
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

      subroutine blom_init
c
c --- ------------------------------------------------------------------
c --- initialize the model
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_swabs, only: iniswa
      use mod_dia
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_bud.h"
#include "common_bud.h"
c
      integer istat,ncid,varid,i,j,k,l,m,n,mm,km
      real q
      logical icrest,fexist
c
      real gettime
      external gettime
c
c --- ------------------------------------------------------------------
c --- Initialize SPMD processing
c --- ------------------------------------------------------------------
c
      call xcspmd
c
c --- ------------------------------------------------------------------
c --- Initialize timing
c --- ------------------------------------------------------------------
c
      call initime
c
c --- print seconds elapsed since startup (should be almost zero)
      if (mnproc.eq.1) then
        write (lp,'(f12.4,a,i8)') gettime(),' Time 0 BLOM starting up'
        call flush(lp)
      endif
c
c --- ------------------------------------------------------------------
c --- Initialize check sum algorithm
c --- ------------------------------------------------------------------
c
      call crcinit
c
c --- ------------------------------------------------------------------
c --- Identify processor and horizontal indexes where detailed
c --- diagnostics are desired
c --- ------------------------------------------------------------------
c
      ptest=0
      if (itest.gt.i0.and.itest.le.i0+ii.and.
     .    jtest.gt.j0.and.jtest.le.j0+jj) then
        write (lp,'(a,i4,a,i4,a,i5)')
     .    ' itest =',itest,', jtest =',jtest,
     .    ' found on processor ',mnproc
        call flush(lp)
        ptest=mnproc
        itest=itest-i0
        jtest=jtest-j0
      endif
c
c --- ------------------------------------------------------------------
c --- Read limits file
c --- ------------------------------------------------------------------
c
      call rdlim
c
c --- ------------------------------------------------------------------
c --- Initialize the geographic environment
c --- ------------------------------------------------------------------
c
      call inigeo
c
c --- ------------------------------------------------------------------
c --- Initialize various arrays
c --- ------------------------------------------------------------------
c
      call inivar
c
c --- ------------------------------------------------------------------
c --- Initialize physical parameterizations
c --- ------------------------------------------------------------------
c
      call iniphy
c
c --- ------------------------------------------------------------------
c --- Initialize forcing
c --- ------------------------------------------------------------------
c
      call inifrc
c
c --- ------------------------------------------------------------------
c --- Define coefficients for equation of state functions
c --- ------------------------------------------------------------------
c
      call inieos
c
c --- ------------------------------------------------------------------
c --- Initialize shortwave radiation absorption
c --- ------------------------------------------------------------------
c
      call iniswa
c
c
#if defined(TRC) && defined(TKE)
c --- ------------------------------------------------------------------
c --- Initialize second order turbulence closure closure
c --- ------------------------------------------------------------------
c
      call initke
c
#endif
c --- ------------------------------------------------------------------
c --- Initialize diagnostic accumulation fields
c --- ------------------------------------------------------------------
c
      call diaini
c
c --- ------------------------------------------------------------------
c --- Set up initial conditions or start from restart file
c --- ------------------------------------------------------------------
c
c --- check whether initial condition file given in namelist is a
c --- restart file
      icrest=.false.
      if (mnproc.eq.1) then
        inquire(file=icfile,exist=fexist)
        if (fexist) then
          istat=nf90_open(icfile,nf90_nowrite,ncid)
          if (istat.eq.nf90_noerr) then
            istat=nf90_inq_varid(ncid,'dp',varid)
            if (istat.eq.nf90_noerr) then
              icrest=.true.
            endif
          endif 
        endif
      endif
      call xcbcst(icrest)
c
      if (nday1+nint(time0).eq.0.and..not.icrest) then
c
c --- ----------------------------------------------------------------
c --- - start from initial conditions derived from climatology
c --- ----------------------------------------------------------------
c
        if (nmonth.ne.1.or.nday.ne.1) then
          if (mnproc.eq.1) then
            write (lp,*)
     .'Warning! date is inconsistent with ini. cond. (Jan 1st)!'
            call flush(lp)
          endif
        endif
c
        call inicon
#ifdef TRC
        call initrc
#endif
c
        delt1=baclin
c
      else !  nday1+nint(time0) > 0 .or. icrest
c
c --- ------------------------------------------------------------------
c --- - start from restart file
c --- ------------------------------------------------------------------
c
        call restart_rd
c
        delt1=baclin+baclin
c
      endif
c
c --- ------------------------------------------------------------------
c --- Set layer thickness at u,v points
c --- ------------------------------------------------------------------
c
      call xctilr(dp, 1,2*kk, 3,3, halo_ps)
c
      n=mod(nstep1,2)+1
c
      do m=n,3-n,3-2*n
        mm=(m-1)*kk
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=-2,jj+2
          do k=1,kk
            do l=1,isp(j)
            do i=max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i,q)
        do j=-1,jj+2
          do k=1,kk
            km=k+mm
            do l=1,isu(j)
            do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
              q=min(p(i,j,kk+1),p(i-1,j,kk+1))
              dpu(i,j,km)=
     .          .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .             +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
            enddo
            enddo
            do l=1,isv(j)
            do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              q=min(p(i,j,kk+1),p(i,j-1,kk+1))
              dpv(i,j,km)=
     .          .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .             +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      nstep=nstep1
      if (mnproc.eq.1.and.expcnf.ne.'cesm') then
        write (lp,'(/2(a,i6),2(a,i9),a/)')
     .    'model starts at day',nday1,', goes to day',nday2,'   (steps',
     .    nstep1,' --',nstep2,')'
          call flush(lp)
      endif
c
      if (cnsvdi) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=pb(i,j,1)*scp2(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(mass0,util1,ips)
      endif
c
c --- ------------------------------------------------------------------
c --- update some halos
c --- ------------------------------------------------------------------
c
      call xctilr(uflx, 1,2*kk, 1,1, halo_uv)
      call xctilr(vflx, 1,2*kk, 1,1, halo_vv)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 2,2, halo_ps)
      call xctilr(pvtrop, 1,2, 1,3, halo_qs)
      call xctilr(pgfxm, 1,2, 1,2, halo_uv)
      call xctilr(xixp, 1,2, 1,2, halo_us)
      call xctilr(xixm, 1,2, 1,2, halo_us)
      call xctilr(pgfym, 1,2, 1,2, halo_vv)
      call xctilr(xiyp, 1,2, 1,2, halo_vs)
      call xctilr(xiym, 1,2, 1,2, halo_vs)
      call xctilr(uml, 1,4, 1,0, halo_uv)
      call xctilr(vml, 1,4, 0,1, halo_vv)
      call xctilr(umlres, 1,2, 1,0, halo_uv)
      call xctilr(vmlres, 1,2, 0,1, halo_vv)
c
c --- with arctic patch, switch xixp and xixm and xiyp and xiym in the
c --- halo region adjacent to the arctic grid intersection
      if (nreg.eq.2.and.nproc.eq.jpr) then
        do j=jj,jj+2
          do i=0,ii+1
            q=xixp(i,j,1)
            xixp(i,j,1)=xixm(i,j,1)
            xixm(i,j,1)=q
            q=xixp(i,j,2)
            xixp(i,j,2)=xixm(i,j,2)
            xixm(i,j,2)=q
          enddo
        enddo
        do i=max(0,itdm/2-i0+1),ii+1
          q=xiyp(i,jj,1)
          xiyp(i,jj,1)=xiym(i,jj,1)
          xiym(i,jj,1)=q
          q=xiyp(i,jj,2)
          xiyp(i,jj,2)=xiym(i,jj,2)
          xiym(i,jj,2)=q
        enddo
        do j=jj+1,jj+2
          do i=0,ii+1
            q=xiyp(i,j,1)
            xiyp(i,j,1)=xiym(i,j,1)
            xiym(i,j,1)=q
            q=xiyp(i,j,2)
            xiyp(i,j,2)=xiym(i,j,2)
            xiym(i,j,2)=q
          enddo
        enddo
      endif
c
      m=3-n
      call cmnfld(m,n,(m-1)*kk,(n-1)*kk,1+(m-1)*kk,1+(n-1)*kk)
c
c --- ------------------------------------------------------------------
c --- Extract reference potential density vector representative of the
c --- dominating ocean domain
c --- ------------------------------------------------------------------
c
      call diasg1
c
c --- print seconds elapsed since last call to system_clock (Time 0)
      if (mnproc.eq.1) then
        write (lp,'(f12.4,a,i8)')
     .    gettime(),' Time 1 Just before main loop'
        call flush(lp)
      endif
c
      return
      end
