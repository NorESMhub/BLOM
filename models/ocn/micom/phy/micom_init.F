      subroutine micom_init
c
c --- ------------------------------------------------------------------
c --- initialize the model
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_bud.h"
#include "common_bud.h"
c
cbud  real vol0j(jdm),vol0
      integer i,j,k,l,m,n,mm,nn,km,nday1,nday2
      real q
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
        write (lp,'(f12.4,a,i8)') gettime(),' Time 0 MICOM starting up'
        call flush(lp)
      endif
c
c --- ------------------------------------------------------------------
c --- Consistency check of preprocessor identifiers
c --- ------------------------------------------------------------------
#if !defined(CCSMCOUPLED) && !defined(CLIM) && !defined(SYN)
      if (mnproc.eq.1) then
        write (lp,'(2a)')
     .    ' Preproc. ident. check: Either CCSMCOUPLED, CLIM ',
     .    ', or SYN must be defined!'
      endif
      call xcstop('(micom)')
             stop '(micom)'
#endif
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
      call rdlim(nday1,nday2)
c
c --- ------------------------------------------------------------------
c --- Initialize the calendar
c --- ------------------------------------------------------------------
c
      call inicln
c
c --- ------------------------------------------------------------------
c --- Set up parameters defining the geographic environment
c --- ------------------------------------------------------------------
c
      call geopar
c
#ifdef ICEDYN
c --- ------------------------------------------------------------------
c --- Initialization of sea ice dynamics
c --- ------------------------------------------------------------------
c
      call iniicedyn(scpx,scpy,corio)
#endif
c
c --- ------------------------------------------------------------------
c --- Set up forcing functions
c --- ------------------------------------------------------------------
c
#ifndef CCSMCOUPLED
      call forfun
#  ifdef CLIM
      call forfun_clim
      call chwbud
#  endif
#  ifdef SYN
      call forfun_syn
#  endif
#endif
c
#ifdef NEST
c --- ------------------------------------------------------------------
c --- Initialize nesting
c --- ------------------------------------------------------------------
c
      call ininst
c
#endif
c --- ------------------------------------------------------------------
c --- Define coefficients for equation of state functions, alphar, and
c --- salmin
c --- ------------------------------------------------------------------
c
      call eosini
c
c --- layer specific volume is defined as  alphar=1/(1+sigmar)
c
      do 14 k=1,kk
      sigmar(k)=sigmar(k)*1.e-3 ! [kg/m^3] -> [g/cm^3]
 14   alphar(k)=1./(1.+sigmar(k))
c
c --- set minimum physical temperature for each isopycnic layer
      call settemmin
c
c --- ------------------------------------------------------------------
c --- Initialize diagnostic accumulation fields
c --- ------------------------------------------------------------------
c
      call diaini2d
      call diaini3d
c
c --- ------------------------------------------------------------------
c --- Set up initial conditions or start from restart file
c --- ------------------------------------------------------------------
c
      if (nstep1.eq.0) then
c
c --- ------------------------------------------------------------------
c --- set up initial conditions
c --- ------------------------------------------------------------------
c
      if (nmonth.ne.1.or.nday.ne.1) then
        if (mnproc.eq.1) then
          write (lp,*)
     .      'Warning! date is inconsistent with ini. cond. (Jan 1st)!'
          call flush(lp)
        endif
      endif
c
      call inicon
#ifdef ICEDYN
      call iniconicedyn
#endif
#ifdef TRC
      call initrc
#endif
c
      delt1=baclin
c
      time0=0.
c
      else !  nstep1 > 0
c
c --- ------------------------------------------------------------------
c --- start from restart file
c --- ------------------------------------------------------------------
c
      call restart_rd(nday1)
c
      delt1=baclin+baclin
c
      end if !  nstep1 > 0  or  = 0
c
c --- ------------------------------------------------------------------
c --- Set layer thickness at u,v points
c --- ------------------------------------------------------------------
c
      call xctilr(dp, 1,2*kk, 3,3, halo_ps)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
      do m=1,2
        mm=(m-1)*kk
c
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO PRIVATE(km)
        do j=-1,jj+2
          do k=1,kk
            km=k+mm
            do l=1,isu(j)
            do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
              dpu(i,j,km)=
     .          max(0.,min(pbu(i,j,m),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .                -min(pbu(i,j,m),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
            enddo
            enddo
            do l=1,isv(j)
            do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              dpv(i,j,km)=
     .          max(0.,min(pbv(i,j,m),.5*(p(i,j  ,k+1)+p(i,j-1,k+1)))
     .                -min(pbv(i,j,m),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      nstep=nstep1
      if (mnproc.eq.1) then
        write (lp,'(/2(a,i6),2(a,i9),a/)')
     .    'model starts at day',nday1,', goes to day',nday2,'   (steps',
     .    nstep1,' --',nstep2,')'
          call flush(lp)
      endif
c
cbud  do j=1,jj
cbud    vol0j(j)=0.
cbud  enddo
c$OMP PARALLEL DO
cbud  do j=1,jj
cbud    do l=1,isp(j)
cbud    do i=ifp(j,l),ilp(j,l)
cbud      vol0j(j)=vol0j(j)+pb(i,j,1)*scp2(i,j)
cbud    enddo
cbud    enddo
cbud  enddo
c$OMP END PARALLEL DO
cbud  vol0=0.
cbud  do j=1,jj
cbud    vol0=vol0+vol0j(j)
cbud  enddo
c
c --- ------------------------------------------------------------------
c --- update some halos
c --- ------------------------------------------------------------------
c
      m=mod(nstep  ,2)+1
      call xctilr(ubflxs, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs, 1,2, 2,2, halo_vv)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)
      call xctilr(pvtrop(1-nbdy,1-nbdy,m), 1,1, 1,3, halo_qs)
      call xctilr(pgfxm(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_uv)
      call xctilr(xixp(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_us)
      call xctilr(xixm(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_us)
      call xctilr(pgfym(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_vv)
      call xctilr(xiyp(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_vs)
      call xctilr(xiym(1-nbdy,1-nbdy,m), 1,1, 1,2, halo_vs)
c
#ifdef ARCTIC
c --- with arctic patch, switch xixp and xixm and xiyp and xiym in the
c --- halo region adjacent to the arctic grid intersection
      if (nproc.eq.jpr) then
        do j=jj+1,jj+2
          do i=0,ii
            q=xixp(i,j,m)
            xixp(i,j,m)=xixm(i,j,m)
            xixm(i,j,m)=q
            q=xiyp(i,j,m)
            xiyp(i,j,m)=xiym(i,j,m)
            xiym(i,j,m)=q
          enddo
        enddo
      endif
#endif
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
