      subroutine mxlayr(m,n,mm,nn,k1m,k1n)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_eos.h"
c
      real, dimension(kdm+1) :: pres
      real, dimension(kdm) :: ttem,ssal,dens,delp
      real sdot,tdp,sdp,thmxm,thmxn,sum1,sum2,value,turgen,
     .     swfbas,swfup,swflo,
     .     alfa,beta,dtemp,dsaln,buoyfl,em,en,thermg,
     .     dpth,ekminv,obuinv,ex,alf1,alf2,cp1,cp3,ape,cc4,spe,
     .     thkold,thknew,ustar3,
     .     pnew,tmxl,smxl,tnew,snew,dpn,tn,sn,z,
     .     a,b,c,d,e,f,thk,q,totem,tosal,
     .     tndcyt,tndcys,corrh,corrs
      integer i,j,k,l,kn,kbase,kbot
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      real, dimension(ntr) :: trdp,trmxl
      integer nt
#endif
c
      real dsigdt,dsigds,sofsig,sig,tsroot
      external dsigdt,dsigds,sofsig,sig,tsroot
c
      real ea1,ea2,em1,em2,em3,em4,em5,athird,onemu,tlim,slim
      parameter (ea1=.6,ea2=.3,em1=.45,em2=2.6,em3=1.9,em4=2.3,em5=.6,
     .           athird=1./3.,onemu=.09806,tlim=5.,slim=1.)
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,3p,f8.3,0p,f8.2,f8.1))
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  entering mxlayr:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),
cdiag.  sigma(itest,jtest,k+nn),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
c$OMP PARALLEL DO PRIVATE(thkold,thknew,ustar3,dpth,ekminv,obuinv,ex,
c$OMP+ alf1,alf2,cp1,cp3,ape,cc4,spe,kn,pnew,ttem,ssal,dens,delp,pres,
c$OMP+ q,tmxl,smxl,dpn,sn,tn,totem,tosal,tndcyt,tndcys,dtemp,
c$OMP+ dsaln,tnew,snew,z,a,b,c,d,e,f,
c$OMP+ thmxm,thmxn,sum1,sum2,value,kbase,tdp,sdp,sdot,turgen,
c$OMP+ swfbas,swfup,swflo,
c$OMP+ alfa,beta,
c$OMP+ buoyfl,em,en,thermg,corrh,corrs,kbot
#ifdef TRC
c$OMP+,ttrc,trdp,trmxl
#endif
c$OMP+ )
      do 1 j=1,jj
      do 1 l=1,isp(j)
      do 1 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- extract single column from 3-d fields
      pres(1)=p(i,j,1)
      do 7 k=1,kk
      kn=k+nn
      ttem(k)=temp(i,j,kn)
      ssal(k)=saln(i,j,kn)
      delp(k)=dp(i,j,kn)
      pres(k+1)=pres(k)+delp(k)
      dens(k)=sigmar(k)
#ifdef TRC
      do nt=1,ntrc
        ttrc(nt,k)=trc(i,j,k,nt)
      enddo
#endif
 7    continue
c
      dens(1)=sigma(i,j,k1n)
c
      thmxm=sigma(i,j,k1m)
      thmxn=sigma(i,j,k1n)
      thkold=delp(1)
c
      if (thermo) then
c
c --- store 'old' t/s column integral in totem/tosal (diagnostic use only)
      totem=ttem(1)*delp(1)
      tosal=ssal(1)*delp(1)
      do k=kk,2,-1
        totem=totem+ttem(k)*delp(k)
        tosal=tosal+ssal(k)*delp(k)
      end do
c
c --- ------------------------------------------------------------------
c --- compute temperature and salinity tendencies due to forcing
c --- ------------------------------------------------------------------
c
c --- if sdot < 0, mixed layer water will be detrained into isopycnic
c --- layer defined in -kbase-. to prevent odd/even time step decoupling
c --- of mixed- layer depth, determine -kbase- from -thmix- at 2
c --- consecutive time levels.  also find deepest physical layer -kbot-
c --- within a prescribed maximum radiation depth
c
      kbase=2
      kbot=1
      do k=2,kk
        if (max(thmxm,thmxn)+sigjmp.ge.dens(k)) kbase=k+1
        if (pres(k).lt.raddep*onem.and.delp(k).gt.onemu)
     .    kbot=k
      enddo
c
c --- adjust temperature profile due to penetrating solar radiation
      swfbas=(1.-redfac)*exp(-pres(kbase)/(betabl*onem))
      swfup=swfbas
      do k=kbase,kbot
        swflo=(1.-redfac)*exp(-min(raddep*onem,pres(k+1))/(betabl*onem))
        if (delp(k).gt.onemu) then
          ttem(k)=ttem(k)
     .           -(swfup-swflo)*sswflx(i,j)*delt1*g/(spcifh*delp(k))
          swfup=swflo
        endif
      enddo
c
      thkold=delp(1)
c
c --- temperature and salinity tendencies for mixed layer due to forcing
      dtemp=-(surflx(i,j)+surrlx(i,j)-(swfbas-swfup)*sswflx(i,j))
     .       *delt1*g/(spcifh*thkold)
      dsaln=-(salflx(i,j)+salrlx(i,j))*delt1*g/thkold
c
c --- buoyfl = buoyancy flux, w_prime_buoyancy_prime_bar (cm**2/sec**3)
c --- note: surface density increases (column is destabilized) if buoyfl > 0
      alfa=-alpha0*dsigdt(ttem(1),ssal(1))
      beta=alpha0*dsigds(ttem(1),ssal(1))
      buoyfl=g*(alfa*(surflx(i,j)-(swfbas-swfup)*sswflx(i,j))/spcifh
     .         -beta*salflx(i,j))
c
c --- ------------------------------------------------------------------
c --- mixed layer entrainment/detrainment
c --- ------------------------------------------------------------------
c
      ustar3=ustar(i,j)**3
c
c --- determine turb.kin.energy generation due to wind stirring (ustar3) and
c --- diabatic forcing (buoyfl).
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- option 1 :   k r a u s  -  t u r n e r    mixed-layer t.k.e. closure
c
ccc   em=0.8*exp(-pres(2)/(50.*onem))   !   hadley centre choice (orig.: 1.25)
ccc   en=0.15                           !   hadley centre choice (orig.: 0.4)
ccc   thermg=.5*((en+1.)*buoyfl+(en-1.)*abs(buoyfl))
ccc   turgen=delt1*(2.*em*g*ustar3/alpha0+thkold*thermg)/alpha0**2
c
c --- find monin-obukhov length in case of receding mixed layer (turgen < 0).
c --- the monin-obukhov length is found by stipulating turgen = 0.
c
ccc   if (turgen.lt.0.) then
ccc     thknew=-2.*em*g*ustar3/min(-epsil,alpha0*thermg)
ccc   else
ccc     thknew=thkold
ccc   end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- option 2 :    g a s p a r   mixed-layer t.k.e. closure
c
      dpth=thkold/onecm
      ekminv=abs(corio(i,j))/max(epsil,ustar(i,j))
      obuinv=-buoyfl/max(epsil,ustar3)
      ex=exp(max(-50.,min(50.,dpth*obuinv)))
      alf1=ea1+ea2*max(1.,2.5*dpth*ekminv)*ex
      alf2=ea1+ea2*ex
      cp1=((1.-em5)*(alf1/alf2)+.5*em4)*athird
      cp3=max(0.,(em4*(em2+em3)-(alf1/alf2)*(em2+em3-em3*em5))*athird)
      ape=cp3*ustar3+cp1*dpth*buoyfl
c
      if(ape.lt.0.) then                                       ! detrainment
      turgen=(g*delt1/alpha0**3)*ape
      thknew=min(thkold,g*cp3/(alpha0*cp1*max(epsil,obuinv)))
c
      else                                                     ! entrainment
      cc4=2.*em4/(em1*em1) * alf1*alf1
      spe=(em2+em3)*ustar3+0.5*dpth*buoyfl
      turgen=(g*delt1/alpha0**3)*(sqrt((.5*ape-cp1*spe)**2
     .        +2.*cc4*ape*spe)-(.5*ape+cp1*spe))/(cc4-cp1)
      thknew=thkold
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- sum1,sum2 are used to evaluate pot.energy changes during entrainment
      sum1=dens(1)*thkold
 86   sum2=dens(1)*thkold**2
c
c --- find pnew in case of mixed layer deepening (turgen > 0).
c --- entrain as many layers as needed to deplete -turgen-
c
      do 85 k=2,kk
      kn=k+nn
      value=min(pres(k+1),(2.*turgen+dens(k)*pres(k)**2-sum2)/
     .                     max(epsil,dens(k)*pres(k)   -sum1))
c --- stop iterating for 'pnew' as soon as pnew < k-th interface pressure
      if (value.lt.pres(k)) value=thknew
c --- substitute 'pnew' for monin-obukhov length if mixed layer is deepening
      if (turgen.ge.0.) thknew=value
c
      sum1=sum1+dens(k)*delp(k)
 85   sum2=sum2+dens(k)*(pres(k+1)**2-pres(k)**2)
c
      else				!  thermo = .false.
c
      surflx(i,j)=0.
      surrlx(i,j)=0.
      salflx(i,j)=0.
      salrlx(i,j)=0.
      thknew=thkold
      dtemp=0.
      dsaln=0.
c
      end if				!  thermo = .true. or .false.
c
cdiag if (i.eq.itest .and. j.eq.jtest .and. turgen.lt.0.)
cdiag.  write (lp,'(i9,2i5,a,f8.2)') nstep,i,j,
cdiag.  '  monin-obukhov length (m):',thknew/onem
c
c --- store (pnew - pold) in 'sdot'.
c --- don't allow mixed layer to get too deep or too shallow.
      sdot=min(pres(kk+1),max(thkmin*onem,thknew))-delp(1)
      pnew=thkold+sdot
c
      tdp=0.
      sdp=0.
#ifdef TRC
      do nt=1,ntrc
        trdp(nt)=0.
      enddo
#endif
c
      do 43 k=2,kk
c --- 'tdp,sdp' will be needed for temp./salin. mixing during entrainment
      q=(min(pnew,pres(k+1))-min(pnew,pres(k)))
      tdp=tdp+ttem(k)*q
      sdp=sdp+ssal(k)*q
#ifdef TRC
      do nt=1,ntrc
        trdp(nt)=trdp(nt)+ttrc(nt,k)*q
      enddo
#endif
c
c --- if sdot > 0, remove water from layers about to be entrained.
      q=delp(k)					! diapyc.flux
      delp(k)=max(pres(k+1),pnew)-max(pres(k),pnew)
      diaflx(i,j,k)=diaflx(i,j,k)+(delp(k)-q)	! diapyc.flux
 43   continue
c
c --- increase temperature in thin layers where water is too cold to fit
c --- into layer k
c 
      corrh=0.
      corrs=0.
      do 45 k=2,kk
      if (delp(k).lt.onemm) then
        q=ttem(k)
        ttem(k)=max(ttem(k),temmin(k))
        corrh=corrh+(ttem(k)-q)*delp(k)
        q=ssal(k)
        ssal(k)=sofsig(dens(k),ttem(k))
        corrs=corrs+(ssal(k)-q)*delp(k)
      end if
 45   continue
c
c --- redistribute temp. and salin. during both de- and entrainment
c
      tmxl=ttem(1)+dtemp
      smxl=ssal(1)+dsaln
#ifdef TRC
      do nt=1,ntrc
        trmxl(nt)=ttrc(nt,1)
      enddo
#endif
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,'(i9,2i5,a,3f7.3,f8.2)')
cdiag.  nstep,i,j,'  t,s,sig,dp after diab.forcing',tmxl,smxl,
cdiag.  sig(tmxl,smxl)*1000.,thkold/onem
c
      if (sdot.ge.0.) then
c
c --- (mixed layer  d e e p e n s)
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,'(i9,2i5,a,f9.3,a)')
cdiag.  nstep,i,j,'  entrain',sdot/onem,' m of water'
c
      q=1./pnew
      tmxl=(tmxl*thkold+tdp)*q
      smxl=(smxl*thkold+sdp)*q
#ifdef TRC
      do nt=1,ntrc
        trmxl(nt)=(trmxl(nt)*thkold+trdp(nt))*q
      enddo
#endif
      delp(1)=pnew
      diaflx(i,j,1)=diaflx(i,j,1)+sdot ! diapyc.flux
c
      else if (sdot.lt.-onecm.and.-alfa*dtemp+beta*dsaln.lt.0.) then ! sdot < 0
c
c --- (mixed layer  r e c e d e s)
c
      k=kbase
      if (k.gt.kk) go to 27
c
      dpn=max(delp(k),0.)
      sn=ssal(k)
      tn=ttem(k)
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,k,
cdiag.  '  t,s,dp before detrainment',tn,sn,dpn/onem
c
c --- distribute last time step's heating and freshwater flux over depth range
c --- 'pnew' (monin-obukhov length). split fossil mixed layer (depth= -sdot=
c --- thkold-pnew) into lower part ('lo') of depth z cooled and detrained into
c --- layer k, and an upper part ('up') heated and/or freshened to match
c --- temperature rise and/or salinity decrease in mixed layer.
c
      dtemp=dtemp*thkold/pnew
      dsaln=dsaln*thkold/pnew
      tnew=ttem(1)+dtemp
      snew=ssal(1)+dsaln
c
c --- (t_lo,s_lo) differ from (tn,sn). main problem now is in maintaining
c --- density in layer k during detrainment. This requires solving for z with
c --- new (t,s) in layer k of the form t=(a*z+b)/(z+d), s=(e*z+f)/(z+d)
c --- substituted into the functional form of potential density.
      a=tnew
      b=(tn*dpn+dtemp*sdot)/onem
      c=1.
      d=dpn/onem
      e=snew
      f=(sn*dpn+dsaln*sdot)/onem
c
      z=tsroot(a,b,c,d,e,f,sig(tn,sn),-sdot/onem,0.)
c
c --- does root fall into appropriate range?
      if (z*onem+sdot.gt.0.) go to 27
      if (z.le.0.005) go to 27
      q=1./(z+d)
      a=(a*z+b)*q
      b=(e*z+f)*q
      if (a.lt.min(tn,temmin(k)).or.
     .    abs(ttem(1)-a).gt.tlim.or.abs(ssal(1)-b).gt.slim) go to 27
c
c --- ready to detrain lowest 'z' meters from mixed layer
c
      ttem(k)=a
      ssal(k)=b
#ifdef TRC
      do nt=1,ntrc
        ttrc(nt,k)=(ttrc(nt,1)*z+ttrc(nt,k)*d)*q
      enddo
#endif
c
      sdot=-z*onem
      delp(1)=thkold+sdot
      delp(k)=dpn   -sdot
      smxl=snew
      tmxl=tnew
      diaflx(i,j,1)=diaflx(i,j,1)+sdot				! diapyc.flux
      diaflx(i,j,k)=diaflx(i,j,k)-sdot				! diapyc.flux
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,k,
cdiag.  '  t,s,dp after  detrainment',ttem(k),ssal(k),
cdiag.  delp(k)/onem
c
      end if				!  sdot > or < 0
c
 27   ttem(1)=tmxl
      ssal(1)=smxl
#ifdef TRC
      do nt=1,ntrc 
        ttrc(nt,1)=trmxl(nt)
      enddo
#endif
c
cdiag if (i.eq.itest.and.j.eq.jtest) write
cdiag.  (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,1,
cdiag.  ' final mixed-layer t,s,dp ',tmxl,smxl,delp(1)/onem
c
c --- compare 'old' with 'new' t/s column integral (diagnostic use only)
c
      tndcyt=-totem-corrh
      tndcys=-tosal-corrs
      do k=kk,1,-1
        tndcyt=tndcyt+ttem(k)*delp(k)
        tndcys=tndcys+ssal(k)*delp(k)
      end do
      q=100.*acurcy*pb_p(i,j)*g
      if (abs((surflx(i,j)+surrlx(i,j))*delt1*g+tndcyt*spcifh).gt.10.*q
     ..or.abs((salflx(i,j)+salrlx(i,j))*delt1*g+tndcys).gt.35.*q) then
      q=1./pb_p(i,j)
        write (lp,'(i9,2i5,3x,a,2f9.5/22x,a,2f9.5)') nstep,i0+i,j0+j,
     .    'saln flux,tndcy:',
     .    q*(salflx(i,j)+salrlx(i,j))*delt1*g       ,q*tndcys,
     .    'heat flux,tndcy:',
     .    q*(surflx(i,j)+surrlx(i,j))*delt1*g/spcifh,q*tndcyt
      end if
c
c --- put single column back into 3-d fields
      do 8 k=1,kk
      kn=k+nn
      temp(i,j,kn)=ttem(k)
      saln(i,j,kn)=ssal(k)
      sigma(i,j,kn)=sig(ttem(k),ssal(k))
      dp(i,j,kn)=delp(k)
#ifdef TRC
      do nt=1,ntrc 
        trc(i,j,k,nt)=ttrc(nt,k)
      enddo
#endif
 8    continue
c
c
 1    continue
c$OMP END PARALLEL DO
c
c --- store 'old' interface pressures in -pu,pv-
c
      do 882 k=2,kk+1
c$OMP PARALLEL DO
      do j=1,jj
c
        do 881 l=1,isu(j)
        do 881 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 881    pu(i,j,k)=min(pbu(i,j,n),.5*(p(i  ,j,k)+p(i-1,j,k)))
c
        do 880 l=1,isv(j)
        do 880 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 880    pv(i,j,k)=min(pbv(i,j,n),.5*(p(i,j  ,k)+p(i,j-1,k)))
      end do
c$OMP END PARALLEL DO
 882  continue
c
c --- store 'new' layer thicknesses in -dpu,dpv-
c
      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
c
c$OMP PARALLEL DO
      do 883 j=-2,jj+2
      do 883 k=1,kk
      do 883 l=1,isp(j)
      do 883 i=max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
 883  p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c$OMP END PARALLEL DO
c
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=-1,jj+2
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            dpu(i,j,kn)=
     .        max(0.,min(pbu(i,j,n),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .              -min(pbu(i,j,n),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            dpv(i,j,kn)=
     .        max(0.,min(pbv(i,j,n),.5*(p(i,j  ,k+1)+p(i,j-1,k+1)))
     .              -min(pbv(i,j,n),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- redistribute momentum in the vertical.
c --- homogenize (u,v) over depth range defined in -util1,util2-
c
c --- thk>0 activates momentum diffusion across mixed-layer interface
      thk=vertmx*onecm*delt1
c
c$OMP PARALLEL DO PRIVATE(delp,kn,q)
      do 97 j=1,jj
c
      do 83 l=1,isu(j)
c
      do 822 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      util1(i,j)=max(dpu(i,j,k1n),pu(i,j,2)+thk)
      uflux(i,j)=0.
 822  util3(i,j)=0.
c
      do 82 k=1,kk
      do 82 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      delp(1)=max(0.,min(util1(i,j),pu(i,j,k+1))
     .              -min(util1(i,j),pu(i,j,k  )))
      uflux(i,j)=uflux(i,j)+u(i,j,k+nn)*delp(1)
 82   util3(i,j)=util3(i,j)            +delp(1)
c
      do 83 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 83   u(i,j,k1n)=uflux(i,j)/util3(i,j)
c
      do 84 l=1,isv(j)
c
      do 844 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      util2(i,j)=max(dpv(i,j,k1n),pv(i,j,2)+thk)
      vflux(i,j)=0.
 844  util4(i,j)=0.
c
      do 80 k=1,kk
      do 80 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      delp(1)=max(0.,min(util2(i,j),pv(i,j,k+1))
     .              -min(util2(i,j),pv(i,j,k  )))
      vflux(i,j)=vflux(i,j)+v(i,j,k+nn)*delp(1)
 80   util4(i,j)=util4(i,j)            +delp(1)
c
      do 84 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 84   v(i,j,k1n)=vflux(i,j)/util4(i,j)
c
      do 97 k=2,kk
      kn=k+nn
c
      do 96 l=1,isu(j)
      do 96 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      pu(i,j,k)=pu(i,j,k-1)+dpu(i,j,kn-1)
      q=max(0.,min(1.,(util1(i,j)-pu(i,j,k))/(dpu(i,j,kn)+epsil)))
 96   u(i,j,kn)=u(i,j,k1n)*q+u(i,j,kn)*(1.-q)
c
      do 97 l=1,isv(j)
      do 97 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      pv(i,j,k)=pv(i,j,k-1)+dpv(i,j,kn-1)
      q=max(0.,min(1.,(util2(i,j)-pv(i,j,k))/(dpv(i,j,kn)+epsil)))
 97   v(i,j,kn)=v(i,j,k1n)*q+v(i,j,kn)*(1.-q)
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  exiting  mxlayr:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),
cdiag.  sigma(itest,jtest,k+nn),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'mxlayr:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(diaflx,ip,kk,'diaflx')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
        call chksummsk(dpu,iu,2*kk,'dpu')
        call chksummsk(dpv,iv,2*kk,'dpv')
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
c
c
c> Revision history:
c>
c> June 1995 - removed restriction  'klist(i,j) .le. kk'
c> June 1995 - added code for setting t/s in massless layers below mix.layer
c> Oct. 1995 - removed bug created while changing klist (June 1995 revision):
c>             'if (k.gt.kk) go to 26' now reads 'if (k.gt.kk) go to 27'
c> May  1997 - changed -sdot- into local array
c> Mar. 1998 - added -th3d-
c> Nov. 1998 - fixed bug in computing tnew,snew in situations where z < 0.005
c> Dec. 1998 - replaced dsaln by (s_up-s1) in definition of 'f'
c> Feb. 1999 - limited 'tofsig' call in loop 45 to cases where saln < salmin
c> Apr. 2000 - search for 3rd deg.polynomial roots done in double precision
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> June 2000 - changed loop hierarchy to have one overarching j/i loop
c> July 2000 - added "min(pres(k+1),..." to calculation of -value- in loop 85
c> Sep. 2000 - changed dimension of -dens- from (1:kdm) to (0:kdm)
