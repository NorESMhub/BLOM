      subroutine mxlayr(m,n,mm,nn,k1m,k1n)
c
c --- -----------------------------------------------------------------
c --- Modify mixed layer depth by enforcing a turbulent kinetic energy
c --- balance and apply surface forcing.
c --- -----------------------------------------------------------------
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
      real, dimension(kdm+1) :: pres,po,pn
      real, dimension(kdm) :: ttem,ssal,delp,dens,densr,uo,un
      real q,ustar3,alfa,beta,bfltot,bflpsw,pmxl,lui,lbi,lei,cus,
     .     cbftot,cbfpsw,crs,rm1,rm2,rm3,rm4,tke,dpmxl,tkeo,dtke,
     .     tfsl,sfsl,dpfsl,pswbas,pswup,pswlo,ttmp,stmp,sigtmp,sigfsl,
     .     tmxl0,smxl0,dpe0,tdps,sdps,tmxl,smxl,dpe,dps,udpn
      integer i,j,k,l,kn,nitr,ko
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      real, dimension(ntr) :: trfsl,trdps
      integer nt
#endif
c
      real sig0,dsigdt0,dsigds0,sig,p_alpha,p_p_alpha
      external sig0,dsigdt0,dsigds0,sig,p_alpha,p_p_alpha
c
      real onemu
      parameter (onemu=.09806)
c
c --- Parameters for Oberhuber (1993) TKE closure:
c ---  ustmin - minium value of ustar used in computing the length
c ---           scales for wind and buoyancy induced mixing [cm/s]
c ---  mldjmp - minimum density jump at the mixed layer base used in the
c ---           computation of potential energy change due to
c ---         - entrainment
      real rm0,kappa,mu,ustmin,mldjmp
      integer maxitr
      parameter (rm0=1.25,kappa=.4,mu=1.,ustmin=.1,mldjmp=1.e-6,
     .           maxitr=20)
c
c --- Parameters for the parameterization of restratification by mixed
c --- layer eddies by Fox-Kemper et al. (2008):
c ---   ce   - efficiency factor []
c ---   rtau - reciprocal of timescale for momentum mixing across mixed
c ---          layer [1/s]
c ---   rlf  - reciprocal of length scale of mixed layer fronts [1/cm]
c ---   ci   - constant that appears when integrating the shape function
c ---          over the mixed layer depth []
      real ce,rtau,rlf,ci
      parameter (ce=.06,rtau=1./(2.*86400.),rlf=1./5.e5,ci=44./63.)
c
c --- ------------------------------------------------------------------
c --- Compute squared lateral buoyancy gradient in the mixed layer and
c --- store it in -util1-.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=g*alpha0*sig0(temp(i,j,k1n),saln(i,j,k1n))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO PRIVATE(q)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          q=(util1(i,j)-util1(i-1,j))*scuxi(i,j)
          util2(i,j)=q*q
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(q)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          q=(util1(i,j)-util1(i,j-1))*scvyi(i,j)
          util3(i,j)=q*q
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if     (ip(i-1,j)+ip(i+1,j).eq.2) then
            util1(i,j)=.5*(util2(i,j)+util2(i+1,j))
          elseif (ip(i-1,j).eq.1) then
            util1(i,j)=util2(i  ,j)
          elseif (ip(i+1,j).eq.1) then
            util1(i,j)=util2(i+1,j)
          else
            util1(i,j)=0.
          endif
          if     (ip(i,j-1)+ip(i,j+1).eq.2) then
            util1(i,j)=util1(i,j)+.5*(util3(i,j)+util3(i,j+1))
          elseif (ip(i,j-1).eq.1) then
            util1(i,j)=util1(i,j)+util3(i,j)
          elseif (ip(i,j+1).eq.1) then
            util1(i,j)=util1(i,j)+util3(i,j+1)
          endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ pres,ttem,ssal,delp,dens,densr,
c$OMP+ q,ustar3,alfa,beta,bfltot,bflpsw,pmxl,lui,lbi,lei,cus,
c$OMP+ cbftot,cbfpsw,crs,rm1,rm2,rm3,rm4,tke,dpmxl,tkeo,dtke,
c$OMP+ tfsl,sfsl,dpfsl,pswbas,pswup,pswlo,ttmp,stmp,sigtmp,sigfsl,
c$OMP+ tmxl0,smxl0,dpe0,tdps,sdps,tmxl,smxl,dpe,dps,
c$OMP+ k,kn,nitr
#ifdef TRC
c$OMP+,ttrc,trfsl,trdps
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Extract single column from 3-d fields.
          pres(1)=p(i,j,1)
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            pres(k+1)=pres(k)+delp(k)
            dens(k)=sigma(i,j,kn)
            densr(k)=sigmar(i,j,k)
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,k)=trc(i,j,k,nt)
            enddo
#endif
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Compute the turbulent kinetic energy (TKE) balance of the
c --- --- mixed layer.
c --- ------------------------------------------------------------------
c
          ustar3=ustar(i,j)**3
c
c --- --- bfltot = total buoyancy flux [cm**2/sec**3]
c --- --- bflpsw = buoyancy flux due to penetrating short-wave radiation
c ---              [cm**2/sec**3]
c --- --- note: surface density increases (column is destabilized) if
c --- --- bfltot > 0
          alfa=-alpha0*dsigdt0(ttem(1),ssal(1))
          beta=alpha0*dsigds0(ttem(1),ssal(1))
          bfltot=g*(alfa*surflx(i,j)/spcifh-beta*salflx(i,j))
          bflpsw=g*alfa*(1.-redfac)*sswflx(i,j)/spcifh
c
          pmxl=pres(2)
          q=alpha0/g
          lui=abs(coriop(i,j))*q/(kappa*max(ustmin,ustar(i,j)))
          if (bfltot.lt.0) then
            lbi=lui
          else
            lbi=abs(coriop(i,j))*q/(mu*max(ustmin,ustar(i,j)))
          endif
          lei=1./(onem*betabl)
          cus=rm0*ustar3
          cbftot=.5*bfltot*q
          cbfpsw=.5*bflpsw*q
          crs=ci*ce*util1(i,j)*q**3
     .        *sqrt(scp2(i,j)/(coriop(i,j)*coriop(i,j)+rtau*rtau))*rlf
c      
          rm1=exp(-lui*pmxl)
          rm2=exp(-lbi*pmxl)
          q=lei*pmxl
          rm3=exp(-q)
          rm4=2./q
          tke=cus*rm1
     .       +((cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))*rm2
     .        -crs*pmxl*pmxl)*pmxl
c      
          if (tke.lt.0.and.pmxl.gt.thkmin*onem) then
c
c --- -----------------------------------------------------------------
c --- ----- If there is a TKE deficit in the mixed layer, reduce the
c --- ----- mixed layer depth until the TKE balance is restored.
c --- -----------------------------------------------------------------
c
            if (pres(2)*lbi.gt.1.) then
              pmxl=1./lbi
              dpmxl=min(pmxl-pres(1),pres(2)-pmxl,onecm)
              pmxl=pmxl-.5*dpmxl
            else
              dpmxl=-onecm
              pmxl=pres(2)+dpmxl
            endif
            tkeo=tke
            nitr=0
            do
              nitr=nitr+1
              rm1=exp(-lui*pmxl)
              rm2=exp(-lbi*pmxl)
              q=lei*max(onecm,pmxl)
              rm3=exp(-q)
              rm4=2./q
              tke=cus*rm1
     .           +((cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))*rm2
     .            -crs*pmxl*pmxl)*pmxl
              if (.not.(nitr.eq.1.and.pres(2)*lbi.gt.1.)) then
                dtke=(tke-tkeo)/dpmxl
                if (abs(dtke)<(abs(tke)+1.e-16)/delp(1)) then
                  if (tke.lt.0.) then
                    dpmxl=.5*(pres(1)-pmxl)
                  else
                    dpmxl=.5*(pres(2)-pmxl)
                  endif
                else
                  dpmxl=max(pres(1)-pmxl,min(pres(2)-pmxl,-tke/dtke))
                endif
              endif
              pmxl=pmxl+dpmxl
              tkeo=tke
              if (abs(dpmxl).lt.onemm.or.nitr.eq.maxitr) exit
            enddo
            if (nitr.eq.maxitr) then
              write (lp,*) 'reached maxitr when detraining',i+i0,j+j0
              write (lp,*) 'dpth=',pres(2)/onem,';'
              write (lp,*) 'pmxl=',pmxl/onem,';'
              write (lp,*) 'corio=',coriop(i,j),';'
              write (lp,*) 'ustar=',ustar(i,j)*1.e-2,';'
              write (lp,*) 'bfltot=',bfltot*1.e-4,';'
              write (lp,*) 'bflpsw=',bflpsw*1.e-4,';'
              write (lp,*) 'bg2=',util1(i,j),';'
              write (lp,*) 'ce=',ce*sqrt(scp2(i,j))*rlf,';'
              write (lp,*)
c             call xchalt('(mxlayr)')
c                    stop '(mxlayr)'
            endif
c
            pmxl=max(thkmin*onem,pmxl)
            tfsl=ttem(1)
            sfsl=ssal(1)
#ifdef TRC
            do nt=1,ntrc
              trfsl(nt)=ttrc(nt,1)
            enddo
#endif
            dpfsl=pres(2)-pmxl
            delp(1)=pmxl-pres(1)
c
c --- ----- Apply forcing to the water column.
            pswbas=(1.-redfac)*exp(-lei*pmxl)
            pswup=pswbas
            pswlo=(1.-redfac)*exp(-lei*min(raddep*onem,pres(2)))
            if (dpfsl.gt.onemu) then
              tfsl=tfsl
     .            -(pswup-pswlo)*sswflx(i,j)*delt1*g/(spcifh*dpfsl)
              pswup=pswlo
            endif
            do k=2,kk
              if (pres(k).lt.raddep*onem.and.delp(k).gt.onemu) then
                pswlo=(1.-redfac)*exp(-lei*min(raddep*onem,pres(k+1)))
                ttem(k)=ttem(k)
     .                 -(pswup-pswlo)*sswflx(i,j)*delt1*g
     .                  /(spcifh*delp(k))
                dens(k)=sig(ttem(k),ssal(k))
                pswup=pswlo
              endif
            enddo
            q=delt1*g/delp(1)
            ttem(1)=ttem(1)
     .             -(surflx(i,j)+surrlx(i,j)-(pswbas-pswup)*sswflx(i,j))
     .              *q/spcifh
            ssal(1)=ssal(1)-(salflx(i,j)+salrlx(i,j))*q
            dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,1)=ttrc(nt,1)-trflx(nt,i,j)*q
            enddo
#endif
c
            if (dpfsl.le.onemu) then
              q=1./(dpfsl+delp(1))
              ttem(1)=(tfsl*dpfsl+ttem(1)*delp(1))*q
              ssal(1)=(sfsl*dpfsl+ssal(1)*delp(1))*q
              dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,1)=(trfsl(nt)*dpfsl+ttrc(nt,1)*delp(1))*q
              enddo
#endif
              delp(1)=dpfsl+delp(1)
            else
c
c --- ------- Place the content of the fossil mixed layer in isopycnic
c --- ------- layers.
              k=min(kk,kfpla(i,j,n))
c
              if (k.eq.2) then
                q=1./(dpfsl+delp(k))
                ttem(k)=(tfsl*dpfsl+ttem(k)*delp(k))*q
                ssal(k)=(sfsl*dpfsl+ssal(k)*delp(k))*q
                dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,k)=(trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                enddo
#endif
                delp(k)=dpfsl+delp(k)
              else
                q=1./(dpfsl+delp(k))
                ttmp=(tfsl*dpfsl+ttem(k)*delp(k))*q
                stmp=(sfsl*dpfsl+ssal(k)*delp(k))*q
                sigtmp=sig(ttmp,stmp)
                sigfsl=sig(tfsl,sfsl)
                if (sigtmp.ge.densr(k)) then
                  if (sigfsl.gt.dens(k).and.
     .                dens(k).le.dens(min(kk,k+1))) then
                    k=k-1
                    q=1./(dpfsl+delp(k))
                    ttem(k)=(tfsl*dpfsl+ttem(k)*delp(k))*q
                    ssal(k)=(sfsl*dpfsl+ssal(k)*delp(k))*q
                    dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=(trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                    enddo
#endif
                    delp(k)=dpfsl+delp(k)
                  else
                    ttem(k)=ttmp
                    ssal(k)=stmp
                    dens(k)=sigtmp
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=(trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                    enddo
#endif
                    delp(k)=dpfsl+delp(k)
                  endif
                else
                  if (delp(k).gt.onemu.and.dens(k).gt.densr(k).and.
     .                sigfsl.lt.densr(k)-1.e-9) then
                    dps=min(dpfsl,delp(k)*(dens(k)-densr(k))
     .                                   /(densr(k)-sigfsl))
                    q=1./(dps+delp(k))
                    ttem(k)=(tfsl*dps+ttem(k)*delp(k))*q
                    ssal(k)=(sfsl*dps+ssal(k)*delp(k))*q
                    dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=(trfsl(nt)*dps+ttrc(nt,k)*delp(k))*q
                    enddo
#endif
                    delp(k)=dps+delp(k)
                    dpfsl=dpfsl-dps
                    if (dpfsl.le.onemu) then
                      q=1./(dpfsl+delp(1))
                      ttem(1)=(tfsl*dpfsl+ttem(1)*delp(1))*q
                      ssal(1)=(sfsl*dpfsl+ssal(1)*delp(1))*q
                      dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                      do nt=1,ntrc
                        ttrc(nt,1)=(trfsl(nt)*dpfsl
     .                             +ttrc(nt,1)*delp(1))*q
                      enddo
#endif
                      delp(1)=dpfsl+delp(1)
                    else
                      k=k-1
                      do while (sigfsl.lt.densr(k))
                        if (k.eq.2) exit
                        k=k-1
                      enddo
                      q=1./(dpfsl+delp(k))
                      ttem(k)=(tfsl*dpfsl+ttem(k)*delp(k))*q
                      ssal(k)=(sfsl*dpfsl+ssal(k)*delp(k))*q
                      dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                      do nt=1,ntrc
                        ttrc(nt,k)=(trfsl(nt)*dpfsl
     .                             +ttrc(nt,k)*delp(k))*q
                      enddo
#endif
                      delp(k)=dpfsl+delp(k)
                    endif
                  else
                    k=k-1
                    do while (sigfsl.lt.densr(k))
                      if (k.eq.2) exit
                      k=k-1
                    enddo
                    q=1./(dpfsl+delp(k))
                    ttem(k)=(tfsl*dpfsl+ttem(k)*delp(k))*q
                    ssal(k)=(sfsl*dpfsl+ssal(k)*delp(k))*q
                    dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=(trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                    enddo
#endif
                    delp(k)=dpfsl+delp(k)
                  endif
                endif
              endif
c
            endif
c
          else
c
            if (tke.lt.0.) then
              pmxl=thkmin*onem
              tdps=ttem(1)*delp(1)
              sdps=ssal(1)*delp(1)
#ifdef TRC
              do nt=1,ntrc
                trdps(nt)=ttrc(nt,1)*delp(1)
              enddo
#endif
              k=kfpla(i,j,n)
              do
                q=min(pmxl,pres(k+1))-pres(k)
                tdps=tdps+ttem(k)*q
                sdps=sdps+ssal(k)*q
#ifdef TRC
                do nt=1,ntrc
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*q
                enddo
#endif
                delp(k)=pres(k+1)-min(pmxl,pres(k+1))
                if (k.eq.kk.or.pres(k+1).gt.pmxl) exit
                k=k+1
              enddo
            else
c
c --- -----------------------------------------------------------------
c --- ------- If there is a TKE surplus in the mixed layer, increase
c --- ------- the mixed layer depth until the TKE balance is restored.
c --- -----------------------------------------------------------------
c
              tmxl0=ttem(1)
              smxl0=ssal(1)
              dpe0=0.
              tkeo=tke
              k=kfpla(i,j,n)
              tdps=ttem(1)*delp(1)
              sdps=ssal(1)*delp(1)
#ifdef TRC
              do nt=1,ntrc
                trdps(nt)=ttrc(nt,1)*delp(1)
              enddo
#endif
              do
                if (k.gt.kk) then
                  exit
                elseif (delp(k).lt.epsil) then
                  k=k+1
                else
                  if (pres(k)*lbi.lt.1..and.pres(k+1)*lbi.gt.1.) then
                    pmxl=1./lbi
                    dpmxl=min(pmxl-pres(k),pres(k+1)-pmxl,onecm)
                    pmxl=pmxl-.5*dpmxl
                  else
                    dpmxl=min(onecm,delp(k))
                    pmxl=pres(k)+dpmxl
                  endif
                  nitr=0
                  do
                    nitr=nitr+1
                    rm1=exp(-lui*pmxl)
                    rm2=exp(-lbi*pmxl)
                    q=lei*pmxl
                    rm3=exp(-q)
                    rm4=2./q
                    tmxl=(tmxl0*(pres(k)-pres(1))
     .                   +ttem(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    smxl=(smxl0*(pres(k)-pres(1))
     .                   +ssal(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    dpe=dpe0
     .                 +max(.5*alpha0*alpha0*mldjmp
     .                        *(pres(k)-pres(1))*(pmxl-pres(k)),
     .                      p_p_alpha(pmxl,pres(1),tmxl,smxl)
     .                     -p_p_alpha(pmxl,pres(k),ttem(k),ssal(k))
     .                     -p_p_alpha(pres(k),pres(1),tmxl0,smxl0)
     .                     -(pres(1)-pres(k))
     .                      *p_alpha(pmxl,pres(k),ttem(k),ssal(k)))
     .                      *alpha0/(delt1*g)
                    tke=cus*rm1
     .                 +((cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))*rm2
     .                  -crs*pmxl*pmxl)*pmxl
     .                 -dpe
                    if (.not.(nitr.eq.1.and.pres(k)*lbi.lt.1..and.
     .                        pres(k+1)*lbi.gt.1.)) then
                      dtke=(tke-tkeo)/dpmxl
                      if (abs(dtke).lt.(abs(tke)+1.e-16)/delp(k)) then
                        if (tke.lt.0.) then
                          dpmxl=.5*(pres(k)-pmxl)
                        else
                          dpmxl=pres(k+1)-pmxl
                        endif
                      else
                        dpmxl=max(pres(k)-pmxl,
     .                            min(pres(k+1)-pmxl,-tke/dtke))
                      endif
                      dpmxl=max(max(thkmin*onem,pres(k))-pmxl,dpmxl)
                    endif
                    pmxl=pmxl+dpmxl
                    tkeo=tke
                    if (abs(dpmxl).lt.onemm.or.nitr.eq.maxitr) exit
                  enddo
                  if (nitr.eq.maxitr) then
                    write (lp,*) 'reached maxitr when entraining',
     .                           i+i0,j+j0
                    write (lp,*) 'dpth=',pres(2)/onem,';'
                    write (lp,*) 'pmxl=',pmxl/onem,';'
                    write (lp,*) 'corio=',coriop(i,j),';'
                    write (lp,*) 'ustar=',ustar(i,j)*1.e-2,';'
                    write (lp,*) 'bfltot=',bfltot*1.e-4,';'
                    write (lp,*) 'bflpsw=',bflpsw*1.e-4,';'
                    write (lp,*) 'bg2=',util1(i,j),';'
                    write (lp,*) 'ce=',ce*sqrt(scp2(i,j))*rlf,';'
                    write (lp,*) 'pres(2)=',pres(2)/onem,';'
                    do kn=kfpla(i,j,n),k
                      write (lp,*) 'pres(',kn+1,')=',pres(kn+1)/onem,';'
                    enddo
                    write (lp,*) 'ttem(1)=',ttem(1),';'
                    do kn=kfpla(i,j,n),k
                      write (lp,*) 'ttem(',kn,')=',ttem(kn),';'
                    enddo
                    write (lp,*) 'ssal(1)=',ssal(1),';'
                    do kn=kfpla(i,j,n),k
                      write (lp,*) 'ssal(',kn,')=',ssal(kn),';'
                    enddo
c                   call xchalt('(mxlayr)')
c                          stop '(mxlayr)'
                  endif
                  if (pmxl.lt.pres(k+1)-epsil.and.nitr.lt.maxitr) then
                    tdps=tdps+ttem(k)*(pmxl-pres(k))
                    sdps=sdps+ssal(k)*(pmxl-pres(k))
#ifdef TRC
                    do nt=1,ntrc
                      trdps(nt)=trdps(nt)+ttrc(nt,k)*(pmxl-pres(k))
                    enddo
#endif
                    pres(k)=pmxl
                    delp(k)=pres(k+1)-pres(k)
                    exit
                  else
                    tdps=tdps+ttem(k)*delp(k)
                    sdps=sdps+ssal(k)*delp(k)
#ifdef TRC
                    do nt=1,ntrc
                      trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                    enddo
#endif
                    pmxl=pres(k+1)
                    tmxl=(tmxl0*(pres(k)-pres(1))
     .                   +ttem(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    smxl=(smxl0*(pres(k)-pres(1))
     .                   +ssal(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    dpe=dpe0
     .                 +max(.5*alpha0*alpha0*mldjmp
     .                        *(pres(k)-pres(1))*(pmxl-pres(k)),
     .                      p_p_alpha(pmxl,pres(1),tmxl,smxl)
     .                     -p_p_alpha(pmxl,pres(k),ttem(k),ssal(k))
     .                     -p_p_alpha(pres(k),pres(1),tmxl0,smxl0)
     .                     -(pres(1)-pres(k))
     .                      *p_alpha(pmxl,pres(k),ttem(k),ssal(k)))
     .                      *alpha0/(delt1*g)
                    dpe0=dpe
                    tmxl0=tmxl
                    smxl0=smxl
                    pres(k)=pmxl
                    delp(k)=0.
                    k=k+1
                  endif
                endif
              enddo
            endif
c
            pres(2)=min(pres(kk+1),pmxl)
            delp(1)=pres(2)-pres(1)
            q=1./delp(1)
            ttem(1)=tdps*q
            ssal(1)=sdps*q
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,1)=trdps(nt)*q
            enddo
#endif
c
c --- ----- Apply forcing to the water column.
            pswbas=(1.-redfac)*exp(-lei*pres(2))
            pswup=pswbas
            do k=2,kk
              if (pres(k).lt.raddep*onem.and.delp(k).gt.onemu) then
                pswlo=(1.-redfac)*exp(-lei*min(raddep*onem,pres(k+1)))
                ttem(k)=ttem(k)
     .                 -(pswup-pswlo)*sswflx(i,j)*delt1*g
     .                  /(spcifh*delp(k))
                dens(k)=sig(ttem(k),ssal(k))
                pswup=pswlo
              endif
            enddo
            q=delt1*g/delp(1)
            ttem(1)=ttem(1)
     .             -(surflx(i,j)+surrlx(i,j)-(pswbas-pswup)*sswflx(i,j))
     .              *q/spcifh
            ssal(1)=ssal(1)-(salflx(i,j)+salrlx(i,j))*q
            dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,1)=ttrc(nt,1)-trflx(nt,i,j)*q
            enddo
#endif
c
          endif
c
c --- --- Define first physical layer.
          k=2
          dps=0.
          do while (delp(k).lt.epsil)
            dps=dps+delp(k)
            delp(k)=0.
            k=k+1
            if (k.gt.kk) exit
          enddo
          if (k.gt.kk) then
            delp(1)=delp(1)+dps
          else
            delp(k)=delp(k)+dps
          endif
          kfpla(i,j,n)=k
c
c --- --- Put single column back into 3-d fields.
          do k=1,kk
            kn=k+nn
            temp(i,j,kn)=ttem(k)
            saln(i,j,kn)=ssal(k)
            sigma(i,j,kn)=dens(k)
            diaflx(i,j,k)=diaflx(i,j,k)+delp(k)-dp(i,j,kn)
            dp(i,j,kn)=delp(k)
#ifdef TRC
            do nt=1,ntrc 
              trc(i,j,k,nt)=ttrc(nt,k)
            enddo
#endif
          enddo
c
        enddo
        enddo
      enddo
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
 881    pu(i,j,k)=min(pu(i,j,kk+1),.5*(p(i  ,j,k)+p(i-1,j,k)))
c
        do 880 l=1,isv(j)
        do 880 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 880    pv(i,j,k)=min(pv(i,j,kk+1),.5*(p(i,j  ,k)+p(i,j-1,k)))
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
c$OMP PARALLEL DO PRIVATE(kn,q)
      do j=-1,jj+2
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i-1,j,k+1))) 
     .              -min(q,.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i,j-1,k+1)))
     .              -min(q,.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kn,uo,po,pn,ko,un,udpn)
      do j=1,jj
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=u(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pu(i,j,k)
            pn(k)=min(pu(i,j,kk+1),.5*(p(i  ,j,k)+p(i-1,j,k)))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            u(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=v(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pv(i,j,k)
            pn(k)=min(pv(i,j,kk+1),.5*(p(i,j  ,k)+p(i,j-1,k)))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            v(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c     do j=1,jj
c       do l=1,isu(j)
c       do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+u(i,j,kn)*dpu(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'mxlayr: u imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c       do l=1,isv(j)
c       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+v(i,j,kn)*dpv(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'mxlayr: v imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c     enddo
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
