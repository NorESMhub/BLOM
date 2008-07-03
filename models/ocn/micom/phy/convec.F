      subroutine convec(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- remove static instabilitites from the water column
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_dia.h"
c
      real, dimension(kdm) :: ttem,ssal,delp,dens
      real llim,ttmp,stmp,dtmp,ttdcy,stdcy,q,aa,bb,cc,dd,ee,ff,
     .     lambda,tplume,splume,f,thtrgt,fq,alpha,beta,thtdcy,dthk1,udpn
      integer i,j,k,l,kn,ko
      logical done
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      real, dimension(ntr) :: trtmp
      real trplum
      integer nt
#endif
c
      real sofsig,dsigdt,dsigds,sig,tsroot
      external sofsig,dsigdt,dsigds,sig,tsroot
c
      real heat0,salt0,heat,salt
c
      real frcday,tlim,slim
      parameter (frcday=90.,tlim=5.,slim=1.)
c
      llim=frcday*86400./delt1
c
c$OMP PARALLEL DO PRIVATE(ttem,ssal,delp,dens,ttmp,stmp,dtmp,ttdcy,
c$OMP+ stdcy,q,aa,bb,cc,dd,ee,ff,lambda,tplume,splume,f,thtrgt,
c$OMP+ fq,alpha,beta,thtdcy,dthk1,done,k,kn,heat0,salt0,heat,salt
#ifdef TRC
c$OMP+,ttrc,trtmp,trplum
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
cdiag     heat0=0.
cdiag     salt0=0.
cdiag     do k=1,kk
cdiag       kn=k+nn
cdiag       heat0=heat0+temp(i,j,kn)*dp(i,j,kn)
cdiag       salt0=salt0+saln(i,j,kn)*dp(i,j,kn)
cdiag     enddo
c
c --- --- copy variables into 1d arrays
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            dens(k)=sigma(i,j,kn)
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,k)=trc(i,j,k,nt)
            enddo
#endif
          enddo
c
c --- ------------------------------------------------------------------
c --- --- remove static instabilities in interior layers
c --- ------------------------------------------------------------------
c
          done=.false.
c
          do while (.not.done)
            done=.true.
            do k=2,kk-1
              if (dens(k).gt.dens(k+1)) then
                done=.false.
                if (delp(k)+delp(k+1).gt.epsil) then
                  if (min(delp(k),delp(k+1)).gt.onem) then
                    write (lp,*) 'static instability',i0+i,j0+j,k
                    write (lp,*) k,delp(k)/onem,ttem(k),ssal(k),dens(k)
                    write (lp,*) k+1,delp(k+1)/onem,ttem(k+1),ssal(k+1),
     .                          dens(k+1)
                  endif
                  q=1./(delp(k)+delp(k+1))
                  ttmp=(ttem(k)*delp(k)+ttem(k+1)*delp(k+1))*q
                  stmp=(ssal(k)*delp(k)+ssal(k+1)*delp(k+1))*q
#ifdef TRC
                  do nt=1,ntrc
                    trtmp(nt)=(ttrc(nt,k  )*delp(k  )
     .                        +ttrc(nt,k+1)*delp(k+1))*q
                  enddo
#endif
                  dtmp=sig(ttmp,stmp)
                  if (abs(dtmp-sigmar(k)).lt.abs(dtmp-sigmar(k+1))) then
                    ttem(k)=ttmp
                    ssal(k)=stmp
                    dens(k)=dtmp
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=trtmp(nt)
                    enddo
#endif
                    delp(k)=delp(k)+delp(k+1)
                    dens(k+1)=sigmar(k+1)
                    ssal(k+1)=sofsig(dens(k+1),ttem(k+1))
                    delp(k+1)=0.
                  else
                    ttem(k+1)=ttmp
                    ssal(k+1)=stmp
                    dens(k+1)=dtmp
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k+1)=trtmp(nt)
                    enddo
#endif
                    delp(k+1)=delp(k)+delp(k+1)
                    dens(k)=sigmar(k)
                    ssal(k)=sofsig(dens(k),ttem(k))
                    delp(k)=0.
                  endif
                else
                  dens(k)=sigmar(k)
                  ssal(k)=sofsig(dens(k),ttem(k))
                  dens(k+1)=sigmar(k+1)
                  ssal(k+1)=sofsig(dens(k+1),ttem(k+1))
                endif
              endif
            enddo
          enddo
c
c --- ------------------------------------------------------------------
c --- --- make sure the mixed layer is lighter than physical isopycnic
c --- --- layers
c --- ------------------------------------------------------------------
c
          if (dens(2)-dens(1)-sigjmp.ge.0.) then
            klist(i,j)=2
          else
c
c --- ----- Mixed layer temperature and salinity tendencies
            ttdcy=ttem(1)-tmxold(i,j,n)
            stdcy=ssal(1)-smxold(i,j,n)
c
            k=2
            do while (dens(min(k+1,kk))-dens(1)-sigjmp.lt.0..and.
     .                k.lt.kk)
              if (delp(k).gt.0.) then
                q=1./(delp(1)+delp(k))
                ttem(1)=(ttem(1)*delp(1)+ttem(k)*delp(k))*q
                ssal(1)=(ssal(1)*delp(1)+ssal(k)*delp(k))*q
                dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,1)=(ttrc(nt,1)*delp(1)+ttrc(nt,k)*delp(k))*q
                enddo
#endif
                convec_dia(i,j)=convec_dia(i,j)+delp(k)
                delp(1)=delp(1)+delp(k)
                delp(k)=0.
                k=k-1
              endif
              k=k+1
            enddo
            klist(i,j)=k+1
c
            alpha=-alpha0*dsigdt(ttem(1),ssal(1))
            beta=alpha0*dsigds(ttem(1),ssal(1))
c
c --- ----- Density tendency of mixed layer
            thtdcy=-alpha*ttdcy+beta*stdcy
            if (thtdcy.lt.epsil.and.k.lt.kk) then
              ttdcy=ttem(k+1)-ttem(1)
              stdcy=ssal(k+1)-ssal(1)
              thtdcy=-alpha*ttdcy+beta*stdcy
            endif
c
            if (k.eq.kk) then
              if (delp(k).gt.0.) then
                q=1./(delp(1)+delp(k))
                ttem(1)=(ttem(1)*delp(1)+ttem(k)*delp(k))*q
                ssal(1)=(ssal(1)*delp(1)+ssal(k)*delp(k))*q
                dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,1)=(ttrc(nt,1)*delp(1)+ttrc(nt,k)*delp(k))*q
                enddo
#endif
                convec_dia(i,j)=convec_dia(i,j)+delp(k)
                delp(1)=delp(1)+delp(k)
                delp(k)=0.
              endif
            elseif (delp(k).gt.0..and.
     .              (delp(k).le.onem.or.thtdcy.lt.epsil)) then
c
c --- ------- If the light layer below the mixed layer is thin or if the
c --- ------- density tendency is not computed to be large enough,
c --- ------- remove layer k by moving the the interface k down
              q=1./(delp(1)+delp(k))
              ttem(1)=(ttem(1)*delp(1)+ttem(k)*delp(k))*q
              ssal(1)=(ssal(1)*delp(1)+ssal(k)*delp(k))*q
              dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,1)=(ttrc(nt,1)*delp(1)+ttrc(nt,k)*delp(k))*q
              enddo
#endif
              convec_dia(i,j)=convec_dia(i,j)+delp(k)
              delp(1)=delp(1)+delp(k)
              delp(k)=0.
c
            elseif (delp(k).gt.onem) then
c
c --- ------- Target density for the mixed layer
              thtrgt=dens(k)-sigjmp-1.e-8
c
c --- ------- Estimate the plume fraction
              f=(dens(1)-thtrgt)/(dens(k+1)-thtrgt)
c
              if (f*delp(1).lt.delp(k).and.
     .            abs(ttdcy)+abs(stdcy).gt.epsil) then
c
c --- --------- Mix the fraction f between mixed layer and layer k
                fq=f*delp(1)/delp(k)
                aa=fq*ttdcy
                bb=ttem(k)+fq*(ttem(1)-ttem(k))
                cc=0.
                dd=1.
                ee=fq*stdcy
                ff=ssal(k)+fq*(ssal(1)-ssal(k))
                lambda=tsroot(aa,bb,cc,dd,ee,ff,dens(k),0.,llim)
                tplume=ttem(1)+lambda*ttdcy
                splume=ssal(1)+lambda*stdcy
                ttmp=fq*tplume+(1.-fq)*ttem(k)
                if (abs(lambda).lt.llim.and.
     .              abs(tplume-ttem(1)).lt.tlim.and.
     .              abs(splume-ssal(1)).lt.slim.and.
     .              ttmp.gt.temmin(k)) then
                  ttem(1)=ttem(1)+f*(ttem(k)-tplume)
                  ssal(1)=ssal(1)+f*(ssal(k)-splume)
                  dens(1)=sig(ttem(1),ssal(1))
                  ttem(k)=ttmp
                  ssal(k)=fq*splume+(1.-fq)*ssal(k)
                  dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                  do nt=1,ntrc
                    trplum=ttrc(nt,1)
                    ttrc(nt,1)=ttrc(nt,1)+f*(ttrc(nt,k)-trplum)
                    ttrc(nt,k)=fq*trplum+(1.-fq)*ttrc(nt,k)
                  enddo
#endif
                endif
              endif
c
              alpha=-alpha0*dsigdt(ttem(1),ssal(1))
              beta=alpha0*dsigds(ttem(1),ssal(1))
c
c --- ------- Update density tendency of mixed layer
              thtdcy=-alpha*ttdcy+beta*stdcy
              if (thtdcy.lt.epsil) then
                ttdcy=ttem(k+1)-ttem(1)
                stdcy=ssal(k+1)-ssal(1)
                thtdcy=-alpha*ttdcy+beta*stdcy
              endif
c
c --- ------- Compute plume fraction
              dthk1=(-alpha*(ttem(k)-ttem(1))+beta*(ssal(k)-ssal(1)))
              aa=delp(1)*(thtdcy*ttem(1)+dthk1*ttdcy)
              bb=delp(k+1)*thtdcy*ttem(k+1)
     .          +delp(1)*(dens(1)-thtrgt)*ttdcy
              cc=delp(1)*thtdcy
              dd=delp(k+1)*thtdcy
              ee=delp(1)*(thtdcy*ssal(1)+dthk1*stdcy)
              ff=delp(k+1)*thtdcy*ssal(k+1)
     .          +delp(1)*(dens(1)-thtrgt)*stdcy
              f=tsroot(aa,bb,cc,dd,ee,ff,dens(k+1),0.,0.)
              lambda=(dens(1)-thtrgt+dthk1*f)/(thtdcy*max(epsil,f))
              tplume=ttem(1)+lambda*ttdcy
              splume=ssal(1)+lambda*stdcy
c
              if (f.le.epsil.or.f.ge.1..or.
     .            abs(lambda).gt.llim.or.abs(tplume-ttem(1)).gt.tlim.or.
     .            abs(splume-ssal(1)).gt.slim) then
                tplume=ttem(k+1)
                splume=ssal(k+1)
                f=(dens(1)-thtrgt)
     .            /(-alpha*(tplume-ttem(k))+beta*(splume-ssal(k)))
              endif
c
              if (f.le.0..or.f.ge.1.) then
c
c --- --------- If the plume fraction is out of range, let the mixed
c --- --------- layer absorb layer k
                q=1./(delp(1)+delp(k))
                ttem(1)=(delp(1)*ttem(1)+delp(k)*ttem(k))*q
                ssal(1)=(delp(1)*ssal(1)+delp(k)*ssal(k))*q
                dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,1)=(delp(1)*ttrc(nt,1)+delp(k)*ttrc(nt,k))*q
                enddo
#endif
                convec_dia(i,j)=convec_dia(i,j)+delp(k)
                delp(1)=delp(1)+delp(k)
                delp(k)=0.
              elseif (f*delp(1).ge.delp(k)) then
c
c --- --------- Limit the plume volume so that it is not greater than
c --- --------- the volume of layer k
                q=1./(delp(k)+delp(k+1))
                aa=delp(k)*ttdcy*q
                bb=(delp(k)*ttem(1)+delp(k+1)*ttem(k+1))*q
                cc=0.
                dd=1.
                ee=delp(k)*stdcy*q
                ff=(delp(k)*ssal(1)+delp(k+1)*ssal(k+1))*q
                lambda=tsroot(aa,bb,cc,dd,ee,ff,dens(k+1),0.,llim)
                tplume=ttem(1)+lambda*ttdcy
                splume=ssal(1)+lambda*stdcy
                ttmp=(delp(k+1)*ttem(k+1)+delp(k)*tplume)*q
                if (abs(lambda).ge.llim.or.
     .              abs(tplume-ttem(1)).ge.tlim.or.
     .              abs(splume-ssal(1)).ge.slim.or.
     .              ttmp.lt.temmin(k+1)) then
                  q=delp(k)/delp(1)
                  ttem(1)=ttem(1)+q*(ttem(k)-ttem(k+1))
                  ssal(1)=ssal(1)+q*(ssal(k)-ssal(k+1))
                  dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                  do nt=1,ntrc
                    ttrc(nt,1)=ttrc(nt,1)+q*(ttrc(nt,k)-ttrc(nt,k+1))
                  enddo
#endif
                  convec_dia(i,j)=convec_dia(i,j)+delp(k)
                  delp(k+1)=delp(k+1)+delp(k)
                  delp(k)=0.
                else
                  f=delp(k)/delp(1)
                  ttem(1)=ttem(1)+f*(ttem(k)-tplume)
                  ssal(1)=ssal(1)+f*(ssal(k)-splume)
                  dens(1)=sig(ttem(1),ssal(1))
                  ttem(k+1)=ttmp
                  ssal(k+1)=(delp(k+1)*ssal(k+1)+delp(k)*splume)*q
                  dens(k+1)=sig(ttem(k+1),ssal(k+1))
#ifdef TRC
                  do nt=1,ntrc
                    trplum=ttrc(nt,1)
                    ttrc(nt,1)=ttrc(nt,1)+f*(ttrc(nt,k)-trplum)
                    ttrc(nt,k+1)=(delp(k+1)*ttrc(nt,k+1)
     .                           +delp(k)*trplum)*q
                  enddo
#endif
                  convec_dia(i,j)=convec_dia(i,j)+delp(k)
                  delp(k+1)=delp(k+1)+delp(k)
                  delp(k)=0.
                endif
              else
                q=1./(f*delp(1)+delp(k+1))
                ttmp=(f*delp(1)*tplume+delp(k+1)*ttem(k+1))*q
                if (ttmp.lt.temmin(k+1)) then
                  q=1./(delp(1)+delp(k))
                  ttem(1)=(ttem(1)*delp(1)+ttem(k)*delp(k))*q
                  ssal(1)=(ssal(1)*delp(1)+ssal(k)*delp(k))*q
                  dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                  do nt=1,ntrc
                    ttrc(nt,1)=(ttrc(nt,1)*delp(1)+ttrc(nt,k)*delp(k))*q
                  enddo
#endif
                  convec_dia(i,j)=convec_dia(i,j)+delp(k)
                  delp(1)=delp(1)+delp(k)
                  delp(k)=0.
                else
                  ttem(1)=ttem(1)+f*(ttem(k)-tplume)
                  ssal(1)=ssal(1)+f*(ssal(k)-splume)
                  dens(1)=sig(ttem(1),ssal(1))
                  ttem(k+1)=ttmp
                  ssal(k+1)=(f*delp(1)*splume+delp(k+1)*ssal(k+1))*q
                  dens(k+1)=sig(ttem(k+1),ssal(k+1))
#ifdef TRC
                  do nt=1,ntrc
                    trplum=ttrc(nt,1)
                    ttrc(nt,1)=ttrc(nt,1)+f*(ttrc(nt,k)-trplum)
                    ttrc(nt,k+1)=(f*delp(1)*trplum
     .                           +delp(k+1)*ttrc(nt,k+1))*q
                  enddo
#endif
                  convec_dia(i,j)=convec_dia(i,j)+f*delp(1)
                  delp(k)=delp(k)-f*delp(1)
                  delp(k+1)=delp(k+1)+f*delp(1)
                  klist(i,j)=k
                  if (sig(ttem(1),ssal(1)).gt.dens(k)-sigjmp) then
                    write(lp,*) 'density miss',i0+i,j0+j,f,
     .                (dens(k)-sigjmp)*1.e3,sig(ttem(1),ssal(1))*1.e3
                  endif
                endif
              endif
c
            endif
c
          endif
c
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
cdiag     heat=0.
cdiag     salt=0.
cdiag     do k=1,kk
cdiag       kn=k+nn
cdiag       heat=heat+temp(i,j,kn)*dp(i,j,kn)
cdiag       salt=salt+saln(i,j,kn)*dp(i,j,kn)
cdiag     enddo
cdiag     if (abs((heat-heat0)/pb(i,j,n)).gt.1.e-11) then
cdiag       write (*,*) 'heat',i,j
cdiag       write (*,*) heat0/pb(i,j,n),(heat-heat0)/pb(i,j,n)
cdiag       write (*,*) salt0/pb(i,j,n),(salt-salt0)/pb(i,j,n)
cdiag       write (*,*)
cdiag       do k=1,kk
cdiag         kn=k+nn
cdiag         write (*,*) k,dp(i,j,kn)/onem,temp(i,j,kn),saln(i,j,kn)
cdiag       enddo
cdiag       write (*,*)
cdiag       stop
cdiag     endif
cdiag     if (abs((salt-salt0)/pb(i,j,n)).gt.1.e-11) then
cdiag       write (*,*) 'salt',i,j
cdiag       write (*,*) heat0/pb(i,j,n),(heat-heat0)/pb(i,j,n)
cdiag       write (*,*) salt0/pb(i,j,n),(salt-salt0)/pb(i,j,n)
cdiag       write (*,*)
cdiag       do k=1,kk
cdiag         kn=k+nn
cdiag         write (*,*) k,dp(i,j,kn)/onem,temp(i,j,kn),saln(i,j,kn)
cdiag       enddo
cdiag       stop
cdiag     endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tmxold(i,j,n)=temp(i,j,k1n)
          smxold(i,j,n)=saln(i,j,k1n)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'convec:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(diaflx,ip,kk,'diaflx')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
        call chksummsk(tmxold,ip,1,'tmxold')
        call chksummsk(smxold,ip,1,'smxold')
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
c
