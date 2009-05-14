      subroutine convec(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Remove static instabilitites from the water column
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
c --- parameters:
c ---   dnseps - minimum density difference between interior layers
c ---            [g/cm^3]
c
      real dnseps
      parameter (dnseps=1.e-12)
c
      real, dimension(kdm) :: ttem,ssal,delp,dens,densr
      real tdps,sdps,dps,ttmp,stmp,dtmp,q
      integer i,j,k,l,kn,kfpl,kfplo,kmix
      logical done
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      real, dimension(ntr) :: trdps
      integer nt
#endif
c
      real sofsig,sig,rho
      external sofsig,sig,rho
c
      integer niter
c
c$OMP PARALLEL DO PRIVATE(ttem,ssal,delp,dens,densr,
c$OMP+ done,k,kn,tdps,sdps,dps,ttmp,stmp,dtmp,q,kfpl,kfplo,kmix,niter
#ifdef TRC
c$OMP+,ttrc,trdps
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Copy variables into 1d arrays
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            dens(k)=sigma(i,j,kn)
            densr(k)=sigmar(i,j,k)
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,k)=trc(i,j,k,nt)
            enddo
#endif
          enddo
c
c         if (mnproc.eq.ptest.and.i.eq.itest.and.j.eq.jtest) then
c           write (lp,*) 'convec init'
c           write (lp,*) 'kfpl',kfpla(i,j,n)
c           do k=1,kk
c             write (lp,'(i3,6e12.5)')
c    .          k,delp(k)/onem,ttem(k),ssal(k),dens(k),densr(k),
c    .          rho(0.,ttem(k),ssal(k))-1.
c           enddo
c           write (lp,*)
c         endif
c
c --- ------------------------------------------------------------------
c --- --- Define first physical interior layer
c --- ------------------------------------------------------------------
c
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
          kfpl=k
          kfplo=kfpla(i,j,n)
          if (kfpl.lt.kfplo) then
            if (kfplo.le.kk) then
              tdps=0.
              sdps=0.
              dps=0.
#ifdef TRC
              do nt=1,ntrc
                trdps(nt)=0.
              enddo
#endif
              do k=kfpl,kfplo
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
#ifdef TRC
                do nt=1,ntrc
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
              enddo
              q=1./dps
              ttmp=tdps*q
              stmp=sdps*q
              dtmp=sig(ttmp,stmp)
              if (dtmp.gt.densr(kfplo)) then
                do k=kfpl,kfplo-1
                  delp(k)=0.
                enddo
                kfpl=kfplo
                ttem(kfpl)=ttmp
                ssal(kfpl)=stmp
                dens(kfpl)=dtmp
                delp(kfpl)=dps
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,kfpl)=trdps(nt)*q
                enddo
#endif
              endif
            else
              tdps=0.
              sdps=0.
              dps=0.
#ifdef TRC
              do nt=1,ntrc
                trdps(nt)=0.
              enddo
#endif
              do k=kfpl,kk
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
#ifdef TRC
                do nt=1,ntrc
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
                delp(k)=0.
              enddo
              q=1./dps
              ttmp=tdps*q
              stmp=sdps*q
              dtmp=sig(ttmp,stmp)
              kfpl=kk
              do while (dtmp.lt.densr(kfpl))
                if (kfpl.eq.2) exit
                kfpl=kfpl-1
              enddo
              ttem(kfpl)=ttmp
              ssal(kfpl)=stmp
              dens(kfpl)=dtmp
              delp(kfpl)=dps
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,kfpl)=trdps(nt)*q
              enddo
#endif
            endif
          endif
c
          if (kfpl.le.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Remove static instabilities
c --- ------------------------------------------------------------------
c
            done=.false.
c
            niter=0
            do while (.not.done)
              niter=niter+1
              if (niter.eq.100) then
                write (lp,*) 'convec: no convergence!',i+i0,j+j0
                exit
c               call xchalt('(convec)')
c                      stop('(convec)')
              endif
c
              done=.true.
c
c --- ------- Remove instabilities between the mixed layer and interior
c --- ------- layers by considering potential density with reference
c --- ------- pressure at the surface
              tdps=ttem(1)*delp(1)
              sdps=ssal(1)*delp(1)
              dps=delp(1)
#ifdef TRC
              do nt=1,ntrc
                trdps(nt)=ttrc(nt,1)*delp(1)
              enddo
#endif
              ttmp=ttem(1)
              stmp=ssal(1)
              k=kfpl
              do while (rho(dps,ttmp,stmp).gt.
     .                  rho(dps,ttem(k),ssal(k)).or.
     .                  delp(k).lt.epsil)
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
                q=1./dps
                ttmp=tdps*q
                stmp=sdps*q
#ifdef TRC
                do nt=1,ntrc
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
                k=k+1
                if (k.gt.kk) exit
              enddo
              kmix=k-1
              if (kmix.ge.kfpl) then
                ttem(1)=ttmp
                ssal(1)=stmp
                dens(1)=sig(ttem(1),ssal(1))
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,1)=trdps(nt)*q
                enddo
#endif
                dps=0.
                do k=kfpl,kmix
                  dps=dps+delp(k)
                  delp(k)=0.
                enddo
                k=kmix
                do while (dens(1).lt.densr(k))
                  if (k.eq.2) exit
                  k=k-1
                enddo
                kfpl=k
                ttem(kfpl)=ttem(1)
                ssal(kfpl)=ssal(1)
                dens(kfpl)=dens(1)
                delp(kfpl)=dps
#ifdef TRC
                do nt=1,ntrc
                  ttrc(nt,kfpl)=ttrc(nt,1)
                enddo
#endif
                do k=kfpl+1,kmix
                  ttem(k)=ttem(1)
                  dens(k)=densr(k)
                  ssal(k)=sofsig(dens(k),ttem(k))
                enddo
              endif
c
c             if (mnproc.eq.ptest.and.i.eq.itest.and.j.eq.jtest) then
c               write (lp,*) 'convec niter',niter
c               write (lp,*) 'kfpl',kfpl
c               do k=1,kk
c                 write (lp,'(i3,6e12.5)')
c    .              k,delp(k)/onem,ttem(k),ssal(k),dens(k),densr(k),
c    .              sig0(ttem(k),ssal(k))
c               enddo
c               write (lp,*)
c             endif
c
c --- ------- Remove instabilities between layers below surface mixed
c --- ------- layers
              do k=kfpl+1,kk-1
                if (min(densr(k+1),dens(k+1))
     .             -max(densr(k  ),dens(k  )).lt.dnseps) then
                  done=.false.
                  if (delp(k).gt.epsil.or.delp(k+1).gt.epsil) then
                    q=1./(delp(k)+delp(k+1))
                    ttmp=(ttem(k)*delp(k)+ttem(k+1)*delp(k+1))*q
                    stmp=(ssal(k)*delp(k)+ssal(k+1)*delp(k+1))*q
                    dtmp=sig(ttmp,stmp)
                    if (abs(dtmp-densr(k)).lt.abs(dtmp-densr(k+1))) then
                      ttem(k)=ttmp
                      ssal(k)=stmp
                      dens(k)=dtmp
#ifdef TRC
                      do nt=1,ntrc
                        ttrc(nt,k)=(ttrc(nt,k  )*delp(k  )
     .                             +ttrc(nt,k+1)*delp(k+1))*q
                      enddo
#endif
                      delp(k)=delp(k)+delp(k+1)
                      dens(k+1)=densr(k+1)
                      ssal(k+1)=sofsig(dens(k+1),ttem(k+1))
                      delp(k+1)=0.
                    else
                      ttem(k+1)=ttmp
                      ssal(k+1)=stmp
                      dens(k+1)=dtmp
#ifdef TRC
                      do nt=1,ntrc
                        ttrc(nt,k+1)=(ttrc(nt,k  )*delp(k  )
     .                             +ttrc(nt,k+1)*delp(k+1))*q
                      enddo
#endif
                      delp(k+1)=delp(k)+delp(k+1)
                      dens(k)=densr(k)
                      ssal(k)=sofsig(dens(k),ttem(k))
                      delp(k)=0.
                    endif
                  else
                    dens(k)=densr(k)
                    ssal(k)=sofsig(dens(k),ttem(k))
                    dens(k+1)=densr(k+1)
                    ssal(k+1)=sofsig(dens(k+1),ttem(k+1))
                  endif
                endif
              enddo
c
            enddo
c
c           if (mnproc.eq.ptest.and.i.eq.itest.and.j.eq.jtest) then
c             write (lp,*) 'convec final'
c             write (lp,*) 'kfpl',kfpl
c             do k=1,kk
c               write (lp,'(i3,6e12.5)')
c    .            k,delp(k)/onem,ttem(k),ssal(k),dens(k),densr(k),
c    .            rho(0.,ttem(k),ssal(k))-1.
c             enddo
c             write (lp,*)
c           endif
c
          endif
c
          kfpla(i,j,n)=kfpl
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
