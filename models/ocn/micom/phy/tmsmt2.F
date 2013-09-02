      subroutine tmsmt2(m,n,mm,nn,k1m,k1n)
c
c --- time smoothing of layer thickness, temperature and salinity
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy) :: aux,pbfaco,pbfacn
      integer i,j,k,l,kn,km,kp
      real pold,pmid,pnew,q
c
#ifdef TKE
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr) :: trbudo
      real trbudn
      integer nt
#  ifdef ATRC
      integer nat
#  endif
#endif
c
c$OMP PARALLEL DO PRIVATE(aux,pbfaco,pbfacn,km,kn,kp,pold,pmid,pnew,q)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          aux(i)=0.
          pbfaco(i)=0.
          pbfacn(i)=0.
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbfaco(i)=pbfaco(i)+dpold(i,j,kn)
            pbfacn(i)=pbfacn(i)+dp(i,j,kn)
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i)=pb(i,j,m)/pbfaco(i)
          pbfacn(i)=pb(i,j,m)/pbfacn(i)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          kn=k+nn
          kp=min(k+1,kk)
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pold=max(0.,dpold(i,j,kn)*pbfaco(i))
            pmid=max(0.,dp(i,j,km))
            pnew=max(0.,dp(i,j,kn)*pbfacn(i))
            dp(i,j,km)=wts1*pmid+wts2*(pold+pnew)
            aux(i)=aux(i)+(dp(i,j,km)-pmid)			! diapyc.flux
            diaflx(i,j,k )=diaflx(i,j,k )+aux(i)		! diapyc.flux
            diaflx(i,j,kp)=diaflx(i,j,kp)-aux(i)		! diapyc.flux
            pold=pold+epsil
            pmid=pmid+epsil
            pnew=pnew+epsil
            temp(i,j,km)=(wts1*pmid*temp(i,j,km)
     .                   +wts2*(pold*told(i,j,k)+pnew*temp(i,j,kn)))
     .                   /(dp(i,j,km)+epsil)
            saln(i,j,km)=(wts1*pmid*saln(i,j,km)
     .                   +wts2*(pold*sold(i,j,k)+pnew*saln(i,j,kn)))
     .                   /(dp(i,j,km)+epsil)
#ifdef TKE
            tke(i,j,km)=(wts1*pmid*tke(i,j,km)
     .                  +wts2*(pold*tke_old(i,j,k)+pnew*tke(i,j,kn)))
     .                  /(dp(i,j,km)+epsil)
#  ifdef GLS
            gls_psi(i,j,km)=(wts1*pmid*gls_psi(i,j,km)
     .                      +wts2*(pold*gls_psi_old(i,j,k)
     .                            +pnew*gls_psi(i,j,kn)))
     .                      /(dp(i,j,km)+epsil)
#  endif
#endif
#ifdef TRC
            if (n.eq.2) then
              q=pmid/(dp(i,j,km)+epsil)
#  ifdef ATRC
              do nt=1,ntr-ntriag-natr
#  else
              do nt=1,ntr-ntriag
#  endif
                trc(i,j,k,nt)=trc(i,j,k,nt)*q
              enddo
            endif
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(dp(1-nbdy,1-nbdy,k1m), 1,kk, 3,3, halo_ps)
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
c$OMP PARALLEL DO PRIVATE(km,q)
      do j=-1,jj+2
        do k=1,kk
          km=k+mm
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,km)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i-1,j,k+1))) 
     .              -min(q,.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,km)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i,j-1,k+1)))
     .              -min(q,.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt2:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(diaflx,ip,kk,'diaflx')
        call chksummsk(dpu,iu,2*kk,'dpu')
        call chksummsk(dpv,iv,2*kk,'dpv')
#ifdef TKE
        call chksummsk(tke,ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(gls_psi,ip,2*kk,'gls_psi')
#  endif
#endif
      endif
c
      return
      end
