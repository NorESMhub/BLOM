      subroutine mod2ice(mo,no,k1m,k1n,ztx,mty,uocbtp,vocbtp,uoc,voc,
     .                   sealv,scuxi,scvyi,hicem,ficem,onecm)
c
      use mod_xc, only: nbdy,idm,jdm,kdm,itdm,jtdm,xctilr,halo_uv,
     .                  halo_vv,ii,jj,ip,iu,iv,xcaget,mnproc,halo_ps
c
      implicit none
c
      integer mo,no,k1m,k1n
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: uoc,voc
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5) ::
     .  uocbtp,vocbtp
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ztx,mty,sealv,scuxi,scvyi,hicem,ficem
      real onecm
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp1,tmp2
      real, dimension(itdm,jtdm) :: tmp2d
      integer i,j
c
c --- ------------------------------------------------------------------------
c --- set input fields to dynamic ice model
c -o-   tauxw           x-component of surf wind stress [N/m**2]
c -o-   tauyw           y-component of surf wind stress [N/m**2]
c -o-   uwat            x-component of ml current velocity [m/s]
c -o-   vwat            y-component of ml current velocity [m/s]
c -o-   gradhx          x-component of ocean surface tilt [dimless]
c -o-   gradhy          y-component of ocean surface tilt [dimless]
c ---
c --- n o t e :  dynamic ocean model is formulated in c-grid, whereas
c ---            dynamic ice model is in b-grid
c --- ------------------------------------------------------------------------
      call xctilr(ztx, 1,1, 1,1, halo_uv)
      call xctilr(mty, 1,1, 1,1, halo_vv)
c$OMP PARALLEL DO
      do j=1,jj
        do i=1,ii
          if (iu(i,j).eq.1) then
            tmp1(i,j)=ztx(i,j)
          else
            tmp1(i,j)=ztx(i-1,j)*iu(i-1,j)+ztx(i+1,j)*iu(i+1,j)
          endif
          if (iv(i,j).eq.1) then
            tmp2(i,j)=mty(i,j)
          else
            tmp2(i,j)=mty(i,j-1)*iv(i,j-1)+mty(i,j+1)*iv(i,j+1)
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2d,tmp1,1)
      if (mnproc.eq.1) then
        call intp_mu2iq(itdm,jtdm,vm,tmp2d,tauxw)
      endif
      call xcaget(tmp2d,tmp2,1)
      if (mnproc.eq.1) then
        call intp_mv2iq(itdm,jtdm,vm,tmp2d,tauyw)
      endif
c
c$OMP PARALLEL DO
      do j=1,jj
        do i=1,ii
          if (iu(i,j).eq.1) then
            tmp1(i,j)=.5*(uocbtp(i,j,mo)+uoc(i,j,k1m)
     .                   +uocbtp(i,j,no)+uoc(i,j,k1n))*1.e-2
          else
            tmp1(i,j)=0.
          endif
          if (iv(i,j).eq.1) then
            tmp2(i,j)=.5*(vocbtp(i,j,mo)+voc(i,j,k1m)
     .                   +vocbtp(i,j,no)+voc(i,j,k1n))*1.e-2
          else
            tmp2(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2d,tmp1,1)
      if (mnproc.eq.1) then
        call intp_mu2iq(itdm,jtdm,vm,tmp2d,uwat2l(1,0,lnew))
      endif
      call xcaget(tmp2d,tmp2,1)
      if (mnproc.eq.1) then
        call intp_mv2iq(itdm,jtdm,vm,tmp2d,vwat2l(1,0,lnew))
      endif
c
      call xctilr(sealv, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
      do j=1,jj
        do i=1,ii
          if (iu(i,j).eq.1) then
            tmp1(i,j)=(sealv(i,j)-sealv(i-1,j))*scuxi(i,j)
          endif
          if (iv(i,j).eq.1) then
            tmp2(i,j)=(sealv(i,j)-sealv(i,j-1))*scvyi(i,j)
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(tmp1, 1,1, 1,1, halo_uv)
      call xctilr(tmp2, 1,1, 1,1, halo_vv)
c$OMP PARALLEL DO
      do j=1,jdm
        do i=1,idm
          if (iu(i,j).eq.0) then
            tmp1(i,j)=tmp1(i-1,j)*iu(i-1,j)+tmp1(i+1,j)*iu(i+1,j)
          endif
          if (iv(i,j).eq.0) then
            tmp2(i,j)=tmp2(i,j-1)*iv(i,j-1)+tmp2(i,j+1)*iv(i,j+1)
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2d,tmp1,1)
      if (mnproc.eq.1) then
        call intp_mu2iq(itdm,jtdm,vm,tmp2d,gradhx2l(1,0,lnew))
      endif
      call xcaget(tmp2d,tmp2,1)
      if (mnproc.eq.1) then
        call intp_mv2iq(itdm,jtdm,vm,tmp2d,gradhy2l(1,0,lnew))
      endif
c
c --- ------------------------------------------------------------------------
c --- set ice fields updated by  -thermf- and advection
c --- ------------------------------------------------------------------------
c
c --- note that the dynamic ice code need mean ice thickness, while the
c --- ocean model operates with thickness of the ice covered fraction
c$OMP PARALLEL DO
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.1) then
            tmp1(i,j)=hicem(i,j)*ficem(i,j)
          else
            tmp1(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2d,tmp1,1)
      if (mnproc.eq.1) then
        call intp_mp2ip(itdm,jtdm,hm,tmp2d,h)
      endif
      call xcaget(tmp2d,ficem,1)
      if (mnproc.eq.1) then
        call intp_mp2ip(itdm,jtdm,hm,tmp2d,a)
      endif
c
      return
      end
