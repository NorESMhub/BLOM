      subroutine iadvdf
c
c --- ------------------------------------------------------------------
c --- Do horizontal advection and diffusion of sea ice scalar fields
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_ice.h"
c
      integer i,j,l
c
c --- obtain volume of ice and snow
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          hicem(i,j)=hicem(i,j)*ficem(i,j)
          hsnwm(i,j)=hsnwm(i,j)*ficem(i,j)
          iagem(i,j)=iagem(i,j)*hicem(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- advect and diffuse
      call iceadv(ficem,uicem,vicem,scuy,scvx,scp2i,scp2,baclin)
      call icedif(ficem,scuy,scvx,scuxi,scvyi,scp2i,baclin,cdiff,cdfac)
c
      call iceadv(hicem,uicem,vicem,scuy,scvx,scp2i,scp2,baclin)
      call icedif(hicem,scuy,scvx,scuxi,scvyi,scp2i,baclin,cdiff,cdfac)
c
      call iceadv(hsnwm,uicem,vicem,scuy,scvx,scp2i,scp2,baclin)
      call icedif(hsnwm,scuy,scvx,scuxi,scvyi,scp2i,baclin,cdiff,cdfac)
c
      call iceadv(iagem,uicem,vicem,scuy,scvx,scp2i,scp2,baclin)
      call icedif(iagem,scuy,scvx,scuxi,scvyi,scp2i,baclin,cdiff,cdfac)
c
c --- get scalar fields from updated volume estimates and limit ice
c --- concentration in a volume conservative manner
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (hicem(i,j).lt.1.e-9) then
            iagem(i,j)=0.
          else
            iagem(i,j)=max(0.,iagem(i,j)/hicem(i,j))
          endif
          if (ficem(i,j).lt.1.e-9) then
            ficem(i,j)=0.
            hicem(i,j)=0.
            hsnwm(i,j)=0.
          else
            ficem(i,j)=max(1.e-4,min(fice_max,ficem(i,j)))
            hicem(i,j)=max(0.,hicem(i,j)/ficem(i,j))
            hsnwm(i,j)=max(0.,hsnwm(i,j)/ficem(i,j))
          endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        call chksummsk(ficem,ip,1,'ficem')
        call chksummsk(hicem,ip,1,'hicem')
        call chksummsk(hsnwm,ip,1,'hsnwm')
        call chksummsk(iagem,ip,1,'iagem')
      endif
c
      return
      end
