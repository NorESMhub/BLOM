      subroutine intp_mp2ip(idm,jdm,ip,fm,fi)
c
c --- ------------------------------------------------------------------
c --- Transfer field from MICOM p-points to sea ice model p-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,ip(0:idm-1,0:jdm-1)
      real fm(idm,jdm),fi(0:idm-1,0:jdm-1)
c
      integer i,j,i1,j1
c
c$OMP PARALLEL DO PRIVATE(j1,i1)
      do j=1,jdm
        j1=j-1
        do i=1,idm
          i1=i-1
          if (ip(i1,j1).eq.1) then
            fi(i1,j1)=fm(i,j)
          else
            fi(i1,j1)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_mq2iq(idm,jdm,fm,fi)
c
c --- ------------------------------------------------------------------
c --- Transfer field from MICOM q-points to sea ice model q-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm
      real fm(idm,jdm),fi(0:idm-1,0:jdm-1)
c
      integer i,j,i1,j1
c
c$OMP PARALLEL DO PRIVATE(j1,i1)
      do j=1,jdm
        j1=j-1
        do i=1,idm
          i1=i-1
          fi(i1,j1)=fm(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_mu2iq(idm,jdm,iq,fm,fi)
c
c --- ------------------------------------------------------------------
c --- Transfer field from MICOM u-points to sea ice model q-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,iq(idm-1,0:jdm-1)
      real fm(idm,jdm),fi(idm-1,0:jdm-1)
c
      integer i,j,i1,j1,jm1
c
c$OMP PARALLEL DO PRIVATE(j1,jm1,i1)
      do j=1,jdm
        j1=j-1
        jm1=mod(j-2+jdm,jdm)+1
        do i=2,idm
          i1=i-1
          if (iq(i1,j1).eq.1) then
            fi(i1,j1)=.5*(fm(i,jm1)+fm(i,j))
          else
            fi(i1,j1)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_mv2iq(idm,jdm,iq,fm,fi)
c
c --- ------------------------------------------------------------------
c --- Transfer field from MICOM v-points to sea ice model q-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,iq(idm-1,0:jdm-1)
      real fm(idm,jdm),fi(idm-1,0:jdm-1)
c
      integer i,j,i1,j1
c
c$OMP PARALLEL DO PRIVATE(j1,i1)
      do j=1,jdm
        j1=j-1
        do i=2,idm
          i1=i-1
          if (iq(i1,j1).eq.1) then
            fi(i1,j1)=.5*(fm(i-1,j)+fm(i,j))
          else
            fi(i1,j1)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_iq2mp(idm,jdm,mp,fi,fm)
c
c --- ------------------------------------------------------------------
c --- Transfer field from sea ice model q-points to MICOM p-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,mp(idm,jdm)
      real fi(idm-1,0:jdm-1),fm(idm,jdm)
c
      integer i,j,i1,j1,j1p1
c
c$OMP PARALLEL DO PRIVATE(j1,j1p1,i1)
      do j=1,jdm
        j1=j-1
        j1p1=mod(j1+1,jdm)
        do i=1,idm
          i1=i-1
          if (mp(i,j).eq.1) then
            fm(i,j)=.25*(fi(i1  ,j1  )+fi(i1+1,j1  )
     .                  +fi(i1  ,j1p1)+fi(i1+1,j1p1))
          else
            fm(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_iq2mu(idm,jdm,mu,fi,fm)
c
c --- ------------------------------------------------------------------
c --- Transfer field from sea ice model q-points to MICOM u-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,mu(idm,jdm)
      real fi(idm-1,0:jdm-1),fm(idm,jdm)
c
      integer i,j,i1,j1,j1p1
c
c$OMP PARALLEL DO PRIVATE(j1,j1p1,i1)
      do j=1,jdm
        j1=j-1
        j1p1=mod(j1+1,jdm)
        do i=1,idm
          i1=i-1
          if (mu(i,j).eq.1) then
            fm(i,j)=.5*(fi(i1,j1)+fi(i1,j1p1))
          else
            fm(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
c
c --- ------------------------------------------------------------------
c
      subroutine intp_iq2mv(idm,jdm,mv,fi,fm)
c
c --- ------------------------------------------------------------------
c --- Transfer field from sea ice model q-points to MICOM v-points
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer idm,jdm,mv(idm,jdm)
      real fi(idm-1,0:jdm-1),fm(idm,jdm)
c
      integer i,j,i1,j1
c
c$OMP PARALLEL DO PRIVATE(j1,i1)
      do j=1,jdm
        j1=j-1
        do i=1,idm
          i1=i-1
          if (mv(i,j).eq.1) then
            fm(i,j)=.5*(fi(i1,j1)+fi(i1+1,j1))
          else
            fm(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
