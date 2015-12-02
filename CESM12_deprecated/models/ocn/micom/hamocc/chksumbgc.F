      subroutine chksumbgc(a,kd,text)
c
c --- compute checksum of bgc field
c
      use mod_xc
c
      implicit none
c
      real, dimension(idm,jdm,kd) :: a
      integer kd
      character text*(*)
c
      real, dimension(itdm,jtdm,kd) :: ag
      real, dimension(itdm,jtdm) :: rmskg
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: al,rmskl
      integer i,j,k
c
      integer crcfast
      external crcfast
c
      do j=1,jj
        do i=1,ii
          rmskl(i,j)=ip(i,j)
        enddo
      enddo
c
      do k=1,kd
        do j=1,jj
          do i=1,ii
            al(i,j)=a(i,j,k)
          enddo
        enddo
        call xcaget(ag(1,1,k),al,1)
      enddo
      call xcaget(rmskg,rmskl,1)
c
      do k=1,kd
        do j=1,jtdm
          do i=1,itdm
            if (rmskg(i,j).lt..5) then
              ag(i,j,k)=0.
            endif
          enddo
        enddo
      enddo
c
      if (mnproc.eq.1) then
        write (lp,'(3a,z8.8)') ' chksum: ',text,': 0x',
     .    crcfast(ag,itdm*jtdm*kd*8)
      endif
c
      return
      end
